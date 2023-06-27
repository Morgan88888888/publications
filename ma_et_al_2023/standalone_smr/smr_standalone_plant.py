#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
This flowsheet is a standalone SMR plant
"""
import xlsxwriter
import os
from collections import OrderedDict
import logging


import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.infeasible import log_infeasible_constraints


from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state as copy_port_values
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util import model_serializer as ms
import idaes.core.util.scaling as iscale
from idaes.models.properties import iapws95
from idaes.models.properties import swco2
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock)



from gt_flue_gas_ideal import FlueGasParameterBlock
from heat_exchanger_1D_cross_flow import HeatExchangerCrossFlow1D
from water_knockout_condenser import WaterKnockoutCondenser

from idaes.models.unit_models import (
    Mixer,
    MomentumMixingType,
    Heater,
    HeatExchanger,
    PressureChanger,
    Turbine,
    StoichiometricReactor,
    Separator,
    SplittingType,
    Translator,
    Valve,
    ValveFunctionType,
    Compressor)

from idaes.models_extra.power_generation.unit_models import (
    Drum1D,
    Downcomer,
    )
from idaes.models_extra.power_generation.unit_models.helm import (
    ValveFunctionType as HelmValveFunctionType,
    HelmTurbineMultistage,
    HelmNtuCondenser,
    HelmMixer,
    MomentumMixingType,
    HelmValve,
    HelmIsentropicCompressor,
    HelmSplitter,
    )
from evaporator import Evaporator, TubeArrangement
from gibbs_reactor_ideal import GibbsReactor
from idaes.models.unit_models.heat_exchanger import \
    delta_temperature_underwood_callback
from idaes.models.unit_models.pressure_changer import \
    ThermodynamicAssumption

from natural_gas_PR import get_prop, get_rxn
import idaes.logger as idaeslog

__author__ = "Jinliang Ma, Maojian Wang"


def build_blue_hydrogen_plant(m):
    # create property packages
    # water property
    m.fs.water_props = iapws95.Iapws95ParameterBlock()
    m.fs.co2_props = swco2.SWCO2ParameterBlock()
    if not hasattr(m.fs, "syn_props"):
        # syngas property
        syn_config = get_prop(
            components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
        m.fs.syn_props = GenericParameterBlock(**syn_config)

    # water gas shift reaction package
    m.fs.wgs_rxn_props = GenericReactionParameterBlock(
        **get_rxn(
            m.fs.syn_props, reactions=["water_gas_shift_rxn"]))

    # build blue hydrogen plant section
    # natural gas preheater 1
    m.fs.ng_preheat1 = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.syn_props,
                          "has_pressure_change": True}})

    # natural gas preheater 2
    m.fs.ng_preheat2 = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.syn_props,
                          "has_pressure_change": True}})

    # mixer before SMR
    m.fs.smr_tube_mix = Mixer(
        **{"inlet_list": ["fuel_inlet", "steam_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    m.fs.smr_shell_mix = Mixer(
        **{"inlet_list": ["fuel_inlet", "air_inlet", "tailgas_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    # tube side smr reactor
    m.fs.smr_tube = GibbsReactor(
        ** {"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    # shell side smr combustor
    m.fs.smr_shell = GibbsReactor(
        **{"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    # flue gas MP cooler as part of HRSG
    m.fs.flue_gas_mp_cooler = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.water_props,
                          "has_pressure_change": True}})

    # flue gas LP cooler as part of HRSG
    m.fs.flue_gas_lp_cooler = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",                
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.water_props,
                          "has_pressure_change": True}})

    # flue gas mp economizer as part of HRSG
    m.fs.flue_gas_mp_econ = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.water_props,
                          "has_pressure_change": True}})

    # postcombustion co2 capture system
    m.fs.co2_post_separator = Separator(
        **{"outlet_list": ["co2_outlet", "fluegas_outlet", "h2o_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

    # syngas cooler, default is counter-current
    m.fs.syngas_cooler = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.water_props,
                          "has_pressure_change": True}})
    
    # water gas shift reactor
    m.fs.wgsr = StoichiometricReactor(
        **{"has_heat_of_reaction": False,
                 "has_heat_transfer": False,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props,
                 "reaction_package": m.fs.wgs_rxn_props})

    # syngas cooler after shift reactor, default is counter-current
    m.fs.wgsr_cooler = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.water_props,
                          "has_pressure_change": True}})

    # lp pump
    m.fs.lp_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # mp pump
    m.fs.mp_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # mp booster pump after mixed with knockout water
    m.fs.mp_booster_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # lp knockout condenser
    m.fs.lp_knockout_condenser = WaterKnockoutCondenser(
        **{"tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.syn_props,
                 "water_property_package": m.fs.water_props})

    # mp knockout condenser
    m.fs.mp_knockout_condenser = WaterKnockoutCondenser(
        **{"tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.syn_props,
                 "water_property_package": m.fs.water_props})

    # cooling water knockout condenser
    m.fs.cw_knockout_condenser = WaterKnockoutCondenser(
        **{"tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.syn_props,
                 "water_property_package": m.fs.water_props})

    # knockout water mixing with mp water, currently ignore pressure drop of knockout water
    m.fs.mp_mix = HelmMixer(
        **{
            "property_package": m.fs.water_props,
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["mp_inlet", "lpkw_inlet", "mpkw_inlet", "cwkw_inlet"]})

    # precombustion co2 capture system
    m.fs.co2_pre_separator = Separator(
        **{"outlet_list": ["co2_outlet", "syngas_outlet", "h2o_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

    m.fs.co2_mix = Mixer(
        **{"inlet_list": ["pre_inlet", "post_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    m.fs.co2_pre_cool = Heater(**{"property_package": m.fs.syn_props})

    m.fs.co2_comp1 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool1 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp2 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool2 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp3 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool3 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp4 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool4 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp5 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool5 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp6 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool6 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp7 = Compressor(**{"property_package": m.fs.co2_props})

    m.fs.co2_inter_cool7 = Heater(**{"property_package": m.fs.co2_props})

    m.fs.co2_comp8 = Compressor(**{"property_package": m.fs.co2_props})

    # h2 separation system, e.g. PSA
    m.fs.h2_separator = Separator(
        **{"outlet_list": ["h2_outlet", "tailgas_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

    m.fs.h2_comp = Compressor(**{"property_package": m.fs.syn_props})

    # splitter for mp steam, mp_outlet is the extra steam for power generation or polygeneration
    m.fs.mp_splitter = HelmSplitter(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "outlet_list": ["smr_outlet", "mp_outlet"]})

    # property package translater from water_props to syn_props
    m.fs.steam_translator = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.water_props,
                 "outlet_property_package": m.fs.syn_props})

    # additional variables
    # constraints for steam_translator from water_props to syn_props
    @m.fs.steam_translator.Constraint(m.fs.time)
    def steam_translator_T(b, t):
        return b.properties_in[t].temperature == b.properties_out[t].temperature

    @m.fs.steam_translator.Constraint(m.fs.time)
    def steam_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.steam_translator.Constraint(m.fs.time)
    def steam_translator_flow(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.steam_translator.Constraint(m.fs.time, m.fs.syn_props.component_list)
    def steam_translator_x(b, t, j):
        if j=='H2O':
            return b.properties_out[t].mole_frac_comp[j] == 1
        else:
            return b.properties_out[t].mole_frac_comp[j] == 0

    # property package translater from water_props to syn_props
    m.fs.co2_translator = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.co2_props})

    # additional variables
    # constraints for co2_translator from syn_props to co2_props
    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_T(b, t):
        return b.properties_in[t].temperature == b.properties_out[t].temperature

    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_flow(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    # water shift reaction conversion factor
    m.fs.wgsr.x_co = pyo.Var(initialize=0.9,doc='conversion factor for reaction based on CO')

    # additional constraints
    @m.fs.Constraint(m.fs.config.time)
    def heat_balance_of_smr(b, t):
        return 1e-8*(b.smr_tube.control_volume.heat[t] + b.smr_shell.control_volume.heat[t]) == 0

    @m.fs.Constraint(m.fs.config.time)
    def outlet_temperature_of_smr(b, t):
        return 0.01*(b.smr_tube.outlet.temperature[t] + 50 - b.smr_shell.outlet.temperature[t]) == 0

    # constraint for water gas shift reaction reaction extent based on conversion
    @m.fs.wgsr.Constraint(m.fs.config.time)
    def reaction_extent(b, t):
        # Assume CO is the limiting species and conversion is 90%
        prp = b.control_volume.properties_in[t]
        stc = -m.fs.wgs_rxn_props.rate_reaction_stoichiometry["water_gas_shift_rxn", "Vap", "CO"]
        return b.rate_reaction_extent[t, "water_gas_shift_rxn"] == prp.flow_mol*prp.mole_frac_comp["CO"]/stc*b.x_co

    # constraints for outlet pressure of smr_shell outlet
    @m.fs.smr_shell_mix.Constraint(m.fs.time)
    def smr_shell_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.air_inlet_state[t].pressure

    # constraints for outlet pressure of smr_shell outlet
    @m.fs.smr_tube_mix.Constraint(m.fs.time)
    def smr_tube_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.fuel_inlet_state[t].pressure

    # constraints for mp mixer related to the knockout water
    @m.fs.mp_mix.Constraint(m.fs.time)
    def mp_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.mp_inlet_state[t].pressure

    # constraints for co2 mixer of pre and post co2 streams
    @m.fs.co2_mix.Constraint(m.fs.time)
    def co2_mix_pressure(b, t):
        return b.mixed_state[t].pressure == 1.4e5

    @m.fs.Expression(m.fs.config.time)
    def pre_capture_reboiler_duty(b, t):
        return (1600*1054/0.453592*b.co2_pre_separator.co2_outlet_state[t].flow_mol
                *b.co2_pre_separator.co2_outlet_state[t].mole_frac_comp['CO2']*0.044/1e6)

    @m.fs.Expression(m.fs.config.time)
    def post_capture_reboiler_duty(b, t):
        return (1217.184*1054/0.453592*b.co2_post_separator.co2_outlet_state[t].flow_mol
                *b.co2_post_separator.co2_outlet_state[t].mole_frac_comp['CO2']*0.044/1e6)

    @m.fs.Expression(m.fs.config.time)
    def total_lp_duty_available(b, t):
        return (b.flue_gas_lp_cooler.tube.properties_out[t].flow_mol*
               (b.flue_gas_lp_cooler.tube.properties_out[t].enth_mol-
                b.flue_gas_lp_cooler.tube.properties_out[t].enth_mol_sat_phase["Liq"])/1e6)


    # constraint for reboiler duty
    @m.fs.Constraint(m.fs.config.time)
    def total_reboiler_duty_constraint(b, t):
        return b.pre_capture_reboiler_duty[t] + b.post_capture_reboiler_duty[t] == b.total_lp_duty_available[t]

    @m.fs.Expression(m.fs.config.time)
    def total_co2_compression_work(b, t):
        return (b.co2_comp1.work_mechanical[t]+
                b.co2_comp2.work_mechanical[t]+
                b.co2_comp3.work_mechanical[t]+
                b.co2_comp4.work_mechanical[t]+
                b.co2_comp5.work_mechanical[t]+
                b.co2_comp6.work_mechanical[t]+
                b.co2_comp7.work_mechanical[t]+
                b.co2_comp8.work_mechanical[t])

    # stream connections
    m.fs.ng_preheat1_to_ng_preheat2 = Arc(
        source=m.fs.ng_preheat1.tube_outlet,
        destination=m.fs.ng_preheat2.tube_inlet)

    m.fs.ng_preheat2_to_smr_tube_mix = Arc(
        source=m.fs.ng_preheat2.tube_outlet,
        destination=m.fs.smr_tube_mix.fuel_inlet)

    m.fs.smr_tube_mix_to_smr_tube = Arc(
        source=m.fs.smr_tube_mix.outlet,
        destination=m.fs.smr_tube.inlet)

    m.fs.smr_shell_mix_to_smr_shell = Arc(
        source=m.fs.smr_shell_mix.outlet,
        destination=m.fs.smr_shell.inlet)

    m.fs.smr_shell_to_flue_gas_mp_cooler = Arc(
        source=m.fs.smr_shell.outlet,
        destination=m.fs.flue_gas_mp_cooler.shell_inlet)

    m.fs.flue_gas_mp_cooler_to_ng_preheat2 = Arc(
        source=m.fs.flue_gas_mp_cooler.shell_outlet,
        destination=m.fs.ng_preheat2.shell_inlet)

    m.fs.ng_preheat2_to_flue_gas_lp_cooler = Arc(
        source=m.fs.ng_preheat2.shell_outlet,
        destination=m.fs.flue_gas_lp_cooler.shell_inlet)

    m.fs.flue_gas_lp_cooler_to_flue_gas_mp_econ = Arc(
        source=m.fs.flue_gas_lp_cooler.shell_outlet,
        destination=m.fs.flue_gas_mp_econ.shell_inlet)

    m.fs.flue_gas_mp_econ_to_ng_preheat1 = Arc(
        source=m.fs.flue_gas_mp_econ.shell_outlet,
        destination=m.fs.ng_preheat1.shell_inlet)

    m.fs.ng_preheat1_to_co2_post_separator = Arc(
        source=m.fs.ng_preheat1.shell_outlet,
        destination=m.fs.co2_post_separator.inlet)

    m.fs.smr_tube_to_syngas_cooler = Arc(
        source=m.fs.smr_tube.outlet,
        destination=m.fs.syngas_cooler.shell_inlet)

    m.fs.syngas_cooler_to_wgsr = Arc(
        source=m.fs.syngas_cooler.shell_outlet,
        destination=m.fs.wgsr.inlet)

    m.fs.wgsr_to_wgsr_cooler = Arc(
        source=m.fs.wgsr.outlet,
        destination=m.fs.wgsr_cooler.shell_inlet)

    m.fs.wgsr_cooler_to_lp_knockout_condenser = Arc(
        source=m.fs.wgsr_cooler.shell_outlet,
        destination=m.fs.lp_knockout_condenser.shell_inlet)

    m.fs.lp_knockout_condenser_to_mp_knockout_condenser = Arc(
        source=m.fs.lp_knockout_condenser.shell_outlet,
        destination=m.fs.mp_knockout_condenser.shell_inlet)

    m.fs.mp_knockout_condenser_to_cw_knockout_condenser = Arc(
        source=m.fs.mp_knockout_condenser.shell_outlet,
        destination=m.fs.cw_knockout_condenser.shell_inlet)

    m.fs.cw_knockout_condenser_to_co2_pre_separator = Arc(
        source=m.fs.cw_knockout_condenser.shell_outlet,
        destination=m.fs.co2_pre_separator.inlet)

    m.fs.co2_pre_separator_to_co2_mix = Arc(
        source=m.fs.co2_pre_separator.co2_outlet,
        destination=m.fs.co2_mix.pre_inlet)

    m.fs.co2_post_separator_to_co2_mix = Arc(
        source=m.fs.co2_post_separator.co2_outlet,
        destination=m.fs.co2_mix.post_inlet)

    m.fs.co2_mix_to_co2_pre_cool = Arc(
        source=m.fs.co2_mix.outlet,
        destination=m.fs.co2_pre_cool.inlet)

    m.fs.co2_pre_cool_to_co2_translator = Arc(
        source=m.fs.co2_pre_cool.outlet,
        destination=m.fs.co2_translator.inlet)

    m.fs.co2_translator_to_co2_comp1 = Arc(
        source=m.fs.co2_translator.outlet,
        destination=m.fs.co2_comp1.inlet)

    m.fs.co2_comp1_to_co2_inter_cool1 = Arc(
        source=m.fs.co2_comp1.outlet,
        destination=m.fs.co2_inter_cool1.inlet)

    m.fs.co2_inter_cool1_to_co2_comp2 = Arc(
        source=m.fs.co2_inter_cool1.outlet,
        destination=m.fs.co2_comp2.inlet)

    m.fs.co2_comp2_to_co2_inter_cool2 = Arc(
        source=m.fs.co2_comp2.outlet,
        destination=m.fs.co2_inter_cool2.inlet)

    m.fs.co2_inter_cool2_to_co2_comp3 = Arc(
        source=m.fs.co2_inter_cool2.outlet,
        destination=m.fs.co2_comp3.inlet)

    m.fs.co2_comp3_to_co2_inter_cool3 = Arc(
        source=m.fs.co2_comp3.outlet,
        destination=m.fs.co2_inter_cool3.inlet)

    m.fs.co2_inter_cool3_to_co2_comp4 = Arc(
        source=m.fs.co2_inter_cool3.outlet,
        destination=m.fs.co2_comp4.inlet)

    m.fs.co2_comp4_to_co2_inter_cool4 = Arc(
        source=m.fs.co2_comp4.outlet,
        destination=m.fs.co2_inter_cool4.inlet)

    m.fs.co2_inter_cool4_to_co2_comp5 = Arc(
        source=m.fs.co2_inter_cool4.outlet,
        destination=m.fs.co2_comp5.inlet)

    m.fs.co2_comp5_to_co2_inter_cool5 = Arc(
        source=m.fs.co2_comp5.outlet,
        destination=m.fs.co2_inter_cool5.inlet)

    m.fs.co2_inter_cool5_to_co2_comp6 = Arc(
        source=m.fs.co2_inter_cool5.outlet,
        destination=m.fs.co2_comp6.inlet)

    m.fs.co2_comp6_to_co2_inter_cool6 = Arc(
        source=m.fs.co2_comp6.outlet,
        destination=m.fs.co2_inter_cool6.inlet)

    m.fs.co2_inter_cool6_to_co2_comp7 = Arc(
        source=m.fs.co2_inter_cool6.outlet,
        destination=m.fs.co2_comp7.inlet)

    m.fs.co2_comp7_to_co2_inter_cool7 = Arc(
        source=m.fs.co2_comp7.outlet,
        destination=m.fs.co2_inter_cool7.inlet)

    m.fs.co2_inter_cool7_to_co2_comp8 = Arc(
        source=m.fs.co2_inter_cool7.outlet,
        destination=m.fs.co2_comp8.inlet)

    m.fs.co2_pre_separator_to_h2_separator = Arc(
        source=m.fs.co2_pre_separator.syngas_outlet,
        destination=m.fs.h2_separator.inlet)

    m.fs.h2_separator_to_h2_comp = Arc(
        source=m.fs.h2_separator.h2_outlet,
        destination=m.fs.h2_comp.inlet)

    m.fs.h2_separator_to_smr_shell_mix = Arc(
        source=m.fs.h2_separator.tailgas_outlet,
        destination=m.fs.smr_shell_mix.tailgas_inlet)

    m.fs.lp_pump_to_lp_knockout_condenser = Arc(
        source=m.fs.lp_pump.outlet,
        destination=m.fs.lp_knockout_condenser.tube_inlet)

    m.fs.lp_knockout_condenser_to_wgsr_cooler = Arc(
        source=m.fs.lp_knockout_condenser.tube_outlet,
        destination=m.fs.wgsr_cooler.tube_inlet)

    m.fs.wgsr_cooler_to_syngas_cooler = Arc(
        source=m.fs.wgsr_cooler.tube_outlet,
        destination=m.fs.syngas_cooler.tube_inlet)

    m.fs.syngas_cooler_to_flue_gas_lp_cooler = Arc(
        source=m.fs.syngas_cooler.tube_outlet,
        destination=m.fs.flue_gas_lp_cooler.tube_inlet)

    m.fs.mp_pump_to_mp_knockout_condenser = Arc(
        source=m.fs.mp_pump.outlet,
        destination=m.fs.mp_knockout_condenser.tube_inlet)

    m.fs.mp_knockout_condenser_to_mp_mix = Arc(
        source=m.fs.mp_knockout_condenser.tube_outlet,
        destination=m.fs.mp_mix.mp_inlet)

    m.fs.mp_knockout_water_to_mp_mix = Arc(
        source=m.fs.mp_knockout_condenser.shell_water_outlet,
        destination=m.fs.mp_mix.mpkw_inlet)

    m.fs.lp_knockout_water_to_mp_mix = Arc(
        source=m.fs.lp_knockout_condenser.shell_water_outlet,
        destination=m.fs.mp_mix.lpkw_inlet)

    m.fs.cw_knockout_water_to_mp_mix = Arc(
        source=m.fs.cw_knockout_condenser.shell_water_outlet,
        destination=m.fs.mp_mix.cwkw_inlet)

    m.fs.mp_mix_to_mp_booster_pump = Arc(
        source=m.fs.mp_mix.outlet,
        destination=m.fs.mp_booster_pump.inlet)

    m.fs.mp_booster_pump_to_flue_gas_mp_econ = Arc(
        source=m.fs.mp_booster_pump.outlet,
        destination=m.fs.flue_gas_mp_econ.tube_inlet)

    m.fs.flue_gas_mp_econ_to_flue_gas_mp_cooler = Arc(
        source=m.fs.flue_gas_mp_econ.tube_outlet,
        destination=m.fs.flue_gas_mp_cooler.tube_inlet)

    m.fs.flue_gas_mp_cooler_to_mp_splitter = Arc(
        source=m.fs.flue_gas_mp_cooler.tube_outlet,
        destination=m.fs.mp_splitter.inlet)

    m.fs.mp_splitter_to_steam_translator = Arc(
        source=m.fs.mp_splitter.smr_outlet,
        destination=m.fs.steam_translator.inlet)

    m.fs.steam_translator_to_smr_tube_mix = Arc(
        source=m.fs.steam_translator.outlet,
        destination=m.fs.smr_tube_mix.steam_inlet)


def set_blue_hydrogen_plant_inputs(m):
    # input data related to unit operations
    # ng preheat1
    m.fs.ng_preheat1.tube.deltaP.fix(-1000)
    m.fs.ng_preheat1.shell.deltaP.fix(-500)
    m.fs.ng_preheat1.area.fix(10000) #7000
    m.fs.ng_preheat1.overall_heat_transfer_coefficient.fix(50)
    # ng preheat2
    m.fs.ng_preheat2.tube.deltaP.fix(-1000)
    m.fs.ng_preheat2.shell.deltaP.fix(-500)
    m.fs.ng_preheat2.area.fix(8000) #5000
    m.fs.ng_preheat2.overall_heat_transfer_coefficient.fix(50)
    # smr reformer tube side
    m.fs.smr_tube.deltaP.fix(-7e4)
    # smr reformer shell side
    m.fs.smr_shell.deltaP.fix(0)
    # syngas cooler
    m.fs.syngas_cooler.tube.deltaP.fix(-10000)
    m.fs.syngas_cooler.shell.deltaP.fix(-30000)
    m.fs.syngas_cooler.area.fix(5000)
    m.fs.syngas_cooler.overall_heat_transfer_coefficient.fix(80)
    # water gas shift reactor
    m.fs.wgsr.x_co.fix(0.96) #0.95
    m.fs.wgsr.deltaP.fix(-3.2e5)
    # water gas shift reactor cooler
    m.fs.wgsr_cooler.tube.deltaP.fix(-10000)
    m.fs.wgsr_cooler.shell.deltaP.fix(-30000)
    m.fs.wgsr_cooler.area.fix(3000)
    m.fs.wgsr_cooler.overall_heat_transfer_coefficient.fix(100)
    # lp pump for returned process water
    m.fs.lp_pump.efficiency_isentropic.fix(0.80)
    m.fs.lp_pump.deltaP.fix(0.35e5) # from 4e5 to 4.2e5
    # mp makeup water pump
    m.fs.mp_pump.efficiency_isentropic.fix(0.80)
    m.fs.mp_pump.deltaP.fix(2.9e6) # from 1e5 to 3.0e6
    # mp booster pump
    m.fs.mp_booster_pump.efficiency_isentropic.fix(0.80)
    m.fs.mp_booster_pump.deltaP.fix(2e5) # from 3.0e6 to 3.2e6
    # lp water knockout condenser
    m.fs.lp_knockout_condenser.deltaP_tube.fix(-10000)
    m.fs.lp_knockout_condenser.deltaP_shell.fix(-10000)
    m.fs.lp_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.lp_knockout_condenser.area.fix(12000) #10000
    # mp water knockout condenser
    m.fs.mp_knockout_condenser.deltaP_tube.fix(-10000)
    m.fs.mp_knockout_condenser.deltaP_shell.fix(-10000)
    m.fs.mp_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.mp_knockout_condenser.area.fix(12000) #10000
    # cooling water knockout condenser
    m.fs.cw_knockout_condenser.deltaP_tube.fix(-2000)
    m.fs.cw_knockout_condenser.deltaP_shell.fix(0)
    m.fs.cw_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.cw_knockout_condenser.area.fix(12000) #10000
    # CO2 pre separator
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'H2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'CO'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'H2O'].fix(0.2)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'CO2'].fix(0.97)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'CH4'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'N2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'O2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'co2_outlet', 'Ar'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'H2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'CO'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'H2O'].fix(0.8)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'CO2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'CH4'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'N2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'O2'].fix(0)
    m.fs.co2_pre_separator.split_fraction[:, 'h2o_outlet', 'Ar'].fix(0)

    # CO2 post separator
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'H2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'CO'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'H2O'].fix(0.005)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'CO2'].fix(0.97)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'CH4'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'N2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'O2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'Ar'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'H2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'CO'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'H2O'].fix(0.98)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'CO2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'CH4'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'N2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'O2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'Ar'].fix(0)
    # co2_pre_cool
    #m.fs.co2_pre_cool.heat_duty.fix(-1e7)
    m.fs.co2_pre_cool.outlet.temperature.fix(303)
    # pre-combustion CO2 compressor1
    m.fs.co2_comp1.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp2.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp3.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp4.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp5.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp6.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp7.efficiency_isentropic.fix(0.85)
    m.fs.co2_comp8.efficiency_isentropic.fix(0.85)
    pout0 = 0.14e6
    pout1 = 0.32e6
    pout2 = 0.7e6
    pout3 = 1.38e6
    pout4 = 2.54e6
    pout5 = 4e6
    pout6 = 6.38e6
    pout7 = 9.91e6
    pout8 = 15.28e6
    m.fs.co2_comp1.ratioP.fix(pout1/pout0)
    m.fs.co2_comp2.ratioP.fix(pout2/pout1)
    m.fs.co2_comp3.ratioP.fix(pout3/pout2)
    m.fs.co2_comp4.ratioP.fix(pout4/pout3)
    m.fs.co2_comp5.ratioP.fix(pout5/pout4)
    m.fs.co2_comp6.ratioP.fix(pout6/pout5)
    m.fs.co2_comp7.ratioP.fix(pout7/pout6)
    m.fs.co2_comp8.ratioP.fix(pout8/pout7)
    # inter cooler
    enth1 = swco2.htpx(T=303*pyo.units.K, P=pout1*pyo.units.Pa)
    enth2 = swco2.htpx(T=303*pyo.units.K, P=pout2*pyo.units.Pa)
    enth3 = swco2.htpx(T=303*pyo.units.K, P=pout3*pyo.units.Pa)
    enth4 = swco2.htpx(T=303*pyo.units.K, P=pout4*pyo.units.Pa)
    enth5 = swco2.htpx(T=303*pyo.units.K, P=pout5*pyo.units.Pa)
    enth6 = swco2.htpx(T=305*pyo.units.K, P=pout6*pyo.units.Pa)
    enth7 = swco2.htpx(T=307*pyo.units.K, P=pout7*pyo.units.Pa)
    m.fs.co2_inter_cool1.outlet.enth_mol.fix(enth1)
    m.fs.co2_inter_cool2.outlet.enth_mol.fix(enth2)
    m.fs.co2_inter_cool3.outlet.enth_mol.fix(enth3)
    m.fs.co2_inter_cool4.outlet.enth_mol.fix(enth4)
    m.fs.co2_inter_cool5.outlet.enth_mol.fix(enth5)
    m.fs.co2_inter_cool6.outlet.enth_mol.fix(enth6)
    m.fs.co2_inter_cool7.outlet.enth_mol.fix(enth7)
    pin_h2 = 2.448e6
    pout_h2 = 6.48e6
    m.fs.h2_comp.ratioP.fix(pout_h2/pin_h2)
    m.fs.h2_comp.efficiency_isentropic.fix(0.85)

    # H2 separator
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'H2'].fix(0.85)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO'].fix(1e-19)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'H2O'].fix(1e-19)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO2'].fix(1e-19)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CH4'].fix(1e-19)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'N2'].fix(0.03)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'O2'].fix(1e-19)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'Ar'].fix(0.05)
    # mp steam splitter
    m.fs.mp_splitter.split_fraction[:,'smr_outlet'].fix(1)
    # flue gas mp cooler
    m.fs.flue_gas_mp_cooler.tube.deltaP.fix(-5000)
    m.fs.flue_gas_mp_cooler.shell.deltaP.fix(-500)
    m.fs.flue_gas_mp_cooler.area.fix(7000)
    m.fs.flue_gas_mp_cooler.overall_heat_transfer_coefficient.fix(80)
    # flue gas lp cooler
    m.fs.flue_gas_lp_cooler.tube.deltaP.fix(-5000)
    m.fs.flue_gas_lp_cooler.shell.deltaP.fix(-500)
    m.fs.flue_gas_lp_cooler.area.fix(1788) #3000
    m.fs.flue_gas_lp_cooler.overall_heat_transfer_coefficient.fix(80)
    # flue gas mp economizer
    m.fs.flue_gas_mp_econ.tube.deltaP.fix(-1000)
    m.fs.flue_gas_mp_econ.shell.deltaP.fix(-500)
    m.fs.flue_gas_mp_econ.area.fix(12000)
    m.fs.flue_gas_mp_econ.overall_heat_transfer_coefficient.fix(80)
    
    # input data or initial guess for streams
    # natural gas feed to ng preheat1
    m.fs.ng_preheat1.tube_inlet.flow_mol.fix(1024.17)  # mol/s
    m.fs.ng_preheat1.tube_inlet.temperature.fix(288)  # K
    m.fs.ng_preheat1.tube_inlet.pressure.fix(2.92e6)  # Pa,
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CH4'].fix(0.974)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO'].fix(1e-19)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO2'].fix(0.01)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2'].fix(1e-19)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2O'].fix(1e-19)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'N2'].fix(0.016)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'O2'].fix(1e-19)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'Ar'].fix(0.000001) # Add 1 ppm to avoid zero Ar in Gibbs 
    # flue gas to ng preheat1 guess
    m.fs.ng_preheat1.shell_inlet.flow_mol[:].value = 6443  # mol/s
    m.fs.ng_preheat1.shell_inlet.temperature[:].value = 414  # K
    m.fs.ng_preheat1.shell_inlet.pressure[:].value = 108000  # Pa,
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CH4'].value = 1e-19
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CO'].value = 1e-19
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CO2'].value = 0.0704
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'H2'].value = 1e-19
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'H2O'].value = 0.1742
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'N2'].value = 0.7135
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'O2'].value = 0.0419
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'Ar'].value = 0.0000001
    # flue gas to ng preheat2 guess
    m.fs.ng_preheat2.shell_inlet.flow_mol[:].value = 6443  # mol/s
    m.fs.ng_preheat2.shell_inlet.temperature[:].value = 574  # K
    m.fs.ng_preheat2.shell_inlet.pressure[:].value = 109500  # Pa,
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CH4'].value = 1e-19
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CO'].value = 1e-19
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CO2'].value = 0.0704
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'H2'].value = 1e-19
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'H2O'].value = 0.1742
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'N2'].value = 0.7135
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'O2'].value = 0.0419
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'Ar'].value = 0.0000001
    # steam feed to tube side of SMR guess
    m.fs.smr_tube_mix.steam_inlet_state[0].flow_mol.value = 2730.6  # mol/s
    m.fs.smr_tube_mix.steam_inlet_state[0].temperature.value = 672  # K
    m.fs.smr_tube_mix.steam_inlet_state[0].pressure.value = 2.9e6  # Pa,
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CH4'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO2'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2O'].value = 1.0
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['N2'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['O2'].value = 1e-19
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['Ar'].value = 1e-19
    # natural gas feed to shell side of SMR
    m.fs.smr_shell_mix.fuel_inlet_state[0].flow_mol.fix(233) #Aspen 185.83  # mol/s
    m.fs.smr_shell_mix.fuel_inlet_state[0].temperature.fix(288)  # K
    m.fs.smr_shell_mix.fuel_inlet_state[0].pressure.fix(1.1e5)  # Pa,
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CH4'].fix(0.974)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO'].fix(1e-19)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO2'].fix(0.01)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2'].fix(1e-19)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2O'].fix(1e-19)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['N2'].fix(0.016)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['O2'].fix(1e-19)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['Ar'].fix(1e-19)
    # air feed to shell side of SMR using flue gas of gas turbine
    # fresh air feed
    m.fs.smr_shell_mix.air_inlet_state[0].flow_mol.fix(5720.1)  # mol/s
    m.fs.smr_shell_mix.air_inlet_state[0].temperature.fix(297)  # K
    m.fs.smr_shell_mix.air_inlet_state[0].pressure.fix(1.1e5)  # Pa,
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CH4'].fix(1e-19)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO'].fix(1e-19)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO2'].fix(0.0003)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2'].fix(1e-19)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2O'].fix(0.0099)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['N2'].fix(0.7732)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['O2'].fix(0.2074)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['Ar'].fix(0.0092)

    # recycled tailgas intial guess
    m.fs.smr_shell_mix.tailgas_inlet.flow_mol[:].value = 767
    m.fs.smr_shell_mix.tailgas_inlet.temperature[:].value = 320
    m.fs.smr_shell_mix.tailgas_inlet.pressure[:].value = 2.448e6
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CH4'].value = 0.284423
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO'].value =  0.033822
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO2'].value = 0.040808
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2'].value = 0.619181
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2O'].value = 1e-19
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'N2'].value = 0.0217634
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'O2'].value = 1e-19
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'Ar'].value = 0.0000013
    # feed water to syngas cooler guess
    m.fs.syngas_cooler.tube_inlet.flow_mol[:].value = 4000  #679792 lb/hr
    m.fs.syngas_cooler.tube_inlet.pressure[:].value = 4e5 #60 psi
    enth_mol_fw = 14506 #iapws95.htpx(T=400*pyo.units.K, P=4e5*pyo.units.Pa) 416.79K x=0.0908
    m.fs.syngas_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # feed water to water gas shift reactor cooler guess
    m.fs.wgsr_cooler.tube_inlet.flow_mol[:].value = 4000  #679792 lb/hr
    m.fs.wgsr_cooler.tube_inlet.pressure[:].value = 4.1e5  #60 psi
    enth_mol_fw3 = iapws95.htpx(T=390*pyo.units.K, P=4.1e5*pyo.units.Pa)
    m.fs.wgsr_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw3
    # water to flue gas mp cooler guess
    m.fs.flue_gas_mp_cooler.tube_inlet.flow_mol[:].value = 2740 #7000
    m.fs.flue_gas_mp_cooler.tube_inlet.pressure[:].value = 3.19e6
    enth_mol_fw = iapws95.htpx(T=400*pyo.units.K, P=3.19e6*pyo.units.Pa)
    m.fs.flue_gas_mp_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # water to flue gas lp cooler guess, actually not needed for initialization
    m.fs.flue_gas_lp_cooler.tube_inlet.flow_mol[:].value = 4000  #679792 lb/hr
    m.fs.flue_gas_lp_cooler.tube_inlet.pressure[:].value = 3.9e5
    enth_mol_fw = 41250 #iapws95.htpx(T=415.85*pyo.units.K, P=3.9e5*pyo.units.Pa) 415.85 K x=0.79
    m.fs.flue_gas_lp_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # lp feed water to lp pump
    m.fs.lp_pump.inlet.flow_mol.fix(4700)  #679792 lb/hr
    m.fs.lp_pump.inlet.pressure.fix(4e5)  #60 psi
    enth_mol_fw2 = iapws95.htpx(T=415*pyo.units.K, P=4e5*pyo.units.Pa)
    m.fs.lp_pump.inlet.enth_mol.fix(enth_mol_fw2)
    # lp knockout condenser heat duty initial guess
    m.fs.lp_knockout_condenser.heat_duty[:].value = 3e7
    # mp makeup water to mp water knockout condenser
    m.fs.mp_pump.inlet.flow_mol.fix(1616.)  #2731 total to SMR
    m.fs.mp_pump.inlet.pressure.fix(1e5)
    enth_mol_fw3 = iapws95.htpx(T=300*pyo.units.K, P=1e5*pyo.units.Pa)
    m.fs.mp_pump.inlet.enth_mol.fix(enth_mol_fw3)
    # mp knockout condenser heat duty initial guess
    m.fs.mp_knockout_condenser.heat_duty[:].value = 1e7
    # cooling water to cooling water knockout condenser
    m.fs.cw_knockout_condenser.tube_inlet.flow_mol.fix(10000) #3000
    m.fs.cw_knockout_condenser.tube_inlet.pressure.fix(1.5e5)
    enth_mol_fw4 = iapws95.htpx(T=293.15*pyo.units.K, P=1.5e5*pyo.units.Pa)
    m.fs.cw_knockout_condenser.tube_inlet.enth_mol.fix(enth_mol_fw4)
    # cooling water knockout condenser heat duty initial guess
    m.fs.cw_knockout_condenser.heat_duty[:].value = 1e7


def scale_blue_hydrogen_plant(m):
    # set syn_props default scaling
    m.fs.syn_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.syn_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.syn_props.set_default_scaling("temperature", 1e-2)
    m.fs.syn_props.set_default_scaling("pressure", 1e-5)
    m.fs.syn_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.syn_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.syn_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.syn_props.set_default_scaling("entr_mol_phase", 1e-1)
    # set lagrange multiplier scaling factor
    iscale.set_scaling_factor(m.fs.smr_tube.lagrange_mult, 1e-4)
    iscale.set_scaling_factor(m.fs.smr_shell.lagrange_mult, 1e-4)

    iscale.set_scaling_factor(m.fs.syngas_cooler.area, 1e-4)
    iscale.set_scaling_factor(m.fs.syngas_cooler.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.syngas_cooler.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.wgsr_cooler.area, 1e-4)
    iscale.set_scaling_factor(m.fs.wgsr_cooler.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.wgsr_cooler.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.mp_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.cw_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_lp_cooler.area, 1e-4)
    iscale.set_scaling_factor(m.fs.flue_gas_lp_cooler.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_lp_cooler.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_cooler.area, 1e-4)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_cooler.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_cooler.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_econ.area, 1e-4)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_econ.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.flue_gas_mp_econ.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ng_preheat1.area, 1e-4)
    iscale.set_scaling_factor(m.fs.ng_preheat1.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ng_preheat1.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ng_preheat2.area, 1e-4)
    iscale.set_scaling_factor(m.fs.ng_preheat2.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ng_preheat2.shell.heat, 1e-7)
    for c in m.fs.smr_tube_mix.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-9)
    for c in m.fs.smr_shell_mix.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)



def initialize_blue_hydrogen_plant(m):
    outlvl = idaeslog.INFO_LOW
    m.fs.ng_preheat1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.ng_preheat2.tube_inlet, m.fs.ng_preheat1.tube_outlet)
    m.fs.ng_preheat2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube_mix.fuel_inlet, m.fs.ng_preheat2.tube_outlet)
    m.fs.smr_tube_mix.initialize(outlvl = outlvl)
    m.fs.smr_shell_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube.inlet, m.fs.smr_tube_mix.outlet)
    m.fs.smr_tube.outlet.temperature.fix(1137.64)
    m.fs.smr_tube.initialize(outlvl = outlvl)
    m.fs.smr_tube.outlet.temperature.unfix()
    copy_port_values(m.fs.smr_shell.inlet, m.fs.smr_shell_mix.outlet)
    m.fs.smr_shell.outlet.temperature.fix(1187.64)
    m.fs.smr_shell.initialize(outlvl = outlvl)
    m.fs.smr_shell.outlet.temperature.unfix()
    copy_port_values(m.fs.syngas_cooler.shell_inlet, m.fs.smr_tube.outlet)
    m.fs.syngas_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.wgsr.inlet, m.fs.syngas_cooler.shell_outlet)
    m.fs.wgsr.initialize(outlvl = outlvl)
    copy_port_values(m.fs.wgsr_cooler.shell_inlet, m.fs.wgsr.outlet)
    m.fs.wgsr_cooler.initialize(outlvl = outlvl)
    m.fs.lp_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.lp_knockout_condenser.tube_inlet, m.fs.lp_pump.outlet)
    copy_port_values(m.fs.lp_knockout_condenser.shell_inlet, m.fs.wgsr_cooler.shell_outlet)
    m.fs.lp_knockout_condenser.initialize(outlvl = outlvl)
    m.fs.mp_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_knockout_condenser.tube_inlet, m.fs.mp_pump.outlet)
    copy_port_values(m.fs.mp_knockout_condenser.shell_inlet, m.fs.lp_knockout_condenser.shell_outlet)
    m.fs.mp_knockout_condenser.initialize(outlvl = outlvl)
    copy_port_values(m.fs.cw_knockout_condenser.shell_inlet, m.fs.mp_knockout_condenser.shell_outlet)
    m.fs.cw_knockout_condenser.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_pre_separator.inlet, m.fs.cw_knockout_condenser.shell_outlet)
    m.fs.co2_pre_separator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.h2_separator.inlet, m.fs.co2_pre_separator.syngas_outlet)
    m.fs.h2_separator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.flue_gas_mp_cooler.shell_inlet, m.fs.smr_shell.outlet)
    m.fs.flue_gas_mp_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.flue_gas_lp_cooler.shell_inlet, m.fs.flue_gas_mp_cooler.shell_outlet)
    copy_port_values(m.fs.flue_gas_lp_cooler.tube_inlet, m.fs.syngas_cooler.tube_outlet)
    m.fs.flue_gas_lp_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_mix.mp_inlet, m.fs.mp_knockout_condenser.tube_outlet)
    copy_port_values(m.fs.mp_mix.mpkw_inlet, m.fs.mp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.lpkw_inlet, m.fs.lp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.cwkw_inlet, m.fs.cw_knockout_condenser.shell_water_outlet)
    m.fs.mp_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_booster_pump.inlet, m.fs.mp_mix.outlet)
    m.fs.mp_booster_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.flue_gas_mp_econ.shell_inlet, m.fs.flue_gas_lp_cooler.shell_outlet)
    copy_port_values(m.fs.flue_gas_mp_econ.tube_inlet, m.fs.mp_booster_pump.outlet)
    m.fs.flue_gas_mp_econ.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_post_separator.inlet, m.fs.ng_preheat1.shell_outlet)
    m.fs.co2_post_separator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_splitter.inlet, m.fs.flue_gas_mp_cooler.tube_outlet)
    m.fs.mp_splitter.initialize(outlvl = outlvl)
    copy_port_values(m.fs.steam_translator.inlet, m.fs.mp_splitter.smr_outlet)
    m.fs.steam_translator.outlet.mole_frac_comp[:,'H2O'].value = 1
    m.fs.steam_translator.initialize(outlvl = outlvl)
    # co2 mixer and compression
    copy_port_values(m.fs.co2_mix.pre_inlet, m.fs.co2_pre_separator.co2_outlet)
    copy_port_values(m.fs.co2_mix.post_inlet, m.fs.co2_post_separator.co2_outlet)
    m.fs.co2_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_pre_cool.inlet, m.fs.co2_mix.outlet)
    m.fs.co2_pre_cool.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_translator.inlet, m.fs.co2_pre_cool.outlet)
    m.fs.co2_translator.initialize(outlvl = outlvl)   
    copy_port_values(m.fs.co2_comp1.inlet, m.fs.co2_translator.outlet)
    m.fs.co2_comp1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool1.inlet, m.fs.co2_comp1.outlet)
    m.fs.co2_inter_cool1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp2.inlet, m.fs.co2_inter_cool1.outlet)
    m.fs.co2_comp2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool2.inlet, m.fs.co2_comp2.outlet)
    m.fs.co2_inter_cool2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp3.inlet, m.fs.co2_inter_cool2.outlet)
    m.fs.co2_comp3.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool3.inlet, m.fs.co2_comp3.outlet)
    m.fs.co2_inter_cool3.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp4.inlet, m.fs.co2_inter_cool3.outlet)
    m.fs.co2_comp4.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool4.inlet, m.fs.co2_comp4.outlet)
    m.fs.co2_inter_cool4.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp5.inlet, m.fs.co2_inter_cool4.outlet)
    m.fs.co2_comp5.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool5.inlet, m.fs.co2_comp5.outlet)
    m.fs.co2_inter_cool5.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp6.inlet, m.fs.co2_inter_cool5.outlet)
    m.fs.co2_comp6.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool6.inlet, m.fs.co2_comp6.outlet)
    m.fs.co2_inter_cool6.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp7.inlet, m.fs.co2_inter_cool6.outlet)
    m.fs.co2_comp7.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_inter_cool7.inlet, m.fs.co2_comp7.outlet)
    m.fs.co2_inter_cool7.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_comp8.inlet, m.fs.co2_inter_cool7.outlet)
    m.fs.co2_comp8.initialize(outlvl = outlvl)



def get_blue_hydrogen_plant_model():
    logging.getLogger('idaes.core.util.scaling').setLevel(logging.ERROR)
    m = pyo.ConcreteModel(name='Standalone SMR plant')
    m.fs = FlowsheetBlock(dynamic = False)
    # solver and options
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 100,
            'bound_push': 1e-16}
    build_blue_hydrogen_plant(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    set_blue_hydrogen_plant_inputs(m)
    scale_blue_hydrogen_plant(m)
    iscale.calculate_scaling_factors(m)
    init_fname = 'smr_standalone_plant.json'
    if os.path.exists(init_fname):
        ms.from_json(m, fname=init_fname, wts=ms.StoreSpec(suffix=False))
    else:
        initialize_blue_hydrogen_plant(m)
        m.fs.smr_shell_mix.air_inlet.fix()
        m.fs.smr_tube_mix.steam_inlet_state[:].flow_mol.fix(2730.6)
        m.fs.mp_pump.inlet.flow_mol.unfix()
        m.fs.total_reboiler_duty_constraint.deactivate()
        solver.solve(m, tee=True)
        m.fs.total_reboiler_duty_constraint.activate()
        m.fs.flue_gas_lp_cooler.area.unfix()
    solver.solve(m, tee=True)
    if not os.path.exists(init_fname):    
        ms.to_json(m, fname=init_fname)
    return m


if __name__ == "__main__":
    m = get_blue_hydrogen_plant_model()
