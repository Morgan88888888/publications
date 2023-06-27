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
This is a flowsheet for NGSC and SMR integrated energy system with CCS
"""
import xlsxwriter
import os
from collections import OrderedDict
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.infeasible import log_infeasible_constraints
from idaes.core.util import model_serializer as ms

# IDAES Imports
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state as copy_port_values
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe

import idaes.core.util.scaling as iscale

from idaes.models.properties import iapws95
from idaes.models.properties import swco2
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock)
from gt_flue_gas_ideal import FlueGasParameterBlock

from water_knockout_condenser import WaterKnockoutCondenser

from idaes.models.unit_models import (
    Mixer,
    MomentumMixingType,
    Heater,
    HeatExchanger,
    PressureChanger,
    Turbine,
    #GibbsReactor,
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

def performance_curves(m, flow_scale=2.0):
    """
    Add head and efficenty curves to the turbine stages in a model.

    Args:
        m: model with expcted structure

    Returns:
        None
    """
    # The flow scale variable less you scale the turbine up or down for differnt
    # nominal flows and differnt full load power outputs.
    m.fs.performance_flow_scale = pyo.Var(initialize=flow_scale)
    m.fs.performance_flow_scale.fix()
    fscale = m.fs.performance_flow_scale

    # Efficiency curves for three stages
    @m.fs.gts1.performance_curve.Constraint(m.fs.config.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            1.4469E-14*f**5 -
            6.3333E-11*f**4 +
            6.6179E-08*f**3 -
            3.1728E-05*f**2 +
            7.7846E-03*f +
            1.0724E-01)
    @m.fs.gts2.performance_curve.Constraint(m.fs.config.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            2.6599E-16*f**5 -
            2.5894E-12*f**4 +
            6.0174E-09*f**3 -
            6.4156E-06*f**2 +
            3.5005E-03*f +
            1.0724E-01)
    @m.fs.gts3.performance_curve.Constraint(m.fs.config.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            5.8407E-18*f**5 -
            1.2203E-13*f**4 +
            6.0863E-10*f**3 -
            1.3927E-06*f**2 +
            1.6310E-03*f +
            1.0724E-01)
    # Head curves for three stages, increased scaling factor by a factor of 2 to reduce outlet pressure
    @m.fs.gts1.performance_curve.Constraint(m.fs.config.time)
    def head_isen_eqn(b, t):
        f = pyo.log(1.3*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -2085.1*f**3 + 38433*f**2 - 150764*f + 422313)
    @m.fs.gts2.performance_curve.Constraint(m.fs.config.time)
    def head_isen_eqn(b, t):
        f = pyo.log(1.3*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1676.3*f**3 + 34916*f**2 - 173801*f + 456957)
    @m.fs.gts3.performance_curve.Constraint(m.fs.config.time)
    def head_isen_eqn(b, t):
        f = pyo.log(1.3*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1373.6*f**3 + 31759*f**2 - 188528*f + 500520)


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
        **{"inlet_list": ["fuel_inlet", "air_inlet", "gt_inlet", "tailgas_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    # tube side smr reactor
    m.fs.smr_tube = GibbsReactor(
        **{"has_heat_transfer": True,
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

    m.fs.fg_splitter = Separator(**{
        "property_package": m.fs.syn_props,
        "outlet_list":["main_outlet", "bypass_outlet"]})


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

    # mixer for 1st stage cooling air
    m.fs.fg_mix = Mixer(**{
        "property_package": m.fs.syn_props,
        "inlet_list":["main_inlet", "bypass_inlet"],
        "momentum_mixing_type": MomentumMixingType.none})
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
        return b.mixed_state[t].pressure == b.gt_inlet_state[t].pressure

    # constraints for outlet pressure of smr_shell outlet
    @m.fs.smr_tube_mix.Constraint(m.fs.time)
    def smr_tube_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.fuel_inlet_state[t].pressure
        
    # constraints for outlet pressure of fg_mix outlet
    @m.fs.fg_mix.Constraint(m.fs.time)
    def fg_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.main_inlet_state[t].pressure

    # constraints for mp mixer related to the knockout water
    @m.fs.mp_mix.Constraint(m.fs.time)
    def mp_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.mp_inlet_state[t].pressure

    # constraints for co2 mixer of pre and post co2 streams
    @m.fs.co2_mix.Constraint(m.fs.time)
    def co2_mix_pressure(b, t):
        return b.mixed_state[t].pressure == 1.4e5

    # constraint for seting smr tube side outlet temperature
    @m.fs.Constraint(m.fs.config.time)
    def outlet_temperature_of_smr_fixed(b, t):
        return 0.001*(b.smr_tube.outlet.temperature[t] - 1155.6) == 0

    # constraint for setting flue gas O2 mole fraction
    @m.fs.Constraint(m.fs.config.time)
    def o2_in_flue_gas(b, t):
        return 100*(b.smr_shell.outlet.mole_frac_comp[t,'O2'] - 0.01740189) == 0  #0.01827295

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

    # constraint for setting flue gas O2 mole fraction
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

    # constraint for heat transfer coefficient as a function of flow rate
    @m.fs.ng_preheat1.Constraint(m.fs.config.time)
    def heat_transfer_coefficient(b, t):
        return b.overall_heat_transfer_coefficient[t] == 50*(b.shell_inlet.flow_mol[t]/8588.1367)

    # constraint for heat transfer coefficient as a function of flow rate
    @m.fs.ng_preheat2.Constraint(m.fs.config.time)
    def heat_transfer_coefficient(b, t):
        return b.overall_heat_transfer_coefficient[t] == 50*(b.shell_inlet.flow_mol[t]/8588.1367)

    # constraint for heat transfer coefficient as a function of flow rate
    @m.fs.flue_gas_mp_cooler.Constraint(m.fs.config.time)
    def heat_transfer_coefficient(b, t):
        return b.overall_heat_transfer_coefficient[t] == 80*(b.shell_inlet.flow_mol[t]/8588.1367)

    # constraint for heat transfer coefficient as a function of flow rate
    @m.fs.flue_gas_mp_econ.Constraint(m.fs.config.time)
    def heat_transfer_coefficient(b, t):
        return b.overall_heat_transfer_coefficient[t] == 80*(b.shell_inlet.flow_mol[t]/8588.1367)

    # stream connections
    m.fs.ng_preheat1_to_ng_preheat2 = Arc(
        source=m.fs.ng_preheat1.tube_outlet,
        destination=m.fs.ng_preheat2.tube_inlet)

    m.fs.ng_preheat2_to_smr_tube_mix = Arc(
        source=m.fs.ng_preheat2.tube_outlet,
        destination=m.fs.smr_tube_mix.fuel_inlet)

    m.fs.air_gas_mix3_to_smr_shell_mix = Arc(
        source=m.fs.air_gas_mix3.outlet,
        destination=m.fs.smr_shell_mix.gt_inlet)

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

    m.fs.ng_preheat2_to_fg_splitter = Arc(
        source=m.fs.ng_preheat2.shell_outlet,
        destination=m.fs.fg_splitter.inlet)

    m.fs.fg_splitter_to_flue_gas_lp_cooler = Arc(
        source=m.fs.fg_splitter.main_outlet,
        destination=m.fs.flue_gas_lp_cooler.shell_inlet)

    m.fs.flue_gas_lp_cooler_to_fg_mix = Arc(
        source=m.fs.flue_gas_lp_cooler.shell_outlet,
        destination=m.fs.fg_mix.main_inlet)

    m.fs.fg_splitter_to_fg_mix = Arc(
        source=m.fs.fg_splitter.bypass_outlet,
        destination=m.fs.fg_mix.bypass_inlet)

    m.fs.fg_mix_to_flue_gas_mp_econ = Arc(
        source=m.fs.fg_mix.outlet,
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


def build_gas_turbine_plant(m):
    '''
    # gas turbine combustion reaction package
    m.fs.gt_rxn_props = GenericReactionParameterBlock(
        default=get_rxn(
            m.fs.syn_props, reactions=["h2_cmb", "co_cmb", "ch4_cmb"]))
    '''

    syn_config = get_prop(
        components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
    m.fs.syn_props = GenericParameterBlock(**syn_config)

    # Compressor inlet guide vane modeled as a valve to control air flow or load
    m.fs.comp_igv = Valve(**{
        "valve_function_callback": ValveFunctionType.linear,
        "property_package": m.fs.syn_props})

    # Air compressor
    m.fs.air_comp = Compressor(**{
        "property_package": m.fs.syn_props,
        "support_isentropic_performance_curves":True})

    # Blade cooling air splitter
    m.fs.air_splitter = Separator(**{
        "property_package": m.fs.syn_props,
        "outlet_list":["comb_outlet", "gts1_outlet", "gts2_outlet", "gts3_outlet"]})

    # 1st cooling air valve for gas turbine stage 1
    m.fs.air_valve1 = Valve(**{
        "valve_function_callback": ValveFunctionType.linear,
        "property_package": m.fs.syn_props})

    # 2nd cooling air valve for gas turbine stage 2
    m.fs.air_valve2 = Valve(**{
        "valve_function_callback": ValveFunctionType.linear,
        "property_package": m.fs.syn_props})

    # 3rd cooling air valve for gas turbine stage 2
    m.fs.air_valve3 = Valve(**{
        "valve_function_callback": ValveFunctionType.linear,
        "property_package": m.fs.syn_props})

    # mixer of fuel and oxidant before gas turbine combustor
    m.fs.gt_mix = Mixer(
        **{"inlet_list": ["fuel_inlet", "air_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    # gas turbine combustor
    m.fs.gt_combustor = GibbsReactor(
        **{"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    # 1st gas turbine stage
    m.fs.gts1 = Turbine(**{
        "property_package": m.fs.syn_props,
        "support_isentropic_performance_curves":True})

    # mixer for 1st stage cooling air
    m.fs.air_gas_mix1 = Mixer(**{
        "property_package": m.fs.syn_props,
        "inlet_list":["gas_inlet", "air_inlet"],
        "momentum_mixing_type": MomentumMixingType.equality})

    # 2nd gas turbine stage
    m.fs.gts2 = Turbine(**{
        "property_package": m.fs.syn_props,
        "support_isentropic_performance_curves":True})

    # mixer for 2nd stage cooling air
    m.fs.air_gas_mix2 = Mixer(**{
        "property_package": m.fs.syn_props,
        "inlet_list":["gas_inlet", "air_inlet"],
        "momentum_mixing_type": MomentumMixingType.equality})

    # 3rd gas turbine stage
    m.fs.gts3 = Turbine(**{
        "property_package": m.fs.syn_props,
        "support_isentropic_performance_curves":True})

    # mixer for 3rd stage cooling air
    m.fs.air_gas_mix3 = Mixer(**{
        "property_package": m.fs.syn_props,
        "inlet_list":["gas_inlet", "air_inlet"],
        "momentum_mixing_type": MomentumMixingType.equality})

    # setup performance curves for three gas turbine stages
    performance_curves(m)

    # additional constraints
    # constraint for the gas turbine mixer outlet pressure
    @m.fs.gt_mix.Constraint(m.fs.time)
    def gt_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.air_inlet_state[t].pressure 

    @m.fs.Expression(m.fs.time)
    def gt_net_power(b, t):
        return b.gts1.work_mechanical[t] + b.gts2.work_mechanical[t] + b.gts3.work_mechanical[t] + b.air_comp.work_mechanical[t]

    # constraint for exit temperature by adjust the air flow rate
    @m.fs.Constraint(m.fs.time)
    def gt_exit_temp_constraint(b, t):
        return b.air_gas_mix3.outlet.temperature[t] == 953.6     

    # constraint for exit pressure by adjusting inlet guide vane opening
    @m.fs.Constraint(m.fs.time)
    def gt_exit_pres_constraint(b, t):
        return b.air_gas_mix3.outlet.pressure[t] == 103046    

    # stream connections
    m.fs.comp_igv_to_air_comp = Arc(
        source=m.fs.comp_igv.outlet,
        destination=m.fs.air_comp.inlet)

    m.fs.air_comp_to_air_splitter = Arc(
        source=m.fs.air_comp.outlet,
        destination=m.fs.air_splitter.inlet)
    
    m.fs.air_splitter_to_gt_mix = Arc(
        source=m.fs.air_splitter.comb_outlet,
        destination=m.fs.gt_mix.air_inlet)

    m.fs.gt_mix_to_gt_combustor = Arc(
        source=m.fs.gt_mix.outlet,
        destination=m.fs.gt_combustor.inlet)

    m.fs.gt_combustor_to_gts1 = Arc(
        source=m.fs.gt_combustor.outlet,
        destination=m.fs.gts1.inlet)

    m.fs.air_splitter_to_air_valve1 = Arc(
        source=m.fs.air_splitter.gts1_outlet,
        destination=m.fs.air_valve1.inlet)

    m.fs.air_valve1_to_air_gas_mix1 = Arc(
        source=m.fs.air_valve1.outlet,
        destination=m.fs.air_gas_mix1.air_inlet)

    m.fs.gts1_to_air_gas_mix1 = Arc(
        source=m.fs.gts1.outlet,
        destination=m.fs.air_gas_mix1.gas_inlet)

    m.fs.air_gas_mix1_to_gts2 = Arc(
        source=m.fs.air_gas_mix1.outlet,
        destination=m.fs.gts2.inlet)

    m.fs.air_splitter_to_air_valve2 = Arc(
        source=m.fs.air_splitter.gts2_outlet,
        destination=m.fs.air_valve2.inlet)

    m.fs.air_valve2_to_air_gas_mix2 = Arc(
        source=m.fs.air_valve2.outlet,
        destination=m.fs.air_gas_mix2.air_inlet)

    m.fs.gts2_to_air_gas_mix2 = Arc(
        source=m.fs.gts2.outlet,
        destination=m.fs.air_gas_mix2.gas_inlet)

    m.fs.air_gas_mix2_to_gts3 = Arc(
        source=m.fs.air_gas_mix2.outlet,
        destination=m.fs.gts3.inlet)

    m.fs.air_splitter_to_air_valve3 = Arc(
        source=m.fs.air_splitter.gts3_outlet,
        destination=m.fs.air_valve3.inlet)

    m.fs.air_valve3_to_air_gas_mix3 = Arc(
        source=m.fs.air_valve3.outlet,
        destination=m.fs.air_gas_mix3.air_inlet)

    m.fs.gts3_to_air_gas_mix3 = Arc(
        source=m.fs.gts3.outlet,
        destination=m.fs.air_gas_mix3.gas_inlet)

def set_blue_hydrogen_plant_inputs(m):
    # input data related to unit operations
    # ng preheat1
    m.fs.ng_preheat1.tube.deltaP.fix(-1000)
    m.fs.ng_preheat1.shell.deltaP.fix(-500)
    m.fs.ng_preheat1.area.fix(7000)#5000
    #m.fs.ng_preheat1.overall_heat_transfer_coefficient.fix(50)
    # ng preheat2
    m.fs.ng_preheat2.tube.deltaP.fix(-1000)
    m.fs.ng_preheat2.shell.deltaP.fix(-500)
    m.fs.ng_preheat2.area.fix(6500)#5000
    #m.fs.ng_preheat2.overall_heat_transfer_coefficient.fix(50)
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
    m.fs.wgsr.x_co.fix(0.96)
    m.fs.wgsr.deltaP.fix(-3.2e5)
    # water gas shift reactor cooler
    m.fs.wgsr_cooler.tube.deltaP.fix(-10000)
    m.fs.wgsr_cooler.shell.deltaP.fix(-30000)
    m.fs.wgsr_cooler.area.fix(3000)
    m.fs.wgsr_cooler.overall_heat_transfer_coefficient.fix(100)
    # lp pump for returned process water
    m.fs.lp_pump.efficiency_isentropic.fix(0.80)
    m.fs.lp_pump.deltaP.fix(0.35e5) # from 1e5 to 4.2e5
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
    m.fs.cw_knockout_condenser.area.fix(12000)
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
    # co2_pre_cool(-1e7)
    m.fs.co2_pre_cool.heat_duty.fix(-1e6)
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
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'H2O'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO2'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CH4'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'N2'].fix(0.03)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'O2'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'Ar'].fix(0.05)
    # mp steam splitter
    m.fs.mp_splitter.split_fraction[:,'smr_outlet'].fix(1)
    # flue gas mp cooler
    m.fs.flue_gas_mp_cooler.tube.deltaP.fix(-5000)
    m.fs.flue_gas_mp_cooler.shell.deltaP.fix(-500)
    m.fs.flue_gas_mp_cooler.area.fix(13000) #10000
    #m.fs.flue_gas_mp_cooler.overall_heat_transfer_coefficient.fix(80)
    # fg_splitter
    m.fs.fg_splitter.split_fraction[:,'main_outlet'].fix(0.7) #set to a lower value at full load
    # flue gas lp cooler
    m.fs.flue_gas_lp_cooler.tube.deltaP.fix(-5000)
    m.fs.flue_gas_lp_cooler.shell.deltaP.fix(-500)
    m.fs.flue_gas_lp_cooler.area.fix(7000)#2974
    m.fs.flue_gas_lp_cooler.overall_heat_transfer_coefficient.fix(80)
    # flue gas mp economizer
    m.fs.flue_gas_mp_econ.tube.deltaP.fix(-1000)
    m.fs.flue_gas_mp_econ.shell.deltaP.fix(-500)
    m.fs.flue_gas_mp_econ.area.fix(14000) #15000
    #m.fs.flue_gas_mp_econ.overall_heat_transfer_coefficient.fix(80)

    
    # input data or initial guess for streams
    # natural gas feed to ng preheat1
    m.fs.ng_preheat1.tube_inlet.flow_mol.fix(1024.17)  # mol/s
    m.fs.ng_preheat1.tube_inlet.temperature.fix(288)  # K
    m.fs.ng_preheat1.tube_inlet.pressure.fix(2.92e6)  # Pa,
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CH4'].fix(0.974)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO'].fix(1e-8)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO2'].fix(0.01)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2'].fix(1e-8)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2O'].fix(1e-8)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'N2'].fix(0.016)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'O2'].fix(1e-8)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'Ar'].fix(0.000001) # Add 1 ppm to avoid zero Ar in Gibbs 
    # flue gas to ng preheat1 guess
    m.fs.ng_preheat1.shell_inlet.flow_mol[:].value = 8602  # mol/s
    m.fs.ng_preheat1.shell_inlet.temperature[:].value = 410 #533.8  # K
    m.fs.ng_preheat1.shell_inlet.pressure[:].value = 101040  # Pa,
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CH4'].value = 1e-8
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CO'].value = 1e-8
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'CO2'].value = 0.076598
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'H2'].value = 1e-8
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'H2O'].value = 0.20841
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'N2'].value = 0.69812
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'O2'].value = 0.01687
    m.fs.ng_preheat1.shell_inlet.mole_frac_comp[:,'Ar'].value = 0.0000001
    # flue gas to ng preheat2 guess
    m.fs.ng_preheat2.shell_inlet.flow_mol[:].value = 8602  # mol/s
    m.fs.ng_preheat2.shell_inlet.temperature[:].value = 700  # K
    m.fs.ng_preheat2.shell_inlet.pressure[:].value = 102540  # Pa,
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CH4'].value = 1e-8
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CO'].value = 1e-8
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'CO2'].value = 0.076598
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'H2'].value = 1e-8
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'H2O'].value = 0.20841
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'N2'].value = 0.69812
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'O2'].value = 0.01687
    m.fs.ng_preheat2.shell_inlet.mole_frac_comp[:,'Ar'].value = 0.0000001
    # steam feed to tube side of SMR guess
    m.fs.smr_tube_mix.steam_inlet_state[0].flow_mol.value = 2730.6  # mol/s
    m.fs.smr_tube_mix.steam_inlet_state[0].temperature.value = 1018.5  # K
    m.fs.smr_tube_mix.steam_inlet_state[0].pressure.value = 3.184e6  # Pa,
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CH4'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO2'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2O'].value = 1.0
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['N2'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['O2'].value = 1e-8
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['Ar'].value = 1e-8
    # natural gas feed to shell side of SMR
    m.fs.smr_shell_mix.fuel_inlet_state[0].flow_mol.fix(53)
    m.fs.smr_shell_mix.fuel_inlet_state[0].temperature.fix(288)  # K
    m.fs.smr_shell_mix.fuel_inlet_state[0].pressure.fix(1.1e5)  # Pa,
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CH4'].fix(0.974)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO'].fix(1e-8)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO2'].fix(0.01)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2'].fix(1e-8)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2O'].fix(1e-8)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['N2'].fix(0.016)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['O2'].fix(1e-8)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['Ar'].fix(1e-8)
    # air feed to shell side of SMR using flue gas of gas turbine
    # fresh air feed
    m.fs.smr_shell_mix.air_inlet_state[0].flow_mol.fix(0.001)  # mol/s
    m.fs.smr_shell_mix.air_inlet_state[0].temperature.fix(298.15)  # K
    m.fs.smr_shell_mix.air_inlet_state[0].pressure.fix(1.1e5)  # Pa,
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CH4'].fix(1e-8)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO'].fix(1e-8)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO2'].fix(0.0003)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2'].fix(1e-8)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2O'].fix(0.0099)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['N2'].fix(0.7732)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['O2'].fix(0.2074)
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['Ar'].fix(0.0092)
    # flue gas from gas turbine
    m.fs.smr_shell_mix.gt_inlet_state[0].flow_mol.value = 8054.6  # mol/s
    m.fs.smr_shell_mix.gt_inlet_state[0].temperature.value = 953.6  # K
    m.fs.smr_shell_mix.gt_inlet_state[0].pressure.value = 1.0304e5  # Pa,
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['CH4'].value = 1e-8
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['CO'].value = 3.5e-6
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['CO2'].value = 0.04471
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['H2'].value = 2.25e-6
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['H2O'].value = 0.10283
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['N2'].value = 0.743438
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['O2'].value = 0.10901
    m.fs.smr_shell_mix.gt_inlet_state[0].mole_frac_comp['Ar'].value = 1e-8
    # recycled tailgas intial guess
    m.fs.smr_shell_mix.tailgas_inlet.flow_mol[:].value = 750.8
    m.fs.smr_shell_mix.tailgas_inlet.temperature[:].value = 341.5
    m.fs.smr_shell_mix.tailgas_inlet.pressure[:].value = 2.448e6
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CH4'].value = 0.25448
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO'].value =  0.035786
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO2'].value = 0.042077
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2'].value = 0.639076
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2O'].value = 0.00675
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'N2'].value = 0.021824
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'O2'].value = 1e-8
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'Ar'].value = 0.0000013
    # feed water to syngas cooler guess
    m.fs.syngas_cooler.tube_inlet.flow_mol[:].value = 4758  #679792 lb/hr
    m.fs.syngas_cooler.tube_inlet.pressure[:].value = 4e5 #60 psi
    enth_mol_fw = 17976
    m.fs.syngas_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # feed water to water gas shift reactor cooler guess
    m.fs.wgsr_cooler.tube_inlet.flow_mol[:].value = 4758  #679792 lb/hr
    m.fs.wgsr_cooler.tube_inlet.pressure[:].value = 4.1e5  #60 psi
    enth_mol_fw3 = 11868 #iapws95.htpx(T=390*pyo.units.K, P=4.1e5*pyo.units.Pa)
    m.fs.wgsr_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw3
    # water to flue gas mp cooler guess
    m.fs.flue_gas_mp_cooler.tube_inlet.flow_mol[:].value = 2774
    m.fs.flue_gas_mp_cooler.tube_inlet.pressure[:].value = 3.189e6
    enth_mol_fw = 16325 #iapws95.htpx(T=484.9*pyo.units.K, P=3.189e6*pyo.units.Pa)
    m.fs.flue_gas_mp_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # water to flue gas lp cooler guess, actually not needed for initialization
    m.fs.flue_gas_lp_cooler.tube_inlet.flow_mol[:].value = 4758  #679792 lb/hr
    m.fs.flue_gas_lp_cooler.tube_inlet.pressure[:].value = 3.9e5
    enth_mol_fw = 43780
    m.fs.flue_gas_lp_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
   # lp feed water to lp pump
    m.fs.lp_pump.inlet.flow_mol.fix(4850)  #679792 lb/hr
    m.fs.lp_pump.inlet.pressure.fix(4e5)  #60 psi
    enth_mol_fw2 = iapws95.htpx(T=415*pyo.units.K, P=4e5*pyo.units.Pa)
    m.fs.lp_pump.inlet.enth_mol.fix(enth_mol_fw2)
    # lp knockout condenser heat duty initial guess
    m.fs.lp_knockout_condenser.heat_duty[:].value = 4.6e7
    # mp makeup water to mp water knockout condenser
    m.fs.mp_pump.inlet.flow_mol.fix(1670)  #2731 total to SMR
    m.fs.mp_pump.inlet.pressure.fix(1e5)
    enth_mol_fw3 = iapws95.htpx(T=300*pyo.units.K, P=1e5*pyo.units.Pa)
    m.fs.mp_pump.inlet.enth_mol.fix(enth_mol_fw3)
    # mp knockout condenser heat duty initial guess
    m.fs.mp_knockout_condenser.heat_duty[:].value = 1e7
    # cooling water to cooling water knockout condenser
    m.fs.cw_knockout_condenser.tube_inlet.flow_mol.fix(10000)
    m.fs.cw_knockout_condenser.tube_inlet.pressure.fix(1.5e5)
    enth_mol_fw4 = iapws95.htpx(T=293.15*pyo.units.K, P=1.5e5*pyo.units.Pa)
    m.fs.cw_knockout_condenser.tube_inlet.enth_mol.fix(enth_mol_fw4)
    # cooling water knockout condenser heat duty initial guess
    m.fs.cw_knockout_condenser.heat_duty[:].value = 4.9e6


def set_gas_turbine_plant_inputs(m):
    # set flow scale of performance curves
    m.fs.performance_flow_scale.fix(2.0/115*250)
    # air to gas turbine
    m.fs.comp_igv.inlet.flow_mol.fix(16720.0*115/250)  # 16500 mol/s
    m.fs.comp_igv.inlet.temperature.fix(298.15)  # K
    m.fs.comp_igv.inlet.pressure.fix(101325)  # Pa,
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CH4'].fix(1e-8)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CO'].fix(1e-8)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CO2'].fix(0.000335)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'H2'].fix(1e-8)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'H2O'].fix(0.015653)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'N2'].fix(0.777811)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'O2'].fix(0.206201)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'Ar'].fix(1e-8)

    m.fs.gt_mix.fuel_inlet.flow_mol.fix(790.0*115/250)
    m.fs.gt_mix.fuel_inlet.temperature.fix(298.15)
    m.fs.gt_mix.fuel_inlet.pressure.fix(2e6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CH4'].fix(0.974)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CO'].fix(1e-8)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CO2'].fix(0.01)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['H2'].fix(1e-8)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['H2O'].fix(1e-8)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['N2'].fix(0.016)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['O2'].fix(1e-8)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['Ar'].fix(1e-8)

    # compressor inlet guide vanes
    m.fs.comp_igv.deltaP.fix(-5000)
    m.fs.comp_igv.valve_opening.fix(0.9)
    m.fs.comp_igv.Cv.unfix()
    #m.fs.comp_igv.pressure_flow_equation.deactivate()
    # air compressor
    m.fs.air_comp.efficiency_isentropic.fix(0.9)
    m.fs.air_comp.ratioP.fix(18.9343)
    # air_splitter
    m.fs.air_splitter.split_fraction[:,'gts1_outlet'].fix(0.04)
    m.fs.air_splitter.split_fraction[:,'gts2_outlet'].fix(0.02)
    m.fs.air_splitter.split_fraction[:,'gts3_outlet'].fix(0.01)
    # split valves
    m.fs.air_valve1.valve_opening.fix(0.7)
    m.fs.air_valve2.valve_opening.fix(0.7)
    m.fs.air_valve3.valve_opening.fix(0.7)
    m.fs.air_valve1.Cv.unfix()
    m.fs.air_valve2.Cv.unfix()
    m.fs.air_valve3.Cv.unfix()
    # gas turbine combustor pressure drop
    m.fs.gt_combustor.deltaP.fix(-100000)
    m.fs.gt_combustor.heat_duty.fix(0)

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

def scale_gas_turbine_plant(m):
    # set gas turbine flow scaling
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts3.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_out[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_out[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
    m.fs.gts3.control_volume.properties_out[0].flow_mol, 1e-5)

    for c in m.fs.gt_mix.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-9)
    for c in m.fs.air_gas_mix1.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)
    for c in m.fs.air_gas_mix2.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)
    for c in m.fs.air_gas_mix3.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)


def initialize_blue_hydrogen_plant(m):
    outlvl = 4 #idaeslog.INFO_LOW
    m.fs.ng_preheat1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.ng_preheat2.tube_inlet, m.fs.ng_preheat1.tube_outlet)
    m.fs.ng_preheat2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube_mix.fuel_inlet, m.fs.ng_preheat2.tube_outlet)
    m.fs.smr_tube_mix.initialize(outlvl = outlvl)
    m.fs.smr_shell_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube.inlet, m.fs.smr_tube_mix.outlet)
    m.fs.smr_tube.outlet.temperature.fix(1155)
    m.fs.smr_tube.initialize(outlvl = outlvl)
    m.fs.smr_tube.outlet.temperature.unfix()
    print('flow_mol=',m.fs.smr_tube.inlet.flow_mol[0].value)
    print('temp in=',m.fs.smr_tube.inlet.temperature[0].value)
    print('heat=',m.fs.smr_tube.heat_duty[0].value)
    print('temp out=',m.fs.smr_tube.outlet.temperature[0].value)
    copy_port_values(m.fs.smr_shell.inlet, m.fs.smr_shell_mix.outlet)
    m.fs.smr_shell.outlet.temperature.fix(1200)
    m.fs.smr_shell.initialize(outlvl = outlvl)
    m.fs.smr_shell.outlet.temperature.unfix()
    print('flow_mol=',m.fs.smr_shell.inlet.flow_mol[0].value)
    print('temp in =',m.fs.smr_shell.inlet.temperature[0].value)
    print('heat=',m.fs.smr_shell.heat_duty[0].value)
    print('temp out=',m.fs.smr_shell.outlet.temperature[0].value)
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
    copy_port_values(m.fs.h2_comp.inlet, m.fs.h2_separator.h2_outlet)
    m.fs.h2_comp.initialize(outlvl = outlvl)    
    copy_port_values(m.fs.flue_gas_mp_cooler.shell_inlet, m.fs.smr_shell.outlet)
    m.fs.flue_gas_mp_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.fg_splitter.inlet, m.fs.flue_gas_mp_cooler.shell_outlet)
    m.fs.fg_splitter.initialize(outlvl = outlvl)
    copy_port_values(m.fs.flue_gas_lp_cooler.shell_inlet, m.fs.fg_splitter.main_outlet)
    copy_port_values(m.fs.flue_gas_lp_cooler.tube_inlet, m.fs.syngas_cooler.tube_outlet)
    m.fs.flue_gas_lp_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.fg_mix.main_inlet, m.fs.flue_gas_lp_cooler.shell_outlet)
    copy_port_values(m.fs.fg_mix.bypass_inlet, m.fs.fg_splitter.bypass_outlet)
    m.fs.fg_mix.initialize(outlvl = outlvl)
    
    copy_port_values(m.fs.mp_mix.mp_inlet, m.fs.mp_knockout_condenser.tube_outlet)
    copy_port_values(m.fs.mp_mix.mpkw_inlet, m.fs.mp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.lpkw_inlet, m.fs.lp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.cwkw_inlet, m.fs.cw_knockout_condenser.shell_water_outlet)
    m.fs.mp_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_booster_pump.inlet, m.fs.mp_mix.outlet)
    m.fs.mp_booster_pump.initialize(outlvl = outlvl)
    
    copy_port_values(m.fs.flue_gas_mp_econ.shell_inlet, m.fs.fg_mix.outlet)
    copy_port_values(m.fs.flue_gas_mp_econ.tube_inlet, m.fs.mp_booster_pump.outlet)
    m.fs.flue_gas_mp_econ.initialize(outlvl = outlvl)
    copy_port_values(m.fs.ng_preheat1.shell_inlet, m.fs.flue_gas_mp_econ.shell_outlet)
    m.fs.ng_preheat1.initialize(outlvl = outlvl)
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
    m.fs.co2_pre_cool.outlet.temperature.fix(303)
    m.fs.co2_pre_cool.heat_duty.unfix()
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


def initialize_gas_turbine_plant(m):
    enth_mol_25 = iapws95.htpx(P=101325*pyo.units.Pa, T=298.15*pyo.units.K,)
    print('standard enthalpy of water=', enth_mol_25)
    outlvl = 4 #idaeslog.INFO_LOW
    m.fs.comp_igv.initialize(outlvl = outlvl)
    copy_port_values(m.fs.air_comp.inlet, m.fs.comp_igv.outlet)
    m.fs.air_comp.initialize(outlvl = outlvl)
    copy_port_values(m.fs.air_splitter.inlet, m.fs.air_comp.outlet)
    m.fs.air_splitter.initialize(outlvl = outlvl)
    copy_port_values(m.fs.gt_mix.air_inlet, m.fs.air_splitter.comb_outlet)
    m.fs.gt_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.gt_combustor.inlet, m.fs.gt_mix.outlet)
    m.fs.gt_combustor.initialize(outlvl = outlvl)
    copy_port_values(m.fs.gts1.inlet, m.fs.gt_combustor.outlet)
    m.fs.gts1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.air_valve1.inlet, m.fs.air_splitter.gts1_outlet)
    m.fs.air_valve1.outlet.pressure.fix(m.fs.gts1.outlet.pressure[0].value)
    m.fs.air_valve1.deltaP.unfix()
    m.fs.air_valve1.initialize(outlvl = outlvl)
    m.fs.air_valve1.outlet.pressure.unfix()
    copy_port_values(m.fs.air_gas_mix1.gas_inlet, m.fs.gts1.outlet)
    copy_port_values(m.fs.air_gas_mix1.air_inlet, m.fs.air_valve1.outlet)
    m.fs.air_gas_mix1.initialize(outlvl = outlvl)
    copy_port_values(m.fs.gts2.inlet, m.fs.air_gas_mix1.outlet)
    m.fs.gts2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.air_valve2.inlet, m.fs.air_splitter.gts2_outlet)
    m.fs.air_valve2.outlet.pressure.fix(m.fs.gts2.outlet.pressure[0].value)
    m.fs.air_valve2.deltaP.unfix()
    m.fs.air_valve2.initialize(outlvl = outlvl)
    m.fs.air_valve2.outlet.pressure.unfix()
    copy_port_values(m.fs.air_gas_mix2.gas_inlet, m.fs.gts2.outlet)
    copy_port_values(m.fs.air_gas_mix2.air_inlet, m.fs.air_valve2.outlet)
    m.fs.air_gas_mix2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.gts3.inlet, m.fs.air_gas_mix2.outlet)
    m.fs.gts3.initialize(outlvl = outlvl)
    copy_port_values(m.fs.air_valve3.inlet, m.fs.air_splitter.gts3_outlet)
    m.fs.air_valve3.outlet.pressure.fix(m.fs.gts3.outlet.pressure[0].value)
    m.fs.air_valve3.deltaP.unfix()
    m.fs.air_valve3.initialize(outlvl = outlvl)
    m.fs.air_valve3.outlet.pressure.unfix()
    copy_port_values(m.fs.air_gas_mix3.gas_inlet, m.fs.gts3.outlet)
    copy_port_values(m.fs.air_gas_mix3.air_inlet, m.fs.air_valve3.outlet)
    m.fs.air_gas_mix3.initialize(outlvl = outlvl)




def get_full_plant_model():
    logging.getLogger('idaes.core.util.scaling').setLevel(logging.ERROR)
    m = pyo.ConcreteModel(name='NGSC and SMR integrated energy system')
    m.fs = FlowsheetBlock(dynamic = False)
    # solver and options
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 120,
            'bound_push': 1e-16
    }
    build_gas_turbine_plant(m)
    build_blue_hydrogen_plant(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    set_gas_turbine_plant_inputs(m)
    set_blue_hydrogen_plant_inputs(m)
    scale_gas_turbine_plant(m)
    scale_blue_hydrogen_plant(m)
    iscale.calculate_scaling_factors(m)
    init_fname = 'ies_ngsc_baseline.json'
    if not os.path.exists(init_fname):
        print('work')
        initialize_gas_turbine_plant(m)
        initialize_blue_hydrogen_plant(m)
    m.fs.smr_tube_mix.steam_inlet_state[:].flow_mol.fix(2730.6)
    m.fs.mp_pump.inlet.flow_mol.unfix()
    m.fs.comp_igv.Cv.fix()
    m.fs.comp_igv.deltaP.unfix()
    m.fs.comp_igv.valve_opening.unfix()
    m.fs.comp_igv.inlet.flow_mol.unfix()
    #m.fs.smr_shell_mix.air_inlet.flow_mol.unfix()
    m.fs.smr_shell_mix.fuel_inlet.flow_mol.unfix()
    m.fs.o2_in_flue_gas.deactivate()
    if os.path.exists(init_fname):
        #print('loading initial conditions')
        ms.from_json(m, fname=init_fname) #,wts=ms.StoreSpec(suffix=False)
        #m.fs.total_reboiler_duty_constraint.deactivate()
        #m.fs.flue_gas_lp_cooler.area.fix()
        m.fs.o2_in_flue_gas.deactivate()
        solver.solve(m, tee=True)
    else:
        m.fs.total_reboiler_duty_constraint.deactivate()
        solver.solve(m, tee=True)
        m.fs.total_reboiler_duty_constraint.activate()
        #m.fs.flue_gas_lp_cooler.area.unfix()
        m.fs.fg_splitter.split_fraction[:,'main_outlet'].unfix()
        solver.solve(m, tee=True)
        ms.to_json(m, fname='ies_ngsc_baseline.json')
       
    return m



if __name__ == "__main__":
    m = get_full_plant_model()
