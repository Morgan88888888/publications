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
This is a flowsheet for co-generation of blue hydrogen and power by SMR and NGCC integration
"""
import xlsxwriter
import os
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.network import Arc

# IDAES Imports
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state as copy_port_values
from idaes.core.util import model_serializer as ms
from idaes.core.util.model_statistics import degrees_of_freedom
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
    #GibbsReactor,
    StoichiometricReactor,
    Separator,
    SplittingType,
    Translator,
    Valve,
    ValveFunctionType,
    Compressor)
from idaes.models_extra.power_generation.unit_models import (
    Drum,
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

from natural_gas_PR import get_prop, get_rxn
#from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop, get_rxn
import idaes.logger as idaeslog

__author__ = "Jinliang Ma, Maojian Wang"

def performance_curves(m, flow_scale=2):
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
        f = pyo.log(1.265*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -2085.1*f**3 + 38433*f**2 - 150764*f + 422313)
    @m.fs.gts2.performance_curve.Constraint(m.fs.config.time)
    def head_isen_eqn(b, t):
        f = pyo.log(1.265*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1676.3*f**3 + 34916*f**2 - 173801*f + 456957)
    @m.fs.gts3.performance_curve.Constraint(m.fs.config.time)
    def head_isen_eqn(b, t):
        f = pyo.log(1.265*
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1373.6*f**3 + 31759*f**2 - 188528*f + 500520)


def build_blue_hydrogen_plant(m):
    # create property packages
    m.fs.flue_props = FlueGasParameterBlock()

    # water property
    m.fs.water_props = iapws95.Iapws95ParameterBlock()
    if not hasattr(m.fs, "syn_props"):
        # syngas property
        syn_config = get_prop(
            components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
        m.fs.syn_props = GenericParameterBlock(**syn_config)

    # water gas shift reaction package
    m.fs.wgs_rxn_props = GenericReactionParameterBlock(
        **get_rxn(
            m.fs.syn_props, reactions=["water_gas_shift_rxn"]))

    # property package translater from syn_props to flue_props
    m.fs.hrsg_translator = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.flue_props})

    # constraints for property translator
    @m.fs.hrsg_translator.Constraint(m.fs.time, m.fs.flue_props.component_list)
    def mol_frac_eqn(b, t, i):
        return b.inlet.flow_mol[t]*b.inlet.mole_frac_comp[t, i] == \
               b.outlet.flow_mol_comp[t, i]

    @m.fs.hrsg_translator.Constraint(m.fs.time)
    def temperature_eqn(b, t):
        return b.outlet.temperature[t] == b.inlet.temperature[t]

    @m.fs.hrsg_translator.Constraint(m.fs.time)
    def pressure_eqn(b, t):
        return b.outlet.pressure[t] == b.inlet.pressure[t]

    # build blue hydrogen plant section
    # natural gas preheater 2
    m.fs.ng_preheat2 = HeatExchanger(
        **{
                 "hot_side_name": "shell",
                 "cold_side_name" : "tube",
                 "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.flue_props,
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
        **{"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    # shell side smr combustor
    m.fs.smr_shell = GibbsReactor(
        **{"has_heat_transfer": True,
                 "has_pressure_change": True,
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

    # reboiler liquid return pump
    m.fs.rb_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    m.fs.rb_splitter = HelmSplitter(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "outlet_list": ["cond_outlet", "mix_outlet"]})

    # mp pump
    m.fs.mp_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # mp booster pump after mixed with knockout water
    m.fs.mp_booster_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # reboiler knockout condenser
    m.fs.rb_knockout_condenser = WaterKnockoutCondenser(
        **{"tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.syn_props,
                 "water_property_package": m.fs.water_props})

    # mp knockout condenser
    m.fs.mp_knockout_condenser = WaterKnockoutCondenser(
        **{"tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.syn_props,
                 "water_property_package": m.fs.water_props})

    # lp knockout condenser
    m.fs.lp_knockout_condenser = WaterKnockoutCondenser(
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
            "inlet_list": ["mp_inlet", "lpkw_inlet", "mpkw_inlet", "cwkw_inlet", "rbkw_inlet"]})

    # precombustion co2 capture system
    m.fs.co2_pre_separator = Separator(
        **{"outlet_list": ["co2_outlet", "syngas_outlet", "h2o_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

    # h2 separation system, e.g. PSA
    m.fs.h2_separator = Separator(
        **{"outlet_list": ["h2_outlet", "tailgas_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

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
            return b.properties_out[t].mole_frac_comp[j] == 0.99999
        else:
            return b.properties_out[t].mole_frac_comp[j] == 1e-6

    # water shift reaction conversion factor
    m.fs.wgsr.x_co = pyo.Var(initialize=0.9, doc='conversion factor for reaction based on CO')

    # water shift reaction conversion factor
    m.fs.deltaT_smr = pyo.Var(initialize=50.0, units=pyo.units.K, doc='SMR exit temperature difference between flue gas and syngas')
    m.fs.deltaT_smr.fix()

    # additional constraints
    @m.fs.Constraint(m.fs.config.time)
    def heat_balance_of_smr(b, t):
        return 1e-8*(b.smr_tube.control_volume.heat[t] + b.smr_shell.control_volume.heat[t]) == 0

    @m.fs.Constraint(m.fs.config.time)
    def outlet_temperature_of_smr(b, t):
        return 0.01*(b.smr_tube.outlet.temperature[t] + b.deltaT_smr - b.smr_shell.outlet.temperature[t]) == 0 #50

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
        return b.mixed_state[t].pressure == b.fuel_inlet_state[t].pressure

    # constraints for outlet pressure of smr_shell outlet
    @m.fs.smr_tube_mix.Constraint(m.fs.time)
    def smr_tube_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.fuel_inlet_state[t].pressure

    # constraints for mp mixer related to the knockout water
    @m.fs.mp_mix.Constraint(m.fs.time)
    def mp_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.cwkw_inlet_state[t].pressure

    # constraint for seting smr tube side outlet temperature
    @m.fs.Constraint(m.fs.config.time)
    def outlet_temperature_of_smr_fixed(b, t):
        return 0.001*(b.smr_tube.outlet.temperature[t] - 1155.0) == 0    

    # stream connections    
    m.fs.ng_preheat2_to_smr_tube_mix = Arc(
        source=m.fs.ng_preheat2.tube_outlet,
        destination=m.fs.smr_tube_mix.fuel_inlet)

    if hasattr(m.fs, "air_gas_mix3"):
        m.fs.air_gas_mix3_to_smr_shell_mix = Arc(
        source=m.fs.air_gas_mix3.outlet,
        destination=m.fs.smr_shell_mix.air_inlet)

    m.fs.smr_tube_mix_to_smr_tube = Arc(
        source=m.fs.smr_tube_mix.outlet,
        destination=m.fs.smr_tube.inlet)

    m.fs.smr_shell_mix_to_smr_shell = Arc(
        source=m.fs.smr_shell_mix.outlet,
        destination=m.fs.smr_shell.inlet)

    m.fs.smr_shell_to_hrsg_translator = Arc(
        source=m.fs.smr_shell.outlet,
        destination=m.fs.hrsg_translator.inlet)

    m.fs.hrsg_translator_to_ng_preheat2 = Arc(
        source=m.fs.hrsg_translator.outlet,
        destination=m.fs.ng_preheat2.shell_inlet)

    m.fs.smr_tube_to_syngas_cooler = Arc(
        source=m.fs.smr_tube.outlet,
        destination=m.fs.syngas_cooler.shell_inlet)

    m.fs.syngas_cooler_to_wgsr = Arc(
        source=m.fs.syngas_cooler.shell_outlet,
        destination=m.fs.wgsr.inlet)

    m.fs.wgsr_to_wgsr_cooler = Arc(
        source=m.fs.wgsr.outlet,
        destination=m.fs.wgsr_cooler.shell_inlet)

    m.fs.wgsr_cooler_to_rb_knockout_condenser = Arc(
        source=m.fs.wgsr_cooler.shell_outlet,
        destination=m.fs.rb_knockout_condenser.shell_inlet)

    m.fs.rb_knockout_condenser_to_mp_knockout_condenser = Arc(
        source=m.fs.rb_knockout_condenser.shell_outlet,
        destination=m.fs.mp_knockout_condenser.shell_inlet)

    m.fs.mp_knockout_condenser_to_lp_knockout_condenser = Arc(
        source=m.fs.mp_knockout_condenser.shell_outlet,
        destination=m.fs.lp_knockout_condenser.shell_inlet)

    m.fs.lp_knockout_condenser_to_cw_knockout_condenser = Arc(
        source=m.fs.lp_knockout_condenser.shell_outlet,
        destination=m.fs.cw_knockout_condenser.shell_inlet)

    m.fs.cw_knockout_condenser_to_co2_pre_separator = Arc(
        source=m.fs.cw_knockout_condenser.shell_outlet,
        destination=m.fs.co2_pre_separator.inlet)

    m.fs.co2_pre_separator_to_h2_separator = Arc(
        source=m.fs.co2_pre_separator.syngas_outlet,
        destination=m.fs.h2_separator.inlet)

    m.fs.h2_separator_to_smr_shell_mix = Arc(
        source=m.fs.h2_separator.tailgas_outlet,
        destination=m.fs.smr_shell_mix.tailgas_inlet)

    m.fs.rb_pump_to_rb_splitter = Arc(
        source=m.fs.rb_pump.outlet,
        destination=m.fs.rb_splitter.inlet)

    m.fs.rb_splitter_to_rb_knockout_condenser = Arc(
        source=m.fs.rb_splitter.cond_outlet,
        destination=m.fs.rb_knockout_condenser.tube_inlet)

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

    m.fs.rb_knockout_water_to_mp_mix = Arc(
        source=m.fs.rb_knockout_condenser.shell_water_outlet,
        destination=m.fs.mp_mix.rbkw_inlet)

    m.fs.mp_mix_to_mp_booster_pump = Arc(
        source=m.fs.mp_mix.outlet,
        destination=m.fs.mp_booster_pump.inlet)

    m.fs.mp_booster_pump_to_wgsr_cooler = Arc(
        source=m.fs.mp_booster_pump.outlet,
        destination=m.fs.wgsr_cooler.tube_inlet)

    m.fs.wgsr_cooler_to_syngas_cooler = Arc(
        source=m.fs.wgsr_cooler.tube_outlet,
        destination=m.fs.syngas_cooler.tube_inlet)

    m.fs.syngas_cooler_to_steam_translator = Arc(
        source=m.fs.syngas_cooler.tube_outlet,
        destination=m.fs.steam_translator.inlet)

    m.fs.steam_translator_to_smr_tube_mix = Arc(
        source=m.fs.steam_translator.outlet,
        destination=m.fs.smr_tube_mix.steam_inlet)


def build_gas_turbine_plant(m):
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
        return b.air_gas_mix3.outlet.temperature[t] == 955.2     

    # constraint for exit pressure by adjusting inlet guide vane opening
    @m.fs.Constraint(m.fs.time)
    def gt_exit_pres_constraint(b, t):
        return b.air_gas_mix3.outlet.pressure[t] == 104446

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


def build_hrsg_steam_turbine_plant(m):
    if not hasattr(m.fs, "flue_props"):
        m.fs.flue_props = FlueGasParameterBlock()
    if not hasattr(m.fs, "water_props"):
        m.fs.water_props = iapws95.Iapws95ParameterBlock()
    if not hasattr(m.fs, "co2_props"):        
        m.fs.co2_props = swco2.SWCO2ParameterBlock()
    # natural gas preheat1
    m.fs.ng_preheat1 = HeatExchanger(
        **{
                "hot_side_name": "shell",
                "cold_side_name" : "tube",
                "delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.flue_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.syn_props,
                          "has_pressure_change": True}})

    # LP economizer
    m.fs.lp_econ = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 5,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # reboiler econ mixer
    m.fs.rb_econ_mixer = HelmMixer(
        **{
            "property_package": m.fs.water_props,
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["econ_inlet", "reboiler_inlet"]})

    @m.fs.rb_econ_mixer.Constraint(m.fs.time)
    def rb_econ_mixer_pressure_constraint(b, t):
        return 1e-6 * b.econ_inlet_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # lp drum mixer
    m.fs.lp_drum_mixer = HelmMixer(
        **{
            "property_package": m.fs.water_props,
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["evap_inlet", "reboiler_inlet"]})

    @m.fs.lp_drum_mixer.Constraint(m.fs.time)
    def lp_drum_mixer_pressure_constraint(b, t):
        return 1e-6 * b.evap_inlet_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # Splitter of heated feed water to LP, IP, and HP streams
    m.fs.lp_ip_hp_split = HelmSplitter(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "outlet_list": ["lp_outlet", "ip_outlet", "hp_outlet"]})

    # LP drum
    m.fs.lp_drum = Drum(
        **{
            "property_package": m.fs.water_props,
            "has_holdup": False,
            "has_heat_transfer": True,
            "has_pressure_change": True,})

    # LP downcomer
    m.fs.lp_downcomer = Downcomer(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "has_holdup": True,
                 "has_heat_transfer": True})    

    # LP evaporator
    m.fs.lp_evap = Evaporator(**{
                 "tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.flue_props,
                 "has_pressure_change": True,
                 "has_holdup": False,
                 "tube_arrangement": TubeArrangement.inLine})

    # Splitter as deareator at LP drum
    m.fs.lp_da_split = HelmSplitter(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "outlet_list": ["lp_outlet", "purge_outlet"]})

    # IP pump
    m.fs.ip_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # HP pump
    m.fs.hp_pump = HelmIsentropicCompressor(
        **{"dynamic": False,
                 "property_package": m.fs.water_props})

    # HP economizer1
    m.fs.hp_econ1 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 3,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # IP economizer1
    m.fs.ip_econ1 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 1,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # HP economizer2
    m.fs.hp_econ2 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 3,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # IP economizer2
    m.fs.ip_econ2 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 1,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # IP drum
    m.fs.ip_drum = Drum(
        **{
            "property_package": m.fs.water_props,
            "has_holdup": False,
            "has_heat_transfer": True,
            "has_pressure_change": True,})


    # IP Downcomer
    m.fs.ip_downcomer = Downcomer(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "has_holdup": True,
                 "has_heat_transfer": True})    

    # IP evaporator
    m.fs.ip_evap = Evaporator(
                 tube_side_property_package=m.fs.water_props,
                 shell_side_property_package=m.fs.flue_props,
                 has_pressure_change=True,
                 has_holdup=False,
                 tube_arrangement=TubeArrangement.inLine)

    # LP superheater
    m.fs.lp_sh = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 4,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # HP economizer3
    m.fs.hp_econ3 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 3,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # IP superheater1
    m.fs.ip_sh1 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 1,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # IP steam mixer
    m.fs.ip_mixer = HelmMixer(
        **{
            "property_package": m.fs.water_props,
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["turb_inlet", "hrsg_inlet"]})

    @m.fs.ip_mixer.Constraint(m.fs.time)
    def ip_mixer_pressure_constraint(b, t):
        return 1e-6 * b.turb_inlet_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # added constraint to control IP pump deltaP
    @m.fs.Constraint(m.fs.time)
    def ip_mixer_pressure_constraint2(b, t):
        return 1e-6 * b.ip_mixer.turb_inlet_state[t].pressure == 1e-6 * 0.95*b.ip_mixer.hrsg_inlet_state[t].pressure

    # HP economizer4
    m.fs.hp_econ4 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 2,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Liq",
                 "has_radiation": False})

    # HP drum
    m.fs.hp_drum = Drum(
        **{
            "property_package": m.fs.water_props,
            "has_holdup": False,
            "has_heat_transfer": True,
            "has_pressure_change": True,})

    # HP Downcomer
    m.fs.hp_downcomer = Downcomer(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "has_holdup": True,
                 "has_heat_transfer": True})    

    # IP evaporator
    m.fs.hp_evap = Evaporator(**{
                 "tube_side_property_package": m.fs.water_props,
                 "shell_side_property_package": m.fs.flue_props,
                 "has_pressure_change": True,
                 "has_holdup": False,
                 "tube_arrangement": TubeArrangement.inLine})

    # HP superheater1, holdup set to False
    m.fs.hp_sh1 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 2,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # IP superheater2
    m.fs.ip_sh2 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 1,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # HP superheater2
    m.fs.hp_sh2 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 2,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # IP superheater3
    m.fs.ip_sh3 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 1,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # HP superheater3
    m.fs.hp_sh3 = HeatExchangerCrossFlow1D(
        **{"tube_side":{"property_package": m.fs.water_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "shell_side":{"property_package": m.fs.flue_props, "has_holdup": False, 
                              "has_pressure_change": True},
                 "finite_elements": 2,
                 "flow_type": "counter_current",
                 "tube_arrangement": "in-line",
                 "tube_side_water_phase": "Vap",
                 "has_radiation": False})

    # Steam turbine
    m.fs.steam_turb = HelmTurbineMultistage(
        **{
            "property_package": m.fs.water_props,
            "num_parallel_inlet_stages": 1,
            "num_hp": 7,  # at full load ave P ratio about 0.8238 with inlet stage
            "num_ip": 10,  # at full load ave P ratio about 0.8264
            "num_lp": 11,  # at full load ave P ratio about 0.7194 with outlet stage
            "hp_disconnect": [7],  # disconected for reheater
            "ip_disconnect": [10],  # disconnected for HRSG LP steam mix
        })

    # Condenser
    m.fs.condenser = HelmNtuCondenser(
        **{
            "shell": {
                "has_pressure_change": False,
                "property_package": m.fs.water_props,
            },
            "tube": {"has_pressure_change": False, "property_package": m.fs.water_props},
        })

    # Condenser hotwell
    m.fs.hotwell = HelmMixer(
        **{
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["condensate_inlet", "makeup_inlet"],
            "property_package": m.fs.water_props,
        })

    @m.fs.hotwell.Constraint(m.fs.time)
    def hotwell_mixer_pressure_constraint(b, t):
        return 1e-6 * b.condensate_inlet_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # Condensate pump
    m.fs.cond_pump = HelmIsentropicCompressor(
        **{"property_package": m.fs.water_props})

    # LP mix between IP outlet stream and HRSG LP stream
    m.fs.lp_mixer = HelmMixer(
        **{
            "property_package": m.fs.water_props,
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["turb_inlet", "hrsg_inlet"],
        })

    @m.fs.lp_mixer.Constraint(m.fs.time)
    def lp_mixer_pressure_constraint(b, t):
        return 1e-6 * b.turb_inlet_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # added constraint to control LP pump deltaP
    @m.fs.Constraint(m.fs.time)
    def lp_mixer_pressure_constraint2(b, t):
        return 1e-6 * b.lp_mixer.turb_inlet_state[t].pressure == 1e-6 * 0.925*b.lp_mixer.hrsg_inlet_state[t].pressure

    m.fs.lp_splitter = HelmSplitter(
        **{"dynamic": False,
                 "property_package": m.fs.water_props,
                 "outlet_list": ["turb_outlet", "steam_outlet"]})

    # postcombustion co2 capture system
    m.fs.co2_post_separator = Separator(
        **{"outlet_list": ["co2_outlet", "fluegas_outlet", "h2o_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.flue_props})

    m.fs.co2_mix = Mixer(
        **{"inlet_list": ["pre_inlet", "post_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.co2_props})

    # constraints for co2 mixer of pre and post co2 streams
    @m.fs.co2_mix.Constraint(m.fs.time)
    def co2_mix_pressure(b, t):
        return b.mixed_state[t].pressure == b.pre_inlet_state[t].pressure

    # property package translater from flue_props to co2_props
    m.fs.co2_translator_flue = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.flue_props,
                 "outlet_property_package": m.fs.co2_props})

    # property package translater from syn_props to co2_props
    m.fs.co2_translator_syn = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.co2_props})

    # additional variables and constraints
    # constraints for co2_translator_flue
    @m.fs.co2_translator_flue.Constraint(m.fs.time)
    def co2_translator_flue_T(b, t):
        return b.properties_out[t].temperature == 303

    @m.fs.co2_translator_flue.Constraint(m.fs.time)
    def co2_translator_flue_P(b, t):
        return b.outlet.pressure[t] == 1.4e5

    @m.fs.co2_translator_flue.Constraint(m.fs.time)
    def co2_translator_flue_flow(b, t):
        return b.inlet.flow_mol_comp[t,"CO2"] == b.outlet.flow_mol[t]

    # constraints for co2_translator_syn
    @m.fs.co2_translator_syn.Constraint(m.fs.time)
    def co2_translator_syn_T(b, t):
        return b.properties_out[t].temperature == 303

    @m.fs.co2_translator_syn.Constraint(m.fs.time)
    def co2_translator_syn_P(b, t):
        return b.outlet.pressure[t] == 1.4e5

    @m.fs.co2_translator_syn.Constraint(m.fs.time)
    def co2_translator_syn_flow(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]


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

    m.fs.h2_comp = Compressor(**{"property_package": m.fs.syn_props})

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

    # reboiler returns saturated steam to the reboiler pump
    @m.fs.Constraint(m.fs.time)
    def rb_pump_inlet_pressure(b, t):
        return 1e-5 * b.rb_pump.inlet.pressure[t] == 1e-5 * b.lp_splitter.steam_outlet_state[t].pressure

    @m.fs.Constraint(m.fs.time)
    def rb_pump_inlet_flow_mol(b, t):
        return 1e-3 * b.rb_pump.inlet.flow_mol[t] == 1e-3 * b.lp_splitter.steam_outlet_state[t].flow_mol

    @m.fs.Constraint(m.fs.time)
    def rb_pump_inlet_enth_mol(b, t):
        return 1e-3 * b.rb_pump.inlet.enth_mol[t] == 1e-3 * b.lp_splitter.steam_outlet_state[t].enth_mol_sat_phase["Liq"]

    @m.fs.Expression(m.fs.time)
    def st_net_power(b, t):
        return (-b.steam_turb.power[t] - b.cond_pump.work_mechanical[t] - b.ip_pump.work_mechanical[t] - b.hp_pump.work_mechanical[t] - b.rb_pump.work_mechanical[t])

    m.fs.duty_per_mol_pre = pyo.Param(initialize=1600*1054/0.453592*0.044,
        units=pyo.units.J/pyo.units.mol,
        doc='required reboiler duty per mol CO2 captured for pre-combustion capture')

    m.fs.duty_per_mol_post = pyo.Param(initialize=1217.184*1054/0.453592*0.044,
        units=pyo.units.J/pyo.units.mol,
        doc='required reboiler duty per mol CO2 captured for post-combustion capture')

    @m.fs.Expression(m.fs.config.time)
    def pre_capture_reboiler_duty(b, t):
        return (b.duty_per_mol_pre*b.co2_pre_separator.co2_outlet_state[t].flow_mol
                *b.co2_pre_separator.co2_outlet_state[t].mole_frac_comp['CO2'])

    @m.fs.Expression(m.fs.config.time)
    def post_capture_reboiler_duty(b, t):
        return (b.duty_per_mol_post*b.co2_post_separator.co2_outlet_state[t].flow_mol_comp['CO2'])

    @m.fs.Expression(m.fs.config.time)
    def total_lp_duty_available(b, t):
        return (b.lp_splitter.steam_outlet_state[t].flow_mol*
               (b.lp_splitter.steam_outlet_state[t].enth_mol-
                b.lp_splitter.steam_outlet_state[t].enth_mol_sat_phase["Liq"]))

    # reqired reboiler duty
    @m.fs.Constraint(m.fs.time)
    def reboiler_duty_constraint(b, t):
        return 1e-8*(b.pre_capture_reboiler_duty[t] + b.post_capture_reboiler_duty[t] - b.total_lp_duty_available[t]) == 0

    # main steam pressure
    m.fs.p_main_steam = pyo.Var(initialize=1.7e7, units=pyo.units.Pa, doc='main steam pressure')
    m.fs.p_main_steam.fix()
    
    # main steam pressure constraint
    @m.fs.Constraint(m.fs.time)
    def main_steam_pressure_constraint(b, t):
        return 1e-7*(b.hp_sh3.tube.properties[t,0].pressure - b.p_main_steam) == 0

    # Stream connections
    # Shell side arcs   
    m.fs.ng_preheat2_to_hp_sh3 = Arc(
        source=m.fs.ng_preheat2.shell_outlet,
        destination=m.fs.hp_sh3.shell_inlet)
    
    m.fs.hp_sh3_to_ip_sh3 = Arc(
        source=m.fs.hp_sh3.shell_outlet,
        destination=m.fs.ip_sh3.shell_inlet)

    m.fs.ip_sh3_to_hp_sh2 = Arc(
        source=m.fs.ip_sh3.shell_outlet,
        destination=m.fs.hp_sh2.shell_inlet)

    m.fs.hp_sh2_to_ip_sh2 = Arc(
        source=m.fs.hp_sh2.shell_outlet,
        destination=m.fs.ip_sh2.shell_inlet)

    m.fs.ip_sh2_to_hp_sh1 = Arc(
        source=m.fs.ip_sh2.shell_outlet,
        destination=m.fs.hp_sh1.shell_inlet)

    m.fs.hp_sh1_to_hp_evap = Arc(
        source=m.fs.hp_sh1.shell_outlet,
        destination=m.fs.hp_evap.shell_inlet)

    m.fs.hp_evap_to_hp_econ4 = Arc(
        source=m.fs.hp_evap.shell_outlet,
        destination=m.fs.hp_econ4.shell_inlet)

    m.fs.hp_econ4_to_ip_sh1 = Arc(
        source=m.fs.hp_econ4.shell_outlet,
        destination=m.fs.ip_sh1.shell_inlet)

    m.fs.ip_sh1_to_hp_econ3 = Arc(
        source=m.fs.ip_sh1.shell_outlet,
        destination=m.fs.hp_econ3.shell_inlet)

    m.fs.hp_econ3_to_lp_sh = Arc(
        source=m.fs.hp_econ3.shell_outlet,
        destination=m.fs.lp_sh.shell_inlet)

    m.fs.lp_sh_to_ip_evap = Arc(
        source=m.fs.lp_sh.shell_outlet,
        destination=m.fs.ip_evap.shell_inlet)

    m.fs.ip_evap_to_ip_econ2 = Arc(
        source=m.fs.ip_evap.shell_outlet,
        destination=m.fs.ip_econ2.shell_inlet)

    m.fs.ip_econ2_to_hp_econ2 = Arc(
        source=m.fs.ip_econ2.shell_outlet,
        destination=m.fs.hp_econ2.shell_inlet)

    m.fs.hp_econ2_to_ip_econ1 = Arc(
        source=m.fs.hp_econ2.shell_outlet,
        destination=m.fs.ip_econ1.shell_inlet)

    m.fs.ip_econ1_to_hp_econ1 = Arc(
        source=m.fs.ip_econ1.shell_outlet,
        destination=m.fs.hp_econ1.shell_inlet)

    m.fs.hp_econ1_to_lp_evap = Arc(
        source=m.fs.hp_econ1.shell_outlet,
        destination=m.fs.lp_evap.shell_inlet)

    m.fs.lp_evap_to_lp_econ = Arc(
        source=m.fs.lp_evap.shell_outlet,
        destination=m.fs.lp_econ.shell_inlet)

    m.fs.lp_econ_to_ng_preheat1 = Arc(
        source=m.fs.lp_econ.shell_outlet,
        destination=m.fs.ng_preheat1.shell_inlet)

    m.fs.ng_preheat1_to_co2_post_separator = Arc(
        source=m.fs.ng_preheat1.shell_outlet,
        destination=m.fs.co2_post_separator.inlet)

    # Tube side arcs
    m.fs.ng_preheat1_to_ng_preheat2 = Arc(
        source=m.fs.ng_preheat1.tube_outlet,
        destination=m.fs.ng_preheat2.tube_inlet)

    m.fs.lp_econ_to_rb_econ_mixer = Arc(
        source=m.fs.lp_econ.tube_outlet,
        destination=m.fs.rb_econ_mixer.econ_inlet)

    m.fs.rb_splitter_to_rb_econ_mixer = Arc(
        source=m.fs.rb_splitter.mix_outlet,
        destination=m.fs.rb_econ_mixer.reboiler_inlet)

    m.fs.rb_econ_mixer_to_lp_ip_hp_split = Arc(
        source=m.fs.rb_econ_mixer.outlet,
        destination=m.fs.lp_ip_hp_split.inlet)

    m.fs.lp_ip_hp_split_to_lp_drum = Arc(
        source=m.fs.lp_ip_hp_split.lp_outlet,
        destination=m.fs.lp_drum.feedwater_inlet)

    m.fs.lp_drum_to_lp_downcomer = Arc(
        source=m.fs.lp_drum.liquid_outlet,
        destination=m.fs.lp_downcomer.inlet)

    m.fs.lp_downcomer_to_lp_evap = Arc(
        source=m.fs.lp_downcomer.outlet,
        destination=m.fs.lp_evap.tube_inlet)

    m.fs.rb_knockout_condenser_to_lp_drum_mixer = Arc(
        source=m.fs.rb_knockout_condenser.tube_outlet,
        destination=m.fs.lp_drum_mixer.reboiler_inlet)

    m.fs.lp_evap_to_lp_drum_mixer = Arc(
        source=m.fs.lp_evap.tube_outlet,
        destination=m.fs.lp_drum_mixer.evap_inlet)

    m.fs.lp_drum_mixer_to_lp_drum = Arc(
        source=m.fs.lp_drum_mixer.outlet,
        destination=m.fs.lp_drum.water_steam_inlet)

    m.fs.lp_drum_to_lp_da_split = Arc(
        source=m.fs.lp_drum.steam_outlet,
        destination=m.fs.lp_da_split.inlet)

    m.fs.lp_ip_hp_split_to_ip_pump = Arc(
        source=m.fs.lp_ip_hp_split.ip_outlet,
        destination=m.fs.ip_pump.inlet)

    m.fs.lp_ip_hp_split_to_hp_pump = Arc(
        source=m.fs.lp_ip_hp_split.hp_outlet,
        destination=m.fs.hp_pump.inlet)

    m.fs.hp_pump_to_hp_econ1 = Arc(
        source=m.fs.hp_pump.outlet,
        destination=m.fs.hp_econ1.tube_inlet)

    m.fs.ip_pump_to_ip_econ1 = Arc(
        source=m.fs.ip_pump.outlet,
        destination=m.fs.ip_econ1.tube_inlet)

    m.fs.hp_econ1_to_hp_econ2 = Arc(
        source=m.fs.hp_econ1.tube_outlet,
        destination=m.fs.hp_econ2.tube_inlet)

    m.fs.ip_econ1_to_ip_econ2 = Arc(
        source=m.fs.ip_econ1.tube_outlet,
        destination=m.fs.ip_econ2.tube_inlet)

    m.fs.ip_econ2_to_ip_drum = Arc(
        source=m.fs.ip_econ2.tube_outlet,
        destination=m.fs.ip_drum.feedwater_inlet)

    m.fs.ip_drum_to_ip_downcomer = Arc(
        source=m.fs.ip_drum.liquid_outlet,
        destination=m.fs.ip_downcomer.inlet)

    m.fs.ip_downcomer_to_ip_evap = Arc(
        source=m.fs.ip_downcomer.outlet,
        destination=m.fs.ip_evap.tube_inlet)

    m.fs.ip_evap_to_ip_drum = Arc(
        source=m.fs.ip_evap.tube_outlet,
        destination=m.fs.ip_drum.water_steam_inlet)

    m.fs.lp_da_split_to_lp_sh = Arc(
        source=m.fs.lp_da_split.lp_outlet,
        destination=m.fs.lp_sh.tube_inlet)

    m.fs.hp_econ2_to_hp_econ3 = Arc(
        source=m.fs.hp_econ2.tube_outlet,
        destination=m.fs.hp_econ3.tube_inlet)

    m.fs.ip_drum_to_ip_sh1 = Arc(
        source=m.fs.ip_drum.steam_outlet,
        destination=m.fs.ip_sh1.tube_inlet)

    m.fs.hp_econ3_to_hp_econ4 = Arc(
        source=m.fs.hp_econ3.tube_outlet,
        destination=m.fs.hp_econ4.tube_inlet)

    m.fs.hp_econ4_to_hp_drum = Arc(
        source=m.fs.hp_econ4.tube_outlet,
        destination=m.fs.hp_drum.feedwater_inlet)

    m.fs.hp_drum_to_hp_downcomer = Arc(
        source=m.fs.hp_drum.liquid_outlet,
        destination=m.fs.hp_downcomer.inlet)

    m.fs.hp_downcomer_to_hp_evap = Arc(
        source=m.fs.hp_downcomer.outlet,
        destination=m.fs.hp_evap.tube_inlet)

    m.fs.hp_evap_to_hp_drum = Arc(
        source=m.fs.hp_evap.tube_outlet,
        destination=m.fs.hp_drum.water_steam_inlet)

    m.fs.hp_drum_to_hp_sh1 = Arc(
        source=m.fs.hp_drum.steam_outlet,
        destination=m.fs.hp_sh1.tube_inlet)

    m.fs.ip_sh1_to_ip_mixer = Arc(
        source=m.fs.ip_sh1.tube_outlet,
        destination=m.fs.ip_mixer.hrsg_inlet)

    m.fs.ip_mixer_to_ip_sh2 = Arc(
        source=m.fs.ip_mixer.outlet,
        destination=m.fs.ip_sh2.tube_inlet)

    m.fs.hp_sh1_to_hp_sh2 = Arc(
        source=m.fs.hp_sh1.tube_outlet,
        destination=m.fs.hp_sh2.tube_inlet)

    m.fs.ip_sh2_to_ip_sh3 = Arc(
        source=m.fs.ip_sh2.tube_outlet,
        destination=m.fs.ip_sh3.tube_inlet)

    m.fs.hp_sh2_to_hp_sh3 = Arc(
        source=m.fs.hp_sh2.tube_outlet,
        destination=m.fs.hp_sh3.tube_inlet)

    # Connection between hrsg and turbine
    m.fs.hp_sh3_to_steam_turb = Arc(
        source=m.fs.hp_sh3.tube_outlet,
        destination=m.fs.steam_turb.inlet_split.inlet)

    m.fs.steam_turb_to_ip_mixer = Arc(
        source=m.fs.steam_turb.hp_stages[7].outlet,
        destination=m.fs.ip_mixer.turb_inlet)

    m.fs.ip_sh3_to_steam_turb = Arc(
        source=m.fs.ip_sh3.tube_outlet,
        destination=m.fs.steam_turb.ip_stages[1].inlet)

    m.fs.steam_turb_to_lp_mixer = Arc(
        source=m.fs.steam_turb.ip_stages[10].outlet,
        destination=m.fs.lp_mixer.turb_inlet)

    m.fs.lp_sh_to_lp_mixer = Arc(
        source=m.fs.lp_sh.tube_outlet,
        destination=m.fs.lp_mixer.hrsg_inlet)

    m.fs.lp_mixer_to_lp_splitter = Arc(
        source=m.fs.lp_mixer.outlet,
        destination=m.fs.lp_splitter.inlet)

    m.fs.lp_splitter_to_steam_turb = Arc(
        source=m.fs.lp_splitter.turb_outlet,
        destination=m.fs.steam_turb.lp_stages[1].inlet)

    m.fs.steam_turb_to_condenser = Arc(
        source=m.fs.steam_turb.outlet_stage.outlet,
        destination=m.fs.condenser.shell_inlet)

    m.fs.condenser_to_hotwell = Arc(
        source=m.fs.condenser.shell_outlet,
        destination=m.fs.hotwell.condensate_inlet)

    m.fs.hotwell_to_cond_pump = Arc(
        source=m.fs.hotwell.outlet,
        destination=m.fs.cond_pump.inlet)

    m.fs.cond_pump_to_lp_knockout_condenser = Arc(
        source=m.fs.cond_pump.outlet,
        destination=m.fs.lp_knockout_condenser.tube_inlet)

    m.fs.lp_knockout_condenser_to_lp_econ = Arc(
        source=m.fs.lp_knockout_condenser.tube_outlet,
        destination=m.fs.lp_econ.tube_inlet)

    # co2 and h2 compression
    m.fs.co2_post_separator_to_co2_translator_flue = Arc(
        source=m.fs.co2_post_separator.co2_outlet,
        destination=m.fs.co2_translator_flue.inlet)

    m.fs.co2_pre_separator_to_co2_translator_syn = Arc(
        source=m.fs.co2_pre_separator.co2_outlet,
        destination=m.fs.co2_translator_syn.inlet)

    m.fs.co2_translator_flue_to_co2_mix = Arc(
        source=m.fs.co2_translator_flue.outlet,
        destination=m.fs.co2_mix.post_inlet)

    m.fs.co2_translator_syn_to_co2_mix = Arc(
        source=m.fs.co2_translator_syn.outlet,
        destination=m.fs.co2_mix.pre_inlet)

    m.fs.co2_mix_to_co2_comp1 = Arc(
        source=m.fs.co2_mix.outlet,
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

    m.fs.h2_separator_to_h2_comp = Arc(
        source=m.fs.h2_separator.h2_outlet,
        destination=m.fs.h2_comp.inlet)


def set_blue_hydrogen_plant_inputs(m):
    # input data related to unit operations
    m.fs.deltaT_smr.fix(50.0)
    # ng preheat
    m.fs.ng_preheat2.tube.deltaP.fix(-40000)
    m.fs.ng_preheat2.shell.deltaP.fix(-100)
    m.fs.ng_preheat2.area.fix(4000)
    m.fs.ng_preheat2.overall_heat_transfer_coefficient.fix(50)
    # smr reformer tube side
    m.fs.smr_tube.deltaP.fix(-1.7e5)
    # smr reformer shell side
    m.fs.smr_shell.deltaP.fix(-500)
    # syngas cooler
    m.fs.syngas_cooler.tube.deltaP.fix(-30000)
    m.fs.syngas_cooler.shell.deltaP.fix(-30000)
    m.fs.syngas_cooler.area.fix(12000)
    m.fs.syngas_cooler.overall_heat_transfer_coefficient.fix(80)
    # water gas shift reactor
    m.fs.wgsr.x_co.fix(0.96)
    m.fs.wgsr.deltaP.fix(-2.9e5)
    # water gas shift reactor cooler
    m.fs.wgsr_cooler.tube.deltaP.fix(-30000)
    m.fs.wgsr_cooler.shell.deltaP.fix(-30000)
    m.fs.wgsr_cooler.area.fix(3000)
    m.fs.wgsr_cooler.overall_heat_transfer_coefficient.fix(100)
    # reboiler water pump
    m.fs.rb_pump.efficiency_isentropic.fix(0.80)
    m.fs.rb_pump.deltaP.fix(0.58e5)
    # rb_splitter
    m.fs.rb_splitter.split_fraction[:, "cond_outlet"].fix(0.15)
    # rb water knockout condenser
    m.fs.rb_knockout_condenser.deltaP_tube.fix(-10000)
    m.fs.rb_knockout_condenser.deltaP_shell.fix(-5000)
    m.fs.rb_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.rb_knockout_condenser.area.fix(7000)
    # mp makeup water pump
    m.fs.mp_pump.efficiency_isentropic.fix(0.80)
    m.fs.mp_pump.deltaP.fix(2.36e6)
    # mp booster pump
    m.fs.mp_booster_pump.efficiency_isentropic.fix(0.80)
    m.fs.mp_booster_pump.deltaP.fix(7.1e5)
    # lp water knockout condenser
    m.fs.lp_knockout_condenser.deltaP_tube.fix(-10000)
    m.fs.lp_knockout_condenser.deltaP_shell.fix(-5000)
    m.fs.lp_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.lp_knockout_condenser.area.fix(1000) # make it very small
    # mp water knockout condenser
    m.fs.mp_knockout_condenser.deltaP_tube.fix(-10000)
    m.fs.mp_knockout_condenser.deltaP_shell.fix(-10000)
    m.fs.mp_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.mp_knockout_condenser.area.fix(8000) #>=7000 causes convergence problem
    # cooling water knockout condenser
    m.fs.cw_knockout_condenser.deltaP_tube.fix(-2000)
    m.fs.cw_knockout_condenser.deltaP_shell.fix(0)
    m.fs.cw_knockout_condenser.overall_heat_transfer_coefficient.fix(80)
    m.fs.cw_knockout_condenser.area.fix(9000)
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
    # H2 separator
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'H2'].fix(0.85)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'H2O'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CO2'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'CH4'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'N2'].fix(0.03)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'O2'].fix(0)
    m.fs.h2_separator.split_fraction[:, 'h2_outlet', 'Ar'].fix(0.05)
    
    # input data or initial guess for streams
    # natural gas feed to ng preheat1
    m.fs.ng_preheat2.tube_inlet.flow_mol.fix(1024.17)  # mol/s
    m.fs.ng_preheat2.tube_inlet.temperature.fix(417)  # K
    m.fs.ng_preheat2.tube_inlet.pressure.fix(3.04e6)  # Pa,
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'CH4'].fix(0.974)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'CO'].fix(1e-6)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'CO2'].fix(0.01)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'H2'].fix(1e-6)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'H2O'].fix(1e-6)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'N2'].fix(0.016)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'O2'].fix(1e-6)
    m.fs.ng_preheat2.tube_inlet.mole_frac_comp[:,'Ar'].fix(0.000001) # Add 1 ppm to avoid zero Ar in Gibbs 
    # flue gas to ng preheat guess
    m.fs.ng_preheat2.shell_inlet.temperature[:].value = 1205  # K
    m.fs.ng_preheat2.shell_inlet.pressure[:].value = 104000  # Pa,
    m.fs.ng_preheat2.shell_inlet.flow_mol_comp[:, "H2O"].value = 2953.6
    m.fs.ng_preheat2.shell_inlet.flow_mol_comp[:, "CO2"].value = 1176.8
    m.fs.ng_preheat2.shell_inlet.flow_mol_comp[:, "N2"].value = 12993
    m.fs.ng_preheat2.shell_inlet.flow_mol_comp[:, "O2"].value = 974.7
    m.fs.ng_preheat2.shell_inlet.flow_mol_comp[:, "Ar"].value = 0.0213
    # steam feed to tube side of SMR guess
    m.fs.smr_tube_mix.steam_inlet_state[0].flow_mol.value = 2730.6  # mol/s
    m.fs.smr_tube_mix.steam_inlet_state[0].temperature.value = 856.6  # K
    m.fs.smr_tube_mix.steam_inlet_state[0].pressure.value = 3.11e6  # Pa,
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CH4'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['CO2'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['H2O'].value = 0.99999
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['N2'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['O2'].value = 1e-6
    m.fs.smr_tube_mix.steam_inlet_state[0].mole_frac_comp['Ar'].value = 1e-6
    # natural gas feed to shell side of SMR
    m.fs.smr_shell_mix.fuel_inlet_state[0].flow_mol.fix(162.3)
    m.fs.smr_shell_mix.fuel_inlet_state[0].temperature.fix(288)  # K
    m.fs.smr_shell_mix.fuel_inlet_state[0].pressure.fix(1.05e5)  # Pa,
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CH4'].fix(0.974)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO'].fix(1e-6)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['CO2'].fix(0.01)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2'].fix(1e-6)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['H2O'].fix(1e-6)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['N2'].fix(0.016)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['O2'].fix(1e-6)
    m.fs.smr_shell_mix.fuel_inlet_state[0].mole_frac_comp['Ar'].fix(1e-6)
    # flue gas from gas turbine, replacing air
    m.fs.smr_shell_mix.air_inlet_state[0].flow_mol.value = 17450  # mol/s
    m.fs.smr_shell_mix.air_inlet_state[0].temperature.value = 955.95  # K
    m.fs.smr_shell_mix.air_inlet_state[0].pressure.value = 1.04446e5  # Pa,
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CH4'].value = 6.66e-8
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO'].value = 6.99e-6
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['CO2'].value = 0.04465
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2'].value = 4.23e-6
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['H2O'].value = 0.1027
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['N2'].value = 0.74348
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['O2'].value = 0.10913
    m.fs.smr_shell_mix.air_inlet_state[0].mole_frac_comp['Ar'].value = 1e-6
    # recycled tailgas intial guess
    m.fs.smr_shell_mix.tailgas_inlet.flow_mol[:].value = 728.9
    m.fs.smr_shell_mix.tailgas_inlet.temperature[:].value = 317
    m.fs.smr_shell_mix.tailgas_inlet.pressure[:].value = 2460000
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CH4'].value = 0.2455
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO'].value =  0.03028
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'CO2'].value = 0.0332
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2'].value = 0.6692
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'H2O'].value = 1e-6
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'N2'].value = 0.02181
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'O2'].value = 1e-6
    m.fs.smr_shell_mix.tailgas_inlet.mole_frac_comp[:,'Ar'].value = 4.89e-6
    # feed water to syngas cooler guess
    m.fs.syngas_cooler.tube_inlet.flow_mol[:].value = 2730.6
    m.fs.syngas_cooler.tube_inlet.pressure[:].value = 3.14e6 
    enth_mol_fw = 26016 #iapws95.htpx(T=509.55*pyo.units.K, P=3.1399e6*pyo.units.Pa) #actually saturated but set to subcooled
    m.fs.syngas_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw
    # feed water to water gas shift reactor cooler guess
    m.fs.wgsr_cooler.tube_inlet.flow_mol[:].value = 2730.6
    m.fs.wgsr_cooler.tube_inlet.pressure[:].value = 3.17e6
    enth_mol_fw3 = 13052 #iapws95.htpx(T=380.9*pyo.units.K, P=3.1699e6*pyo.units.Pa)
    m.fs.wgsr_cooler.tube_inlet.enth_mol[:].value = enth_mol_fw3
    # reboiler saturated liquid return
    m.fs.rb_pump.inlet.flow_mol[:].value = 6000 #6000
    m.fs.rb_pump.inlet.pressure[:].value = 2.6e5
    enth_mol_rb = iapws95.htpx(T=395*pyo.units.K, P=2.6e5*pyo.units.Pa)
    m.fs.rb_pump.inlet.enth_mol[:].value = enth_mol_rb
    m.fs.rb_knockout_condenser.heat_duty[:].value = 2.3e7
    # lp knockout condenser heat duty initial guess
    m.fs.lp_knockout_condenser.tube_inlet.flow_mol.fix(1501)
    m.fs.lp_knockout_condenser.tube_inlet.pressure.fix(3.14e5)
    enth_mol_fw2 = 2175 #iapws95.htpx(T=300*pyo.units.K, P=3.14e5*pyo.units.Pa)
    m.fs.lp_knockout_condenser.tube_inlet.enth_mol.fix(enth_mol_fw2)
    m.fs.lp_knockout_condenser.heat_duty[:].value = 2.7e7
    # mp makeup water to mp water knockout condenser
    m.fs.mp_pump.inlet.flow_mol.fix(1625.5) #2731 total to SMR after mixed with condensate
    m.fs.mp_pump.inlet.pressure.fix(1e5)
    enth_mol_fw3 = iapws95.htpx(T=300*pyo.units.K, P=1e5*pyo.units.Pa)
    m.fs.mp_pump.inlet.enth_mol.fix(enth_mol_fw3)
    # mp knockout condenser heat duty initial guess
    m.fs.mp_knockout_condenser.heat_duty[:].value = 2.5e7
    # cooling water to cooling water knockout condenser
    m.fs.cw_knockout_condenser.tube_inlet.flow_mol.fix(10000)
    m.fs.cw_knockout_condenser.tube_inlet.pressure.fix(1.5e5)
    enth_mol_fw4 = iapws95.htpx(T=293.15*pyo.units.K, P=1.5e5*pyo.units.Pa)
    m.fs.cw_knockout_condenser.tube_inlet.enth_mol.fix(enth_mol_fw4)
    # cooling water knockout condenser heat duty initial guess
    m.fs.cw_knockout_condenser.heat_duty[:].value = 1e7

def set_gas_turbine_plant_inputs(m):
    # set flow scale of performance curves
    m.fs.performance_flow_scale.fix(2/1.0)
    # air to gas turbine
    m.fs.comp_igv.inlet.flow_mol.fix(16850*1.0)
    m.fs.comp_igv.inlet.temperature.fix(298.15)
    m.fs.comp_igv.inlet.pressure.fix(101325)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CH4'].fix(1e-6)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CO'].fix(1e-6)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'CO2'].fix(0.000335)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'H2'].fix(1e-6)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'H2O'].fix(0.015653)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'N2'].fix(0.777811)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'O2'].fix(0.206201)
    m.fs.comp_igv.inlet.mole_frac_comp[:,'Ar'].fix(1e-6)

    m.fs.gt_mix.fuel_inlet.flow_mol.fix(795*1.0)
    m.fs.gt_mix.fuel_inlet.temperature.fix(298.15)
    m.fs.gt_mix.fuel_inlet.pressure.fix(2e6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CH4'].fix(0.974)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CO'].fix(1e-6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['CO2'].fix(0.01)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['H2'].fix(1e-6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['H2O'].fix(1e-6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['N2'].fix(0.016)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['O2'].fix(1e-6)
    m.fs.gt_mix.fuel_inlet_state[0].mole_frac_comp['Ar'].fix(1e-6)

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

def set_hrsg_steam_turbine_plant_inputs(m):
    # ng_preheat1
    m.fs.ng_preheat1.tube.deltaP.fix(-1000)
    m.fs.ng_preheat1.shell.deltaP.fix(-100)
    m.fs.ng_preheat1.area.fix(6000)
    m.fs.ng_preheat1.overall_heat_transfer_coefficient.fix(50)
    m.fs.ng_preheat1.tube_inlet.flow_mol.fix(1024.17)  # mol/s
    m.fs.ng_preheat1.tube_inlet.temperature.fix(288)  # K
    m.fs.ng_preheat1.tube_inlet.pressure.fix(3.041e6)  # Pa,
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CH4'].fix(0.974)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO'].fix(1e-6)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'CO2'].fix(0.01)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2'].fix(1e-6)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'H2O'].fix(1e-6)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'N2'].fix(0.016)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'O2'].fix(1e-6)
    m.fs.ng_preheat1.tube_inlet.mole_frac_comp[:,'Ar'].fix(0.000001) # Add 1 ppm to avoid zero Ar in Gibbs 
    
    # hrsg and steam turbine inputs
    m.fs.lp_econ.tube_inlet.flow_mol[:].value = 1461
    m.fs.lp_econ.tube_inlet.pressure[:].value = 3.1e5
    enth_mol_lp_econ_in = iapws95.htpx(T=400*pyo.units.K, P=3.1e5*pyo.units.Pa)
    m.fs.lp_econ.tube_inlet.enth_mol[:].value = enth_mol_lp_econ_in

    m.fs.ip_mixer.turb_inlet.flow_mol[:].value = 4976.2
    m.fs.ip_mixer.turb_inlet.pressure[:].value = 4578539
    enth_mol_cold_ip_in = iapws95.htpx(T=630.2*pyo.units.K, P=4578539*pyo.units.Pa)
    m.fs.ip_mixer.turb_inlet.enth_mol[:].value = enth_mol_cold_ip_in
    
    # flue gas to lp econ
    m.fs.lp_econ.shell_inlet.temperature[:].value = 438.3
    m.fs.lp_econ.shell_inlet.pressure[:].value = 100265  # Pa,
    m.fs.lp_econ.shell_inlet.flow_mol_comp[:, "H2O"].value = 2921
    m.fs.lp_econ.shell_inlet.flow_mol_comp[:, "CO2"].value = 1164.2
    m.fs.lp_econ.shell_inlet.flow_mol_comp[:, "N2"].value = 12290.6
    m.fs.lp_econ.shell_inlet.flow_mol_comp[:, "O2"].value = 811.83
    m.fs.lp_econ.shell_inlet.flow_mol_comp[:, "Ar"].value = 0.02045
    #m.fs.lp_econ.shell_inlet.unfix()

    m.fs.lp_econ.di_tube.fix(0.042)
    m.fs.lp_econ.thickness_tube.fix(0.004)
    m.fs.lp_econ.pitch_x.fix(0.1)
    m.fs.lp_econ.pitch_y.fix(0.1)
    m.fs.lp_econ.length_tube_seg.fix(10)
    m.fs.lp_econ.nseg_tube.fix(5)
    m.fs.lp_econ.ncol_tube.fix(100)
    m.fs.lp_econ.nrow_inlet.fix(2)
    m.fs.lp_econ.delta_elevation.fix(0)
    m.fs.lp_econ.therm_cond_wall = 43.0
    m.fs.lp_econ.rfouling_tube = 0.0001
    m.fs.lp_econ.rfouling_shell = 0.0001
    m.fs.lp_econ.fcorrection_htc_tube.fix(1)
    m.fs.lp_econ.fcorrection_htc_shell.fix(1)
    m.fs.lp_econ.fcorrection_dp_tube.fix(1)
    m.fs.lp_econ.fcorrection_dp_shell.fix(1)

    m.fs.lp_ip_hp_split.split_fraction[:, "lp_outlet"].fix(0.09)
    m.fs.lp_ip_hp_split.split_fraction[:, "ip_outlet"].fix(0.13)

    m.fs.lp_drum.drum_diameter.fix(1)
    m.fs.lp_drum.drum_length.fix(8.0)
    m.fs.lp_drum.drum_level.fix(0.5)
    m.fs.lp_drum.number_downcomers.fix(2)
    m.fs.lp_drum.downcomer_diameter.fix(0.5)
    m.fs.lp_drum.heat_duty[:].fix(0.0)

    m.fs.lp_downcomer.diameter.fix(0.15)
    m.fs.lp_downcomer.height.fix(10)
    m.fs.lp_downcomer.number_downcomers.fix(2)
    m.fs.lp_downcomer.heat_duty[:].fix(0.0)
    
    m.fs.lp_evap.tube_diameter_inner.fix(0.042)
    m.fs.lp_evap.tube_thickness.fix(0.004)
    m.fs.lp_evap.pitch_x.fix(0.1)
    m.fs.lp_evap.pitch_y.fix(0.1)
    m.fs.lp_evap.height.fix(10)
    m.fs.lp_evap.number_tube_rows.fix(100)
    m.fs.lp_evap.number_tube_cols.fix(100)
    m.fs.lp_evap.fcorrection_dp_tube.fix(3.0)
    m.fs.lp_evap.fcorrection_dp_shell.fix(1.0)

    m.fs.lp_da_split.split_fraction[:, "purge_outlet"].fix(0.01)

    m.fs.ip_pump.efficiency_isentropic.fix(0.80)
    m.fs.ip_pump.deltaP.fix(4.2e+006)

    m.fs.hp_pump.efficiency_isentropic.fix(0.80)
    m.fs.hp_pump.deltaP.fix(1.8e7)

    m.fs.hp_econ1.di_tube.fix(0.042)
    m.fs.hp_econ1.thickness_tube.fix(0.004)
    m.fs.hp_econ1.pitch_x.fix(0.1)
    m.fs.hp_econ1.pitch_y.fix(0.1)
    m.fs.hp_econ1.length_tube_seg.fix(10)
    m.fs.hp_econ1.nseg_tube.fix(15)
    m.fs.hp_econ1.ncol_tube.fix(100)
    m.fs.hp_econ1.nrow_inlet.fix(1)
    m.fs.hp_econ1.delta_elevation.fix(0)
    m.fs.hp_econ1.therm_cond_wall = 43.0
    m.fs.hp_econ1.rfouling_tube = 0.0001
    m.fs.hp_econ1.rfouling_shell = 0.0001
    m.fs.hp_econ1.fcorrection_htc_tube.fix(1)
    m.fs.hp_econ1.fcorrection_htc_shell.fix(1)
    m.fs.hp_econ1.fcorrection_dp_tube.fix(1)
    m.fs.hp_econ1.fcorrection_dp_shell.fix(1)

    m.fs.ip_econ1.di_tube.fix(0.042)
    m.fs.ip_econ1.thickness_tube.fix(0.004)
    m.fs.ip_econ1.pitch_x.fix(0.1)
    m.fs.ip_econ1.pitch_y.fix(0.1)
    m.fs.ip_econ1.length_tube_seg.fix(10)
    m.fs.ip_econ1.nseg_tube.fix(4)#5
    m.fs.ip_econ1.ncol_tube.fix(100)
    m.fs.ip_econ1.nrow_inlet.fix(1)
    m.fs.ip_econ1.delta_elevation.fix(0)
    m.fs.ip_econ1.therm_cond_wall = 43.0
    m.fs.ip_econ1.rfouling_tube = 0.0001
    m.fs.ip_econ1.rfouling_shell = 0.0001
    m.fs.ip_econ1.fcorrection_htc_tube.fix(1)
    m.fs.ip_econ1.fcorrection_htc_shell.fix(1)
    m.fs.ip_econ1.fcorrection_dp_tube.fix(1)
    m.fs.ip_econ1.fcorrection_dp_shell.fix(1)

    m.fs.hp_econ2.di_tube.fix(0.042)
    m.fs.hp_econ2.thickness_tube.fix(0.004)
    m.fs.hp_econ2.pitch_x.fix(0.1)
    m.fs.hp_econ2.pitch_y.fix(0.1)
    m.fs.hp_econ2.length_tube_seg.fix(10)
    m.fs.hp_econ2.nseg_tube.fix(20)
    m.fs.hp_econ2.ncol_tube.fix(100)
    m.fs.hp_econ2.nrow_inlet.fix(1)
    m.fs.hp_econ2.delta_elevation.fix(0)
    m.fs.hp_econ2.therm_cond_wall = 43.0
    m.fs.hp_econ2.rfouling_tube = 0.0001
    m.fs.hp_econ2.rfouling_shell = 0.0001
    m.fs.hp_econ2.fcorrection_htc_tube.fix(1)
    m.fs.hp_econ2.fcorrection_htc_shell.fix(1)
    m.fs.hp_econ2.fcorrection_dp_tube.fix(1)
    m.fs.hp_econ2.fcorrection_dp_shell.fix(1)

    m.fs.ip_econ2.di_tube.fix(0.042)
    m.fs.ip_econ2.thickness_tube.fix(0.004)
    m.fs.ip_econ2.pitch_x.fix(0.1)
    m.fs.ip_econ2.pitch_y.fix(0.1)
    m.fs.ip_econ2.length_tube_seg.fix(10)
    m.fs.ip_econ2.nseg_tube.fix(4)
    m.fs.ip_econ2.ncol_tube.fix(100)
    m.fs.ip_econ2.nrow_inlet.fix(1)
    m.fs.ip_econ2.delta_elevation.fix(0)
    m.fs.ip_econ2.therm_cond_wall = 43.0
    m.fs.ip_econ2.rfouling_tube = 0.0001
    m.fs.ip_econ2.rfouling_shell = 0.0001
    m.fs.ip_econ2.fcorrection_htc_tube.fix(1)
    m.fs.ip_econ2.fcorrection_htc_shell.fix(1)
    m.fs.ip_econ2.fcorrection_dp_tube.fix(1)
    m.fs.ip_econ2.fcorrection_dp_shell.fix(1)

    m.fs.ip_drum.drum_diameter.fix(1)
    m.fs.ip_drum.drum_length.fix(8.0)
    m.fs.ip_drum.drum_level.fix(0.5)
    m.fs.ip_drum.number_downcomers.fix(2)
    m.fs.ip_drum.downcomer_diameter.fix(0.5)
    m.fs.ip_drum.heat_duty[:].fix(0.0)

    m.fs.ip_downcomer.diameter.fix(0.15)
    m.fs.ip_downcomer.height.fix(10)
    m.fs.ip_downcomer.number_downcomers.fix(2)
    m.fs.ip_downcomer.heat_duty[:].fix(0.0)
    
    m.fs.ip_evap.tube_diameter_inner.fix(0.042)
    m.fs.ip_evap.tube_thickness.fix(0.004)
    m.fs.ip_evap.pitch_x.fix(0.1)
    m.fs.ip_evap.pitch_y.fix(0.1)
    m.fs.ip_evap.height.fix(10)
    m.fs.ip_evap.number_tube_rows.fix(50)
    m.fs.ip_evap.number_tube_cols.fix(100)
    m.fs.ip_evap.fcorrection_dp_tube.fix(3.0)
    m.fs.ip_evap.fcorrection_dp_shell.fix(1.0)

    m.fs.lp_sh.di_tube.fix(0.042)
    m.fs.lp_sh.thickness_tube.fix(0.004)
    m.fs.lp_sh.pitch_x.fix(0.1)
    m.fs.lp_sh.pitch_y.fix(0.1)
    m.fs.lp_sh.length_tube_seg.fix(10)
    m.fs.lp_sh.nseg_tube.fix(8)#5
    m.fs.lp_sh.ncol_tube.fix(100)
    m.fs.lp_sh.nrow_inlet.fix(4)
    m.fs.lp_sh.delta_elevation.fix(0)
    m.fs.lp_sh.therm_cond_wall = 43.0
    m.fs.lp_sh.rfouling_tube = 0.0001
    m.fs.lp_sh.rfouling_shell = 0.0001
    m.fs.lp_sh.fcorrection_htc_tube.fix(1)
    m.fs.lp_sh.fcorrection_htc_shell.fix(1)
    m.fs.lp_sh.fcorrection_dp_tube.fix(1)
    m.fs.lp_sh.fcorrection_dp_shell.fix(1)

    m.fs.hp_econ3.di_tube.fix(0.042)
    m.fs.hp_econ3.thickness_tube.fix(0.004)
    m.fs.hp_econ3.pitch_x.fix(0.1)
    m.fs.hp_econ3.pitch_y.fix(0.1)
    m.fs.hp_econ3.length_tube_seg.fix(10)
    m.fs.hp_econ3.nseg_tube.fix(20)
    m.fs.hp_econ3.ncol_tube.fix(100)
    m.fs.hp_econ3.nrow_inlet.fix(1)
    m.fs.hp_econ3.delta_elevation.fix(0)
    m.fs.hp_econ3.therm_cond_wall = 43.0
    m.fs.hp_econ3.rfouling_tube = 0.0001
    m.fs.hp_econ3.rfouling_shell = 0.0001
    m.fs.hp_econ3.fcorrection_htc_tube.fix(1)
    m.fs.hp_econ3.fcorrection_htc_shell.fix(1)
    m.fs.hp_econ3.fcorrection_dp_tube.fix(1)
    m.fs.hp_econ3.fcorrection_dp_shell.fix(1)

    m.fs.ip_sh1.di_tube.fix(0.042)
    m.fs.ip_sh1.thickness_tube.fix(0.004)
    m.fs.ip_sh1.pitch_x.fix(0.1)
    m.fs.ip_sh1.pitch_y.fix(0.1)
    m.fs.ip_sh1.length_tube_seg.fix(10)
    m.fs.ip_sh1.nseg_tube.fix(2)
    m.fs.ip_sh1.ncol_tube.fix(100)
    m.fs.ip_sh1.nrow_inlet.fix(1)
    m.fs.ip_sh1.delta_elevation.fix(0)
    m.fs.ip_sh1.therm_cond_wall = 43.0
    m.fs.ip_sh1.rfouling_tube = 0.0001
    m.fs.ip_sh1.rfouling_shell = 0.0001
    m.fs.ip_sh1.fcorrection_htc_tube.fix(1)
    m.fs.ip_sh1.fcorrection_htc_shell.fix(1)
    m.fs.ip_sh1.fcorrection_dp_tube.fix(1)
    m.fs.ip_sh1.fcorrection_dp_shell.fix(1)

    m.fs.hp_econ4.di_tube.fix(0.042)
    m.fs.hp_econ4.thickness_tube.fix(0.004)
    m.fs.hp_econ4.pitch_x.fix(0.1)
    m.fs.hp_econ4.pitch_y.fix(0.1)
    m.fs.hp_econ4.length_tube_seg.fix(10)
    m.fs.hp_econ4.nseg_tube.fix(10)
    m.fs.hp_econ4.ncol_tube.fix(100)
    m.fs.hp_econ4.nrow_inlet.fix(1)
    m.fs.hp_econ4.delta_elevation.fix(0)
    m.fs.hp_econ4.therm_cond_wall = 43.0
    m.fs.hp_econ4.rfouling_tube = 0.0001
    m.fs.hp_econ4.rfouling_shell = 0.0001
    m.fs.hp_econ4.fcorrection_htc_tube.fix(1)
    m.fs.hp_econ4.fcorrection_htc_shell.fix(1)
    m.fs.hp_econ4.fcorrection_dp_tube.fix(1)
    m.fs.hp_econ4.fcorrection_dp_shell.fix(1)

    m.fs.hp_drum.drum_diameter.fix(2)
    m.fs.hp_drum.drum_length.fix(8.0)
    m.fs.hp_drum.drum_level.fix(1)
    m.fs.hp_drum.number_downcomers.fix(4) # use large number of downcomer to reduce velocity head
    m.fs.hp_drum.downcomer_diameter.fix(1) # use a larger diameter to reduce velocity head
    m.fs.hp_drum.heat_duty[:].fix(0.0)

    m.fs.hp_downcomer.diameter.fix(0.3)
    m.fs.hp_downcomer.height.fix(10)
    m.fs.hp_downcomer.number_downcomers.fix(2)
    m.fs.hp_downcomer.heat_duty[:].fix(0.0)
    
    m.fs.hp_evap.tube_diameter_inner.fix(0.042)
    m.fs.hp_evap.tube_thickness.fix(0.004)
    m.fs.hp_evap.pitch_x.fix(0.1)
    m.fs.hp_evap.pitch_y.fix(0.1)
    m.fs.hp_evap.height.fix(10)
    m.fs.hp_evap.number_tube_rows.fix(45)
    m.fs.hp_evap.number_tube_cols.fix(100)
    m.fs.hp_evap.fcorrection_dp_tube.fix(1.0)
    m.fs.hp_evap.fcorrection_dp_shell.fix(1.0)

    m.fs.hp_sh1.di_tube.fix(0.042)
    m.fs.hp_sh1.thickness_tube.fix(0.004)
    m.fs.hp_sh1.pitch_x.fix(0.1)
    m.fs.hp_sh1.pitch_y.fix(0.1)
    m.fs.hp_sh1.length_tube_seg.fix(10)
    m.fs.hp_sh1.nseg_tube.fix(12)
    m.fs.hp_sh1.ncol_tube.fix(100)
    m.fs.hp_sh1.nrow_inlet.fix(1)
    m.fs.hp_sh1.delta_elevation.fix(0)
    m.fs.hp_sh1.therm_cond_wall = 43.0
    m.fs.hp_sh1.rfouling_tube = 0.0001
    m.fs.hp_sh1.rfouling_shell = 0.0001
    m.fs.hp_sh1.fcorrection_htc_tube.fix(1)
    m.fs.hp_sh1.fcorrection_htc_shell.fix(1)
    m.fs.hp_sh1.fcorrection_dp_tube.fix(1)
    m.fs.hp_sh1.fcorrection_dp_shell.fix(1)

    m.fs.ip_sh2.di_tube.fix(0.042)
    m.fs.ip_sh2.thickness_tube.fix(0.004)
    m.fs.ip_sh2.pitch_x.fix(0.1)
    m.fs.ip_sh2.pitch_y.fix(0.1)
    m.fs.ip_sh2.length_tube_seg.fix(10)
    m.fs.ip_sh2.nseg_tube.fix(6)
    m.fs.ip_sh2.ncol_tube.fix(100)
    m.fs.ip_sh2.nrow_inlet.fix(1)
    m.fs.ip_sh2.delta_elevation.fix(0)
    m.fs.ip_sh2.therm_cond_wall = 43.0
    m.fs.ip_sh2.rfouling_tube = 0.0001
    m.fs.ip_sh2.rfouling_shell = 0.0001
    m.fs.ip_sh2.fcorrection_htc_tube.fix(1)
    m.fs.ip_sh2.fcorrection_htc_shell.fix(1)
    m.fs.ip_sh2.fcorrection_dp_tube.fix(1)
    m.fs.ip_sh2.fcorrection_dp_shell.fix(1)

    m.fs.hp_sh2.di_tube.fix(0.042)
    m.fs.hp_sh2.thickness_tube.fix(0.004)
    m.fs.hp_sh2.pitch_x.fix(0.1)
    m.fs.hp_sh2.pitch_y.fix(0.1)
    m.fs.hp_sh2.length_tube_seg.fix(10)
    m.fs.hp_sh2.nseg_tube.fix(4)
    m.fs.hp_sh2.ncol_tube.fix(100)
    m.fs.hp_sh2.nrow_inlet.fix(1)
    m.fs.hp_sh2.delta_elevation.fix(0)
    m.fs.hp_sh2.therm_cond_wall = 43.0
    m.fs.hp_sh2.rfouling_tube = 0.0001
    m.fs.hp_sh2.rfouling_shell = 0.0001
    m.fs.hp_sh2.fcorrection_htc_tube.fix(1)
    m.fs.hp_sh2.fcorrection_htc_shell.fix(1)
    m.fs.hp_sh2.fcorrection_dp_tube.fix(1)
    m.fs.hp_sh2.fcorrection_dp_shell.fix(1)

    m.fs.ip_sh3.di_tube.fix(0.042)
    m.fs.ip_sh3.thickness_tube.fix(0.004)
    m.fs.ip_sh3.pitch_x.fix(0.1)
    m.fs.ip_sh3.pitch_y.fix(0.1)
    m.fs.ip_sh3.length_tube_seg.fix(10)
    m.fs.ip_sh3.nseg_tube.fix(6)#5
    m.fs.ip_sh3.ncol_tube.fix(100)
    m.fs.ip_sh3.nrow_inlet.fix(1)
    m.fs.ip_sh3.delta_elevation.fix(0)
    m.fs.ip_sh3.therm_cond_wall = 43.0
    m.fs.ip_sh3.rfouling_tube = 0.0001
    m.fs.ip_sh3.rfouling_shell = 0.0001
    m.fs.ip_sh3.fcorrection_htc_tube.fix(1)
    m.fs.ip_sh3.fcorrection_htc_shell.fix(1)
    m.fs.ip_sh3.fcorrection_dp_tube.fix(1)
    m.fs.ip_sh3.fcorrection_dp_shell.fix(1)

    m.fs.hp_sh3.di_tube.fix(0.042)
    m.fs.hp_sh3.thickness_tube.fix(0.004)
    m.fs.hp_sh3.pitch_x.fix(0.1)
    m.fs.hp_sh3.pitch_y.fix(0.1)
    m.fs.hp_sh3.length_tube_seg.fix(10)
    m.fs.hp_sh3.nseg_tube.fix(4)
    m.fs.hp_sh3.ncol_tube.fix(100)
    m.fs.hp_sh3.nrow_inlet.fix(1)
    m.fs.hp_sh3.delta_elevation.fix(0)
    m.fs.hp_sh3.therm_cond_wall = 43.0
    m.fs.hp_sh3.rfouling_tube = 0.0001
    m.fs.hp_sh3.rfouling_shell = 0.0001
    m.fs.hp_sh3.fcorrection_htc_tube.fix(1)
    m.fs.hp_sh3.fcorrection_htc_shell.fix(1)
    m.fs.hp_sh3.fcorrection_dp_tube.fix(1)
    m.fs.hp_sh3.fcorrection_dp_shell.fix(1)

    for i, s in m.fs.steam_turb.inlet_stage.items():
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
        s.ratioP[0] = 0.979  # use low pressure drop in inlet stage
        m.fs.steam_turb.throttle_valve[i].Cv.fix(0.002) #0.01
        m.fs.steam_turb.throttle_valve[i].valve_opening.fix(0.85)

    m.fs.steam_turb.inlet_mix.use_equal_pressure_constraint()
    m.fs.steam_turb.inlet_stage[1].flow_coeff.fix(0.001)
    m.fs.steam_turb.inlet_stage[1].efficiency_mech.fix(0.99)
    m.fs.steam_turb.inlet_stage[1].eff_nozzle.fix(0.93)
    m.fs.steam_turb.inlet_stage[1].blade_reaction.fix(0.9)

    for i, s in m.fs.steam_turb.hp_stages.items():
        s.ratioP[:] = 0.8238
        s.efficiency_isentropic[:] = 0.89
        s.efficiency_mech.fix(0.99)
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
    for i, s in m.fs.steam_turb.ip_stages.items():
        s.ratioP[:] = 0.79 #0.8264
        s.efficiency_isentropic[:] = 0.89
        s.efficiency_mech.fix(0.99)
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
    for i, s in m.fs.steam_turb.lp_stages.items():
        s.ratioP[:] = 0.79 #0.75
        s.efficiency_isentropic[:] = 0.89
        s.efficiency_mech.fix(0.99)
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)

    m.fs.steam_turb.outlet_stage.flow_coeff.fix(0.025)
    m.fs.steam_turb.outlet_stage.efficiency_mech.fix(0.99)
    m.fs.steam_turb.outlet_stage.eff_dry.fix(0.93)
    m.fs.steam_turb.outlet_stage.design_exhaust_flow_vol.fix(3000)

    m.fs.condenser.tube_inlet.flow_mol.fix(1.9e5)
    m.fs.condenser.tube_inlet.enth_mol.fix(1800)
    m.fs.condenser.tube_inlet.pressure.fix(5e5)
    m.fs.condenser.area.fix(7000)
    m.fs.condenser.overall_heat_transfer_coefficient.fix(3100)

    m.fs.hotwell.makeup_inlet.flow_mol[:].value = 12 #4118.6
    m.fs.hotwell.makeup_inlet.enth_mol.fix(1800)
    m.fs.hotwell.makeup_inlet.pressure.fix(101325)

    m.fs.cond_pump.efficiency_isentropic.fix(0.80)
    m.fs.cond_pump.deltaP.fix(4.0e5)
    
    m.fs.lp_splitter.split_fraction[:, "turb_outlet"].fix(0.22385)

    # CO2 post separator
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'H2O'].fix(0.005)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'CO2'].fix(0.97)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'N2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'O2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'co2_outlet', 'Ar'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'H2O'].fix(0.98)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'CO2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'N2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'O2'].fix(0)
    m.fs.co2_post_separator.split_fraction[:, 'h2o_outlet', 'Ar'].fix(0)

    # CO2 compression
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
    # H2 compression
    pin_h2 = 2.448e6
    pout_h2 = 6.48e6
    m.fs.h2_comp.ratioP.fix(pout_h2/pin_h2)
    m.fs.h2_comp.efficiency_isentropic.fix(0.85)


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
    iscale.set_scaling_factor(m.fs.rb_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.mp_knockout_condenser.tube.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.cw_knockout_condenser.tube.heat, 1e-7)
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


def scale_hrsg_steam_turbine_plant(m):
    iscale.set_scaling_factor(m.fs.lp_econ.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.lp_econ.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.lp_econ.tube.area, 10)
    iscale.set_scaling_factor(m.fs.lp_econ.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.lp_econ.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.lp_econ.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_econ.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_econ.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_econ.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.lp_drum.control_volume.volume, 1)
    iscale.set_scaling_factor(m.fs.lp_drum.control_volume.heat, 1e-3)
    #iscale.set_scaling_factor(m.fs.lp_drum.control_volume.energy_holdup, 1e-9)
    #iscale.set_scaling_factor(m.fs.lp_drum.control_volume.material_holdup, 1e-5)
    iscale.set_scaling_factor(m.fs.lp_drum.deltaP_gravity, 1e-3)
    iscale.set_scaling_factor(m.fs.lp_drum.deltaP_contraction, 1e-3)

    iscale.set_scaling_factor(m.fs.lp_evap.tube_cv.volume, 1e-1)
    iscale.set_scaling_factor(m.fs.lp_evap.shell_cv.volume, 1e-2)
    iscale.set_scaling_factor(m.fs.lp_evap.shell_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_evap.tube_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_evap.heat_flux_tube, 1e-3)
    iscale.set_scaling_factor(m.fs.lp_evap.heat_flux_shell, 1e-3)
    iscale.set_scaling_factor(m.fs.lp_evap.energy_holdup_metal, 1e-6)
    iscale.set_scaling_factor(m.fs.lp_evap.ratio_density, 1e-2)

    iscale.set_scaling_factor(m.fs.hp_econ1.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_econ1.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hp_econ1.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_econ1.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.hp_econ1.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_econ1.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ1.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ1.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ1.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.ip_econ1.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_econ1.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ1.tube.area, 10)
    iscale.set_scaling_factor(m.fs.ip_econ1.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.ip_econ1.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.ip_econ1.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_econ1.shell.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ1.tube.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ1.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_econ2.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_econ2.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hp_econ2.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_econ2.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.hp_econ2.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_econ2.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ2.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ2.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ2.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.ip_econ2.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_econ2.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ2.tube.area, 10)
    iscale.set_scaling_factor(m.fs.ip_econ2.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.ip_econ2.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.ip_econ2.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_econ2.shell.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ2.tube.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_econ2.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.ip_drum.control_volume.volume, 1)
    iscale.set_scaling_factor(m.fs.ip_drum.control_volume.heat, 1e-3)
    #iscale.set_scaling_factor(m.fs.ip_drum.control_volume.energy_holdup, 1e-9)
    #iscale.set_scaling_factor(m.fs.ip_drum.control_volume.material_holdup, 1e-5)
    iscale.set_scaling_factor(m.fs.ip_drum.deltaP_gravity, 1e-3)
    iscale.set_scaling_factor(m.fs.ip_drum.deltaP_contraction, 1e-3)
    
    iscale.set_scaling_factor(m.fs.ip_evap.tube_cv.volume, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_evap.shell_cv.volume, 1e-2)
    iscale.set_scaling_factor(m.fs.ip_evap.shell_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_evap.tube_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_evap.heat_flux_tube, 1e-3)
    iscale.set_scaling_factor(m.fs.ip_evap.heat_flux_shell, 1e-3)
    iscale.set_scaling_factor(m.fs.ip_evap.energy_holdup_metal, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_evap.ratio_density, 1e-1)

    iscale.set_scaling_factor(m.fs.lp_sh.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.lp_sh.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.lp_sh.tube.area, 10)
    iscale.set_scaling_factor(m.fs.lp_sh.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.lp_sh.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.lp_sh.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_sh.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_sh.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.lp_sh.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_econ3.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_econ3.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ3.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_econ3.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_econ3.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_econ3.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ3.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ3.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ3.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.ip_sh1.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_sh1.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh1.tube.area, 10)
    iscale.set_scaling_factor(m.fs.ip_sh1.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.ip_sh1.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.ip_sh1.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh1.shell.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh1.tube.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh1.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_econ4.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_econ4.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ4.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_econ4.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_econ4.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_econ4.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_econ4.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ4.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_econ4.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_drum.control_volume.volume, 1)
    iscale.set_scaling_factor(m.fs.hp_drum.control_volume.heat, 1e-3)
    #iscale.set_scaling_factor(m.fs.hp_drum.control_volume.energy_holdup, 1e-9)
    #iscale.set_scaling_factor(m.fs.hp_drum.control_volume.material_holdup, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_drum.deltaP_gravity, 1e-3)
    iscale.set_scaling_factor(m.fs.hp_drum.deltaP_contraction, 1e-3)

    iscale.set_scaling_factor(m.fs.hp_evap.tube_cv.volume, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_evap.shell_cv.volume, 1e-2)
    iscale.set_scaling_factor(m.fs.hp_evap.shell_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_evap.tube_cv.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_evap.heat_flux_tube, 1e-4)
    iscale.set_scaling_factor(m.fs.hp_evap.heat_flux_shell, 1e-4)
    iscale.set_scaling_factor(m.fs.hp_evap.energy_holdup_metal, 1e-6)
    iscale.set_scaling_factor(m.fs.hp_evap.ratio_density, 1)

    iscale.set_scaling_factor(m.fs.hp_sh1.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_sh1.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh1.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_sh1.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_sh1.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_sh1.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_sh1.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh1.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh1.heat_holdup, 1e-8)
    for t, c in m.fs.hp_sh1.v_tube_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3)
    for t, c in m.fs.hp_sh1.deltaP_tube_uturn_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3)

    iscale.set_scaling_factor(m.fs.ip_sh2.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_sh2.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh2.tube.area, 10)
    iscale.set_scaling_factor(m.fs.ip_sh2.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.ip_sh2.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.ip_sh2.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh2.shell.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh2.tube.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh2.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_sh2.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_sh2.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh2.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_sh2.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_sh2.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_sh2.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh2.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh2.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh2.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.ip_sh3.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.ip_sh3.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh3.tube.area, 10)
    iscale.set_scaling_factor(m.fs.ip_sh3.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.ip_sh3.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.ip_sh3.tube._enthalpy_flow, 1e-7)
    iscale.set_scaling_factor(m.fs.ip_sh3.shell.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh3.tube.enthalpy_flow_dx, 1e-6)
    iscale.set_scaling_factor(m.fs.ip_sh3.heat_holdup, 1e-8)

    iscale.set_scaling_factor(m.fs.hp_sh3.shell.area, 1e-1)
    iscale.set_scaling_factor(m.fs.hp_sh3.shell.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh3.tube.area, 10)
    iscale.set_scaling_factor(m.fs.hp_sh3.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.hp_sh3.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_sh3.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(m.fs.hp_sh3.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh3.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(m.fs.hp_sh3.heat_holdup, 1e-8)
    
    iscale.set_scaling_factor(m.fs.hp_pump.inlet.pressure, 1e-7)

    iscale.set_scaling_factor(m.fs.steam_turb.outlet_stage.control_volume.work, 1e-6)

    iscale.set_scaling_factor(m.fs.condenser.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.condenser.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.cond_pump.control_volume.work, 1e-6)
    
    iscale.set_scaling_factor(m.fs.ng_preheat1.area, 1e-3)
    iscale.set_scaling_factor(m.fs.ng_preheat1.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ng_preheat1.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ng_preheat1.overall_heat_transfer_coefficient, 1e-2)

    iscale.set_scaling_factor(
        m.fs.h2_comp.control_volume.properties_in[0].flow_mol, 1e-3)
    iscale.set_scaling_factor(
        m.fs.h2_comp.control_volume.properties_out[0].flow_mol, 1e-3)

def initialize_blue_hydrogen_plant(m):
    outlvl = 4 #idaeslog.INFO_LOW
    m.fs.ng_preheat2.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube_mix.fuel_inlet, m.fs.ng_preheat2.tube_outlet)
    m.fs.smr_tube_mix.initialize(outlvl = outlvl)
    m.fs.smr_shell_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.smr_tube.inlet, m.fs.smr_tube_mix.outlet)
    m.fs.smr_tube.outlet.temperature.fix(1152)
    m.fs.smr_tube.initialize(outlvl = outlvl)
    m.fs.smr_tube.outlet.temperature.unfix()
    copy_port_values(m.fs.smr_shell.inlet, m.fs.smr_shell_mix.outlet)
    m.fs.smr_shell.outlet.temperature.fix(1202)
    m.fs.smr_shell.initialize(outlvl = outlvl)
    m.fs.smr_shell.outlet.temperature.unfix()
    m.fs.hrsg_translator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.syngas_cooler.shell_inlet, m.fs.smr_tube.outlet)
    m.fs.syngas_cooler.initialize(outlvl = outlvl)
    copy_port_values(m.fs.wgsr.inlet, m.fs.syngas_cooler.shell_outlet)
    m.fs.wgsr.initialize(outlvl = outlvl)
    copy_port_values(m.fs.wgsr_cooler.shell_inlet, m.fs.wgsr.outlet)
    m.fs.wgsr_cooler.initialize(outlvl = outlvl)
    m.fs.rb_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.rb_splitter.inlet, m.fs.rb_pump.outlet)
    split_save = m.fs.rb_splitter.split_fraction[0, "cond_outlet"].value
    m.fs.rb_splitter.split_fraction[:, "cond_outlet"].fix(0.4)
    m.fs.rb_splitter.initialize(outlvl = outlvl)
    m.fs.rb_splitter.split_fraction[:, "cond_outlet"].fix(split_save)
    copy_port_values(m.fs.rb_knockout_condenser.tube_inlet, m.fs.rb_splitter.cond_outlet)
    copy_port_values(m.fs.rb_knockout_condenser.shell_inlet, m.fs.wgsr_cooler.shell_outlet)
    m.fs.rb_knockout_condenser.initialize(outlvl = outlvl)    
    m.fs.mp_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_knockout_condenser.tube_inlet, m.fs.mp_pump.outlet)
    copy_port_values(m.fs.mp_knockout_condenser.shell_inlet, m.fs.rb_knockout_condenser.shell_outlet)
    m.fs.mp_knockout_condenser.initialize(outlvl = outlvl)
    # inlet condition connected to steam turbine
    copy_port_values(m.fs.lp_knockout_condenser.shell_inlet, m.fs.mp_knockout_condenser.shell_outlet)
    m.fs.lp_knockout_condenser.initialize(outlvl = outlvl)
    copy_port_values(m.fs.cw_knockout_condenser.shell_inlet, m.fs.lp_knockout_condenser.shell_outlet)
    m.fs.cw_knockout_condenser.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_pre_separator.inlet, m.fs.cw_knockout_condenser.shell_outlet)
    m.fs.co2_pre_separator.initialize(outlvl = outlvl)    
    copy_port_values(m.fs.h2_separator.inlet, m.fs.co2_pre_separator.syngas_outlet)
    m.fs.h2_separator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_mix.mp_inlet, m.fs.mp_knockout_condenser.tube_outlet)
    copy_port_values(m.fs.mp_mix.mpkw_inlet, m.fs.mp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.lpkw_inlet, m.fs.lp_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.cwkw_inlet, m.fs.cw_knockout_condenser.shell_water_outlet)
    copy_port_values(m.fs.mp_mix.rbkw_inlet, m.fs.rb_knockout_condenser.shell_water_outlet)
    m.fs.mp_mix.initialize(outlvl = outlvl)
    copy_port_values(m.fs.mp_booster_pump.inlet, m.fs.mp_mix.outlet)
    m.fs.mp_booster_pump.initialize(outlvl = outlvl)
    copy_port_values(m.fs.steam_translator.inlet, m.fs.syngas_cooler.tube_outlet)
    m.fs.steam_translator.outlet.mole_frac_comp[:,'H2O'].value = 0.99999
    m.fs.steam_translator.initialize(outlvl = outlvl)


def initialize_gas_turbine_plant(m):
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
    m.fs.gts1.ratioP[0] = 0.39
    m.fs.gts1.efficiency_isentropic[0] = 0.9
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


def initialize_hrsg_steam_turbine_plant(m):
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 200,
            'bound_push': 1e-16
    }
    outlvl = 4
    copy_port_values(m.fs.hp_sh3.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_sh3.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_sh2.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_sh2.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_sh1.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_evap.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_econ4.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_sh1.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_econ3.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.lp_sh.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_evap.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_econ2.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_econ2.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ip_econ1.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.hp_econ1.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.lp_evap.shell_inlet, m.fs.lp_econ.shell_inlet)
    copy_port_values(m.fs.ng_preheat1.shell_inlet, m.fs.lp_econ.shell_inlet)

    # initialize lp_econ
    copy_port_values(m.fs.lp_econ.tube_inlet, m.fs.lp_knockout_condenser.tube_outlet)
    m.fs.lp_econ.shell_inlet.temperature[:].value = 438.3
    m.fs.lp_econ.shell_inlet.pressure[:].value = 100265
    m.fs.lp_econ.initialize(outlvl = outlvl)

    # initialize ng_preheat1
    copy_port_values(m.fs.ng_preheat1.shell_inlet, m.fs.lp_econ.shell_outlet)
    m.fs.ng_preheat1.initialize(outlvl = outlvl)

    # initialize econ_mixer
    copy_port_values(m.fs.rb_econ_mixer.econ_inlet, m.fs.lp_econ.tube_outlet)
    copy_port_values(m.fs.rb_econ_mixer.reboiler_inlet, m.fs.rb_splitter.mix_outlet)
    m.fs.rb_econ_mixer.initialize(outlvl = outlvl)

    # initialize lp_ip_hp_split
    copy_port_values(m.fs.lp_ip_hp_split.inlet, m.fs.rb_econ_mixer.outlet)
    m.fs.lp_ip_hp_split.initialize(outlvl = outlvl)

    # initialize lp_evap water circuit
    copy_port_values(m.fs.lp_drum.feedwater_inlet, m.fs.lp_ip_hp_split.lp_outlet)
    # Estimate water/steam mixture from evaporator
    m.fs.lp_drum.water_steam_inlet.flow_mol[:].value = 34627
    m.fs.lp_drum.water_steam_inlet.pressure[:].value = pyo.value(m.fs.lp_ip_hp_split.lp_outlet.pressure[0])
    pres_lp_drum = pyo.value(m.fs.lp_ip_hp_split.lp_outlet.pressure[0])
    enth_mol_lp_drum = iapws95.htpx(P=pres_lp_drum*pyo.units.Pa, x=0.023663)
    m.fs.lp_drum.water_steam_inlet.enth_mol[:].value = enth_mol_lp_drum
    m.fs.lp_drum.initialize(outlvl = outlvl)
    copy_port_values(m.fs.lp_downcomer.inlet, m.fs.lp_drum.liquid_outlet)
    m.fs.lp_downcomer.initialize(outlvl = outlvl)
    copy_port_values(m.fs.lp_evap.tube_inlet, m.fs.lp_downcomer.outlet)
    m.fs.lp_evap.shell_inlet.temperature[:].value = 498.03
    m.fs.lp_evap.shell_inlet.pressure[:].value = 100790
    m.fs.lp_evap.initialize(outlvl = outlvl)

    # initialize lp_drum_mixer
    copy_port_values(m.fs.lp_drum_mixer.evap_inlet, m.fs.lp_evap.tube_outlet)
    copy_port_values(m.fs.lp_drum_mixer.reboiler_inlet, m.fs.rb_knockout_condenser.tube_outlet)
    m.fs.lp_drum_mixer.initialize(outlvl = outlvl)

    # initialize lp_da_split
    copy_port_values(m.fs.lp_da_split.inlet, m.fs.lp_drum.steam_outlet)
    m.fs.lp_da_split.initialize(outlvl = outlvl)

    # initialize ip_pump
    copy_port_values(m.fs.ip_pump.inlet, m.fs.lp_ip_hp_split.ip_outlet)
    m.fs.ip_pump.initialize(outlvl = outlvl)

    # initialize ip_pump
    copy_port_values(m.fs.hp_pump.inlet, m.fs.lp_ip_hp_split.hp_outlet)
    m.fs.hp_pump.initialize(outlvl = outlvl)

    # initialize hp_econ1
    copy_port_values(m.fs.hp_econ1.tube_inlet, m.fs.hp_pump.outlet)
    m.fs.hp_econ1.shell_inlet.temperature[:].value = 522.75
    m.fs.hp_econ1.shell_inlet.pressure[:].value = 100924
    m.fs.hp_econ1.initialize(outlvl = outlvl)

    # initialize ip_econ1
    copy_port_values(m.fs.ip_econ1.tube_inlet, m.fs.ip_pump.outlet)
    m.fs.ip_econ1.shell_inlet.temperature[:].value = 533.57
    m.fs.ip_econ1.shell_inlet.pressure[:].value = 100971
    m.fs.ip_econ1.initialize(outlvl = outlvl)

    # initialize hp_econ2
    copy_port_values(m.fs.hp_econ2.tube_inlet, m.fs.hp_econ1.tube_outlet)
    m.fs.hp_econ2.shell_inlet.temperature[:].value = 566.9
    m.fs.hp_econ2.shell_inlet.pressure[:].value = 101164
    m.fs.hp_econ2.initialize(outlvl = outlvl)

    # initialize ip_econ2
    copy_port_values(m.fs.ip_econ2.tube_inlet, m.fs.ip_econ1.tube_outlet)
    m.fs.ip_econ2.shell_inlet.temperature[:].value = 573.77
    m.fs.ip_econ2.shell_inlet.pressure[:].value = 101205
    m.fs.ip_econ2.initialize(outlvl = outlvl)

    # initialize ip_evap water circuit
    copy_port_values(m.fs.ip_drum.feedwater_inlet, m.fs.ip_econ2.tube_outlet)
    # Estimate water/steam mixture from evaporator
    m.fs.ip_drum.water_steam_inlet.flow_mol[:].value = 21954
    m.fs.ip_drum.water_steam_inlet.pressure[:].value = pyo.value(m.fs.ip_econ2.tube_outlet.pressure[0])
    pres_ip_drum = pyo.value(m.fs.ip_econ2.tube_outlet.pressure[0])
    enth_mol_ip_drum = iapws95.htpx(P=pres_ip_drum*pyo.units.Pa, x=0.04843)
    m.fs.ip_drum.water_steam_inlet.enth_mol[:].value = enth_mol_ip_drum
    m.fs.ip_drum.initialize(outlvl = outlvl)
    copy_port_values(m.fs.ip_downcomer.inlet, m.fs.ip_drum.liquid_outlet)
    m.fs.ip_downcomer.initialize(outlvl = outlvl)
    copy_port_values(m.fs.ip_evap.tube_inlet, m.fs.ip_downcomer.outlet)
    m.fs.ip_evap.shell_inlet.temperature[:].value = 633.2
    m.fs.ip_evap.shell_inlet.pressure[:].value = 101659
    m.fs.ip_evap.initialize(outlvl = outlvl)

    # initialize lp_sh
    copy_port_values(m.fs.lp_sh.tube_inlet, m.fs.lp_da_split.lp_outlet)
    m.fs.lp_sh.shell_inlet.temperature[:].value = 645.7
    m.fs.lp_sh.shell_inlet.pressure[:].value = 101890
    m.fs.lp_sh.initialize(outlvl = outlvl)

    # initialize hp_econ3
    copy_port_values(m.fs.hp_econ3.tube_inlet, m.fs.hp_econ2.tube_outlet)
    m.fs.hp_econ3.shell_inlet.temperature[:].value = 706.35
    m.fs.hp_econ3.shell_inlet.pressure[:].value = 102129
    m.fs.hp_econ3.initialize(outlvl = outlvl)

    # initialize ip_sh1
    copy_port_values(m.fs.ip_sh1.tube_inlet, m.fs.ip_drum.steam_outlet)
    m.fs.ip_sh1.shell_inlet.temperature[:].value = 713
    m.fs.ip_sh1.shell_inlet.pressure[:].value = 102154.6
    m.fs.ip_sh1.initialize(outlvl = outlvl)

    # initialize ip_mixer
    copy_port_values(m.fs.ip_mixer.hrsg_inlet, m.fs.ip_sh1.tube_outlet)
    m.fs.ip_mixer.initialize(outlvl = outlvl)

    # initialize hp_econ4
    copy_port_values(m.fs.hp_econ4.tube_inlet, m.fs.hp_econ3.tube_outlet)
    m.fs.hp_econ4.shell_inlet.temperature[:].value = 742.6
    m.fs.hp_econ4.shell_inlet.pressure[:].value = 102285
    m.fs.hp_econ4.initialize(outlvl = outlvl)

    # initialize hp_evap water circuit
    copy_port_values(m.fs.hp_drum.feedwater_inlet, m.fs.hp_econ4.tube_outlet)
    # Estimate water/steam mixture from evaporator
    m.fs.hp_drum.water_steam_inlet.flow_mol[:].value = 71383.5
    m.fs.hp_drum.water_steam_inlet.pressure[:].value = pyo.value(m.fs.hp_econ4.tube_outlet.pressure[0])
    pres_hp_drum = pyo.value(m.fs.hp_econ4.tube_outlet.pressure[0])
    enth_mol_hp_drum = iapws95.htpx(P=pres_hp_drum*pyo.units.Pa, x=0.06971)
    m.fs.hp_drum.water_steam_inlet.enth_mol[:].value = enth_mol_hp_drum
    m.fs.hp_drum.initialize(outlvl = outlvl)
    copy_port_values(m.fs.hp_downcomer.inlet, m.fs.hp_drum.liquid_outlet)
    m.fs.hp_downcomer.initialize(outlvl = outlvl)
    copy_port_values(m.fs.hp_evap.tube_inlet, m.fs.hp_downcomer.outlet)
    m.fs.hp_evap.shell_inlet.temperature[:].value = 893
    m.fs.hp_evap.shell_inlet.pressure[:].value = 102859
    m.fs.hp_evap.initialize(outlvl = outlvl)

    # initialize hp_sh1
    copy_port_values(m.fs.hp_sh1.tube_inlet, m.fs.hp_drum.steam_outlet)
    m.fs.hp_sh1.shell_inlet.temperature[:].value = 964
    m.fs.hp_sh1.shell_inlet.pressure[:].value = 103060
    m.fs.hp_sh1.initialize(outlvl = outlvl)

    # initialize ip_sh2
    copy_port_values(m.fs.ip_sh2.tube_inlet, m.fs.ip_mixer.outlet)
    m.fs.ip_sh2.shell_inlet.temperature[:].value = 1015.9
    m.fs.ip_sh2.shell_inlet.pressure[:].value = 103184
    m.fs.ip_sh2.initialize(outlvl = outlvl)

    # initialize hp_sh2
    copy_port_values(m.fs.hp_sh2.tube_inlet, m.fs.hp_sh1.tube_outlet)
    m.fs.hp_sh2.shell_inlet.temperature[:].value = 1051.1
    m.fs.hp_sh2.shell_inlet.pressure[:].value = 103280
    m.fs.hp_sh2.initialize(outlvl = outlvl)

    # initialize ip_sh3
    copy_port_values(m.fs.ip_sh3.tube_inlet, m.fs.ip_sh2.tube_outlet)
    m.fs.ip_sh3.shell_inlet.temperature[:].value = 1092.3
    m.fs.ip_sh3.shell_inlet.pressure[:].value = 103397
    m.fs.ip_sh3.initialize(outlvl = outlvl)

    # initialize hp_sh3
    copy_port_values(m.fs.hp_sh3.tube_inlet, m.fs.hp_sh2.tube_outlet)
    m.fs.hp_sh3.shell_inlet.temperature[:].value = 1128.6
    m.fs.hp_sh3.shell_inlet.pressure[:].value = 103500
    m.fs.hp_sh3.initialize(outlvl = outlvl)

    # initialize lp_mixer
    copy_port_values(m.fs.lp_mixer.hrsg_inlet, m.fs.lp_sh.tube_outlet)
    m.fs.lp_mixer.turb_inlet.flow_mol[:].value = 6039.5
    m.fs.lp_mixer.turb_inlet.pressure[:].value = 364725
    enth_mol_ip_outlet = iapws95.htpx(P=364725*pyo.units.Pa, T=524.36*pyo.units.K,)
    m.fs.lp_mixer.turb_inlet.enth_mol[:].value = enth_mol_ip_outlet
    m.fs.lp_mixer.initialize(outlvl = outlvl)

    # initialize lp_splitter
    copy_port_values(m.fs.lp_splitter.inlet, m.fs.lp_mixer.outlet)
    m.fs.lp_splitter.initialize(outlvl = outlvl)

    # initialize steam_turb
    copy_port_values(m.fs.steam_turb.inlet_split.inlet, m.fs.hp_sh3.tube_outlet)
    copy_port_values(m.fs.steam_turb.ip_stages[1].inlet, m.fs.ip_sh3.tube_outlet)
    copy_port_values(m.fs.steam_turb.lp_stages[1].inlet, m.fs.lp_splitter.turb_outlet)
    m.fs.steam_turb.inlet_split.inlet.flow_mol[:].value = 4976.2
    m.fs.steam_turb.inlet_split.inlet.enth_mol[:].value = 62560.8
    m.fs.steam_turb.inlet_split.inlet.pressure[:].value = 18172562
    # need to specify the back pressure for initialization
    m.fs.steam_turb.outlet_stage.control_volume.properties_out[:].pressure.value = 5226.7
    m.fs.steam_turb.outlet_stage.control_volume.properties_out[:].pressure.fix()
    m.fs.steam_turb.initialize(calculate_inlet_cf=True,
                               calculate_outlet_cf=True,
                               copy_disconneted_flow=False,
                               copy_disconneted_pressure=False,
                               flow_iterate=1,
                               outlvl = outlvl)
    m.fs.steam_turb.outlet_stage.control_volume.properties_out[:].pressure.unfix()
    m.fs.steam_turb.outlet_stage.flow_coeff.fix(0.025)
    m.fs.steam_turb.inlet_stage[1].flow_coeff.fix(0.001)

    # initialize condenser
    copy_port_values(m.fs.condenser.shell_inlet, m.fs.steam_turb.outlet_stage.outlet)
    m.fs.condenser.initialize(unfix='pressure', outlvl = outlvl)

    # initialize hotwell mixer
    copy_port_values(m.fs.hotwell.condensate_inlet, m.fs.condenser.shell_outlet)
    m.fs.hotwell.initialize(outlvl = outlvl)

    # initialize cond_pump
    copy_port_values(m.fs.cond_pump.inlet, m.fs.hotwell.outlet)
    m.fs.cond_pump.initialize(outlvl = outlvl)

    copy_port_values(m.fs.co2_post_separator.inlet, m.fs.ng_preheat1.shell_outlet)
    m.fs.co2_post_separator.initialize(outlvl = outlvl)

    # co2 mixer and compression
    copy_port_values(m.fs.co2_translator_syn.inlet, m.fs.co2_pre_separator.co2_outlet)
    m.fs.co2_translator_syn.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_translator_flue.inlet, m.fs.co2_post_separator.co2_outlet)
    m.fs.co2_translator_flue.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_mix.pre_inlet, m.fs.co2_translator_syn.outlet)
    copy_port_values(m.fs.co2_mix.post_inlet, m.fs.co2_translator_flue.outlet)
    m.fs.co2_mix.initialize(outlvl = outlvl)   
    copy_port_values(m.fs.co2_comp1.inlet, m.fs.co2_mix.outlet)
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
    copy_port_values(m.fs.h2_comp.inlet, m.fs.h2_separator.h2_outlet)
    m.fs.h2_comp.inlet.fix()
    solver.solve(m.fs.h2_comp, tee=True)
    #m.fs.h2_comp.initialize(outlvl = outlvl)
    m.fs.h2_comp.inlet.unfix()



def get_full_plant_model():
    logging.getLogger('idaes.core.util.scaling').setLevel(logging.ERROR)
    # create model and flowsheet
    m = pyo.ConcreteModel(name='NGCC plus SMR IES')
    m.fs = FlowsheetBlock(dynamic=False)
    # solver and options
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 200,
            'bound_push': 1e-16
    }
    init_fname = 'ies_ngcc_baseline.json'
    if os.path.exists(init_fname):
        build_gas_turbine_plant(m)
        build_blue_hydrogen_plant(m)
        build_hrsg_steam_turbine_plant(m)
        pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
        set_gas_turbine_plant_inputs(m)
        set_blue_hydrogen_plant_inputs(m)
        set_hrsg_steam_turbine_plant_inputs(m)
        scale_gas_turbine_plant(m)
        scale_blue_hydrogen_plant(m)
        scale_hrsg_steam_turbine_plant(m)
        iscale.calculate_scaling_factors(m)
        ms.from_json(m, fname=init_fname, wts=ms.StoreSpec(suffix=False))
        solver.solve(m, tee=True)
        return m
    build_gas_turbine_plant(m)
    build_blue_hydrogen_plant(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    set_gas_turbine_plant_inputs(m)
    set_blue_hydrogen_plant_inputs(m)
    scale_gas_turbine_plant(m)
    scale_blue_hydrogen_plant(m)
    iscale.calculate_scaling_factors(m)
    initialize_gas_turbine_plant(m)
    initialize_blue_hydrogen_plant(m)
    m.fs.comp_igv.Cv.fix()
    m.fs.comp_igv.deltaP.unfix()
    m.fs.comp_igv.valve_opening.unfix()
    m.fs.comp_igv.inlet.flow_mol.unfix()
    m.fs.lp_knockout_condenser.tube_inlet.fix()
    m.fs.steam_translator.inlet.flow_mol.fix(2730.6)
    m.fs.rb_pump.inlet.flow_mol.fix()
    m.fs.rb_pump.inlet.enth_mol.fix()
    m.fs.rb_pump.inlet.pressure.fix()
    m.fs.mp_pump.inlet.flow_mol.unfix()
    m.fs.smr_shell_mix.fuel_inlet.flow_mol.unfix()
    solver.solve(m, tee=True)
    # add HRSG and steam cycle to the plant
    build_hrsg_steam_turbine_plant(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    set_hrsg_steam_turbine_plant_inputs(m)
    scale_hrsg_steam_turbine_plant(m)
    iscale.calculate_scaling_factors(m)
    initialize_hrsg_steam_turbine_plant(m)
    m.fs.lp_ip_hp_split.split_fraction[:, "lp_outlet"].unfix()
    m.fs.lp_ip_hp_split.split_fraction[:, "ip_outlet"].unfix()
    m.fs.lp_knockout_condenser.tube_inlet.unfix()
    m.fs.hp_pump.deltaP.unfix()
    m.fs.ip_pump.deltaP.unfix()
    m.fs.cond_pump.deltaP.unfix()
    m.fs.ng_preheat2.tube_inlet.unfix()
    m.fs.rb_pump.inlet.flow_mol.unfix()
    m.fs.rb_pump.inlet.pressure.unfix()
    m.fs.rb_pump.inlet.enth_mol.unfix()
    m.fs.p_main_steam.fix(17277250)
    m.fs.steam_turb.throttle_valve[1].valve_opening.unfix()
    m.fs.reboiler_duty_constraint.deactivate()
    solver.solve(m, tee=True)
    m.fs.reboiler_duty_constraint.activate()
    m.fs.lp_splitter.split_fraction[:, "turb_outlet"].unfix()
    solver.solve(m, tee=True)
    return m

if __name__ == "__main__":
    m = get_full_plant_model()
