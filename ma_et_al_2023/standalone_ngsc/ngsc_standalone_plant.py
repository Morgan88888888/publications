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
This is a flowsheet for standalone NGSC plant
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

from idaes.models.properties import swco2
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock)

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

from gibbs_reactor_ideal import GibbsReactor
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop
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




def build_gas_turbine_plant(m):
    syn_config = get_prop(
        components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
    m.fs.syn_props = GenericParameterBlock(**syn_config)
    m.fs.co2_props = swco2.SWCO2ParameterBlock()
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
        
    # postcombustion co2 capture system
    m.fs.co2_post_separator = Separator(
        **{"outlet_list": ["co2_outlet", "fluegas_outlet", "h2o_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.syn_props})

    # property package translater from water_props to syn_props
    m.fs.co2_translator = Translator(
        **{"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.co2_props})

    # additional variables
    # constraints for co2_translator from syn_props to co2_props
    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_T(b, t):
        return b.properties_out[t].temperature == 303

    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_P(b, t):
        return b.outlet.pressure[t] == 1.4e5

    @m.fs.co2_translator.Constraint(m.fs.time)
    def co2_translator_flow(b, t):
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

    m.fs.duty_per_mol_post = pyo.Param(initialize=1217.184*1054/0.453592*0.044,
        units=pyo.units.J/pyo.units.mol,
        doc='required reboiler duty per mol CO2 captured for post-combustion capture')

    @m.fs.Expression(m.fs.config.time)
    def post_capture_reboiler_duty(b, t):
        return (b.duty_per_mol_post*b.co2_post_separator.co2_outlet_state[t].flow_mol_comp['CO2'])

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

    m.fs.air_gas_mix3_to_co2_post_separator = Arc(
        source=m.fs.air_gas_mix3.outlet,
        destination=m.fs.co2_post_separator.inlet)

    m.fs.co2_post_separator_to_co2_translator = Arc(
        source=m.fs.co2_post_separator.co2_outlet,
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


def set_gas_turbine_plant_inputs(m):
    # set flow scale of performance curves
    m.fs.performance_flow_scale.fix(2.0/115*250)
    # air to gas turbine
    m.fs.comp_igv.inlet.flow_mol.fix(16720.0*115/250)
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

    m.fs.gt_mix.fuel_inlet.flow_mol.fix(790.0*115/250)
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

def scale_gas_turbine_plant(m):
    # set gas turbine flow scaling
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_in[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_in[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts3.control_volume.properties_in[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_out[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_out[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts3.control_volume.properties_out[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts1.properties_isentropic[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts2.properties_isentropic[0].flow_mol, 1e-4)
    iscale.set_scaling_factor(
        m.fs.gts3.properties_isentropic[0].flow_mol, 1e-4)

    
    for c in m.fs.gt_mix.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-9)
    for c in m.fs.air_gas_mix1.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)
    for c in m.fs.air_gas_mix2.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)
    for c in m.fs.air_gas_mix3.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-7)



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
    m.fs.gts1.ratioP[0] = 0.7
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
    # co2 separation and compression
    copy_port_values(m.fs.co2_post_separator.inlet, m.fs.air_gas_mix3.outlet)
    m.fs.co2_post_separator.initialize(outlvl = outlvl)
    copy_port_values(m.fs.co2_translator.inlet, m.fs.co2_post_separator.co2_outlet)
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



def get_ngsc_plant_model():
    logging.getLogger('idaes.core.util.scaling').setLevel(logging.ERROR)
    # create model and flowsheet
    m = pyo.ConcreteModel(name='NGSC plant')
    m.fs = FlowsheetBlock(dynamic=False)
    # solver and options
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 100,
            'bound_push': 1e-16
    }
    build_gas_turbine_plant(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    set_gas_turbine_plant_inputs(m)
    scale_gas_turbine_plant(m)
    iscale.calculate_scaling_factors(m)
    init_fname = 'ngsc_baseline.json'
    if os.path.exists(init_fname):

        ms.from_json(m, fname=init_fname, wts=ms.StoreSpec(suffix=False))
        solver.solve(m, tee=True)
    else:
        initialize_gas_turbine_plant(m)
        m.fs.comp_igv.Cv.fix()
        m.fs.comp_igv.deltaP.unfix()
        m.fs.comp_igv.valve_opening.unfix()
        m.fs.comp_igv.inlet.flow_mol.unfix()
        solver.solve(m, tee=True)
        ms.to_json(m, fname='ngsc_baseline.json')
    return m


if __name__ == "__main__":
    m = get_ngsc_plant_model()
