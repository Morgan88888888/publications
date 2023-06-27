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
Standard IDAES Gibbs reactor model.
"""
# Import Pyomo libraries
from pyomo.environ import Constraint, Param, Reals, Reference, Var, log, exp, value
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import constants
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Jinliang Ma, Andrew Lee"


@declare_process_block_class("GibbsReactor")
class GibbsReactorData(UnitModelBlockData):
    """
    Standard Gibbs Reactor Unit Model Class

    This model assume all possible reactions reach equilibrium such that the
    system partial molar Gibbs free energy is minimized.
    Since some species mole flow rate might be very small,
    the natural log of the species molar flow rate is used.
    Instead of specifying the system Gibbs free energy as an objective
    function, the equations for zero partial derivatives of the grand function
    with Lagrangian multiple terms with repect to product species mole flow
    rates and the multiples are specified as constraints.
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Gibbs reactors do not support dynamic models, thus this must be
False."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag",
        doc="""Gibbs reactors do not have defined volume, thus this must be
False."""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.none,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of material balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))


    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(GibbsReactorData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(**{
                "dynamic": self.config.dynamic,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args})

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Set references to balance terms at unit level
        if (self.config.has_heat_transfer is True and
                self.config.energy_balance_type != EnergyBalanceType.none):
            self.heat_duty = Reference(self.control_volume.heat[:])
        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != MomentumBalanceType.none):
            self.deltaP = Reference(self.control_volume.deltaP[:])

        # Add performance equations
        # Add Lagrangian multiplier variables
        e_units = self.config.property_package.get_metadata(
            ).get_derived_units("energy_mole")
        self.lagrange_mult = Var(self.flowsheet().config.time,
                                 self.config.property_package.element_list,
                                 domain=Reals,
                                 initialize=1e5,
                                 doc="Lagrangian multipliers",
                                 units=e_units)

        # Use Lagrangian multiple method to derive equations for Out_Fi
        # Use RT*lagrange as the Lagrangian multiple such that lagrange is in
        # a similar order of magnitude as log(Yi)

        self.element_mole = Var(self.flowsheet().config.time,
                                 self.config.property_package.element_list,
                                 domain=Reals,
                                 initialize=2,
                                 doc="Specific mole of element")

        self.ln_specific_mole = Var(self.flowsheet().config.time,
                                 self.config.property_package.component_list,
                                 domain=Reals,
                                 initialize=0,
                                 doc="specific mole of species")

        self.sum_specific_mole = Var(self.flowsheet().config.time, initialize=35, doc="sum_of_specific_moles")


        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.element_list,
                         doc="specific mole of element")
        def element_mole_in_eqn(b, t, j):
            return b.element_mole[t,j]*b.control_volume.properties_in[t].mw == \
            sum(b.config.property_package.get_component(c).config.elemental_composition[j]*\
            b.control_volume.properties_in[t].mole_frac_comp[c] for c in b.config.property_package.component_list \
            if j in b.config.property_package.get_component(c).config.elemental_composition)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.element_list,
                         doc="specific mole of element at outlet")
        def element_mole_out_eqn(b, t, j):
            return b.element_mole[t,j] == \
            sum(b.config.property_package.get_component(c).config.elemental_composition[j]*\
            exp(b.ln_specific_mole[t,c]) for c in b.config.property_package.component_list \
            if j in b.config.property_package.get_component(c).config.elemental_composition)

        @self.Constraint(self.flowsheet().config.time,
                         doc="sum_of_specific_mole")
        def sum_specific_mole_eqn(b, t):
            return b.sum_specific_mole[t] == sum(exp(b.ln_specific_mole[t,c]) for c in b.config.property_package.component_list)

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="component molar flow")
        def flow_mol_comp_eqn(b, t, j):
            return b.control_volume.properties_out[t].flow_mol_comp[j] == exp(b.ln_specific_mole[t,j])*b.control_volume.properties_in[t].flow_mass

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.component_list,
                         doc="Gibbs energy minimisation constraint")
        def gibbs_minimization(b, tim, j):
            temp = b.control_volume.properties_out[tim].temperature
            t = temp/1000
            p = b.control_volume.properties_out[tim].pressure
            pk = b.config.property_package
            r_gas = constants.Constants.gas_constant
            h = 1000*(pk.get_component(j).cp_mol_ig_comp_coeff_A*t +
                    pk.get_component(j).cp_mol_ig_comp_coeff_B*t**2 / 2 +
                    pk.get_component(j).cp_mol_ig_comp_coeff_C*t**3 / 3 +
                    pk.get_component(j).cp_mol_ig_comp_coeff_D*t**4 / 4 -
                    pk.get_component(j).cp_mol_ig_comp_coeff_E/t +
                    pk.get_component(j).cp_mol_ig_comp_coeff_F)
            s = (pk.get_component(j).cp_mol_ig_comp_coeff_A*log(t) +
                    pk.get_component(j).cp_mol_ig_comp_coeff_B*t +
                    pk.get_component(j).cp_mol_ig_comp_coeff_C*t**2 / 2 +
                    pk.get_component(j).cp_mol_ig_comp_coeff_D*t**3 / 3 -
                    pk.get_component(j).cp_mol_ig_comp_coeff_E/t**2 / 2 +
                    pk.get_component(j).cp_mol_ig_comp_coeff_G)
            return 0 == 1e-4*(h-temp*s + r_gas*temp*(log(p) - log(101325) +
                    b.ln_specific_mole[tim,j]-log(b.sum_specific_mole[tim])) + 
                    sum(b.lagrange_mult[tim, e] *
                    pk.get_component(j).config.elemental_composition[e]
                    for e in pk.get_component(j).config.elemental_composition))


    def initialize_build(blk, state_args=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
        Gibbs reactor initialization routine.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) for the control_volume to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        '''
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)
        opt.options = {
                  "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 15,
            }
        flags = blk.control_volume.initialize(outlvl=outlvl+1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args)
        init_log.info_high("Initialization Step 1 Complete.")
        # Fix elemental composition based on inlet and guess outlet specific moles
        for t in blk.flowsheet().config.time:
            for j in blk.config.property_package.element_list:
                blk.element_mole[t, j].fix(
                value(sum(blk.config.property_package.get_component(c).config.elemental_composition[j]*\
                blk.control_volume.properties_in[0].mole_frac_comp[c] for c in blk.config.property_package.component_list \
                if j in blk.config.property_package.get_component(c).config.elemental_composition)/blk.control_volume.properties_in[0].mw))
                if blk.element_mole[t, j].value <= 0:
                    init_log.warning("Element {} is not in the inlet stream.",j)
            # Find if the reactor is fuel rich or fuel lean
            cov_fuel = 0
            cov_oxidizer = 0
            for j in blk.config.property_package.element_list:
                if j=='C':
                    cov_fuel += blk.element_mole[0,j].value*4
                if j=='H':
                    cov_fuel += blk.element_mole[0,j].value
                if j=='O':
                    cov_oxidizer += blk.element_mole[0,j].value*2
                if j=='S':
                    cov_fuel += blk.element_mole[0,j].value*4
            cov_diff = cov_fuel - cov_oxidizer
            for i in blk.config.property_package.component_list:
                cov_f = 0
                cov_o = 0
                ele_comp = blk.config.property_package.get_component(i).config.elemental_composition
                for j in ele_comp:
                    if j=='C':
                        cov_f += ele_comp[j]*4
                    if j=='H':
                        cov_f += ele_comp[j]
                    if j=='O':
                        cov_o += ele_comp[j]*2
                    if j=='S':
                        cov_f += ele_comp[j]*4
                cov_d = cov_f - cov_o
                if cov_d==0: # neither fuel nor oxidizer species
                    if len(ele_comp)>1:   # typical combustion product such as CO2, SO2, H2O, and H2S
                        blk.ln_specific_mole[t,i].value = 1 #0
                    else:   # typical inert single element species such as Ar and N2
                        for e in ele_comp:
                            blk.ln_specific_mole[t,i].value = log(blk.element_mole[0,e].value/ele_comp[e])
                elif cov_d > 0: # fuel species 
                    if cov_diff > 0:   # fuel rich
                        blk.ln_specific_mole[t,i].value = 2
                        if 'C' in ele_comp and 'H' in ele_comp:   # hydrocarbons such as CH4 and C2H6 are usually not stable at high temperature
                            blk.ln_specific_mole[t,i].value = -5
                    else:   # fuel lean
                        blk.ln_specific_mole[t,i].value = -5
                else: # oxidizer species
                    if cov_diff > 0:   # fuel rich
                        blk.ln_specific_mole[t,i].value = -10
                    else:   # fuel lean
                        blk.ln_specific_mole[t,i].value = 5
                #print ('ln specific mole =', i, blk.ln_specific_mole[t,i].value)
        blk.element_mole_in_eqn.deactivate()
        tout_fixed = blk.control_volume.properties_out[0].temperature.is_fixed()
        blk.heat_duty.unfix()
        if (not tout_fixed):
            print('unfixed temp')
            blk.control_volume.properties_out[:].temperature.fix(1000)
            heat_duty_save = blk.heat_duty[0].value
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 Complete.")            
        blk.element_mole.unfix()
        blk.element_mole_in_eqn.activate()
        if (not tout_fixed):
            blk.heat_duty.fix(heat_duty_save)
            blk.control_volume.properties_out[:].temperature.unfix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 Complete.")

        blk.control_volume.release_state(flags, outlvl+1)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for v in self.lagrange_mult.values():
            if iscale.get_scaling_factor(v, warning=True) is None:
                iscale.set_scaling_factor(v, 1e-5)

