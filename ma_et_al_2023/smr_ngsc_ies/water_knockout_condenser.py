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
Water knockout condenser model

main equations:

* Heat removal by cooling water
* Calculate saturated liquid water being knocked out

"""

# Import Pyomo libraries
from pyomo.environ import SolverFactory, value, Var, Param, ExternalFunction,\
     PositiveReals, Reference, Constraint, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.models.properties import iapws95
from idaes.core.util.functions import functions_lib

#M
from idaes.models.unit_models import HeatExchangerFlowPattern

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog


# Additional import for the unit operation
from idaes.core.util.constants import Constants as const


__author__ = "J. Ma"
__version__ = "1.0.0"

@declare_process_block_class("WaterKnockoutCondenser")
class WaterKnockoutCondenserData(UnitModelBlockData):
    """
    WaterKnockoutCondenser Unit Class
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([useDefault, True, False]),
        default=useDefault,
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentPhase,
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
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
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
        default=True,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("tube_side_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("tube_side_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("shell_side_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("shell_side_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("water_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("water_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare(
        "flow_pattern",
        ConfigValue(
            default=HeatExchangerFlowPattern.countercurrent,
            domain=In(HeatExchangerFlowPattern),
            description="Heat exchanger flow pattern",
            doc="""Heat exchanger flow pattern,
**default** - HeatExchangerFlowPattern.countercurrent.
**Valid values:** {
**HeatExchangerFlowPattern.countercurrent** - countercurrent flow,
**HeatExchangerFlowPattern.cocurrent** - cocurrent flow,
**HeatExchangerFlowPattern.crossflow** - cross flow, factor times
countercurrent temperature difference.}""",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation)


        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(WaterKnockoutCondenserData,self).build()

        # Build Control Volumes
        self.tube = ControlVolume0DBlock(**{
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.tube_side_property_package,
                "property_package_args": self.config.tube_side_property_package_args})

        self.tube.add_state_blocks(has_phase_equilibrium=False)

        self.tube.add_material_balances(
            balance_type=self.config.material_balance_type)

        self.tube.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.tube.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        self.shell = ControlVolume0DBlock(**{
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.shell_side_property_package,
                "property_package_args": self.config.shell_side_property_package_args})

        self.shell.add_state_blocks(has_phase_equilibrium=False)

        self.shell.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        
        self.shell_water_out = self.config.water_property_package.\
            build_state_block(
                self.flowsheet().config.time,
                **self.config.water_property_package_args)
        

        # Add Ports
        self.add_inlet_port(name="tube_inlet", block=self.tube)
        self.add_outlet_port(name="tube_outlet", block=self.tube)
        self.add_inlet_port(name="shell_inlet", block=self.shell)
        self.add_outlet_port(name="shell_outlet", block=self.shell)
        self.add_port("shell_water_outlet", self.shell_water_out)

        # Set references to balance terms at unit level
        if (self.config.has_heat_transfer is True and
                self.config.energy_balance_type != EnergyBalanceType.none):
            self.heat_duty = Reference(self.tube.heat)

        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != 'none'):
            self.deltaP_tube = Reference(self.tube.deltaP)
            self.deltaP_shell = Reference(self.shell.deltaP)

        # Construct performance equations
        self._make_performance()

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """
        s1_metadata = self.config.tube_side_property_package.get_metadata()

        q_units = s1_metadata.get_derived_units("power")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        temp_units = s1_metadata.get_derived_units("temperature")
        pres_units = s1_metadata.get_derived_units("pressure")

        
        self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().config.time,
            domain=PositiveReals,
            initialize=100.0,
            doc="Overall heat transfer coefficient",
            units=u_units
        )
        self.area = Var(
            domain=PositiveReals,
            initialize=1000.0,
            doc="Heat exchange area",
            units=a_units
        )
        self.delta_temperature_in = Var(
            self.flowsheet().config.time,
            initialize=10.0,
            doc="Temperature difference at the hot inlet end",
            units=temp_units
        )
        self.delta_temperature_out = Var(
            self.flowsheet().config.time,
            initialize=10.1,
            doc="Temperature difference at the hot outlet end",
            units=temp_units
        )

        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_in_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_in[t]
                    == b.shell.properties_in[t].temperature
                    - b.tube.properties_in[t].temperature
                )
            else:
                return (
                    b.delta_temperature_in[t]
                    == b.shell.properties_in[t].temperature
                    - b.tube.properties_out[t].temperature
                )

        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_out_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (
                    b.delta_temperature_out[t]
                    == b.shell.properties_out[t].temperature
                    - b.tube.properties_out[t].temperature
                )
            else:
                return (
                    b.delta_temperature_out[t]
                    == b.shell.properties_out[t].temperature
                    - b.tube.properties_in[t].temperature
                )

        @self.Constraint(self.flowsheet().config.time,
                         doc="Water outlet pressure")
        def shell_water_outlet_pressure_eqn(b, t):
            return (b.shell_water_out[t].pressure*1e-5 
                == (b.shell.properties_in[t].pressure
                    + b.deltaP_shell[t]) * 1e-5)

        self.cbrt = ExternalFunction(
            library=functions_lib(),
            function="cbrt",
            arg_units=[temp_units])

        @self.Expression(self.flowsheet().config.time, doc="underwood temperature difference")
        def delta_temperature(b, t):
            dT1 = b.delta_temperature_in
            dT2 = b.delta_temperature_out
            temp_units = pyunits.get_units(dT1)
            return ((b.cbrt(dT1[t]) + b.cbrt(dT2[t])) / 2.0) ** 3 * temp_units

        @self.Constraint(self.flowsheet().config.time)
        def heat_transfer_equation(b, t):
            return pyunits.convert(b.heat_duty[t], to_units=q_units) == (
                    b.overall_heat_transfer_coefficient[t] * b.area * b.delta_temperature[t])

        @self.Constraint(self.flowsheet().config.time)
        def shell_outlet_temperature_equation(b, t):
            return b.shell.properties_out[t].temperature == b.shell_water_out[t].temperature

        @self.Constraint(self.flowsheet().config.time)
        def shell_energy_balance_equation(b, t):
            return (b.shell.properties_in[t].enth_mol*b.shell.properties_in[t].flow_mol
                    == b.shell.properties_out[t].enth_mol*b.shell.properties_out[t].flow_mol
                    + b.shell_water_out[t].flow_mol
                    *(b.shell_water_out[t].enth_mol - 1890.165 - 285840)  # consider difference in reference enthalpy values
                    + b.tube.heat[t]
                    )

        # Note that the current version does not consider the case without water condensation or not enough cooling
        # It may calculate negative liquid water flow at the water outlet port
        @self.Constraint(self.flowsheet().config.time)
        def water_in_shell_gas_out_equation(b, t):
            return 1e-4*b.shell_water_out[t].pressure_sat == 1e-4*b.shell.properties_out[t].pressure * b.shell.properties_out[t].mole_frac_comp['H2O']
        
        @self.Constraint(self.flowsheet().config.time, self.config.shell_side_property_package.component_list)
        def shell_species_mass_balance_equation(b, t, j):
            if j=='H2O':
                return Constraint.Skip
            else:
                return (b.shell.properties_in[t].flow_mol_comp[j]
                        == b.shell.properties_out[t].flow_mol_comp[j])

        @self.Constraint(self.flowsheet().config.time)
        def shell_overall_mass_balance_equation(b, t):
            return b.shell.properties_in[t].flow_mol == b.shell.properties_out[t].flow_mol + b.shell_water_out[t].flow_mol


    def set_initial_condition(self):
        ''' Initialization of dynamic accumulation terms '''
        if self.config.dynamic is True:
            self.tube.material_accumulation[:, :, :].value = 0
            self.tube.energy_accumulation[:, :].value = 0
            self.tube.material_accumulation[0, :, :].fix(0)
            self.tube.energy_accumulation[0, :].fix(0)

    def initialize(blk, state_args_tube=None, state_args_shell=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
        Water knockout condenser initialization routine.

        Keyword Arguments:
            state_args_tube : a dict of arguments to be passed to the property
                           package(s) for the control_volume of tube side to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            state_args_shell : a dict of arguments to be passed to the property
                           package(s) for the control_volume of shell side to
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
        flags_tube = blk.tube.initialize(outlvl=outlvl+1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args_tube)

        flags_shell = blk.shell.initialize(outlvl=outlvl+1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args_shell)

        temp_water = value(blk.tube.properties_in[0].temperature + 10)
        pres_water = value(blk.shell.properties_in[0].pressure - 10000)
        enth_mol_water = iapws95.htpx(P=pres_water*pyunits.Pa, T=temp_water*pyunits.K,)
        state_args = {"flow_mol": value(blk.shell.properties_in[0].flow_mol*0.8*blk.shell.properties_in[0].mole_frac_comp['H2O']),
                      "flow_enth": enth_mol_water,
                      "pressure": pres_water}
                      
        blk.shell_water_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args
        )

        dof = degrees_of_freedom(blk)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        blk.tube.release_state(flags_tube, outlvl+1)
        blk.shell.release_state(flags_shell, outlvl+1)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.heat_transfer_equation.items():
            s = iscale.get_scaling_factor(
                self.tube.heat[t], default=1e-7, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.shell_energy_balance_equation.items():
            s = iscale.get_scaling_factor(
                self.tube.heat[t], default=1e-7, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
