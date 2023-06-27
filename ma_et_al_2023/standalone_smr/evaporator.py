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
Convective evaporator model

main equations:

* Two-phase flow and heat transfer in evaporator riser tubes


"""
from enum import Enum
# Import Pyomo libraries
from pyomo.environ import SolverFactory, value, Var, Param, \
     sqrt, log10, PositiveReals, Reference, units as pyunits
from pyomo.dae import DerivativeVar
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
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog


# Additional import for the unit operation
from idaes.core.util.constants import Constants as const


__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "1.0.0"

class TubeArrangement(Enum):
    inLine = 0
    staggered = 1

@declare_process_block_class("Evaporator")
class EvaporatorData(UnitModelBlockData):
    """
    Evaporator Unit Class
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
    CONFIG.declare("tube_arrangement", ConfigValue(
        default=TubeArrangement.inLine,
        domain=In(TubeArrangement),
        description='tube configuration',
        doc='Tube arrangement could be in-line and staggered'))

    def build(self):
        """
        Begin building model (pre-DAE transformation)


        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(EvaporatorData,self).build()

        # Build Control Volumes
        self.tube_cv = ControlVolume0DBlock(default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.tube_side_property_package,
                "property_package_args": self.config.tube_side_property_package_args})

        self.shell_cv = ControlVolume0DBlock(default={
                "dynamic": self.config.dynamic,
                "has_holdup": self.config.has_holdup,
                "property_package": self.config.shell_side_property_package,
                "property_package_args": self.config.shell_side_property_package_args})
                
        self.tube_cv.add_geometry()
        self.shell_cv.add_geometry()
        # This model requires the IAPWS95 property package with the mixed phase
        # option, therefore, phase equilibrium calculations are handled by
        # the property package.
        self.tube_cv.add_state_blocks(has_phase_equilibrium=False)

        self.tube_cv.add_material_balances(
            balance_type=self.config.material_balance_type)

        self.tube_cv.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.tube_cv.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        self.shell_cv.add_state_blocks(has_phase_equilibrium=False)

        self.shell_cv.add_material_balances(
            balance_type=self.config.material_balance_type)

        self.shell_cv.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.shell_cv.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=True)

        # Add Ports
        self.add_inlet_port(name="tube_inlet", block=self.tube_cv)
        self.add_inlet_port(name="shell_inlet", block=self.shell_cv)
        self.add_outlet_port(name="tube_outlet", block=self.tube_cv)
        self.add_outlet_port(name="shell_outlet", block=self.shell_cv)

        # Add object references
        self.tube_volume = Reference(self.tube_cv.volume)
        self.shell_volume = Reference(self.shell_cv.volume)

        # Set references to balance terms at unit level
        if (self.config.has_heat_transfer is True and
                self.config.energy_balance_type != EnergyBalanceType.none):
            self.tube_heat_duty = Reference(self.tube_cv.heat)
            self.shell_heat_duty = Reference(self.shell_cv.heat)

        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != 'none'):
            self.tube_deltaP = Reference(self.tube_cv.deltaP)
            self.shell_deltaP = Reference(self.shell_cv.deltaP)

        # Set Unit Geometry and Holdup Volume
        self._set_geometry()

        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary, and link to holdup volume
        """
        units_meta = self.config.tube_side_property_package.get_metadata()
        # Number of tube rows
        self.number_tube_rows = Var(
                initialize=4,
                doc="Number of tube rows")
        # Number of tube cols
        self.number_tube_cols = Var(
                initialize=4,
                doc="Number of tube columns")
        @self.Expression(doc="Total number of tubes")
        def number_tubes(b):
            return b.number_tube_rows*b.number_tube_cols
        # Height of tube bundle from bottom to top for the evaporator
        self.height = Var(
                initialize=5.0,
                doc="Height of tube bundle",
                units=units_meta.get_derived_units("length"))
        # Inner diameter of tubes
        self.tube_diameter_inner = Var(
                initialize=0.05,
                doc="Inside diameter of tube",
                units=units_meta.get_derived_units("length"))
        # Total cross section area of fluid flow
        @self.Expression(doc="Cross section area of fluid")
        def area_cross_fluid_total(b):
            return 0.25*const.pi*b.tube_diameter_inner**2*b.number_tubes
        # Tube thickness
        self.tube_thickness = Var(
                initialize=0.005,
                doc="Thickness of tube",
                units=units_meta.get_derived_units("length"))
        # Outer diameter tube
        @self.Expression(doc="Outside diameter of tube")
        def tube_diameter_outer(b):
            return b.tube_diameter_inner + b.tube_thickness*2
        # Pitch in flow direction between two neighboring tube rows
        self.pitch_x = Var(
                initialize=0.05,
                doc="Pitch in flow direction",
                units=units_meta.get_derived_units("length"))
        # Pitch between two neighboring tube columns perpendicular to flow
        self.pitch_y = Var(
                initialize=0.05,
                doc="Pitch between tube columns",
                units=units_meta.get_derived_units("length"))
        # Cross section area of a single tube
        @self.Expression(doc="Cross section area of tube metal")
        def area_cross_metal(b):
            return 0.25*const.pi*(b.tube_diameter_outer**2-b.tube_diameter_inner**2)
        # flow area of shell side
        @self.Expression(doc="Cross section area of tube metal")
        def area_flow_shell(b):
            return b.height*(b.pitch_y-b.tube_diameter_outer)*b.number_tube_cols
        # pitch x to outer diameter ratio
        @self.Expression(doc="Pitch x to outer diameter ratio")
        def pitch_x_to_do(b):
            return b.pitch_x/b.tube_diameter_outer
        # pitch y to outer diameter ratio
        @self.Expression(doc="Pitch y to outer diameter ratio")
        def pitch_y_to_do(b):
            return b.pitch_y/b.tube_diameter_outer
        # Tube volume constraint
        @self.Constraint(self.flowsheet().config.time,
                         doc="fluid volume of all tubes")
        def tube_volume_eqn(b, t):
            return b.tube_volume[t] == b.area_cross_fluid_total*b.height
        # Tube volume constraint
        @self.Constraint(self.flowsheet().config.time,
                         doc="fluid volume of all tubes")
        def shell_volume_eqn(b, t):
            return b.shell_volume[t] == b.height*b.number_tubes\
            *(b.pitch_x*b.pitch_y - 0.25*const.pi*b.tube_diameter_outer**2)


    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """
        units_meta = self.config.tube_side_property_package.get_metadata()
        # Correction factor for pressure drop on tube side
        self.fcorrection_dp_tube = Var(
                initialize=1.2,
                within=PositiveReals,
                doc='correction factor for pressure drop due to acceleration'
                'and unsmooth tube applied to friction term')
        # Correction factor for pressure drop on shell side
        self.fcorrection_dp_shell = Var(
                initialize=1.0,
                within=PositiveReals,
                doc='correction factor for pressure drop on shell side')
        # Thermal conductivity of metal
        self.therm_cond_metal = Param(
                initialize=43.0,
                mutable=True,
                doc='Thermal conductivity of tube metal',
                units=units_meta.get_derived_units("thermal_conductivity"))
        # Heat capacity of metal
        self.cp_metal = Param(
                initialize=500.0,
                mutable=True,
                doc='Heat capacity of tube metal',
                units=units_meta.get_derived_units("heat_capacity_mass"))
        # Density of metal
        self.dens_metal = Param(
                initialize=7800.0,
                mutable=True,
                doc='Density of tube metal',
                units=units_meta.get_derived_units("density_mass"))

        # Heat conduction resistance of half of metal wall thickness
        # based on interface perimeter
        @self.Expression(doc="half metal layer conduction resistance")
        def half_resistance_metal(b):
            return b.tube_thickness/2/b.therm_cond_metal

        # Add performance variables
        # Tube inner wall temperature
        self.temp_tube_inner = Var(
                self.flowsheet().config.time,
                initialize=400.0,
                doc='Temperature of tube inner wall',
                units=units_meta.get_derived_units("temperature"))
        # Tube center point wall temperature
        self.temp_tube_center = Var(
                self.flowsheet().config.time,
                initialize=450.0,
                doc='Temperature of tube center wall',
                units=units_meta.get_derived_units("temperature"))
        # Tube outer wall temperature
        self.temp_tube_outer = Var(
                self.flowsheet().config.time,
                initialize=600.0,
                doc='Temperature of tube outer wall',
                units=units_meta.get_derived_units("temperature"))
        # Energy holdup of tube metal of a single tube per tube length
        self.energy_holdup_metal = Var(
                self.flowsheet().config.time,
                initialize=1e6,
                doc='Energy holdup of tube metal per length',
                units=(units_meta.get_derived_units("energy") *
                       units_meta.get_derived_units("length")**-1))

        # Energy accumulation for metal
        if self.config.dynamic is True:
            self.energy_accumulation_metal = DerivativeVar(
                self.energy_holdup_metal,
                wrt=self.flowsheet().config.time,
                doc='Energy accumulation of tube metal')

        def energy_accumulation_term_metal(b, t):
            return b.energy_accumulation_metal[t] if b.config.dynamic else 0

        # Velocity of liquid only
        self.velocity_liquid = Var(
                self.flowsheet().config.time,
                initialize=3.0,
                doc='Velocity of liquid only',
                units=units_meta.get_derived_units("velocity"))

        # Reynolds number based on liquid only flow
        self.N_Re = Var(
                self.flowsheet().config.time,
                initialize=1.0e6,
                doc='Reynolds number')

        # Prandtl number of liquid phase
        self.N_Pr = Var(
                self.flowsheet().config.time,
                initialize=2.0,
                doc='Reynolds number')

        # Darcy friction factor
        self.friction_factor_darcy = Var(
                self.flowsheet().config.time,
                initialize=0.01,
                doc='Darcy friction factor')

        # Vapor fraction at inlet
        self.vapor_fraction = Var(
                self.flowsheet().config.time,
                initialize=0.0,
                doc='Vapor fractoin of vapor-liquid mixture')

        # Liquid fraction at inlet
        self.liquid_fraction = Var(
                self.flowsheet().config.time,
                initialize=1.0,
                doc='Liquid fractoin of vapor-liquid mixture')

        # Density ratio of liquid to vapor
        self.ratio_density = Var(
                self.flowsheet().config.time,
                initialize=2.0,
                doc='Liquid to vapor density ratio')

        # Void fraction at inlet
        self.void_fraction = Var(
                self.flowsheet().config.time,
                initialize=0.0,
                doc='void fraction at inlet')

        # Exponent n for gamma at inlet, typical range in (0.75,0.8294)
        self.n_exp = Var(
                self.flowsheet().config.time,
                initialize=0.82,
                doc='exponent for gamma at inlet')

        # Gamma for velocity slip at inlet
        self.gamma = Var(
                self.flowsheet().config.time,
                initialize=3.0,
                doc='gamma for velocity slip at inlet')

        # Mass flux
        self.mass_flux = Var(
                self.flowsheet().config.time,
                initialize=1000.0,
                doc='mass flux',
                units=units_meta.get_derived_units("flux_mass"))

        # Reduced pressure
        self.reduced_pressure = Var(
                self.flowsheet().config.time,
                initialize=0.85,
                doc='reduced pressure')

        # Two-phase correction factor
        self.phi_correction = Var(
                self.flowsheet().config.time,
                initialize=1.01,
                doc='Two-phase flow correction factor')

        # Convective heat transfer coefficient on shell side,
        self.hconv_shell = Var(
                self.flowsheet().config.time,
                initialize=50.0,
                doc='Convective heat transfer coefficient on shell side',
                units=units_meta.get_derived_units(
                    "heat_transfer_coefficient"))

        # Convective heat transfer coefficient on tube side,
        # typically in range (1000, 5e5)
        self.hconv_tube = Var(
                self.flowsheet().config.time,
                initialize=10000.0,
                doc='Convective heat transfer coefficient on tube side',
                units=units_meta.get_derived_units(
                    "heat_transfer_coefficient"))

        # Convective heat transfer coefficient for liquid only,
        # typically in range (1000.0, 1e5)
        self.hconv_liquid = Var(
                self.flowsheet().config.time,
                initialize=20000.0,
                doc='Convective heat transfer coefficient of liquid only',
                units=units_meta.get_derived_units(
                    "heat_transfer_coefficient"))

        # Pool boiling heat transfer coefficient, typically in range (1e4, 5e5)
        self.hpool = Var(
                self.flowsheet().config.time,
                initialize=1e5,
                doc='Pool boiling heat transfer coefficient',
                units=units_meta.get_derived_units(
                    "heat_transfer_coefficient"))

        # Boiling number, typical range in (1e-7, 5e-4) in original formula.
        # we define here as boiling_number_scaled == 1e6*boiling_number
        self.boiling_number_scaled = Var(
                self.flowsheet().config.time,
                initialize=1,
                doc='Scaled boiling number')

        # Enhancement factor, typical range in (1.0, 3.0)
        self.enhancement_factor = Var(
                self.flowsheet().config.time,
                initialize=1.3,
                doc='Enhancement factor')

        # Suppression factor, typical range in (0.005, 1.0)
        self.suppression_factor = Var(
                self.flowsheet().config.time,
                initialize=0.03,
                doc='Suppression factor')

        # Convective heat flux to fluid on tube side
        self.heat_flux_tube = Var(
                self.flowsheet().config.time,
                initialize=7e4,
                doc='Convective heat flux to fluid',
                units=units_meta.get_derived_units("flux_energy"))

        # Convective heat flux from flue gas on shell side
        self.heat_flux_shell = Var(
                self.flowsheet().config.time,
                initialize=100000.0,
                doc='Shell side heat flux',
                units=units_meta.get_derived_units("flux_energy"))

        # Pressure change due to friction
        self.deltaP_friction_tube = Var(
                self.flowsheet().config.time,
                initialize=-1000.0,
                doc='Pressure change due to friction on tube side',
                units=units_meta.get_derived_units("pressure"))

        # Pressure change due to gravity
        self.deltaP_gravity_tube = Var(
                self.flowsheet().config.time,
                initialize=-1000.0,
                doc='Pressure change due to gravity on tube side',
                units=units_meta.get_derived_units("pressure"))

        # Equation to calculate tube outside wall temperature
        @self.Constraint(self.flowsheet().config.time,
                         doc="tube outside wall temperature")
        def outside_wall_temperature_eqn(b, t):
            return b.heat_flux_shell[t] * b.half_resistance_metal == \
                   (b.temp_tube_outer[t] - b.temp_tube_center[t])

        # Equation to calculate tube inside wall temperature
        @self.Constraint(self.flowsheet().config.time,
                         doc="tube inside wall temperature")
        def inside_wall_temperature_eqn(b, t):
            return b.heat_flux_tube[t] * b.half_resistance_metal == \
               (b.temp_tube_center[t] - b.temp_tube_inner[t])

        # Equation to calculate heat flux at inside tube boundary
        @self.Constraint(self.flowsheet().config.time,
                         doc="convective heat flux at tube inside boundary")
        def heat_flux_tube_eqn(b, t):
            return b.heat_flux_tube[t] == b.hconv_tube[t] * (b.temp_tube_inner[t] \
                   - (b.tube_cv.properties_in[t].temperature \
                   +  b.tube_cv.properties_out[t].temperature)/2 )

        # Equation to calculate heat flux at outside tube boundary
        @self.Constraint(self.flowsheet().config.time,
                         doc="convective heat flux at tube outside boundary")
        def heat_flux_shell_eqn(b, t):
            return b.heat_flux_shell[t] == b.hconv_shell[t] * (\
                   ( b.shell_cv.properties_in[t].temperature \
                   + b.shell_cv.properties_out[t].temperature)/2 \
                   - b.temp_tube_outer[t])

        # Equation to calculate energy holdup for tube metal
        # per tube length for a single tube
        @self.Constraint(self.flowsheet().config.time,
                         doc="energy holdup for metal")
        def energy_holdup_metal_eqn(b, t):
            return b.energy_holdup_metal[t] == b.temp_tube_center[t] \
                * b.cp_metal * b.dens_metal * b.area_cross_metal

        # Energy balance for metal
        @self.Constraint(self.flowsheet().config.time,
                         doc="energy balance for metal")
        def energy_balance_metal_eqn(b, t):
            return energy_accumulation_term_metal(b, t) == const.pi*(\
                   b.heat_flux_shell[t]*b.tube_diameter_outer - \
                   b.heat_flux_tube[t]*b.tube_diameter_inner)

        # Equations for calculate pressure drop
        # and convective heat transfer coefficient for 2-phase flow
        # Equation to calculate average liquid to vapor density ratio
        @self.Constraint(self.flowsheet().config.time,
                         doc="liquid to vapor density ratio")
        def ratio_density_eqn(b, t):
            '''
            return 2e-5*b.ratio_density[t] \
                * (b.tube_cv.properties_in[t].dens_mol_phase["Vap"] \
                + b.tube_cv.properties_out[t].dens_mol_phase["Vap"]) == \
                2e-5*(b.tube_cv.properties_in[t].dens_mol_phase["Liq"] \
                   + b.tube_cv.properties_out[t].dens_mol_phase["Liq"])
            '''
            return 2e-5*b.ratio_density[t] \
                * b.tube_cv.properties_out[t].dens_mol_phase["Vap"] == \
                2e-5*b.tube_cv.properties_out[t].dens_mol_phase["Liq"]

        # Equation for calculating velocity if the flow is liquid only
        # use inlet density since liquid density does not change much
        @self.Constraint(self.flowsheet().config.time,
                         doc="Vecolity of fluid if liquid only")
        def velocity_lo_eqn(b, t):
            return 1e-4*b.velocity_liquid[t]*b.area_cross_fluid_total * \
                   b.tube_cv.properties_in[t].dens_mol_phase["Liq"] \
                   == 1e-4*b.tube_cv.properties_in[t].flow_mol

        # Equation for calculating Reynolds number if liquid only
        # use inlet property since liquid phase property does not change much
        @self.Constraint(self.flowsheet().config.time,
                         doc="Reynolds number if liquid only")
        def Reynolds_number_eqn(b, t):
            return b.N_Re[t] * \
                   b.tube_cv.properties_in[t].visc_d_phase["Liq"] == \
                   b.tube_diameter_inner * b.velocity_liquid[t] * \
                   b.tube_cv.properties_in[t].dens_mass_phase["Liq"]

        # Friction factor depending on laminar or turbulent flow,
        # usually always turbulent (>1187.385)
        @self.Constraint(self.flowsheet().config.time,
                         doc="Darcy friction factor")
        def friction_factor_darcy_eqn(b, t):
            return b.friction_factor_darcy[t]*b.N_Re[t]**0.25/0.3164 == 1.0

        # Vapor fractoin equation, use half of outlet as average since inlet is zero
        # add 1e-5 such that vapor fraction is always positive
        @self.Constraint(self.flowsheet().config.time,
                         doc="Average vapor fractoin")
        def vapor_fraction_eqn(b, t):
            return 100*b.vapor_fraction[t] == \
                100*(b.tube_cv.properties_out[t].vapor_frac/2+ 1e-5)

        # n-exponent equation for inlet
        @self.Constraint(self.flowsheet().config.time, doc="n-exponent")
        def n_exp_eqn(b, t):
            return (0.001*(0.8294 - b.n_exp[t]) *
                    b.tube_cv.properties_in[t].pressure ==
                    8.0478*units_meta.get_derived_units("pressure"))

        # Gamma equation for inlet
        @self.Constraint(self.flowsheet().config.time, doc="Gamma at inlet")
        def gamma_eqn(b, t):
            return b.gamma[t] == b.ratio_density[t]**b.n_exp[t]

        # void faction at inlet equation
        @self.Constraint(self.flowsheet().config.time,
                         doc="Void fractoin at inlet")
        def void_fraction_eqn(b, t):
            return b.void_fraction[t] * (1.0 + b.vapor_fraction[t]
                                         * (b.gamma[t] - 1.0)) == \
                b.vapor_fraction[t] * b.gamma[t]

        # Two-phase flow correction factor equation
        @self.Constraint(self.flowsheet().config.time, doc="Correction factor")
        def correction_factor_eqn(b, t):
            return (b.phi_correction[t] - 0.027*b.liquid_fraction[t])**2 == \
                (0.973*b.liquid_fraction[t] + b.vapor_fraction[t]
                 * b.ratio_density[t]) * \
                (0.973*b.liquid_fraction[t] + b.vapor_fraction[t])

        # Pressure change equation due to friction,
        # -1/2*density*velocity^2*fD/diameter*length*phi^2
        @self.Constraint(self.flowsheet().config.time,
                         doc="pressure change due to friction")
        def pressure_change_friction_eqn(b, t):
            return 0.01 * b.deltaP_friction_tube[t] * b.tube_diameter_inner == \
                -0.01 * b.fcorrection_dp_tube * 0.5 \
                * b.tube_cv.properties_in[t].dens_mass_phase["Liq"] * \
                b.velocity_liquid[t]**2 * b.friction_factor_darcy[t] \
                * b.height * b.phi_correction[t]**2

        # Pressure change equation due to gravity,
        # density_mixture*gravity*height, use vapor density as average density
        @self.Constraint(self.flowsheet().config.time,
                         doc="pressure change due to gravity")
        def pressure_change_gravity_eqn(b, t):
            return 1e-4 * b.deltaP_gravity_tube[t] == -1e-4 \
                * const.acceleration_gravity * b.height \
                * (b.tube_cv.properties_out[t].dens_mass_phase["Vap"] \
                   * b.void_fraction[t] +
                   b.tube_cv.properties_out[t].dens_mass_phase["Liq"]
                   * (1.0 - b.void_fraction[t]))

        # Mass flux of vapor-liquid mixture
        # (density*velocity or mass_flow/area)
        @self.Constraint(self.flowsheet().config.time, doc="mass flux")
        def mass_flux_eqn(b, t):
            return b.mass_flux[t] * b.area_cross_fluid_total == \
                   b.tube_cv.properties_in[t].flow_mol * \
                   b.tube_cv.properties_in[0].mw

        # Liquid fraction at inlet
        @self.Constraint(self.flowsheet().config.time, doc="liquid fraction")
        def liquid_fraction_eqn(b, t):
            return b.liquid_fraction[t] + b.vapor_fraction[t] == 1.0

        # Total pressure change equation on tube side
        @self.Constraint(self.flowsheet().config.time, doc="pressure drop")
        def pressure_change_total_tube_eqn(b, t):
            return b.tube_deltaP[t] == b.deltaP_friction_tube[t] + b.deltaP_gravity_tube[t]

        # Total heat added to tube control_volume
        @self.Constraint(self.flowsheet().config.time,
                         doc="total heat added to fluid control_volume on tube side")
        def tube_heat_eqn(b, t):
            return b.tube_heat_duty[t] == b.number_tubes * b.heat_flux_tube[t] \
                * b.height * b.tube_diameter_inner * const.pi

        # Total heat added to shell control_volume
        @self.Constraint(self.flowsheet().config.time,
                         doc="total heat added to fluid control_volume on shell side")
        def shell_heat_eqn(b, t):
            return b.shell_heat_duty[t] == -b.number_tubes * b.heat_flux_shell[t] \
                * b.height * b.tube_diameter_outer * const.pi

        # Reduced pressure
        @self.Constraint(self.flowsheet().config.time, doc="reduced pressure")
        def reduced_pressure_eqn(b, t):
            return b.reduced_pressure[t] \
                * self.config.tube_side_property_package.pressure_crit == \
                b.tube_cv.properties_in[t].pressure

        # Prandtl number of liquid
        @self.Constraint(self.flowsheet().config.time,
                         doc="liquid Prandtl number")
        def N_Pr_eqn(b, t):
            return b.N_Pr[t] \
                * b.tube_cv.properties_in[t].therm_cond_phase["Liq"] \
                * b.tube_cv.properties_in[0].mw == \
                b.tube_cv.properties_in[t].cp_mol_phase["Liq"] * \
                b.tube_cv.properties_in[t].visc_d_phase["Liq"]

        # Forced convection heat transfer coefficient for liquid only
        @self.Constraint(self.flowsheet().config.time,
                         doc="forced convection heat transfer "
                         "coefficient for liquid only")
        def hconv_lo_eqn(b, t):
            return b.hconv_liquid[t] * b.tube_diameter_inner == \
                0.023 * b.N_Re[t]**0.8 * b.N_Pr[t]**0.4 * \
                b.tube_cv.properties_in[t].therm_cond_phase["Liq"]

        # Pool boiling heat transfer coefficient
        @self.Constraint(self.flowsheet().config.time,
                         doc="pool boiling heat transfer coefficient")
        def hpool_eqn(b, t):
            return (1e-4*b.hpool[t]*sqrt(pyunits.convert(
                        b.tube_cv.properties_in[0].mw,
                        to_units=pyunits.g/pyunits.mol)) *
                    (-log10(b.reduced_pressure[t]))**(0.55) ==
                    1e-4*55.0*b.reduced_pressure[t]**0.12 *
                    b.heat_flux_tube[t]**0.67)

        # Boiling number scaled by a factor of 1e6
        @self.Constraint(self.flowsheet().config.time, doc="boiling number")
        def boiling_number_eqn(b, t):
            return 1e-10*b.boiling_number_scaled[t] \
                * b.tube_cv.properties_in[t].dh_vap_mol \
                * b.mass_flux[t] == \
                b.heat_flux_tube[t]*b.tube_cv.properties_in[0].mw*1e-4

        # Enhancement factor
        # due to low contribution the reciprocal Martinalli parameter,
        # it can be removed from the enhancement factor equation
        @self.Constraint(self.flowsheet().config.time,
                         doc="Forced convection enhancement factor")
        def enhancement_factor_eqn(b, t):
            return b.enhancement_factor[t] == 1.0 + 24000.0 \
                    * (b.boiling_number_scaled[t]/1e6)**1.16

        # Suppression factor
        @self.Constraint(self.flowsheet().config.time,
                         doc="Pool boiler suppression factor")
        def suppression_factor_eqn(b, t):
            return b.suppression_factor[t] \
                * (1.0 + 1.15e-6*b.enhancement_factor[t]**2
                   * b.N_Re[t]**1.17) == 1.0

        @self.Constraint(self.flowsheet().config.time,
                         doc="convective heat transfer coefficient")
        def hconv_tube_eqn(b, t):
            return 1e-3*b.hconv_tube[t] == 1e-3 * b.hconv_liquid[t] \
                * b.enhancement_factor[t] + 1e-3 * b.hpool[t] \
                * b.suppression_factor[t]

        # Shell side heat transfer coefficient and pressure drop
        if self.config.tube_arrangement == TubeArrangement.inLine:
            self.f_arrangement = Param(initialize=0.788,
                                       doc="In-line tube arrangement factor")
        elif self.config.tube_arrangement == TubeArrangement.staggered:
            self.f_arrangement = Param(initialize=1.0,
                                       doc="Staggered tube arrangement factor")
        else:
            raise Exception('tube arrangement type not supported')
        # Velocity on shell side
        self.v_shell = Var(self.flowsheet().time,
                           initialize=10.0,
                           doc="Velocity on shell side - m/s",
                           units=pyunits.m/pyunits.s)

        # Reynalds number on shell side
        self.N_Re_shell = Var(self.flowsheet().time,
                              initialize=10000.0,
                              doc="Reynolds number on shell side")

        # Friction factor on shell side
        self.friction_factor_shell = Var(self.flowsheet().time,
                                         initialize=0.1,
                                         doc='Friction factor on shell side')

        # Prandtl number on shell side
        self.N_Pr_shell = Var(self.flowsheet().time,
                              initialize=0.8,
                              doc="Prandtl number on shell side")

        # Nusselt number on shell side
        self.N_Nu_shell = Var(self.flowsheet().time,
                              initialize=10.0,
                              doc="Nusselts number on shell side")

        # Velocity equation on shell side
        @self.Constraint(self.flowsheet().time, doc="Velocity on shell side")
        def v_shell_eqn(b, t):
            return b.v_shell[t] * \
                b.shell_cv.properties_in[t].dens_mol_phase["Vap"] * \
                b.area_flow_shell == b.shell_cv.properties_in[t].flow_mol

        # Reynolds number
        @self.Constraint(self.flowsheet().time,
                         doc="Reynolds number equation on shell side")
        def N_Re_shell_eqn(b, t):
            return b.N_Re_shell[t] * b.shell_cv.properties_in[t].visc_d == \
                b.tube_diameter_outer * b.v_shell[t] \
                * b.shell_cv.properties_in[t].dens_mol_phase["Vap"] \
                * b.shell_cv.properties_in[t].mw

        # Friction factor on shell side
        if self.config.tube_arrangement == TubeArrangement.inLine:
            @self.Constraint(self.flowsheet().time,
                             doc="In-line friction factor on shell side")
            def friction_factor_shell_eqn(b, t):
                return b.friction_factor_shell[t] \
                    * b.N_Re_shell[t]**0.15 == \
                    (0.044 + 0.08 * b.pitch_x_to_do
                     / (b.pitch_y_to_do - 1.0)**(0.43 + 1.13
                                                 / b.pitch_x_to_do)
                     ) * b.fcorrection_dp_shell
        elif self.config.tube_arrangement == TubeArrangement.staggered:
            @self.Constraint(self.flowsheet().time,
                             doc="Staggered friction factor on shell side")
            def friction_factor_shell_eqn(b, t):
                return b.friction_factor_shell[t] \
                    * b.N_Re_shell[t]**0.16 == \
                    (0.25 + 0.118 / (b.pitch_y_to_do - 1.0)**1.08) \
                    * b.fcorrection_dp_shell
        else:
            raise Exception('tube arrangement type not supported')

        # Pressure drop on shell side
        @self.Constraint(self.flowsheet().time,
                         doc="Pressure change on shell side")
        def shell_deltaP_eqn(b, t):
            return (
                b.shell_deltaP[t] ==
                -1.4 * b.friction_factor_shell[t] * b.number_tube_rows *
                b.shell_cv.properties_in[t].dens_mol_phase["Vap"] *
                b.shell_cv.properties_in[t].mw *
                b.v_shell[t]**2)

        # Prandtl number
        @self.Constraint(self.flowsheet().time,
                         doc="Prandtl number equation on shell side")
        def N_Pr_shell_eqn(b, t):
            return b.N_Pr_shell[t] * b.shell_cv.properties_in[t].therm_cond \
                * b.shell_cv.properties_in[t].mw == \
                b.shell_cv.properties_in[t].cp_mol * \
                b.shell_cv.properties_in[t].visc_d

        # Nusselt number, currently assume Re>300
        @self.Constraint(self.flowsheet().time,
                         doc="Nusselts number equation on shell side")
        def N_Nu_shell_eqn(b, t):
            return b.N_Nu_shell[t] == b.f_arrangement * 0.33 \
                * b.N_Re_shell[t]**0.6 * b.N_Pr_shell[t]**0.333333

        # Convective heat transfer coefficient on shell side due to convection
        @self.Constraint(self.flowsheet().time,
                         doc="Convective heat transfer coefficient equation"
                         "on shell side due to convection")
        def hconv_shell_eqn(b, t):
            return b.hconv_shell[t] * b.tube_diameter_outer  == \
                b.N_Nu_shell[t] * b.shell_cv.properties_in[t].therm_cond



    def set_initial_condition(self):
        ''' Initialization of dynamic accumulation terms '''

        if self.config.dynamic is True:
            self.tube_cv.material_accumulation[:, :, :].value = 0
            self.tube_cv.energy_accumulation[:, :].value = 0
            self.tube_cv.material_accumulation[0, :, :].fix(0)
            self.tube_cv.energy_accumulation[0, :].fix(0)
            self.energy_accumulation_metal[:].value = 0
            self.energy_accumulation_metal[0].fix(0)

    def initialize(blk, state_args_tube=None, state_args_shell=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
        Evaporator initialization routine.

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
        flags_tube = blk.tube_cv.initialize(outlvl=outlvl+1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args_tube)
        flags_shell = blk.shell_cv.initialize(outlvl=outlvl+1,
                                              optarg=optarg,
                                              solver=solver,
                                              state_args=state_args_shell)
        init_log.info_high("Initialization Step 1 Complete.")
        # Fix outlet enthalpy and pressure for tube side
        for t in blk.flowsheet().config.time:
            blk.tube_cv.properties_out[t].enth_mol.fix(
                value(blk.tube_cv.properties_in[t].enth_mol) + 100)
            blk.tube_cv.properties_out[t].pressure.fix(
                value(blk.tube_cv.properties_in[t].pressure) - 1.0)
        blk.tube_heat_eqn.deactivate()
        blk.pressure_change_total_tube_eqn.deactivate()
        # Fix outlet temperature and pressure for shell side
        for t in blk.flowsheet().config.time:
            blk.shell_cv.properties_out[t].temperature.fix(
                value(blk.shell_cv.properties_in[t].temperature) - 1)
            blk.shell_cv.properties_out[t].pressure.fix(
                value(blk.shell_cv.properties_in[t].pressure) - 1.0)
        # Estimate wall temperatures
        blk.temp_tube_inner[:].value = value(blk.tube_cv.properties_in[0].temperature) + 1
        blk.temp_tube_center[:].value = value(blk.tube_cv.properties_in[0].temperature) + 2
        blk.temp_tube_outer[:].value = value(blk.tube_cv.properties_in[0].temperature) + 3
        blk.hconv_shell.fix(100)
        blk.hconv_shell_eqn.deactivate()
        blk.shell_heat_eqn.deactivate()
        blk.shell_deltaP_eqn.deactivate()
        # Solve with fixed outlet conditions
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )
        # Unfix outlet enthalpy and pressure
        for t in blk.flowsheet().config.time:
            blk.tube_cv.properties_out[t].enth_mol.unfix()
            blk.tube_cv.properties_out[t].pressure.unfix()
        blk.tube_heat_eqn.activate()
        blk.pressure_change_total_tube_eqn.activate()
        for t in blk.flowsheet().config.time:
            blk.shell_cv.properties_out[t].temperature.unfix()
            blk.shell_cv.properties_out[t].pressure.unfix()
        blk.hconv_shell.unfix()
        blk.hconv_shell_eqn.activate()
        blk.shell_heat_eqn.activate()
        blk.shell_deltaP_eqn.activate()
        # Solve with all constraits activated
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 3 {}.".format(idaeslog.condition(res))
            )

        blk.tube_cv.release_state(flags_tube, outlvl+1)
        blk.shell_cv.release_state(flags_shell, outlvl+1)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.heat_flux_tube_eqn.items():
            s = iscale.get_scaling_factor(
                self.heat_flux_tube[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.heat_flux_shell_eqn.items():
            s = iscale.get_scaling_factor(
                self.heat_flux_shell[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.tube_heat_eqn.items():
            s = iscale.get_scaling_factor(
                self.tube_heat_duty[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.shell_heat_eqn.items():
            s = iscale.get_scaling_factor(
                self.shell_heat_duty[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.energy_holdup_metal_eqn.items():
            s = iscale.get_scaling_factor(
                self.energy_holdup_metal[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)

