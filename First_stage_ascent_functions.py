import math
import numpy as np
import pygmo as pg
from matplotlib import pyplot as plt
import First_Mass_Geometry
import First_Propulsion
import First_Aerodynamics

import tudatpy
from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.math import interpolators


"""
First the take-off is modelled using a simple euler integration according to: m dv/dt = T - D - mu*(m*g - L) 
The propagation is terminated once Lift off velocity is reached, defined as in W = L = 0.5*rho*V^2_LOF*S*C_L_LOF
The final conditions are initial conditions of the created vehicle.
"""

###########################################################################
# Define function for EOM during take-off roll #######################################################
###########################################################################
def take_off(CL, CD, rho, V, m, F_T , S_ref):
    g_0 = 9.81
    mu_f = 0.03 # 0.02 - 0.05
    L = CL * 0.5 * rho * V**2 * S_ref
    D = CD * 0.5 * rho * V**2 * S_ref
    dV_dt = (F_T - D - mu_f * (m*g_0 - L))/m
    return dV_dt



###########################################################################
# Define function for vehicle set-up #######################################################
###########################################################################
def vehicle_configuration(t_b_r, t_b_ab, rocket_engine_input, rocket_engine_design, upper_design):

    rocket_data_take_off = First_Propulsion.rocket_engine_model_2(rocket_engine_input, 0, rocket_engine_design)
    AB_engine_data_take_off = {}
    geometry_output = First_Mass_Geometry.first_stage_geometry(t_b_r,
                                                               t_b_ab,
                                                               rocket_data_take_off,
                                                               AB_engine_data_take_off,
                                                               upper_design,
                                                               False)

    vehicle_geometry = geometry_output["vehicle_geometry_1"]
    tank_geometry = geometry_output["tank_geometry_1"]
    vehicle_mass_output = First_Mass_Geometry.first_stage_mass(vehicle_geometry,
                                                               tank_geometry,
                                                               rocket_data_take_off,
                                                               AB_engine_data_take_off,
                                                               None,
                                                               t_b_r,
                                                               t_b_ab,
                                                               upper_design["m_upper"])

    vehicle = {"Gross_mass": vehicle_mass_output["Gross_mass"],
               "Propellant_mass": vehicle_mass_output["Propellant_mass"],
               "Dry_mass": vehicle_mass_output["Gross_mass"] - vehicle_mass_output["Propellant_mass"],
               "Fuel_mass": vehicle_mass_output["Fuel_mass"], "Ox_mass": vehicle_mass_output["Ox_mass"],
               "Wing_surface_area": vehicle_geometry["surface_area"], "p_e": rocket_data_take_off["p_e"],
               "U_e": rocket_data_take_off["U_e"], "A_e": rocket_data_take_off["A_e_1"]}
    return vehicle

###########################################################################
# Define function to get take_off_initial_state #######################################################
###########################################################################
def get_initial_take_off_state():
    V_0 = 0
    t_0 = 0
    x_0 = 0
    aoa_0 = 0
    gamma_0 = 0

    T = 288.15
    rho = 1.225
    aoa_LOF = 15    # deg

    return t_0, V_0, x_0, aoa_0, gamma_0, T, rho, aoa_LOF

###########################################################################
# Define function for plotting take-off #######################################################
###########################################################################
def plot_take_off_results(t_list, V_list, x_list, mass_list, thrust_list):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize = (20,17))

    ax1.plot(t_list, V_list)
    ax1.set_ylabel("Velocity [m/s]")

    ax2.plot(t_list, mass_list)
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("mass [kg]")

    ax3.plot(t_list, x_list)
    ax3.set_ylabel("Take-off distance [m]")

    ax4.plot(t_list, thrust_list)
    ax4.set_ylabel("Thrust [N]")

    for ax in fig.get_axes():
        ax.set_xlabel("Time [s]")
        ax.grid()
    plt.show()

"""From here on, the created functions are those for the ascent trajectory up until separation"""

###########################################################################
# Define function to get initial state for the ascent in inertial frame###############################################
###########################################################################
def get_initial_ascent_state(V_LOF, t_LOF, simulation_start_epoch,
                             bodies) -> np.ndarray:
    # Set initial time
    ascent_start_epoch = simulation_start_epoch + t_LOF
    # Set initial mass

    # Set initial spherical elements
    V = V_LOF

    radius = spice_interface.get_average_radius("Earth")
    gamma = 0
    latitude = np.deg2rad(51.951)
    longitude = np.deg2rad(4.4344)
    heading = np.deg2rad(360-14.758)


    # Transform spherical elements to body-fixed cartesian
    initial_cartesian_body_fixed_state = element_conversion.spherical_to_cartesian_elementwise(
        radius, latitude, longitude, V, gamma, heading
    )
    # Get rotationa emphermerides of the earth
    earth_rotation_model = bodies.get_body("Earth").rotation_model

    # Transform body-fixed cartesian state to global inertial frame
    initial_state_inertial_coordinates = environment.transform_to_inertial_orientation(
        initial_cartesian_body_fixed_state,
        ascent_start_epoch,
        earth_rotation_model)
    return  initial_state_inertial_coordinates

class CustomGuidanceModel:
    """"
    Custom guidance model which provides:
        - a thrust magnitude and Isp function as input for the custom engine model
        - aerodynamic angles for the rotation model
        - aerodynamic coefficients for the vehicle aerodynamic interface based on current angle of attack and Mach number
    The aerodynamic angles are all zero except for the angle of attack which is computed based on interpolated control
    nodes containing flight path angles, following the control law: aoa = Kc * (gamme_desired - gamma_current)
    """

    def __init__(self, bodies, engine_input, engine_design, K_c, initial_time, control_node_input):
        self.bodies = bodies
        self.vehicle = bodies.get_body("Mk-III")
        self.earth = bodies.get_body("Earth")

        environment_setup.add_flight_conditions(bodies, "Mk-III", "Earth")
        self.vehicle_flight_conditions = bodies.get_body("Mk-III").flight_conditions
        print(self.vehicle_flight_conditions.altitude)

        self.gain_factor = K_c
        self.initial_time = initial_time
        self.control_node_input = control_node_input
        self.time_interval = self.control_node_input[0]

        self.engine_input = engine_input
        self.engine_design = engine_design

        time_now = initial_time
        self.flight_path_angle_dict = {}
        self.flight_path_angle_derivative_dict = list()

        self.current_time = float("NaN")


        for i in range(len(self.control_node_input)-1):
            # Store time as key, control node flight path angle as value
            self.flight_path_angle_dict[time_now] = self.control_node_input[i + 1]
            # If using hermite spline interpolator, compute derivatives
            # Derivatives are computed by taking previous and next node, divided by 2 times time interval, except for
            # first and last node where derivative is set to 0
            if i == 0 or i + 1 == (len(self.control_node_input) - 1):
                self.flight_path_angle_derivative_dict.append(0.0)
            else:
                self.flight_path_angle_derivative_dict.append(
                    (control_node_input[i+2] - control_node_input[i]) / (2.0 * self.time_interval)
                )
            time_now += self.time_interval
        # Add final node 1 day after the final and set flight path angle to latest value and its derivative to 0
        self.flight_path_angle_dict[time_now + 86400.0] = control_node_input[-1]
        self.flight_path_angle_derivative_dict.append(0.0)

        # Create interpolator settings
        cubic_interpolator_settings = interpolators.hermite_spline_interpolation(
            boundary_interpolation=interpolators.use_boundary_value
        )
        linear_interpolator_settings = interpolators.linear_interpolation(
            boundary_interpolation=interpolators.use_boundary_value
        )

        # Create interpolator
        self.flight_path_cubic_interpolator = interpolators.create_one_dimensional_scalar_interpolator(
            self.flight_path_angle_dict, cubic_interpolator_settings, self.flight_path_angle_derivative_dict
        )

        self.flight_path_linear_interpolator = interpolators.create_one_dimensional_scalar_interpolator(
            self.flight_path_angle_dict, linear_interpolator_settings
        )

    def get_thrust_magnitude(self, current_time: float):
        self.update_guidance(current_time)
        print(current_time, self.current_time)
        return self.F_magnitude

    def get_Isp(self, current_time: float):
        self.update_guidance(current_time)

        return self.I_sp

    def get_aerodynamic_angles(self, current_time: float):
        self.update_guidance(current_time)

        return np.array([self.angle_of_attack, 0, 0])

    def get_aerodynamic_coefficients(self, current_time: float):
        self.update_guidance(current_time)

        return np.array([self.CD, 0, self.CL])

    def update_guidance(self, current_time: float):
        # if (math.isnan(current_time)):
        #     print("dis_not_work")
        #     self.current_time = float("NaN")
        #
        # elif (current_time != self.current_time):

        # Interpolate flight path angle between nodes at time
        desired_flight_path_angle = self.flight_path_cubic_interpolator.interpolate(current_time)
        # Let TUDAT know to update all values
        self.vehicle.flight_conditions.update_conditions(current_time)

        # Obtain flight path angle
        current_flight_path_calculator = self.vehicle.flight_conditions.aerodynamic_angle_calculator
        current_flight_path_angle = current_flight_path_calculator.get_angle(environment.flight_path_angle)
        # Compute current angle of attack using a parametric control law and return aerodynamic angles of attack
        self.angle_of_attack = self.gain_factor * (desired_flight_path_angle - current_flight_path_angle)

        # Thrust magnitude function which determines the thrust based on the engine design and the current time density
        self.F_magnitude = First_Propulsion.rocket_engine_model_2(
            self.engine_input,
            self.bodies.get("Earth").flight_conditions.density,
            self.engine_design)["F_net_1"] * self.engine_input["N_eng"]

        # Specific impulse function based on current thrust and the (constant) mass flow (and number of engines)
        g_0 = 9.80665
        self.I_sp = self.F_magnitude / (g_0 * self.engine_design["m_dot_1"] * self.engine_input["N_eng"])

        # Get current aerodynamic force coefficients
        current_mach = self.earth.flight_conditions.mach_number
        Mach_extrapolation_min = 2
        angle_of_attack_spacing = 20
        Mach_spacing = 20
        self.CL, self.CD = First_Aerodynamics.get_aerodynamic_data(Mach_spacing,
                                                                    angle_of_attack_spacing,
                                                                    current_mach,
                                                                    self.angle_of_attack,
                                                                    Mach_extrapolation_min)
        self.current_time = current_time
            #return self.F_magnitude, self.I_sp, np.array([self.angle_of_attack, 0, 0]), np.array([self.CD, 0, self.CL])

def get_propagator_settings(bodies,
                            simulation_start_epoch,
                            vehicle_initial_mass,
                            termination_settings,
                            dependent_variables_to_save,
                            mass_flow_settings,
                            initial_state):
    """Class which returns propagator settings for translational motion and mass based on:
    - Bodies: Earth and Mk-III
    - Mass flow settings:
        Mass flow rate from rocket_engine_design_variables
    - Simulation start epoch: Initial start epoch + t_takeoff
    - Termination settings: To be defined
    - Dependent variables to save: Position, altitude, velocity, flight path angle, angle of attack mass
    - Propagator: cowell
    - Initial state: array containing the inertial coordinates of the vehicle
    """

    # Define bodies to be propagated and their central body of propagation
    bodies_to_propagate = ["Mk-III"]
    central_bodies = ["Earth"]


    # Define accelerations acting on the vehicle: point mass gravity, aerodynamic forces, thrust
    thrust_acceleration_settings = propagation_setup.acceleration.thrust_from_engine("MainEngine")

    acceleration_settings_on_vehicle = {
        "Earth": [propagation_setup.acceleration.point_mass_gravity(),
                  propagation_setup.acceleration.aerodynamic()],
        "Mk-III": [thrust_acceleration_settings]
    }

    acceleration_settings = {"Mk-III": acceleration_settings_on_vehicle}
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        acceleration_settings,
        bodies_to_propagate,
        central_bodies
    )

    # Create fixed step runge kutta 4 integrator
    fixed_step_size = 1.0
    fixed_integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step_size(
        initial_time_step = fixed_step_size,
        coefficient_set= propagation_setup.integrator.rk_4
    )
    # Create variable stepsize RK4(5) integrator
    current_phase_start_time = simulation_start_epoch
    initial_time_step = 0.1
    minimum_step_size = 0.01
    maximum_step_size = 100

    #validation_settings = propagation_setup.integrator.step_size_validation(0.01, 100)
    integrator_settings_variable = propagation_setup.integrator.runge_kutta_variable_step_size(
        current_phase_start_time,
        initial_time_step,
        propagation_setup.integrator.rkf_45,
        minimum_step_size,
        maximum_step_size,
        10*10**(-6), 10*10**(-6)
    )

    # Create translational propagator settings
    translational_propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        simulation_start_epoch,
        integrator_settings_variable,
        termination_settings,
        output_variables=dependent_variables_to_save
    )

    # Create mass rate model
    mass_rate_settings_on_vehicle = {"Mk-III": [propagation_setup.mass_rate.from_thrust()]}#[propagation_setup.mass_rate.custom_mass_rate(mass_flow_settings)]}
    mass_rate_models = propagation_setup.create_mass_rate_models(bodies,
                                                                 mass_rate_settings_on_vehicle,
                                                                 acceleration_models)
    mass_propagator_settings = propagation_setup.propagator.mass(bodies_to_propagate,
                                                                 mass_rate_models,
                                                                 np.array([vehicle_initial_mass]),
                                                                 simulation_start_epoch,
                                                                 integrator_settings_variable,
                                                                 termination_settings
                                                                 )
    propagator_settings_list = [translational_propagator_settings,
                                mass_propagator_settings]

    propagator_settings = propagation_setup.propagator.multitype(propagator_settings_list,
                                                                 integrator_settings_variable,
                                                                 simulation_start_epoch,
                                                                 termination_settings,
                                                                 dependent_variables_to_save)
    return propagator_settings

def get_termination_settings(simulation_start_epoch: float,
                             maximum_duration: float,
                             termination_altitude: float,
                             vehicle_dry_mass: float) \
        -> tudatpy.kernel.numerical_simulation.propagation_setup.propagator.PropagationTerminationSettings:
    """
    # TODO: Review  upper termination altitude and mass termination setting
    :param simulation_start_epoch: start of ascent propagation = initial start epoch + takeoff duration
    :param maximum_duration: maximum duration of the simulation [s]
    :param termination_altitude: maximum altitude [m]
    :param vehicle_dry_mass: Dry mass of the vehicle (when all propellant has been expelled
    :return: hybrid_termination_settings
    """
    # Create PropagationTerminationSettings objects
    # Time
    time_termination_settings = propagation_setup.propagator.time_termination(
        simulation_start_epoch + maximum_duration,
        terminate_exactly_on_final_condition = False
    )
    # Altitude
    upper_altitude_termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings = propagation_setup.dependent_variable.altitude("Mk-III", "Earth"),
        limit_value = termination_altitude,
        use_as_lower_limit = False,
        terminate_exactly_on_final_condition = False
    )
    lower_altitude_termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings = propagation_setup.dependent_variable.altitude("Mk-III", "Earth"),
        limit_value = 0.0,
        use_as_lower_limit = True,
        terminate_exactly_on_final_condition = False

    )
    # Vehicle mass
    mass_termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings = propagation_setup.dependent_variable.body_mass("Mk-III"),
        limit_value = vehicle_dry_mass, # TODO: Change this to leave mass for return boost or something
        use_as_lower_limit = True,
        terminate_exactly_on_final_condition = False
    )

    termination_settings_list = [time_termination_settings,
                                 upper_altitude_termination_settings,
                                 lower_altitude_termination_settings,
                                 mass_termination_settings]
    hybrid_termination_settings = propagation_setup.propagator.hybrid_termination(termination_settings_list,
                                                                                  fulfill_single_condition = True)

    return hybrid_termination_settings

def get_dependent_variable_save_settings() -> list:
    """"
    Retrieves dependent variables to save

    Altitude w.r.t. earth
    Velocity w.r.t. earth
    flight path angle of the vehicle
    angle of attack of the vehicle
    downrange of the vehicle
    Mach number
    """
    dependent_variables_to_save = [propagation_setup.dependent_variable.altitude("Mk-III", "Earth"),
                                   propagation_setup.dependent_variable.mach_number("Mk-III", "Earth"),
                                   propagation_setup.dependent_variable.relative_speed("Mk-III", "Earth"),
                                   propagation_setup.dependent_variable.flight_path_angle("Mk-III", "Earth"),
                                   propagation_setup.dependent_variable.body_mass("Mk-III"),
                                   propagation_setup.dependent_variable.relative_position("Mk-III", "Earth"),
                                   propagation_setup.dependent_variable.angle_of_attack("Mk-III", "Earth")]

    return dependent_variables_to_save


