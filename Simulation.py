###########################################################################
# IMPORT STATEMENTS #######################################################
###########################################################################

from math import *
import numpy as np
import pygmo as pg
from matplotlib import pyplot as plt
import First_Aerodynamics
import First_Mass_Geometry
import First_Propulsion
from First_stage_ascent_functions import take_off
from First_stage_ascent_functions import vehicle_configuration
from First_stage_ascent_functions import get_initial_take_off_state
from First_stage_ascent_functions import plot_take_off_results
from First_stage_ascent_functions import get_initial_ascent_state
from First_stage_ascent_functions import CustomGuidanceModel
from First_stage_ascent_functions import get_termination_settings
from First_stage_ascent_functions import get_dependent_variable_save_settings
from First_stage_ascent_functions import get_propagator_settings

import tudatpy
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import environment

spice_interface.load_standard_kernels()

"""
First the take-off is modelled using a simple euler integration according to: m dv/dt = T - D - mu*(m*g - L) 
The propagation is terminated once Lift off velocity is reached, defined as in W = L = 0.5*rho*V^2_LOF*S*C_L_LOF
The final conditions are initial conditions of the created vehicle.
"""
###########################################################################
# Set design variables and vehicle configuration #######################################################
###########################################################################
t_b_r = 200
t_b_ab = 0
print(1)
# Set-up vehicle design
rocket_engine_input = {"p_c": 10*10**6 , "T_c": 2869.9225 , "O/F": 8, "gamma_c":1.1443, "Molar_mass": 22.3894*10**(-3), "zeta_F": 0.95 ,
                "F_net_SL": 110*10**3, "N_eng": 2}
rocket_engine_design = {"m_dot_1": 42.5, "D_e_1": 0.3677}

upper_design = {"l_upper": 5.54, "m_upper": 6008} # {"l_upper": 5.86, "m_upper": 7002} {"l_upper": 4.93, "m_upper": 3687}
vehicle_design = vehicle_configuration(t_b_r,
                                       t_b_ab,
                                       rocket_engine_input,
                                       rocket_engine_design,
                                       upper_design)
mass, propellant_mass, vehicle_dry_mass, fuel_mass, ox_mass, S_ref, A_e, p_e, U_e  = vehicle_design["Gross_mass"],\
                                                                                     vehicle_design["Propellant_mass"],\
                                                                                     vehicle_design["Dry_mass"],\
                                                                                     vehicle_design["Fuel_mass"],\
                                                                                     vehicle_design["Ox_mass"],\
                                                                                     vehicle_design["Wing_surface_area"],\
                                                                                     vehicle_design["A_e"],\
                                                                                     vehicle_design["p_e"],\
                                                                                     vehicle_design["U_e"]

N_eng = rocket_engine_input["N_eng"]
F_TO = rocket_engine_input["F_net_SL"] * N_eng
m_dot = rocket_engine_design["m_dot_1"] * N_eng
###########################################################################
# Set propagator settings for take off: dt #######################################################
###########################################################################
dt = 0.5    # m/s

###########################################################################
# Propagate take-off roll #######################################################
###########################################################################
# Set take-off conditions
LOF = False
g_0 = 9.81

t, V, x, aoa, gamma, T, rho, aoa_LOF = get_initial_take_off_state()

# Create lists to store data
t_list_TO = [t]
V_list_TO = [V]
x_list_TO = [x]
mass_list_TO = [mass]
Mach_list_TO = [0]
Thrust_list_TO = [F_TO]

# Propagate until lift-off speed is reached at which aoa is increased to 15 degrees
while LOF == False:
    Mach = V / np.sqrt(T * 287 * 1.4)
    CL, CD = First_Aerodynamics.get_aerodynamic_data(20,20, Mach, aoa, 2)
    dV_dt = take_off(CL, CD, rho, V, mass, F_TO, S_ref)

    CL_LOF = First_Aerodynamics.get_aerodynamic_data(20,20, Mach, aoa_LOF,2)[0]
    V_LOF = np.sqrt((mass * g_0)/(CL_LOF * 0.5 * rho * S_ref))

    if V >= V_LOF:
        LOF = True
        break
    mass = mass - dt * m_dot
    x = x + V * dt
    V = V + dV_dt * dt
    t = t + dt

    t_list_TO.append(t)
    mass_list_TO.append(mass)
    x_list_TO.append(x)
    V_list_TO.append(V)
    Mach_list_TO.append(Mach)
    Thrust_list_TO.append(F_TO)

V_TO = V_list_TO[-1]
mass_TO = mass_list_TO[-1]
aoa_TO = aoa_LOF
t_TO = t_list_TO[-1]

# First_stage_ascent_functions.plot_take_off_results(t_list_TO, V_list_TO, x_list_TO, mass_list_TO, Thrust_list_TO)

""""From here one, the ascent simulation is coded"""

###########################################################################
# Create environment ######################################################
###########################################################################
bodies_to_create = ["Earth"]
global_frame_origin = "Earth"
global_fame_orientation = 'J2000'
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create,
    global_frame_origin,
    global_fame_orientation
)
#body_settings.get("Earth").atmosphere_settings = environment_setup.atmosphere.us76()
bodies = environment_setup.create_system_of_bodies(body_settings)

# Create Mk-III vehicle object
bodies.create_empty_body("Mk-III")
# Set vehicle mass
bodies.get_body("Mk-III").mass = mass_TO

###########################################################################
# Set initial conditions based on take-off performance #######################################################
###########################################################################
simulation_start_epoch = 0.0    # s
vehicle_initial_state_inertial = get_initial_ascent_state(V_TO, t_TO, simulation_start_epoch, bodies)
vehicle_mass = mass_TO
aoa = aoa_LOF
simulation_start_epoch = simulation_start_epoch + t_TO
###########################################################################
# Control node input and simulation settings ######################################################
###########################################################################
K_c = 1.5
t_b_r_remaining = t_b_r - t_TO
N_control_nodes = 5
control_nodes = [t_b_r_remaining/N_control_nodes, 0, 45*np.pi/180, 50*np.pi/180, 65*np.pi/180, 60*np.pi/180]
maximum_duration = 3600 # s
termination_altitude = 100.0*10**3  # m

###########################################################################
# Create thrust settings acting on the vehicle ######################
###########################################################################

# Create Guidance model model for both thrust magnitude (always acting along body axis) and aerodynamic angles
Guidance_model = CustomGuidanceModel(bodies,
                                     rocket_engine_input,
                                     rocket_engine_design,
                                     K_c,
                                     simulation_start_epoch,
                                     control_nodes)

thrust_magnitude_settings = (
    propagation_setup.thrust.custom_thrust_magnitude(
        lambda time: Guidance_model.get_thrust_magnitude, lambda time: Guidance_model.get_Isp))


# # Create engine model using magnitude settings and thrust direction function with input time
environment_setup.add_engine_model("Mk-III", "MainEngine", thrust_magnitude_settings, bodies, np.array([1, 0, 0]))


###########################################################################
# Create aerodynamic acceleration settings acting on the vehicle ######################
###########################################################################

# Get aerodynamic coefficients as function of "time" (time determines current angle of attack and Mach number)
aerodynamic_coefficients = Guidance_model.get_aerodynamic_coefficients
# Define aerodynamic coefficient settings and add aerodynamic interface
aerodynamic_coefficient_settings = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
    lambda time: aerodynamic_coefficients,
    S_ref,
    independent_variable_names = [environment.AerodynamicCoefficientsIndependentVariables.time_dependent])
environment_setup.add_aerodynamic_coefficient_interface(bodies, "Mk-III", aerodynamic_coefficient_settings)

# Get aerodynamic angles (only angle of attack is non-zero) at time t.
aerodynamic_angles = Guidance_model.get_aerodynamic_angles
# Define settings for rotation model from inertial to body frame?
rotation_model_settings = environment_setup.rotation_model.aerodynamic_angle_based(
    "Earth", "J2000", "Vehicle_fixed", lambda time: aerodynamic_angles)

environment_setup.add_rotation_model(bodies, "Mk-III", rotation_model_settings)


###########################################################################
# CREATE PROPAGATOR SETTINGS ##############################################
###########################################################################
# Get termination settings
termination_settings = get_termination_settings(simulation_start_epoch,
                                                maximum_duration,
                                                termination_altitude,
                                                vehicle_dry_mass )
# Get dependent variables to save
dependent_variables_to_save = get_dependent_variable_save_settings()

# Define propagator settings
propagator_settings = get_propagator_settings(bodies,
                                              simulation_start_epoch,
                                              vehicle_mass,
                                              termination_settings,
                                              dependent_variables_to_save,
                                              lambda time: rocket_engine_design["m_dot_1"]*rocket_engine_input["N_eng"],
                                              vehicle_initial_state_inertial)
print(bodies.get_body("Mk-III"))
# Create dynamics_simulator

dynamics_simulator = numerical_simulation.create_dynamics_simulator(bodies, propagator_settings)

# Output simulation / process data
state_history = dynamics_simulator.state_history
dependent_variable_history = dynamics_simulator.dependent_variable_history

print(state_history)














