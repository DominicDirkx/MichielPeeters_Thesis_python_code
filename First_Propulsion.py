from math import *
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import interp2d

def Kerckhove_function(gamma_c):
    return np.sqrt(gamma_c * ((1 + gamma_c)/2)**((1+gamma_c)/(1-gamma_c)))

def FindPratio(gamma, eps):
    Guess = 0.1
    for _ in range(10):
        Guess = (Kerckhove_function(gamma)**2 / (eps**2 * (2 * gamma / (gamma - 1)) *
                                                 (1 - Guess**((gamma - 1) / gamma))))**(gamma / 2)
    return Guess


# engine_input = [100*10**5, 2869.9225, 0, 8, 22.3894*10**(-3), 1.1443, 0.92]
'''Rocket engine model'''
# General workings of this model is as follows: First rocket_engine_model_1 and find_rocket_engine_performance are
# used to find a rocket engine design (mass flow and exit diameter) which is capable of delivering 110 kN of sea level
# thrust. Then this mass flow is assumed constant during the duration of the rocket operations and the thrust level
# and Isp of the rocket engine is computed at all input p_a using rocket_engine_model_2
def rocket_engine_model(engine_input: dict, ambient_pressure):
    p_c = engine_input["p_c"]
    T_c = engine_input["T_c"]
    OF = engine_input["O/F"]
    gamma_c = engine_input["gamma_c"]
    M_m = engine_input["Molar_mass"]
    zeta_F = engine_input["zeta_F"]
    F_net_SL = engine_input["F_net_SL"]
    N_eng = engine_input["N_eng"]
    I_sp_vac = engine_input["Impulse"]
    p_a = 101325 # 110 kN given for sea level thrust
    R = 8.31446
    g_0 = 9.81  # m/s2




    return {"F_net_SL": F_net_SL, "F_net_1": F_net, "I_sp_1": I_sp,"O/F": OF, "m_dot_1": m_dot, "C_F_1": C_F,
                   "c_star_1": c_star, "D_e_1": D_e ,"A_t_1": A_t, "A_e_1": A_e, "E_ratio_1": A_e/A_t  ,"pe_pc_1":  pe_pc, "l_eng_1": l_eng,
                   "N_eng_1": N_eng}

def rocket_engine_model_1(engine_input:dict, engine_variables):
    # Unpack input
    p_c = engine_input["p_c"]
    T_c = engine_input["T_c"]
    OF = engine_input["O/F"]
    gamma_c = engine_input["gamma_c"]
    M_m = engine_input["Molar_mass"]
    zeta_F = engine_input["zeta_F"]
    F_net_SL = engine_input["F_net_SL"]
    N_eng = engine_input["N_eng"]

    m_dot = engine_variables
    p_a = 101325 # 110 kN given for sea level thrust

    R = 8.31446
    g_0 = 9.81  # m/s2

    # Compute Vandenkerckhove function
    Kerckhove = Kerckhove_function(gamma_c)

    # Compute c-star and throat area based on c-star
    c_star = 1/Kerckhove * np.sqrt(R * T_c/M_m)
    A_t = c_star*m_dot/p_c

    # Compute gross thrust required to provide the given net thrust, accounting for efficiencies
    F_t = F_net_SL / zeta_F
    # Compute thrust coefficient
    C_F = F_t/(p_c*A_t)

    # Create list for pressure ratios which provide the correct thrust coefficient values
    p_ratio = []
    C_F_long_list = []
    # Loop through several possible pressure ratios
    for x in np.arange(0.00001, 0.15, 0.00001):
        # Compute thrust coefficient using alternative formulation
        C_F_long = Kerckhove * np.sqrt(2*gamma_c/(gamma_c - 1) * (1 - x**((gamma_c - 1)/gamma_c))) + \
                   (x - p_a/p_c) * Kerckhove/np.sqrt(2*gamma_c/(gamma_c - 1) * x**(2/gamma_c) *
                                                                 (1 - x**((gamma_c - 1)/gamma_c)))
        # Check difference with known thrust coefficient, if almost equal, pressure ratio is stored as solution
        if abs(C_F - C_F_long) < 0.001 and x < 0.5:
            p_ratio.append(x)
            C_F_long_list.append(C_F_long)
    # Take minimum pressure ratio as solution (Large pressure ratios give undesired solutions
    if len(p_ratio) != 0:
        pe_pc = min(p_ratio)
        i = p_ratio.index(min(p_ratio))
        C_F_long = C_F_long_list[i]


        # Compute nozzle exit area and diameter
        A_e = A_t * Kerckhove/np.sqrt(2*gamma_c/(gamma_c - 1) * pe_pc**(2/gamma_c) *
                                                  (1 - pe_pc**((gamma_c - 1)/gamma_c)))
        D_e = np.sqrt(4 * A_e/np.pi)

        I_sp = C_F * c_star / g_0
        m_dot_fuel = m_dot / (1 + OF)
        m_dot_ox = OF / (1+OF) * m_dot
        rocket_design = {"F_net_SL": F_net_SL,"I_sp_SL": I_sp,"pe_pc_1": pe_pc,"D_e_1": D_e,"m_dot_1": m_dot, "O/F": OF,
                         "m_dot_ox_1": m_dot_ox, "m_dot_f_1": m_dot_fuel, "C_F_1": C_F,"C_F_long": C_F_long,"c_star_1": c_star,
                         "A_e_1": A_e, "A_t_1": A_t, "E_ratio_1": A_e/A_t,"N_eng_1": N_eng}

        return rocket_design, F_t, I_sp
        #f"Thrust:  {F_t} N, I_sp: {I_sp} s, Pressure ratio: {pe_pc}, D_e: {D_e} m, m_dot: {m_dot} kg/s," \
        #f" m_dot_fuel: {m_dot_fuel} kg/s, m_dot_ox: {m_dot_ox} kg/s, C_F: {C_F}, C_F_long: {C_F_long}"
    else:
        error = f"The required thrust cannot be produced by a realistic engine for {m_dot} kg/s"
        return

def find_rocket_engine_performance(engine_input):
    for m_dot in np.arange(30, 50, 0.1):
        if rocket_engine_model_1(engine_input, m_dot) != None:
            return {"m_dot_1": m_dot, "D_e_1": rocket_engine_model_1(engine_input, m_dot)[0]["D_e_1"], "F_t":
                    rocket_engine_model_1(engine_input, m_dot)[1],"I_sp": rocket_engine_model_1(engine_input, m_dot)[2]}

# print(find_rocket_engine_performance({"p_c": 10*10**6 , "T_c": 2869.9225 , "O/F": 8, "gamma_c":1.1443,
#                                       "Molar_mass": 22.3894*10**(-3), "zeta_F": 0.95 , "F_net_SL": 110*10**3, "N_eng": 2}))

def rocket_engine_model_2(engine_input:dict, environment_input, design_variables):
    # Define input parameters for rocket engine model

    # Unpack the engine input values
    #p_c, T_c, OF, M_m, gamma_c, zeta_F, F_net_SL, N_eng = engine_input
    p_c = engine_input["p_c"]
    T_c = engine_input["T_c"]
    OF = engine_input["O/F"]
    gamma_c = engine_input["gamma_c"]
    zeta_F = engine_input["zeta_F"]
    F_net_SL = engine_input["F_net_SL"]
    #I_sp = engine_input["I_sp"]
    N_eng = engine_input["N_eng"]
    M_m = engine_input["Molar_mass"]

    # Make sure quality factor = 0.95 --> Validation result
    zeta_F = 0.95
    # Unpack the engine design variables
    m_dot = design_variables["m_dot_1"]
    D_e = design_variables["D_e_1"]

    # unpack environment variables
    p_a = environment_input

    R = 8.31446
    g_0 = 9.81  # m/s2

    # Compute engine performance using ideal rocket theory
    # Compute the Vandenkerckhove function
    Kerckhove = Kerckhove_function(gamma_c)

    # Compute nozzle exit area
    A_e = np.pi/4 * D_e**2

    # Compute characteristic velocity
    c_star = 1/Kerckhove * np.sqrt(R * T_c/M_m)

    # Compute throat area, area ratio and obtain pressure ratio
    A_t = c_star * m_dot/p_c
    eps = A_e/A_t
    pe_pc = FindPratio(gamma_c, eps)
    U_e = np.sqrt(2 * gamma_c/(gamma_c-1) * R/M_m * T_c * (1 - pe_pc**((gamma_c-1)/gamma_c)))
    p_e = pe_pc * p_c

    # Compute thrust coefficient
    C_F = Kerckhove * np.sqrt((2 * gamma_c)/(gamma_c - 1) * (1 - (pe_pc)**((gamma_c - 1)/gamma_c))) +\
          (pe_pc - p_a/p_c) * eps

    # Compute specific impulse, thrust and corrected thrust.
    I_sp = C_F * c_star/ g_0
    F = m_dot * I_sp * g_0
    F_net = F*zeta_F



    # TODO Compute nozzle length and more detailed engine size in terms of diameter, nozzle length, etc
    l_eng = 0.1362 * F_net_SL**0.2279

    return {"F_net_SL": F_net_SL, "F_net_1": F_net, "I_sp_1": I_sp,"O/F": OF, "m_dot_1": m_dot, "C_F_1": C_F,
            "c_star_1": c_star, "D_e_1": D_e ,"A_t_1": A_t, "A_e_1": A_e, "E_ratio_1": A_e/A_t  ,"pe_pc_1":  pe_pc, "l_eng_1": l_eng,
            "N_eng_1": N_eng, "U_e": U_e, "p_e": p_e }


# print(rocket_engine_model_2(engine_input, environment))


# TODO: Add engine drag such that the installed thrust can be computed
'''Ramjet engine model'''
# Define input parameters for the airbreathing engine model
def ramjet_engine_model(environment_input:dict, ramjet_design:dict):
    p_a = environment_input["p_a"]
    T_a = environment_input["T_a"]
    V =  environment_input["V"]
    aoa = environment_input["AoA"]

    T_04 = ramjet_design["T_04"]
    A_i = ramjet_design["A_i"]
    r_d_max = ramjet_design["r_d_max"]
    r_c = ramjet_design["r_c"]
    eta_b = ramjet_design["eta_b"]
    N_ab_eng = ramjet_design["N_ab_eng"]

    zeta_F = 0.925
    zeta_m_dot = 1.38
    g_0 = 9.81
    gamma_air = 1.4
    Cp_air = 1005
    LHV_kerosene = 43*10**6
    R_air = Cp_air * (gamma_air - 1)/gamma_air

    M_a = V / np.sqrt(gamma_air * R_air * T_a)
    a_a = np.sqrt(gamma_air * R_air * T_a)

    T_0a = T_a * (1 + (gamma_air - 1)/2 * M_a**2)

    p_e = p_a


    # Compute inlet efficiency
    if M_a <= 1:
        eta_r = 1
    elif M_a > 1 and M_a <= 5:
        eta_r = 1 - 0.075*(M_a - 1)**1.35
    else:
        eta_r = 800/(M_a**4 + 935)

    r_d = r_d_max * eta_r

    m = 1 + (gamma_air - 1)/2 * M_a**2 * (r_d * r_c * r_n * (p_a/p_e)**((gamma_air-1)/gamma_air))
    M_e = 2 / (gamma_air - 1) * (m - 1)

    f = (Cp_air * T_04 - Cp_air * T_0a) / (eta_b * LHV_kerosene - Cp_air * T_04)

    u_e = np.sqrt((2 * gamma_air * R_air * T_04 * (m - 1))/((gamma_air - 1) * m))

    m_dot_a = p_a / (R_air * T_a) * A_i * V

    F_t = m_dot_a * ((1 + f)*u_e - V) #+ A_e * (p_e - p_a)
    F_t = F_t * zeta_F

    # Account for inlet and nozzle drag to obtain installed thrust.
    if M_a < 1:
        phi_inlet = 0.02
        phi_nozzle = 0.03
    else:
        A1_A0 = 1.2
        phi_inlet = ((A1_A0 - 1)*(1 - 1/M_a * (2/(gamma_air + 1) * (gamma_air - 1)/(gamma_air + 1) * M_a**2)**0.5))/(
        F_t / (m_dot_a * a_a ) * gamma_air * M_a)
        phi_nozzle = 0.03

    T_t = F_t - F_t*phi_inlet - F_t*phi_nozzle
    m_dot_f = f * m_dot_a
    m_dot_f = m_dot_f * zeta_m_dot

    TSFC = m_dot_f/F_t # Uninstalled
    I_sp = F_t/(m_dot_f * g_0) # Uninstalled

    return {"F_t_ab":T_t, "I_sp_ab": I_sp, "TSFC_ab": TSFC,"m_dot_a": m_dot_a,"m_dot_f_ab": m_dot_f,"T_04": T_04,
            "N_ab_eng": N_ab_eng}

# def ramjet_engine_model_Tc_assumed(vehicle_input, ramjet_variables, efficiencies):
#     # Input should be atmospheric pressure and temperature at which the vehicle is flying and velocity and angle of attack
#     pa, Ta, V, aoa = vehicle_input
#
#     # Ramjet variables are intake area and fuel mass flow
#     T0_4, A_inlet = ramjet_variables
#
#     # Get efficiencies, currently only Pi_d_max
#     pi_d_max, eta_b = efficiencies
#
#     # Define flow parameters assuming ideal and calorically perfect gas
#     g_0 = 9.81
#     gamma_air = 1.4
#     Cp_air = 1005 # J/kg*K
#     LHV_kerosene = 43*10**6 # J/kg
#     R_air = Cp_air * (gamma_air - 1)/gamma_air
#
#     # Account for vehicle aoa, assuming the engine inlet is perpendicular to the longitudinal body axis
#     Va = V*np.cos(np.deg2rad(aoa))
#
#     # Compute flight mach number, ambient stagnation temperature and stagnation temperature after inlet
#     M0 = Va/np.sqrt(Ta*gamma_air*R_air)
#     a0 = np.sqrt(Ta*gamma_air*R_air)
#
#     T0_a = Ta * (1 + (gamma_air - 1)/2 * M0**2)
#     T0_2 = T0_a
#
#     # Free stream temperature ratio
#     tau_r = 1 + (gamma_air -1)/2 * M0**2
#     tau_l = T0_4/Ta
#     # Temperature ratio at combustor
#     tau_b = tau_l/tau_r
#
#     # Compute pressure ratios
#     if M0 <= 1:
#         eta_r = 1
#     elif M0 > 1 and M0 <= 5:
#         eta_r = 1 - 0.075*(M0 - 1)**1.35
#     else:
#         eta_r = 800/(M0**4 + 935)
#
#     # Compute actual diffuser pressure ratio
#     pi_d = pi_d_max * eta_r
#
#     # Compute V9/a0
#     V9_a0 = np.sqrt(2*tau_l/(gamma_air - 1) * (1 - 1/(tau_r * pi_d**((gamma_air-1)/gamma_air))))
#     V9 = V9_a0 * a0
#     T_9 = Ta * tau_b / (pi_d**((gamma_air - 1)/gamma_air))
#     rho_9 = pa / (T_9 * R_air)
#     # Compute specific thrust
#     F_m0 = a0 * (V9_a0 - M0)
#
#
#     # Compute mass flow of the air through the intake
#     m_dot_air = A_inlet * Va * pa/(R_air*Ta)
#
#     # Compute thrust
#     F_t = F_m0 * m_dot_air
#
#     # Compute the fuel/air mass flow ratio f
#     # f = Cp_air * (T0_4 - T0_2) / (LHV_kerosene - Cp_air * T0_4)
#     # f = Cp_air * Ta * (tau_l - tau_r)/LHV_kerosene
#     f = Cp_air * (tau_l - tau_r)/(eta_b * LHV_kerosene/Ta - Cp_air * tau_l)
#
#     # Compute fuel mass flow
#     m_dot_f = f * m_dot_air
#
#     # Compute exit area
#     A_exit = (m_dot_air + m_dot_f)/(rho_9 * V9)
#
#
#     # Compute thrust specific fuel consumption and specific impulse
#     TSFC = m_dot_f/F_t
#     I_sp = F_t/(m_dot_f * g_0)
#
#     # Compute efficiencies
#     eta_t = (tau_b - tau_b/(tau_r*pi_d**((gamma_air-1)/gamma_air)) - 1 + 1/tau_r)/(tau_b - 1)
#
#     eta_p = 2/(np.sqrt(2*tau_l/(gamma_air-1) * (1 - 1/(tau_r*pi_d**((gamma_air-1)/gamma_air)))) + 1)
#
#     eta_o = eta_p*eta_t
#
#     # TODO Compute engine length and the cross-section profile
#     # TODO !!! ADD output dictionaries to several models! e.g. prop_to_mass{A_ratio: [], F_rocket: [], ...} !!!
#
#     return {"F_t_ab":F_t, "I_sp_ab": I_sp, "TSFC_ab": TSFC,"m_dot_a": m_dot_air,"m_dot_f_ab": m_dot_f,"T_04": T0_4}


'''Turbojet engine model'''
def turbojet_model(turbojet_design: dict, Mach ,altitude, afterburner: bool):
    N_ab_eng = turbojet_design["N_ab_eng"]
    filename = turbojet_design["filename"]
    sheet_1 = turbojet_design["sheet_no_afterburner"]
    sheet_2 = turbojet_design["sheet_afterburner"]
    Afterburner = afterburner



    # Extract data from excel sheet using pandas
    data_no_AB = pd.read_excel(filename, sheet_name=sheet_1)
    data_AB = pd.read_excel(filename, sheet_name=sheet_2)

    # Obtain separate data into desired dataframes
    mass_flow_TBF_df = data_no_AB.iloc[:,1:13]
    thrust_TBF_df = data_no_AB.iloc[:,15:27]
    mass_flow_AB_TBF_df = data_AB.iloc[:,1:13]
    thrust_AB_TBF_df = data_AB.iloc[:,15:27]

    # Transform dataframes to numpy matrices
    mdot_TBF_matrix = mass_flow_TBF_df.to_numpy()
    thrust_TBF_matrix = thrust_TBF_df.to_numpy()
    mdot_AB_TBF_matrix = mass_flow_AB_TBF_df.to_numpy()
    thrust_AB_TBF_matrix = thrust_AB_TBF_df.to_numpy()
    Mach_number_range = np.arange(0, 2.4, 0.2)
    altitude_range = np.arange(0, 20400, 400)

    # Interpolate between Mach numbers and altitude with spacing of 0.01 Mach and 100m
    m_dot_TBF_interp = interp2d(Mach_number_range, altitude_range, mdot_TBF_matrix)
    thrust_TBF_interp = interp2d(Mach_number_range, altitude_range, thrust_TBF_matrix)
    mdot_AB_TBF_interp = interp2d(Mach_number_range, altitude_range, mdot_AB_TBF_matrix)
    thrust_AB_TBF_interp = interp2d(Mach_number_range, altitude_range, thrust_AB_TBF_matrix)
    interp_Mach = np.arange(0,2.21, 0.01)
    interp_altitude = np.arange(0,20100, 100)
    interpolated_mdot = m_dot_TBF_interp(interp_Mach, interp_altitude)
    interpolated_thrust = thrust_TBF_interp(interp_Mach, interp_altitude)
    interpolated_AB_mdot = mdot_AB_TBF_interp(interp_Mach, interp_altitude)
    interpolated_AB_thrust = thrust_AB_TBF_interp(interp_Mach, interp_altitude)

    # Get array element closest to input Mach and altitude combination
    i = 0
    while abs(interp_Mach[i+1] - Mach) <= abs(interp_Mach[i] - Mach):
        i = i + 1
    closest_Mach_number = interp_Mach[i]
    j = 0
    while abs(interp_altitude[j + 1] - altitude) <= abs(interp_altitude[j] - altitude):
        j = j + 1
    closest_altitude = interp_altitude[j]

    # Compute installed thrust
    if Mach < 1:
        phi_inlet = 0.02
        phi_nozzle = 0.03
    else:
        phi_inlet = 0.015
        phi_nozzle = 0.03

    if Afterburner == True:
        F_t_ab =  interpolated_AB_thrust[j][i] * 1000
        T_t_ab = F_t_ab - F_t_ab*phi_inlet - F_t_ab*phi_nozzle
        m_dot_f_ab = interpolated_AB_mdot[j][i]
        I_sp = F_t_ab/(m_dot_f_ab * 9.81)
        TSFC = m_dot_f_ab/F_t_ab

    else:
        F_t_ab = interpolated_thrust[j][i] * 1000
        T_t_ab = F_t_ab - F_t_ab*phi_inlet - F_t_ab*phi_nozzle
        m_dot_f_ab = interpolated_mdot[j][i]
        I_sp = F_t_ab/(m_dot_f_ab * 9.81)
        TSFC = m_dot_f_ab/F_t_ab
    print(interpolated_thrust)
    # for j in range(len(interp_Mach)):
    #     plt.plot(interp_altitude, interpolated_AB_thrust[:,j], label = interp_Mach[j])
    # plt.xlabel("Altitude")
    # plt.legend()
    # plt.grid()
    # plt.show()
    # for j in range(len(interp_Mach)):
    #     plt.plot(interp_altitude, interpolated_AB_mdot[:,j], label = interp_Mach[j])
    # plt.grid()
    # plt.xlabel("Altitude")
    # plt.legend()
    # plt.show()

    for j in range(len(Mach_number_range)):
        plt.plot(altitude_range, thrust_TBF_matrix[:,j], label = Mach_number_range[j])
    plt.legend()
    plt.xlabel("Altitude [m]")
    plt.ylabel("Thrust [kN]")
    plt.title("Thrust vs Altitude")
    plt.grid()
    plt.show()

    for j in range(len(Mach_number_range)):
        plt.plot(altitude_range, thrust_AB_TBF_matrix[:,j], label = Mach_number_range[j])
    plt.grid()
    plt.legend()
    plt.xlabel("Altitude [m]")
    plt.ylabel("Thrust [kN]")
    plt.title("Thrust vs Altitude with Afterburner ")
    plt.show()

    for i in range(len(altitude_range)):
        plt.plot(Mach_number_range, thrust_AB_TBF_matrix[i], label = altitude_range[i])
    plt.grid()
    plt.legend()
    plt.xlabel("Mach [-]")
    plt.ylabel("Thrust [kN]")
    plt.title("Thrust vs Mach number with Afterburner ")
    plt.show()

    for i in range(len(altitude_range)):
        plt.plot(Mach_number_range, thrust_TBF_matrix[i], label = altitude_range[i])
    plt.grid()
    plt.legend()
    plt.xlabel("Mach [-]")
    plt.ylabel("Thrust [kN]")
    plt.title("Thrust vs Mach number without Afterburner ")
    plt.show()

    return {"F_t_ab": T_t_ab,"I_sp_ab": I_sp, "TSFC_ab": TSFC,"m_dot_a": 90.45,"m_dot_f_ab": m_dot_f_ab,"T_04": 1755,
            "N_ab_eng": N_ab_eng}




#turbojet_design = {"filename": "Input_TJ_data.xlsx", "sheet_no_afterburner": "Constant_TIT",
#                   "sheet_afterburner": "Constant_TIT_AB", "N_ab_eng": 2}
#print(turbojet_model(turbojet_design, 1.5, 12000, False))





















#
