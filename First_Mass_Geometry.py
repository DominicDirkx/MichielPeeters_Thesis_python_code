from math import *
import numpy as np

# TODO: Change fuel tank to partial storage inside the wings

Scaling_factor = 5
''''Geometry model'''
# First stage has same geometrical shape as the mk-II
def first_stage_geometry(t_b_rocket, t_b_ab, rocket_data: dict, AB_engine_data: dict, upper_design: dict, AB_engine:bool):
    # Define rocket propellant and tank material parameters

    m_dot_rocket = rocket_data["m_dot_1"]
    F_rocket = rocket_data["F_net_SL"]
    OF = rocket_data["O/F"]
    N_eng = rocket_data["N_eng_1"]

    if AB_engine == True:
        m_dot_ab = AB_engine_data["m_dot_f_ab"]
        N_ab_eng = AB_engine_data["N_ab_eng"]
    else:
        m_dot_ab = 0
        N_ab_eng = 0
        t_b_ab = 0

    l_upper = upper_design["l_upper"]


    rho_ox = 1.43*10**3  # 98% purity
    rho_f = 0.81*10**3
    sigma_ult = 603 * 10**6  # Yield/ultimate strength of the tank material
    FOSU = 1.7  # Factor of Safety for ultimate strength = KP*KM * FOSU = 1.05*1.15*1.375
    MEOP = 5*10**5  # Maximum expected operating pressure
    t_min = 1 * 10**(-3)  # Minimum tank material thickness

    # Compute propellant tank volume
    V_ox_1 = OF/(OF + 1) * m_dot_rocket * N_eng * t_b_rocket / rho_ox
    V_tank_ox_1 = 1.05 * V_ox_1

    V_f_1 = (m_dot_rocket* N_eng * t_b_rocket)/((OF + 1) * rho_f) + m_dot_ab * N_ab_eng * t_b_ab / rho_f
    V_tank_f_1 = 1.05 * V_f_1

    l_engine_1 = 0.1362 * F_rocket**0.2279

    # Compute cylindrical segment length based on propellant volume (for common bulkhead tanks), upper stage length and engine length
    # l_cyl_guess = 20
    # length_found = False
    # D_t_f_1_spherical = (V_tank_f_1 * 6/np.pi)**(1/3)
    # while length_found == False:
    #     print(3, V_tank_f_1, V_tank_ox_1, l_upper, l_engine_1, l_cyl_guess, D_t_f_1_spherical, l_cyl_guess/10-0.1)
    #     #if D_t_f_1_spherical < (l_cyl_guess/10 - 0.1):
    #     l_cyl = l_upper + l_engine_1 +\
    #             (V_tank_ox_1 - np.pi/6 * (l_cyl_guess/10 - 0.1)**3)/(np.pi/4 * (l_cyl_guess/10 - 0.1)**2) +\
    #             (l_cyl_guess/10 - 0.1) + D_t_f_1_spherical
    #     #else:
    #     #    l_cyl = l_upper + l_engine_1 +\
    #     #            (V_tank_ox_1 - np.pi/6 * (l_cyl_guess/10 - 0.1)**3)/(np.pi/4 * (l_cyl_guess/10 - 0.1)**2) +\
    #     #            (l_cyl_guess/10 - 0.1) + (V_tank_f_1)/(np.pi/4 * (l_cyl_guess/10 - 0.1)**2)
    #
    #     if abs(l_cyl_guess - l_cyl) < 0.01:
    #         length_found = True
    #     else:
    #         l_cyl_guess = l_cyl

    # Set cylindrical and fuselage lengths to 20 and 25m respectively (minimum size because of upper stage design that is frozen)

    l_cyl = 20
    l_fus = 5/4 * l_cyl

    if l_fus < 25:
        SF = 5
        l_fus = 25
        l_cyl = 4/5 * l_fus
    else:
        SF = l_fus/5
    l_nose = 1/5 * l_fus
    D_fus = 0.4 * SF

    # Compute propellant tank sizes, spherical if possible, if not then cylindrical (and common bulkhead)
    D_sphere_f = 2 * (V_tank_f_1/(4 * np.pi/3))**(1/3)
    D_sphere_ox = 2 * (V_tank_ox_1/(4 * np.pi/3))**(1/3)

    if D_sphere_ox <= D_fus - 0.1 and D_sphere_f <= D_fus - 0.1:
        D_t_ox_1 = D_sphere_ox
        D_t_f_1 = D_sphere_f
        l_tank_ox_1 = D_t_ox_1
        l_tank_f_1 = D_t_f_1

        l_cyl_ox_1 = 0
        l_cyl_f_1 = 0
        l_tanks = l_tank_ox_1 + l_tank_f_1

        t_tank_ox = MEOP * D_t_ox_1 / (4 * sigma_ult) * FOSU
        t_tank_f = MEOP * D_t_f_1 / (4 * sigma_ult) * FOSU
        t_cyl_ox = t_tank_ox
        t_cyl_f = t_tank_f

    elif D_sphere_ox > D_fus - 0.1 and D_sphere_f <= D_fus - 0.1:
        D_t_f_1 = D_sphere_f
        l_tank_f_1 = D_t_f_1
        l_cyl_f_1 = 0
        t_tank_f = MEOP * D_t_f_1 / (4 * sigma_ult) * FOSU
        t_cyl_f = t_tank_f

        D_t_ox_1 = D_fus - 0.1
        l_cyl_ox_1 = (V_tank_ox_1 - np.pi/6 * D_t_ox_1**3)/(np.pi/4 * D_t_ox_1**2)
        l_tank_ox_1 = l_cyl_ox_1 + D_t_ox_1
        l_tanks = l_tank_f_1 + l_tank_ox_1

        t_tank_ox = MEOP * D_t_ox_1 / (4 * sigma_ult) * FOSU
        t_cyl_ox = 2 * t_tank_ox

    else: #D_sphere_ox and D_sphere_f > 1.8:
        D_t_ox_1 = D_fus - 0.1
        D_t_f_1 = D_fus - 0.1

        l_tank_ox_1 = (V_tank_ox_1 - np.pi/6 * D_t_ox_1**3)/(np.pi/4 * D_t_ox_1**2) + D_t_ox_1
        l_cyl_ox_1 = l_tank_ox_1 - D_t_ox_1

        l_tank_f_1 = V_tank_f_1 / (np.pi/4 * D_t_f_1**2) + D_t_f_1/2
        l_cyl_f_1 = l_tank_f_1 - D_t_f_1/2
        l_tanks = l_tank_f_1 + l_tank_ox_1 - D_t_f_1/2

        t_tank_ox = MEOP * D_t_ox_1 / (4 * sigma_ult) * FOSU
        t_tank_f = MEOP * D_t_f_1 / (4 * sigma_ult) * FOSU
        t_cyl_ox = 2 * t_tank_ox
        t_cyl_f = 2 * t_tank_f

    # Compute propellant tank dimensions
    # t_tank_ox = MEOP * D_t_ox_1 / (4 * sigma_ult) * FOSU  # Spherical cap tank wall thickness
    # t_tank_f = MEOP * D_t_f_1 / (4 * sigma_ult) * FOSU
    # t_cyl_ox = 2 * t_tank_ox
    # t_cyl_f = 2 * t_tank_f

    # Add condition for minimum tank material thickness
    if t_tank_ox < t_min:
        t_tank_ox = t_min
        if t_cyl_ox < t_min:
            t_cyl_ox = 2 * t_min
    if t_tank_f < t_min:
        t_tank_f = t_min
        if t_cyl_f < t_min:
            t_cyl_f = 2 * t_min

    # Define wing planform parameters
    Lambda_LE = 50.0  # degrees
    Lambda_c4 = 41.7  # degrees
    Lambda_c2 = 30.53
    Lambda_TE = -0.685 # degrees
    lambda_w = 0.2

    b = 2.4 * SF
    S = 2.6 * SF**2
    c_r = 2 * S / ((1 + lambda_w) * b)
    c_t = c_r * lambda_w
    c_a = 1.1 * SF
    A = b**2 / S
    tc = 0.1

    # Define vertical wing parameters
    Lambda_LE_v = 43
    Lambda_c4_v = 39.16
    Lambda_TE_v = 24
    lambda_v = 0.5

    b_v = 0.65 * SF
    S_v = 0.3 * SF**2
    c_r_v = 2 * S_v / ((1 + lambda_v) * b_v)
    c_t_v = c_r_v * lambda_v
    A_v = b_v**2/S_v

    tank_geometry = {"V_tank_ox_1": V_tank_ox_1, "V_tank_f_1": V_tank_f_1, "l_tank_ox_1": l_tank_ox_1,
                     "l_cyl_ox_1": l_cyl_ox_1, "l_tank_f_1": l_tank_f_1, "l_cyl_f_1": l_cyl_f_1,"t_tank_ox": t_tank_ox,
                     "t_tank_f": t_tank_f,"t_cyl_ox": t_cyl_ox, "t_cyl_f": t_cyl_f, "D_ox_tank": D_t_ox_1,
                     "D_f_tank": D_t_f_1}
    vehicle_geometry = {"Lambda_LE": Lambda_LE, "Lambda_c2": Lambda_c2, "Lambda_c4": Lambda_c4, "Lambda_TE": Lambda_TE ,
                        "wing_taper": lambda_w, "span": b, "surface_area": S, "root_chord": c_r, "tip_chord": c_t,
                        "aerodynamic_chord": c_a, "aspect_ratio": A, "t/c": tc, "Lambda_LE_v": Lambda_LE_v,
                        "Lambda_c4_v": Lambda_c4_v, "Lambda_TE_v": Lambda_TE_v, "tail_taper": lambda_v, "span_v": b_v,
                        "surface_area_v": S_v, "root_chord_v": c_r_v, "tip_chord_v": c_t_v, "aspect_ratio_v": A_v,
                        "l_fuselage": l_fus, "D_fuselage": D_fus, "l_nose": l_nose}

    return {"vehicle_geometry_1": vehicle_geometry, "tank_geometry_1": tank_geometry, "l_engine_1": l_engine_1}


'''Mass model'''
def first_stage_mass(vehicle_geometry: dict, tank_geometry: dict, rocket_data, AB_data, AB_engine_type,
                     t_b_rocket, t_b_ab, M_upper):
    # Unit transformations
    lb_kg = 0.4535924
    ft_m = 0.3048
    g0 = 9.81

    # Tank material density
    rho_t_m = 1550  # For Carbon-epoxy composite UD prepeg and QI lay-up


    # Define modifying factor
    mf = 1.0 # See HASA by NASA for graph
    # Unpack variables and change m, m2 to ft, ft2
    l_fus = vehicle_geometry["l_fuselage"]/ft_m
    l_n = vehicle_geometry["l_nose"]/ft_m
    D_fus = vehicle_geometry["D_fuselage"]/ft_m
    S_ref = vehicle_geometry["surface_area"]/ft_m**2
    b = vehicle_geometry["span"]/ft_m
    A = vehicle_geometry["aspect_ratio"]
    lambda_w = np.deg2rad(vehicle_geometry["wing_taper"])
    Lambda_c2 = np.deg2rad(vehicle_geometry["Lambda_c2"])
    Lambda_LE = np.deg2rad(vehicle_geometry["Lambda_LE"])
    Lambda_TE = np.deg2rad(vehicle_geometry["Lambda_TE"]) # TODO: check sign
    tc = vehicle_geometry["t/c"]
    S_wfv = vehicle_geometry["surface_area_v"]/ft_m**2
    c_r = vehicle_geometry["root_chord"]/ft_m

    E_ratio = rocket_data["E_ratio_1"]
    F_rock = rocket_data["F_net_SL"] / (g0 * lb_kg)  # lbf
    m_dot_rocket = rocket_data["m_dot_1"]/lb_kg
    OF = rocket_data["O/F"]
    N_eng = rocket_data["N_eng_1"]

    if AB_engine_type != None:
        m_dot_a = AB_data["m_dot_a"]/lb_kg
        m_dot_f = AB_data["m_dot_f_ab"]/lb_kg
        F_ab = AB_data["F_t_ab"] / (g0 * lb_kg)  # lbf
        N_ab_eng = AB_data["N_ab_eng"]
    else:
        m_dot_a = 0
        m_dot_f = 0
        F_ab = 0
        N_ab_eng = 0

    t_tank_ox = tank_geometry["t_tank_ox"]
    t_tank_f = tank_geometry["t_tank_f"]
    t_cyl_ox = tank_geometry["t_cyl_ox"]
    t_cyl_f = tank_geometry["t_cyl_f"]
    D_t_ox_1 = tank_geometry["D_ox_tank"]
    D_t_f_1 = tank_geometry["D_f_tank"]
    l_tank_ox_1 = tank_geometry["l_tank_ox_1"]
    l_cyl_ox_1 = tank_geometry["l_cyl_ox_1"]
    l_tank_f_1 = tank_geometry["l_tank_f_1"]
    l_cyl_f_1 = tank_geometry["l_cyl_f_1"]


    W_ins =  5.29  # Specific weight of tps insulation kg/m2 (Titanium from Philippe's work)
    W_ins = W_ins / 4.88243  # lb/ft2

    # Compute propellant masses
    if AB_engine_type != None:
        W_fuel = t_b_rocket * m_dot_rocket * N_eng / (1 + OF) + t_b_ab * m_dot_f * N_ab_eng
    else:
        W_fuel = t_b_rocket * m_dot_rocket * N_eng / (1 + OF)
    M_fuel = W_fuel * lb_kg
    W_ox = t_b_rocket * m_dot_rocket * OF * N_eng / (1 + OF)
    M_ox = W_ox * lb_kg
    #
    print(M_ox)

    W_prop = W_fuel + W_ox
    M_propellant = M_fuel + M_ox

    # Compute body wetted area = Cylindrical area + nose cone area
    # S_w_fus = (c_r + (c_r - np.tan(Lambda_LE) * D_fus/2 - np.tan(Lambda_TE) * D_fus/2))/2 * D_fus/2
    S_b_tot = (l_fus - l_n) * np.pi * D_fus + 2/3 * np.pi * D_fus * np.sqrt(125/36 * D_fus**2) + np.pi/18 * D_fus**2
    #  + 2 * (S_ref - S_w_fus) + 2 * S_wfv +\
    S_tb = 0.5 * S_b_tot

    ULF = 8.15 # max long acceleration = 100m/s2, lateral acceleration 40m/s2 ULF = 1.5*Limit load 40/g0 * 1.5
    Q_max = 90*10**3 # Max dynamic pressure 57.5 - 90 kPa (Van kesteren)
    Q_max = (Q_max * ft_m**2)/(g0 * lb_kg) # lb/ft2

    # Initiate iteration to find the gross take-off weight based on an initial estimate
    W_gtot = 23*10**3 / lb_kg # lb
    diff = 101
    while diff > 50:
        guess = W_gtot
        # Compute structural weight
        # Body weight
        sigma = (l_fus * ULF/D_fus)**0.15 * Q_max**0.16 * S_b_tot**1.05
        W_b = 0.341 * mf * sigma  # lbs

        # Wing weight
        W_empty = (guess - W_prop) # lbs
        W_w = 0.2958 * mf * ((W_empty * ULF/1000)**0.52 * S_ref**0.7 * A**0.47 * ((1 + lambda_w)/tc)**0.4 *
                             (0.3 + 0.7/np.cos(Lambda_c2)))**1.017  # lbs

        # Vertical tail fin weight
        W_fv = 5 * S_wfv**1.09 # lbs

        # TPS weight
        W_tps = W_ins * (S_tb + S_ref)  # lbs

        # Landing gear weight
        W_lg = 0.00916 * guess**1.124  # lbs

        # Thrust structure weight
        W_thr_r = 0.0025 * F_rock * N_eng   # lbs

        # Propulsion system component weights
        W_rocket = 0.00766 * F_rock * N_eng + 0.00033 * F_rock * N_eng * E_ratio**0.5 + 130 * N_eng
        # TODO: Change this to fixed value of F100-PW-229
        if AB_engine_type == 'turbojet':
            W_ab = N_ab_eng * 3036 # lbs
            #W_ab = N_ab_eng * (m_dot_a * 133.3 - 16600)/4
            W_thr_ab = 0.00625 * F_ab * N_ab_eng + 69

        elif AB_engine_type == "ramjet":
            W_ab = N_ab_eng * 0.01 * F_ab
            W_thr_ab = 0.00625 * F_ab * N_ab_eng + 69

        else:
            W_ab = 0
            W_thr_ab = 0

        W_thr = W_thr_r + W_thr_ab
        W_eng = W_rocket + W_ab

        # Compute tank weight based on shell mass and component correction factor
        V_shell_ox =  4 * np.pi * (D_t_ox_1/2)**2 * t_tank_ox + np.pi * D_t_ox_1 * l_cyl_ox_1 * t_cyl_ox
        V_shell_f = 4/2 * np.pi * (D_t_f_1/2)**2 * t_tank_f + np.pi * D_t_f_1 * l_cyl_f_1 * t_cyl_f
        M_shell_tank_ox = V_shell_ox * rho_t_m
        M_shell_tank_f = V_shell_f * rho_t_m
        M_shell_tank = M_shell_tank_f + M_shell_tank_ox
        K_shell = 1.85  # Shell mass correction factor TRP reader 1.2-1.3 and 2.2-2.5 (Averaged) but 3.5 for composite?
        M_tank = M_shell_tank * K_shell  # kg
        W_tank = M_tank / lb_kg  # lbs


        # Subsystem component weights
        # Hydraulic component weight
        psi = ((S_ref + S_wfv) * Q_max/1000)**0.334 * (l_fus + b)**0.5
        W_hydr = 2.64 * psi
        # Avionics weight
        W_av = 66.37 * guess**0.361
        # Electrical system weight
        theta = guess**0.5 * l_fus**0.25
        W_el = 1.167 * theta

        # Total structural weight
        W_str = W_b + W_w + W_fv + W_tps + W_lg + W_thr  # lbs
        M_str = W_str * lb_kg  # kg

        # Total propulsion weight
        W_prop = W_eng + W_tank  # lbs
        M_prop = W_prop * lb_kg  # kg

        # Total subsystem weight
        W_sub = W_hydr + W_el + W_av  # lbs
        M_sub = W_sub * lb_kg  # kg

        # Vehicle gross weight
        W_upper = M_upper / lb_kg
        W_gtot = W_str + W_prop + W_sub + W_fuel + W_ox + W_upper
        M_gtot = W_gtot * lb_kg

        diff = abs(W_gtot - guess)

    Mass_budget_first = {'Gross_mass': M_gtot, 'Structural_mass': M_str, 'Propulsion_mass': M_prop,
                         'Subsystem_mass': M_sub, 'Propellant_mass': M_propellant, 'Fuel_mass': M_fuel, 'Ox_mass': M_ox}
    return Mass_budget_first