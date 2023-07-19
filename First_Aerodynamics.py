import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator

def fit_function(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

# TODO: Add least square to find optimal starting point
def get_aerodynamic_data(Mach_spacing, AoA_spacing, M, aoa, extrapolation_minimum):
    # Read the data from the text file
    filename = 'trimmed_CL_CD.txt' #'CL_CD_coefficients_Mk-III.txt'
    with open(filename, 'r') as file:
        data = file.readlines()

    # Process the data and extract the coefficients
    angle_of_attack = []
    mach_numbers = []
    lift_coefficients = []
    drag_coefficients = []

    is_angle_of_attack = False

    for line in data[3:]:
        line = line.strip().split('\t')
        if line[0] != '':
            if line[0] == 'AoA [deg]':
                is_angle_of_attack = True
                mach_numbers = list(map(float, line[1:]))
            else:
                try:
                    if is_angle_of_attack:
                        angle_of_attack.append(float(line[0]))
                        drag_coefficients.append(list(map(float, line[1:])))
                    else:
                        lift_coefficients.append(list(map(float, line[1:])))
                except ValueError:
                    continue
        else:
            continue

    # Create numpy arrays for the coefficients
    lift_matrix = np.array(lift_coefficients)
    drag_matrix = np.array(drag_coefficients)

    # Display the matrices
    # print("Matrix for Lift Coefficients (CL):")
    # print(lift_matrix)
    # print("\nMatrix for Drag Coefficients (CD):")
    # print(drag_matrix)
    # print(mach_numbers)
    # print(angle_of_attack)
    mach_numbers = np.array(mach_numbers)
    angle_of_attack = np.array(angle_of_attack)

    lift_interp = RegularGridInterpolator((angle_of_attack, mach_numbers), lift_matrix)
    drag_interp = RegularGridInterpolator((angle_of_attack, mach_numbers), drag_matrix)
    desired_mach_grid, desired_aoa_grid = np.meshgrid(
        np.linspace(min(mach_numbers), max(mach_numbers), Mach_spacing),
        np.linspace(min(angle_of_attack), max(angle_of_attack), AoA_spacing), indexing="ij"
    )

    desired_mach_numbers = np.linspace(min(mach_numbers), max(mach_numbers), Mach_spacing)
    desired_angles_of_attack = np.linspace(min(angle_of_attack), max(angle_of_attack), AoA_spacing)

    interpolated_lift = lift_interp((desired_aoa_grid, desired_mach_grid))
    interpolated_drag = drag_interp((desired_aoa_grid, desired_mach_grid))
    interpolated_lift = interpolated_lift.T
    interpolated_drag = interpolated_drag.T
    # lift_interp = RectBivariateSpline(mach_numbers, angle_of_attack, lift_matrix)
    # drag_interp = RectBivariateSpline(mach_numbers, angle_of_attack, drag_matrix)
    # interpolated_lift = lift_interp.ev(desired_mach_numbers, desired_angles_of_attack)
    # interpolated_drag = drag_interp.ev(desired_mach_numbers, desired_angles_of_attack)


    # print(interpolated_lift)
    # print(interpolated_drag)
    # Extrapolation
    # Loop to cover all rows in the matrix and interpolate all rows

    Mach_numbers_extrap = np.append(desired_mach_numbers, np.linspace(4.0, 6.0, Mach_spacing))
    new_interpolated_lift = np.zeros((len(desired_angles_of_attack), len(Mach_numbers_extrap)))
    new_interpolated_drag = np.zeros((len(desired_angles_of_attack), len(Mach_numbers_extrap)))

    # Find index of Mach number greater than 1.5
    for i in range(len(desired_mach_numbers)):
        if desired_mach_numbers[i] > extrapolation_minimum:
            #print(i, desired_mach_numbers[i])
            index_extrapolation = i
            break

    for i in range(len(desired_angles_of_attack)):
        p_lift, _ = curve_fit(fit_function, desired_mach_numbers[index_extrapolation:], interpolated_lift[i,index_extrapolation:])
        p_drag, _ = curve_fit(fit_function, desired_mach_numbers[index_extrapolation:], interpolated_drag[i,index_extrapolation:])
        Mach_extrapolation = np.linspace(4.0, 6.0, Mach_spacing)
        lift_extrapolation = fit_function(Mach_extrapolation, *p_lift)
        drag_extrapolation = fit_function(Mach_extrapolation, *p_drag)

        new_interpolated_lift[i, :len(desired_mach_numbers)] = interpolated_lift[i, :]
        new_interpolated_lift[i, len(desired_mach_numbers):] = lift_extrapolation
        new_interpolated_drag[i, :len(desired_mach_numbers)] = interpolated_drag[i, :]
        new_interpolated_drag[i, len(desired_mach_numbers):] = drag_extrapolation

    interpolated_lift = new_interpolated_lift
    interpolated_drag = new_interpolated_drag


    # Get index of the Mach number closest to input Mach number
    i = 0
    #print(len(Mach_numbers_extrap))
    while abs(Mach_numbers_extrap[i + 1] - M) <= abs(Mach_numbers_extrap[i] - M):
        i = i + 1
    Mach_number = Mach_numbers_extrap[i]
    # Get index j of the AoA closest to the input AoA
    j = 0
    while abs(desired_angles_of_attack[j + 1] - aoa) < abs(desired_angles_of_attack[j] - aoa):
        j = j + 1
    angle_of_attack = desired_angles_of_attack[j]

    # print(interpolated_lift)
    # print(interpolated_drag)

    CL = interpolated_lift[j][i]
    CD = interpolated_drag[j][i]

    # for k in range(len(desired_angles_of_attack)):
    #     plt.plot(Mach_numbers_extrap, interpolated_lift[k, :], label= f'Lift {desired_angles_of_attack[k]} deg')
    #     plt.xlabel('Mach')
    #     plt.ylabel('CL')
    #     plt.rc('xtick', labelsize=12)
    #     plt.rc('ytick', labelsize=12)
    #     plt.rc('axes', labelsize=12)
    # plt.grid()
    # #plt.legend()
    # plt.show()
    # for k in range(len(desired_angles_of_attack)):
    #     plt.plot(Mach_numbers_extrap, interpolated_drag[k,:], label = f'Drag {desired_angles_of_attack[k]} deg')
    #     plt.xlabel('Mach')
    #     plt.ylabel('CD')
    # plt.grid()
    # #plt.legend()
    # plt.show()
    #
    # for k in range(len(Mach_numbers_extrap)):
    #     plt.plot(desired_angles_of_attack, interpolated_lift[:,k], label = f'CL: Mach {Mach_numbers_extrap[k]}')
    #     plt.xlabel('AoA')
    #     plt.ylabel('CL')
    # plt.grid()
    # #plt.legend()
    # plt.show()
    # for k in range(len(Mach_numbers_extrap)):
    #     plt.plot(desired_angles_of_attack, interpolated_drag[:, k], label = f'CD: Mach {Mach_numbers_extrap[k]}')
    #     plt.xlabel('AoA')
    #     plt.ylabel('CD')
    # plt.grid()
    # #plt.legend()
    # plt.show()


    return  CL, CD #, Mach_number, angle_of_attack,
    # return Mach_number, angle_of_attack, desired_mach_numbers, desired_angles_of_attack #interpolated_lift, interpolated_drag,

# print(get_aerodynamic_data(10, 10, 0.1, 3.6, 2.0))
