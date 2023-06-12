from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from tc_python import *

"""
This program create a single equilibrium calculation from a "Fe" system, loop every element "Cr", "Co", "Co", "C", "N"
while changing their concentration, then calculate the density and plot the result as a 3D surface.                           
"""


def plot_3d(list_of_x, list_of_y, list_of_z, xlabel, ylabel, zlabel, title):
    """
    Plot a 3d figure using matplotlib given data and labels on the three axes.
    """
    fig = plt.figure()
    fig.suptitle(title, fontsize=14, fontweight='bold')
    ax = fig.gca(projection='3d')
    z = np.empty([len(list_of_x), len(list_of_y)])
    k = 0
    for index_x, x in enumerate(list_of_x):
        for index_y, y in enumerate(list_of_y):
            z[index_x, index_y] = list_of_z[k]
            k = k + 1

    xx, yy = np.meshgrid(list_of_x, list_of_y, indexing='ij')
    ax.plot_surface(xx, yy, z, cmap=cm.coolwarm, linewidth=1, antialiased=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.show()


def list_stable_phases(calc_result):
    """
    List the stable phases and their amount on screen.

    Args:
        calc_result: Calculation result from an equilibrium calculation
        output_string : a string list
    """
    stable_phases = calc_result.get_stable_phases()
    phase_amounts = {} # initialize an empty dictionary to store the phase and corresponding phase amounts
    output_string_list = []
    output_string_list_2 = []
    for phase in stable_phases:
        amount = calc_result.get_value_of('NP(' + phase + ')')
        phase_amounts[phase] = amount
        output_string_list.append(phase + " = {0:.4f}".format(amount))
        output_string_list_2.append("{0:.4f}".format(amount))
    output = ', '.join(output_string_list)
    output_2 = ", ".join(output_string_list_2)
    # print (output)
    return phase_amounts, output, output_2




with TCPython() as start:
    # create and configure a single equilibrium calculation
    calculation = (
        start
            .set_cache_folder(os.path.basename(__file__) + "_cache")
            .select_database_and_elements("TCFE9", ["Fe", "Cr", "Co", "C", "N"])
            .get_system()
            .with_single_equilibrium_calculation()
            .set_condition(ThermodynamicQuantity.temperature(), 300) # room temperature
            .disable_global_minimization()
            # .with_matrix_phase(MatrixPhase("BCC_A2")
            #                    .set_grain_radius(1.e-4)
            #                     .add_precipitate_phase(PrecipitatePhase("M23C6")
            #                                         .set_interfacial_energy(0.252)
            #                                         .set_nucleation_at_grain_boundaries()
            #                                         )
            #                    )
    )
    # in total: 19*11*6 = 1254 combinations
    list_of_x_Cr = np.linspace(10e-2, 145e-3, 19)*100 # (start, stop, number of points). Cr: 19 levels (10-14wt%)
    list_of_x_Co = np.linspace(0, 5e-2, 11)*100 # Co: 11 levels (0-5 wt%)
    list_of_x_C_N = np.linspace(15e-4, 4e-3, 6)*100 # C and N has the same percentage. C/N: 6 levels (0.15-0.4 wt%)
    list_of_density = []
    list_of_index = []
    list_of_element_Cr = []
    list_of_element_Co = []
    list_of_element_C_N = []
    index = 1
    for x_Cr in list_of_x_Cr:
        for x_Co in list_of_x_Co:
            for x_C_N in list_of_x_C_N:
                calc_result = (calculation
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Cr"), x_Cr/100)
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Co"), x_Co/100)
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("C"), x_C_N/100)
                               .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("N"), x_C_N/100)
                            .calculate()
                            )


                mass = calc_result.get_value_of('BM')
                volume = calc_result.get_value_of('VM')
                density = 1e-3 * mass / volume
                list_of_density.append(density)
                list_of_index.append(index)
                list_of_element_Cr.append(x_Cr)
                list_of_element_Co.append(x_Co)
                list_of_element_C_N.append(x_C_N)
                
                output_string_list = []
                output_string_1 = f"Index: {index}" + ", X(Cr)={0:.2f}".format(x_Cr) + " , X(Co)={0:.2f}".format(x_Co) + " , X(C/N)={0:.2f}".format(x_C_N) + ", Density = {0:.4f}".format(density) + "[kg/m3]"  
                output_string_list.append(output_string_1)
                phase_amounts, phase_string, phase_string_2 = list_stable_phases(calc_result)
                output_string_list.append(phase_string)

                output_string = ', '.join(output_string_list)
                print(output_string)

                save_file = open("single_equalibrium.txt", mode = "a")
                save_file.write(output_string + "\n")

                string_2 = f"{index}" + ", {0:.2f}".format(x_Cr) + ", {0:.2f}".format(x_Co) + ", {0:.4f}".format(x_C_N) + ", {0:.4f}".format(density)
                output_string_2 = ", ".join([string_2, phase_string_2])
                save_file_2 = open("single_equalibrium_numerical_result.txt", mode = "a")
                save_file_2.write(output_string_2 + "\n")

                index +=1

    # single_equalibrium_data = np.asarray(list(zip(list_of_index, list_of_element_Cr, list_of_element_Co, list_of_element_C_N, list_of_density)))
    # np.savetxt("single_equalibrium.txt", single_equalibrium_data) #fmt='%.f'  
    # plot_3d(list_of_x_Al, list_of_x_Cr, list_of_density, 'X(Al)', 'X(Cr)', 'Density [kg/m3]',
    #         "Density for Ni-Al-Cr alloy at 800K")
