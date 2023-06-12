from tc_python import *
import matplotlib.pyplot as plt
import numpy as np

"""
This program simulates the kinetics of precipitation of both stable and metastable carbides from ferrite phase.
It demonstrates that metastable carbides (cementite and M7C3) may first emerge and then disappear and the stable phase
(M23C6) prevails.
"""
# ---- Plot result function ------#
def plot_result(x_time, y_number_density, index, precipitate, fig_extension="PNG", resolution=300):
    # save_path = os.path.join(os.path.basename(__file__) + "_cache_", "Figures")
    save_path = os.path.join(".", "preciptation_figures", precipitate)
    path = os.path.join(save_path, str(index) + "_" + precipitate + "." + fig_extension)
    fig, ax = plt.subplots(1)
    fig.suptitle('Precipitation', fontsize=14, fontweight='bold')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Number Density')
    # ax.semilogx(time_1, volume_fraction_1, 'b-', label="Volume fraction of CEMENTITE (Grain boundaries)")
    # ax.semilogx(time_2, volume_fraction_2, 'r-', label="Volume fraction of M7C3 (Grain boundaries)")
    ax.plot(x_time, y_number_density, 'g-', label=("number density of " +  precipitate))
    ax.legend()
    # plt.show()
    plt.savefig(path, format=fig_extension, dpi=resolution)
    plt.close("all")
    # plt.clf()
    # plt.cla()


##############----------------------updated implementation -----------------------------################
with TCPython() as start:
    # create and configure a single equilibrium calculation
    calculation = ( start
                   .set_cache_folder(os.path.basename(__file__) + "_cache")
                   .select_thermodynamic_and_kinetic_databases_with_elements("TCFE9", "MOBFE5", ["Fe", "Cr", "Co", "C", "N"])
                   .get_system()
                   .with_isothermal_precipitation_calculation()
                   .set_composition_unit(CompositionUnit.MASS_FRACTION)
                   .with_numerical_parameters(NumericalParameters()
                                                .set_max_time_step(10)
                                             )
                   
                   .with_matrix_phase(MatrixPhase("BCC_A2")
                                    #  .add_precipitate_phase(PrecipitatePhase("M23C6")
                                    #                         .set_interfacial_energy_estimation_prefactor(1.0)
                                    #                         .set_nucleation_in_bulk()
                                    #                        )
                                     .add_precipitate_phase(PrecipitatePhase("HCP_A3#2")
                                                            .set_interfacial_energy_estimation_prefactor(1.0)
                                                            .set_nucleation_in_bulk()
                                                           )
                                     )
                   .set_temperature(763.15) # 490 degree
                   .set_simulation_time(600)
                #    .calculate()
                   )
    
    # in total: 19*11*6 = 1254 combinations
    list_of_x_Cr = np.linspace(10e-2, 145e-3, 19) # (start, stop, number of points). Cr: 19 levels (10-14 wt%)
    list_of_x_Co = np.linspace(0, 5e-2, 11) # Co: 11 levels (0-5 wt%)
    list_of_x_C_N = np.linspace(5e-4, 2e-3, 6) # C and N has the same percentage. C/N: 6 levels (0.1-0.4 wt%) --> 0.05 - 0.2, 
    # list_of_density_M23C6 = []
    list_of_density_HCP_A3 = []
    index = 1
    for x_Cr in list_of_x_Cr:
        for x_Co in list_of_x_Co:
            for x_C_N in list_of_x_C_N:
                sim_results = (calculation
                            #    .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Cr"), x_Cr)
                            #    .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("Co"), x_Co)
                            #    .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("C"), x_C_N)
                            #    .set_condition(ThermodynamicQuantity.mole_fraction_of_a_component("N"), x_C_N)
                            .set_composition("Cr", x_Cr)
                            .set_composition("Co", x_Co)
                            .set_composition("C", x_C_N)
                            .set_composition("N", x_C_N)
                            .calculate()
                            )
                # time_1, number_density_M23C6 = sim_results.get_number_density_of("M23C6")
                time_2, number_density_HCP_A3 = sim_results.get_number_density_of("HCP_A3#2")
          
                # list_of_density_M23C6.append(number_density_M23C6[-1])
                # list_of_density_HCP_A3.append(number_density_HCP_A3[-1])
                output_string = f"Index: {index}" + ", X(Cr)={0:.4f}".format(x_Cr) + " , X(Co)={0:.4f}".format(x_Co) + " , X(C/N)={0:.6f}".format(x_C_N) +  ", Maximum number density of HCP_A3#2 = {0:.4f}".format(max(number_density_HCP_A3)) + "[m-3]" + ", Precipitation speed of HCP_A3#2 = {0:.4f}".format(max(number_density_HCP_A3)/time_2[number_density_HCP_A3.index(max(number_density_HCP_A3))]) + "[m-3 s-1]"
                # output_string = f"Index: {index}" + ", X(Cr)={0:.4f}".format(x_Cr) + " , X(Co)={0:.4f}".format(x_Co) + " , X(C/N)={0:.6f}".format(x_C_N) + ", Number density of M23C6 at 600s = {0:.4f}".format(number_density_M23C6[-1]) + "[m-3]" +  ", Number density of HCP_A3#2 at 600s = {0:.4f}".format(number_density_HCP_A3[-1]) + "[m-3]"
                print(output_string)
                save_file = open("precipitation_data_HCPA3.txt", mode = "a")
                save_file.write(output_string + "\n")
                string_2 = f"{index}" + ", {0:.4f}".format(x_Cr) + ", {0:.4f}".format(x_Co) + ", {0:.6f}".format(x_C_N) + ", {0:.4f}".format(max(number_density_HCP_A3))+ ", {0:.4f}".format(max(number_density_HCP_A3)/time_2[number_density_HCP_A3.index(max(number_density_HCP_A3))])
                # string_2 = f"{index}" + ", {0:.4f}".format(x_Cr) + ", {0:.4f}".format(x_Co) + ", {0:.6f}".format(x_C_N) + ", {0:.4f}".format(number_density_M23C6[-1]) + ", {0:.4f}".format(number_density_HCP_A3[-1])
                save_file_2 = open("precipitation_data_numerical_HCPA3.txt", mode = "a")
                save_file_2.write(string_2 + "\n")

                # # save the plot
                # plot_result(time_1, number_density_M23C6, index, "M23C6")
                plot_result(time_2, number_density_HCP_A3, index, "HCP_A3")
                index += 1

    # np.savetxt("precipitation_data_2(M23C6).txt", list_of_density_M23C6, fmt='%.4f')  
    # np.savetxt("precipitation_data_HCPA3_2(HCP_A3).txt", list_of_density_M23C6, fmt='%.4f') 

##############----------------------original implementation -----------------------------################
# with TCPython():
#     sim_results = (SetUp()
#                    .set_cache_folder(os.path.basename(__file__) + "_cache")
#                    .select_thermodynamic_and_kinetic_databases_with_elements("FEDEMO", "MFEDEMO", ["Fe", "C", "Cr"])
#                    .get_system()
#                    .with_isothermal_precipitation_calculation()
#                    .set_composition_unit(CompositionUnit.MASS_PERCENT)
#                    .set_composition("C", 0.1)
#                    .set_composition("Cr", 12)
#                    .with_matrix_phase(MatrixPhase("BCC_A2")
#                                      .set_grain_radius(1.e-4)
#                                     #  .add_precipitate_phase(PrecipitatePhase("CEMENTITE")
#                                     #                        .set_interfacial_energy(0.167)
#                                     #                        .set_nucleation_at_grain_boundaries()
#                                     #                         )
#                                     #  .add_precipitate_phase(PrecipitatePhase("M7C3")
#                                     #                        .set_interfacial_energy(0.282)
#                                     #                        .set_nucleation_at_grain_boundaries()
#                                     #                         )
#                                     .add_precipitate_phase(PrecipitatePhase("M23C6")
#                                                             .set_interfacial_energy(0.252)
#                                                             .set_nucleation_at_grain_boundaries()
#                                                            )
#                                      )
#                    .set_temperature(1053)
#                    .set_simulation_time(300)
#                    .calculate()
#                    )

#     # time_1, volume_fraction_1 = sim_results.get_volume_fraction_of("CEMENTITE")
#     # time_2, volume_fraction_2 = sim_results.get_volume_fraction_of("M7C3")
#     time_1, number_density_1 = sim_results.get_number_density_of("M23C6")

