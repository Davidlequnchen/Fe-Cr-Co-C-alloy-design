import numpy as np
import matplotlib.pyplot as plt
from tc_python import *

"""
Shows the basic usage of Scheil-calculations in TC-Python and mixing them with equilibrium calculations.
"""



def hot_cracking_susceptibility(x_mole_fraction_list, y_temperature_list):
    '''
    calculate the HCS based on scheil curve. 
    x_mole_fraction_list: mole fraction of solid
    y_temperature_list: the temperature values from scheil calculation
    '''
    # interperate y value from x value: y_interp np.interp(x_val, x, y)
    y_interp_09 = round(np.interp(0.9, x_mole_fraction_list, y_temperature_list), 4)
    y_interp_10 = round(np.interp(1.0, x_mole_fraction_list, y_temperature_list), 4)
    y_interp_04 = round(np.interp(0.4, x_mole_fraction_list, y_temperature_list), 4)
    y_interp_00 = round(np.interp(0, x_mole_fraction_list, y_temperature_list), 4)

    t_09 = (y_interp_00 - y_interp_09)/100
    t_10 = (y_interp_00 - y_interp_10)/100
    t_04 = (y_interp_00 - y_interp_04)/100 
    hcs = (t_09 - t_10)/(t_04 - t_09)
    return hcs

def growth_restriction_factor(x_mole_fraction_list, y_temperature_list):
    """
    calculate the GRF based on scheil curve. The GRF is the slope of the curve approach to zero
    """
    y_interp_005 = round(np.interp(0.05, x_mole_fraction_list, y_temperature_list), 4)
    y_interp_00 = round(np.interp(0, x_mole_fraction_list, y_temperature_list), 4)
    grf = (y_interp_005 - y_interp_00)/0.05
    return abs(grf)


def plot_scheil_curve(scheil_curve, index, fig_extension="PNG", resolution=150):
    save_path = os.path.join(".", "scheil_curve_figures")
    path = os.path.join(save_path, str(index) + "." + fig_extension)

    fig, ax = plt.subplots(1)
    for label in scheil_curve:
        section = scheil_curve[label]
        x = section.x  
        y = np.array(section.y) - 273.15
        ax.plot(x, y,  label=label)

    ax.set_xlabel("Mole fraction of all solid phases [-]")
    ax.set_ylabel("Temperature [\N{DEGREE SIGN} C]")

    ax.legend(loc="lower left")
    ax.set_title("Solidification of Cr23C6")

    plt.savefig(path, format=fig_extension, dpi=resolution)
    plt.close("all")


##############----------------------updated implementation -----------------------------################
database = "TCFE9" 
dependent_element = "Fe" 
elements = ["Fe", "Cr", "Co", "C", "N"]

with TCPython() as session:
    system = (session.
              set_cache_folder(os.path.basename(__file__) + "_cache").
              select_database_and_elements(database, [dependent_element] + elements).
              get_system_for_scheil_calculations())

    scheil_calculation = (system.
                          with_scheil_calculation().
                          set_composition_unit(CompositionUnit.MASS_PERCENT)
                        #   .set_composition("Cr", 14)
                        #   .set_composition("Co", 2)
                        #   .set_composition("C", 0.2)
                        #   .set_composition("N", 0.2)
                        )


    # in total: 19*11*6 = 1254 combinations
    list_of_x_Cr = np.linspace(10e-2, 145e-3, 19)*100 # (start, stop, number of points). Cr: 19 levels (10-14wt%)
    list_of_x_Co = np.linspace(0, 5e-2, 11)*100 # Co: 11 levels (0-5 wt%)
    list_of_x_C_N = np.linspace(15e-4, 4e-3, 6)*100 # C and N has the same percentage. C/N: 6 levels (0.15-0.4 wt%)
    list_of_time = []
    list_of_density = []
    index = 1
    for x_Cr in list_of_x_Cr:
        for x_Co in list_of_x_Co:
            for x_C_N in list_of_x_C_N:
                solidification_results = (scheil_calculation
                                        .set_composition("Cr", x_Cr)
                                        .set_composition("Co", x_Co)
                                        .set_composition("C", x_C_N)
                                        .set_composition("N", x_C_N)
                                        .calculate()
                                         )
                # --- solidification curve (mole fraction solid phases vs. T) including the equilibrium ------
                scheil_curve = solidification_results.get_values_grouped_by_stable_phases_of(
                    ScheilQuantity.mole_fraction_of_all_solid_phases(),
                    ScheilQuantity.temperature())

                temp_min = 1e6
                temp_max = -1e6
                x_mole_fraction_list = []
                y_temperature_list = []
                for label in scheil_curve:
                    section = scheil_curve[label]
                    temp_min = min(np.min(section.y), temp_min)
                    temp_max = max(np.max(section.y), temp_max)

                    x_mole_fraction_list.extend(section.x)
                    y_temperature_list.extend(section.y)
                
                # plot and save figure
                plot_scheil_curve(scheil_curve, index)

                hcs = hot_cracking_susceptibility(x_mole_fraction_list, y_temperature_list )
                grf = growth_restriction_factor(x_mole_fraction_list, y_temperature_list)
                output_string =  f"Index: {index}" + ", X(Cr)={0:.2f}".format(x_Cr) + " , X(Co)={0:.2f}".format(x_Co) + " , X(C/N)={0:.4f}".format(x_C_N) + " , Hot cracking susceptibility (HCS) = {0:.4f}".format(hcs) + ", Growth restriction factor (GRF)= {0:.4f}".format(grf) 
                print(output_string)
                save_file = open("scheil_curve_calculation.txt", mode = "a")
                save_file.write(output_string + "\n")

                string_2 = f"{index}" + ", {0:.2f}".format(x_Cr) + ", {0:.2f}".format(x_Co) + ", {0:.4f}".format(x_C_N) + ", {0:.4f}".format(hcs) + ", {0:.4f}".format(grf) 
                save_file_2 = open("scheil_curve_calculation_numerical_results.txt", mode = "a")
                save_file_2.write(string_2 + "\n")
                index += 1



##############----------------------original implementation -----------------------------################
'''
with TCPython() as session:
    system = (session.
              set_cache_folder(os.path.basename(__file__) + "_cache").
              select_database_and_elements(database, [dependent_element] + elements).
              get_system_for_scheil_calculations())

    scheil_calculation = (system.
                          with_scheil_calculation().
                          set_composition_unit(CompositionUnit.MASS_PERCENT).
                          set_composition("Si", wt_pct_si))

    solidification = scheil_calculation.calculate()

    # 1. Plot the solidification curve (mole fraction solid phases vs. T) including the equilibrium
    scheil_curve = solidification.get_values_grouped_by_stable_phases_of(
        ScheilQuantity.mole_fraction_of_all_solid_phases(),
        ScheilQuantity.temperature())

    temp_min = 1e6
    temp_max = -1e6
    fig, (ax_1, ax_2) = plt.subplots(2)
    for label in scheil_curve:
        section = scheil_curve[label]
        temp_min = min(np.min(section.y), temp_min)
        temp_max = max(np.max(section.y), temp_max)
        ax_1.plot(section.x, np.array(section.y) - 273.15, label=label)

    # calculate the equilibrium solidification line (starting at the liquidus temperature)
    prop_calculation = system.with_property_diagram_calculation()
    result = (prop_calculation.
              with_axis(CalculationAxis(ThermodynamicQuantity.temperature()).
                        set_min(temp_min).
                        set_max(temp_max)).
              set_condition("W(Si)", wt_pct_si / 100).
              calculate())

    temp_eq_frac, liquid_eq_frac = result.get_values_of(ThermodynamicQuantity.temperature(),
                                                        ThermodynamicQuantity.mole_fraction_of_a_phase("LIQUID"))
    solid_eq_frac = 1 - np.array(liquid_eq_frac)
    temp_eq_frac = np.array(temp_eq_frac) - 273.15
    valid_indices = np.where(solid_eq_frac < 1.0)  # cutting off all values with 100% solid

    ax_1.plot(solid_eq_frac[valid_indices], temp_eq_frac[valid_indices], '--', label="Equilibrium")
    ax_1.set_xlabel("Mole fraction of all solid phases [-]")
    ax_1.set_ylabel("Temperature [\N{DEGREE SIGN} C]")

    ax_1.legend(loc="lower left")
    ax_1.set_title("Solidification of AlSi1")

    # 2. Plot the mole fraction of the solid phases separated for each
    groups = \
        solidification.get_values_grouped_by_stable_phases_of(ScheilQuantity.mole_fraction_of_a_solid_phase(ALL_PHASES),
                                                              ScheilQuantity.temperature())

    for group in groups.values():
        ax_2.plot(group.x, np.array(group.y) - 273.15, label=group.label)

    ax_2.set_xlabel("Mole fraction of each solid phase [-]")
    ax_2.set_ylabel("Temperature  [\N{DEGREE SIGN} C]")
    ax_2.legend(loc="lower center")

    fig.set_size_inches(5, 7)
    plt.show()
'''