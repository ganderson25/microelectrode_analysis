# Microelectrode Pol Curve Data Analysis from condensed ec-lab data file
# Author: Grace Anderson 6/2/22 ECG Group

from pathlib import Path
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import condensed_analysis_functions as fn
import water_content_functions as wc


def main():
    ## INPUTS ##
    reaction = "HER"
    electrodeDiameter = 50
    path = Path("C:/Users/Grace/Documents/Data/data_processing/solid_state/")
    n_lines = 100
    membrane_identifier = "Electrolyte"
    EW_identifier = "Elec_Conc"

    potential_header = "Ewe/V"
    current_header = "<I>/mA"
    rh_header = "RH"


    experimental_conditions = fn.experimental_conditions_to_dictionary(os.listdir(path / reaction)[0])
    membrane_chemistry = experimental_conditions[membrane_identifier]
    membrane_EW = experimental_conditions[EW_identifier]

    # Dictionaries
    water_uptake_details = {
        "membrane": membrane_chemistry + "_" + str(membrane_EW),
        "path": path / "water_uptake",
        "polynomial_degree_fit": 3,
    }

    # Units: scan rate: V/s, a: uC/cm^2
    ecsa_details = {
        "path": path / "ecsa",
        "starting_string": "Ewe/V",
        "scan_rate": 0.100,
        "a": 210,
        "full_cycle_loop": 2,
        "min_V_cap_calc": 0.2,
        "max_V_cap_calc": 1.0,
    }

    # Program - Data Analysis
    df = fn.get_df_from_file(path / reaction, n_lines)

    geometric_area = fn.electrode_geometric_area_cm(electrodeDiameter)
    df["I_geo_mA/cm^2"] = df["<I>/mA"] / geometric_area

    ecsa = fn.ecsa_calculation(ecsa_details)
    df["I_ecsa_mA/cm^2"] = df["<I>/mA"] / ecsa

    water_content_for_rh = wc.calculate_water_content(
        df[rh_header].to_numpy(), water_uptake_details
    )
    df["Water_Content_mol_H2O/mol_SO3-"] = water_content_for_rh


    ### 1.0 Polarization Curve Plot

    # variable you'd like to plot ('trial', 'temperature', 'date', 'electrolyteState', 'electrolyte', 'concentrationOrEquivalentWeight', 'relativeHumidity', 'gasConcentration', 'reaction', 'electrode')
    pol_curve_plot_variable = "RH"
    palette = sns.color_palette("Blues", df[pol_curve_plot_variable].nunique())
    df.sort_values([pol_curve_plot_variable])
    df = df.reset_index(drop=True)

    fig1, ax1 = plt.subplots()

    sns.lineplot(
        x="Applied_Potential",
        y="I_ecsa_mA/cm^2",
        data=df,
        hue=pol_curve_plot_variable,
        palette=palette,
        marker="o",
        legend="full",
        ci=None,
    )


    ax1.set_ylabel(r"$I \ [mA \ cm^{-2}]$")
    ax1.set_xlabel(r"$Overpotential [V \ vs. \ RHE]$")
    # ax.set_title(r'$75 \ \mu m \ Ir$')
    # fig1.savefig('Pt_HOR_RHs.png', dpi=1000 )

    ###2.0 Polarization Curve with Error Bars

    sns.set_style("ticks")
    fig2, ax2 = plt.subplots()
    sns.lineplot(
        x="Applied_Potential",
        y="I_ecsa_mA/cm^2",
        data=df,
        hue=pol_curve_plot_variable,
        palette=palette,
        marker="o",
        legend="full",
        ci=95,
        err_style="bars",
        err_kws={"capsize": 5},
    )

    ax2.set_ylabel(r"$I \ [mA \ cm^{-2}]$")
    ax2.set_xlabel(r"$Overpotential [V \ vs. \ RHE]$")
    # ax2.set_xticks([-0.02, -0.015, -0.01, -0.005, 0])
    # fig2.savefig("Pt_HER_nafion.png", dpi=1000)

    if reaction == "OER":
        ax2.set_xlim(1.4, 1.65)


    # Tafel Plot
    df["log(I)"] = fn.calculate_logarithm(df[current_header])
    average_df = fn.average_trials_df(df, experimental_conditions)


    fig3, ax3 = plt.subplots()
    sns.lineplot(
        x="log(I)",
        y="Applied_Potential",
        data=average_df,
        hue=pol_curve_plot_variable,
        marker="o",
        palette=palette,
        legend="full",
        ci=None,
    )

    ax3.set_ylabel(r"$Overpotential [V]$")
    ax3.set_xlabel(r"$log |I| $")

    # Linear regression fit to Tafel plot for each rh
    tafel_calcs = []
    bv_fit = []

    for rh in average_df["RH"].unique():
        average_rh_df = fn.sort_df_abs_potential(average_df[average_df["RH"] == rh])
        best_fit_tafel = fn.linear_tafel_fit(average_rh_df)
        best_fit_tafel["RH"] = rh
        best_fit_tafel["Water_Content"] = fn.calculate_water_content(water_uptake_details, rh)
        tafel_calcs.append(pd.DataFrame.from_dict(best_fit_tafel))

    # Semi-Emirical Fitting Bulter-Volmer with Least Squares Fit

        # Experimental Data Arrays
        num_points_fit = 4
        fit_starting_index = best_fit_tafel["i"]
        fit_end_index = fit_starting_index + num_points_fit 

        bv_inputs = {
            "n": average_rh_df["Applied_Potential"].to_numpy()[fit_starting_index:fit_end_index],
            "i_n": average_rh_df["I_geo_mA/cm^2"].to_numpy()[fit_starting_index:fit_end_index],
            # x0 is initial guess for proton activity
            "x0": np.array([0.5]),
            # bounds is bounds for value for proton activity
            "bounds": ([0.1], [100.0]),
            "C": 0.015,
            "beta": 0.7,
        }

        bv_fit_result = fn.bv_calculation(bv_inputs, rh)

        bv_fit.append(pd.DataFrame([bv_fit_result]))

    tafel_fit_results = pd.concat(tafel_calcs)
    bv_least_squares_fit = pd.concat(bv_fit)

    # bv fit visualization
    display_fit_df = df[["Applied_Potential", "I_geo_mA/cm^2", "RH"]]
    display_fit_df["I_n_fit_mA/cm^2"] = fn.calc_i_n_fit(
        bv_least_squares_fit,
        display_fit_df["RH"],
        bv_inputs,
        display_fit_df["Applied_Potential"],
    )

    fig4, ax4 = plt.subplots()
    sns.lineplot(
        x="Applied_Potential",
        y="I_geo_mA/cm^2",
        data=display_fit_df,
        hue=pol_curve_plot_variable,
        linestyle="",
        marker="o",
        legend="full",
        ci=None,
    )
    sns.lineplot(
        x="Applied_Potential",
        y="I_n_fit_mA/cm^2",
        data=display_fit_df,
        hue=pol_curve_plot_variable,
        markers=False,
        legend="full",
        ci=None,
    )
    ax4.set_xlabel(r"$Overpotential [V]$")
    ax4.set_ylabel(r"I mA/cm2")

    # Export processed data
    merge_tafel_result = tafel_fit_results.drop(columns = ["RH"])
    fit_results = pd.concat([merge_tafel_result, bv_least_squares_fit], axis=1)
    processed_filename = reaction + '_' + membrane_chemistry + '_' + membrane_EW + '.xlsx' 
    fit_results.to_excel(path / "processed" / processed_filename, index=False)

    print(bv_fit_result)
    print(bv_least_squares_fit)
    plt.show()

    
    return tafel_fit_results, df
    

tafel, df = main()

print("tada")