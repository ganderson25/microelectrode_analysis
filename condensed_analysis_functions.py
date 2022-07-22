import pandas as pd
import os
import math
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.optimize import least_squares
import scipy.integrate

# Inputs
potential_header = "Ewe/V"
current_header = "<I>/mA"
rh_header = "RH"


# Data extraction
def get_df_from_file(folderpath, n_lines):
    dfs = []
    fileList = os.listdir(folderpath)
    for filename in fileList:
        Ns = 0
        Ewe = 0
        path = os.path.join(folderpath, filename)
        with open(path, "r") as f:
            loop_indices = []
            while True:
                line = f.readline()

                # Collect experimental settings and details for reading the data in the file
                if Ns == False:
                    if "Ns" in line:
                        Ns = line.split()[1:]
                if "Ei (V)" in line:
                    Ewe = line.split()[2:]
                    Ewe = [float(x) for x in Ewe]
                if Ns and Ewe:
                    E_dict = dict(zip(Ns, Ewe))
                if "Number of loops" in line:
                    num_loops = line.split()[-1]
                if "from point number" in line:
                    end_index = line.split()[-1]
                    loop_indices.append(end_index)

                # 'Ewe/V' key to start reading data
                if "Ewe/V" in line:
                    columns = line.split("\t")[:-1]
                    full_df = pd.read_csv(f, delimiter="\t", names=columns)
                    break

            # Add trial number to df
            trial_num = get_trial_list(full_df, num_loops, loop_indices)
            full_df["Trial"] = trial_num

            # Add set applied potential to df
            applied_potential = corresponding_potential(full_df, E_dict)
            full_df["Applied_Potential"] = applied_potential

            # Get averaged values for each applied potential to plot polarization curve
            i = 0
            avg_df = []
            while i < int(num_loops):
                j = 0
                i += 1
                while j < len(Ns):
                    filtered_df = full_df[
                        (full_df["Ns"] == j) & (full_df["Trial"] == i)
                    ]
                    filtered_df = filtered_df.iloc[-n_lines:]
                    avg_filtered_df = filtered_df.mean().to_frame().transpose()
                    avg_df.append(avg_filtered_df)
                    j += 1

            averaged_df = pd.concat(avg_df)

        # Add experimental conditions (such as T, RH, etc) to df
        exp_condition_to_value = experimental_conditions_to_dictionary(filename)
        for condition in exp_condition_to_value.keys():
            averaged_df[condition] = exp_condition_to_value[condition]

        dfs.append(averaged_df)

    df = pd.concat(dfs)
    return df


def get_trial_list(df, num_loops, loop_indices):
    i = 0
    trial_list = []
    while i < len(df.index):
        j = 0
        while j < int(num_loops):
            if i < int(loop_indices[j]):
                trial = 1 + j
                break
            j += 1
        trial_list.append(trial)
        i += 1
    return trial_list


def corresponding_potential(df, dictionary):
    i = 0
    applied_potential = []
    while i < len(df.index):
        Ns_local = str(df.at[i, "Ns"])
        potential = dictionary[Ns_local]
        applied_potential.append(potential)
        i += 1
    return applied_potential


def experimental_conditions_to_dictionary(filename):
    experimental_inputs = extract_experimental_conditions(filename)
    condition_to_value = dict()
    for element in experimental_inputs:
        condition = element.split("=")[0]
        value = element.split("=")[1]
        condition_to_value[condition] = value
    return condition_to_value


def extract_experimental_conditions(filename):
    experimental_inputs = []
    filename_details = filename.split()
    for item in filename_details:
        if "=" in item:
            experimental_inputs.append(item)
    return experimental_inputs


# ECSA Functions
def ecsa_calculation(ecsa_details):
    ecsa_df = read_file_ecsa(ecsa_details["path"], ecsa_details["starting_string"])
    hupd_area = calculate_hupd_area(ecsa_df, ecsa_details)
    ecsa = hupd_area / (ecsa_details["scan_rate"] * ecsa_details["a"])
    return ecsa


def read_file_ecsa(path, ecsa_starting_string):
    file_list = os.listdir(path)
    for file in file_list:
        filepath = os.path.join(path, file)
        with open(filepath, "r") as f:
            columns = []
            while True:
                line = f.readline()
                if ecsa_starting_string in line:
                    columns = line.split("\t")[:-1]
                    print(f"cols in file: {columns}")
                    df = pd.read_csv(f, delimiter="\t", names=columns)
                    break
    return df


def calculate_hupd_area(df, ecsa_details):
    positive_cycle_df = get_positive_half_of_cycle(df, ecsa_details["full_cycle_loop"])
    cap_height = calculate_local_minimum_current(
        positive_cycle_df,
        ecsa_details["min_V_cap_calc"],
        ecsa_details["max_V_cap_calc"],
    )
    integrate_data = get_integrate_data(positive_cycle_df, cap_height)
    hupd_area = get_hupd_area(integrate_data, cap_height)
    return hupd_area


def get_positive_half_of_cycle(df, full_cycle_num):
    cycle_df = df[df["cycle number"] == full_cycle_num]
    cyc_pos_df = cycle_df[cycle_df["<I>/mA"] > 0]
    # sns.lineplot(x='Ewe/V', y='<I>/mA', data = cyc_pos_df, marker='o', ci=None)
    return cyc_pos_df


def calculate_local_minimum_current(df, left_bound, right_bound):
    df_between_peaks = df[(df["Ewe/V"] > left_bound) & (df["Ewe/V"] < right_bound)]
    cap_height = df_between_peaks["<I>/mA"].min()
    return cap_height


def get_integrate_data(df, capacitance_height_current):
    cutoff_V = get_V_for_I_value(df, capacitance_height_current)
    units = get_header_units(current_header)
    integrate_data = transform_to_uA(df[df["Ewe/V"] < cutoff_V], units)
    return integrate_data


def get_V_for_I_value(df, current):
    row_with_current_value = df[df["<I>/mA"] == current]
    V = row_with_current_value["Ewe/V"].iloc[0]
    return V


def get_header_units(header_string):
    loops = len(header_string)
    i = 0
    while i < loops:
        if header_string[i] == "/":
            units = header_string[i + 1 :]
        i += 1
    units = str(units)
    return units


def transform_to_uA(df, units):
    if units == "mA":
        df["<I>/uA"] = df["<I>/mA"] * 1000
    else:
        print("Unit conversion incorrect for given units of current")
    return df


def get_hupd_area(df, capacitance_thickness):
    total_curve_area = trapezoid_integration(df, potential_header, current_header)
    capacitance_area = calculate_capacitance_area(df, capacitance_thickness)
    hupd_area = total_curve_area - capacitance_area
    return hupd_area


def trapezoid_integration(df, x_column_header, y_column_header):
    x = df[x_column_header]
    y = df[y_column_header]
    area = scipy.integrate.trapezoid(y, x)
    return area


def calculate_capacitance_area(df, capacitance_thickness):
    V_range = df[potential_header].max() - df[potential_header].min()
    capacitance_area = V_range * capacitance_thickness
    return capacitance_area


def calculate_relative_roughness(ecsa, diameter):
    geometric_electrode_area = electrode_geometric_area_cm(diameter)
    roughness = ecsa / geometric_electrode_area
    return roughness


def electrode_geometric_area_cm(diameter_microns):
    diameter_cm = diameter_microns / 10000
    geometric_area = calculate_circle_area(diameter_cm)
    return geometric_area


def calculate_circle_area(diameter):
    circle_area = (math.pi * (diameter) ** 2) / 4
    return circle_area


# Water Content Functions
def calculate_water_content_for_rh_list(numpy_rh_series, water_uptake_details):
    water_content_list = []
    for rh in numpy_rh_series:
        water_content = calculate_water_content(water_uptake_details, rh)
        water_content_list.append(water_content)
    return water_content_list

def calculate_water_content(water_uptake_details, rh):
    water_content_function = fit_water_uptake(water_uptake_details)
    water_content = water_content_function(get_float_RH(rh))
    return water_content

def fit_water_uptake(uptake_dict):
    df = get_water_uptake_df(uptake_dict["path"], uptake_dict["membrane"])
    x, y = get_water_uptake_arrays(df)
    fit_parameters = np.polyfit(x, y, uptake_dict["polynomial_degree_fit"])
    fit_function = np.poly1d(fit_parameters)
    return fit_function


def get_water_uptake_df(water_content_path, membrane):
    filenames = os.listdir(water_content_path)
    for file in filenames:
        if membrane in file:
            water_content_path = os.path.join(water_content_path, file)
            df = pd.read_excel(water_content_path)
    return df


def get_water_uptake_arrays(df):
    columns = list(df.columns)
    x = df[columns[0]].to_numpy()
    y = df[columns[1]].to_numpy()
    return x, y


def get_float_RH(rh):
    rh_val = rh.replace("%", "")
    rh_float = float(rh_val)
    return rh_float


# Tafel Calculations
def average_trials_df(df, conditions):
    conditions_keys = list(conditions.keys())
    conditions_keys.append("Applied_Potential")
    tafel_df = df.groupby(conditions_keys, as_index=False).mean()
    return tafel_df


def calculate_logarithm(df_column):
    log_of_column = []
    for item in df_column:
        log_of_column.append(np.log10(np.abs(item)))
    return log_of_column


def sort_df_abs_potential(df):
    df["abs_n/V"] = np.abs(df["Applied_Potential"])
    df = df.sort_values(["abs_n/V"])
    df.reset_index(drop=True, inplace=True)
    return df


def linear_tafel_fit(df, points_fit=4):
    tafel_compare = {
        "r_squared": 0,
        "tafel_slope": 0,
        "I_0": 0,
        "i" : 0
    }
    n_rows = len(df)
    if n_rows < points_fit:
        print("Error num points less than number of points used in tafel fit")
        return None
    else:
        num_it = n_rows - points_fit
        i = 0
        while i <= num_it:
            sub_df = df[i : (i + points_fit)]
            tafel_calculated = calc_tafel_parameters(
                sub_df["log(I)"], sub_df["abs_n/V"], i
            )
            if tafel_calculated["r_squared"] > tafel_compare["r_squared"]:
                tafel_compare = tafel_calculated.copy()
            i += 1
        return tafel_compare


def calc_tafel_parameters(df_column_x, df_column_y, index):
    X = df_column_x.to_numpy().reshape(-1, 1)
    Y = df_column_y
    reg = LinearRegression().fit(X, Y)  # perform linear regression
    y_prediction = reg.predict(X)

    tafel_calculations = {
        "r_squared": reg.score(X, Y),
        "tafel_slope": calculate_tafel_slope(X, y_prediction),
        "I_0": calculate_I_0_tafel(X, y_prediction),
        "i" : index
    }
    return tafel_calculations


def calculate_tafel_slope(X, Y):
    slope = float((Y[2] - Y[0]) / (X[2] - X[0]))  # V/dec
    taf_slope = (np.abs(slope)) * 1000  # mV/dec
    return taf_slope


def calculate_I_0_tafel(X, Y):
    slope = float((Y[2] - Y[0]) / (X[2] - X[0]))  # V/dec
    intercept = Y[0] - slope * X[0]
    log_I_0 = -intercept / slope
    I_0 = 10 ** log_I_0
    return I_0


# Bulter Volmer Calculations
def bv_calculation(bv_inputs, rh):
    res_1 = least_squares(
        bv_residuals,
        bv_inputs["x0"],
        max_nfev=1e10,
        bounds=bv_inputs["bounds"],
        args=(bv_inputs["n"], bv_inputs["i_n"], bv_inputs["C"], bv_inputs["beta"]),
    )

    fit_parameters = res_1.x
    cost_fit = res_1.cost

    bv_fit_parameters = {
        "RH": rh,
        "proton_activity": fit_parameters[0],
        "error": cost_fit,
    }
    return bv_fit_parameters


def bv_residuals(x, n, i_n, C, beta):
    # n = list of overpotentials, i_n = list of current density values, x = empirical fit value of proton activity
    R = 8.31446  # Universal gas constant [J/mol*K]
    F = 96485.3365  # Faraday's constant [C/mol]
    T = 303.15  # Temperature [K]

    j = 0
    res = []
    while j < len(n):
        calc_i = butler_volmer(C, x[0], beta, n[j])
        res.append(i_n[j] - calc_i)
        j += 1
    return res


def butler_volmer(C, a_h, beta, n):
    R = 8.31446  # Universal gas constant [J/mol*K]
    F = 96485.3365  # Faraday's constant [C/mol]
    T = 303.15  # Temperature [K]

    i_0 = C * a_h ** (1 - beta)
    anodic = np.exp(((1 - beta) * F * n) / (R * T))
    cathodic = np.exp((-(1 + beta) * F * n) / (R * T))

    calc_i = i_0 * (anodic - cathodic)

    return calc_i


def calc_i_n_fit(fit_df, rh_series, bv_inputs, n_series):
    i_n_fit_list = []
    i = 0
    while i < len(rh_series):
        mask = fit_df["RH"] == rh_series[i]
        a_h_fit = fit_df[mask]["proton_activity"].item()
        i_n_fit = butler_volmer(bv_inputs["C"], a_h_fit, bv_inputs["beta"], n_series[i])
        i_n_fit_list.append(i_n_fit)
        i += 1
    return i_n_fit_list
