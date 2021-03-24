import numpy as np
import os
import pandas as pd
import pdb


def plot_scale_extremes(sd_input, min_lat, max_lat, min_lon, max_lon, depths_to_plot, selected_components,
                        viscosity_input, b_input, python_base_path, run_vec, coordinate_system_input):
    runs = run_vec
    extremes_matrix = np.zeros((len(runs), 2 * len(selected_components)))
    new_scale_limits = np.zeros((1, len(extremes_matrix[0])))
    for run_index in range(len(runs)):
        run_folder = os.path.join(python_base_path, 'run_' + str(runs[run_index]))
        file_list = []
        for root, dirs, files in os.walk(run_folder, topdown=False):
            for name in files:
                if coordinate_system_input == 0:
                    if name == "Complete_file_EARTH.csv":
                        file_list.append(os.path.join(root, name))
                else:
                    if name == "Geographical_complete_file_EARTH.csv":
                        file_list.append(os.path.join(root, name))
        full_stress_files = []
        full_defo_files = []
        for i in range(len(file_list)):
            if 'stress' in file_list[i]:
                full_stress_files.append(file_list[i])
            else:
                full_defo_files.append(file_list[i])
        scale_limits = np.array([])
        min_depth = depths_to_plot[0]
        max_depth = depths_to_plot[1]
        if sd_input == 0 and viscosity_input == 0:
            stress_extremes_matrix = np.zeros((len(full_stress_files), len(selected_components) * 2))
            for j in range(len(full_stress_files)):
                stress_matrix = pd.read_csv(full_stress_files[j], delimiter=",", skiprows=1)
                selected_columns = []
                for k in range(len(selected_components)):
                    if selected_components[k] == 'Mises':
                        selected_columns.append(1)
                    elif selected_components[k] == 'S11':
                        selected_columns.append(2)
                    elif selected_components[k] == 'S22':
                        selected_columns.append(3)
                    elif selected_components[k] == 'S33':
                        selected_columns.append(4)
                    elif selected_components[k] == 'S12':
                        selected_columns.append(5)
                    elif selected_components[k] == 'S13':
                        selected_columns.append(6)
                    elif selected_components[k] == 'S23':
                        selected_columns.append(7)
                    else:
                        pass
                depth = stress_matrix.iloc[:, -3] / 1e3
                lat = stress_matrix.iloc[:, -2]
                lon = stress_matrix.iloc[:, -1]
                depth_condition = (depth > min_depth) & (depth < max_depth)
                lat_condition = (lat > min_lat) & (lat < max_lat)
                lon_condition = (lon > min_lon) & (lon < max_lon)
                condition_matrix = stress_matrix[depth_condition & lat_condition & lon_condition]
                for column in range(len(selected_columns)):
                    variable = condition_matrix.iloc[:, selected_columns[column]]
                    variable_extremes = [np.min(variable), np.max(variable)]
                    stress_extremes_matrix[j, 2*column:2*column+2] = variable_extremes
            minima = np.min(stress_extremes_matrix, axis=0)
            maxima = np.max(stress_extremes_matrix, axis=0)
            scale_limits = np.hstack((minima[0::2], maxima[1::2]))
        elif sd_input == 1 and viscosity_input == 0:
            defo_extremes_matrix = np.zeros((len(full_defo_files), len(selected_components) * 2))
            for k in range(len(full_defo_files)):
                defo_matrix = pd.read_csv(full_stress_files[k], delimiter=",", skiprows=1)
                selected_columns = []
                for l in range(len(selected_components)):
                    if selected_components[l] == 'Magnitude':
                        selected_columns.append(1)
                    elif selected_components[l] == 'U1':
                        selected_columns.append(2)
                    elif selected_components[l] == 'U2':
                        selected_columns.append(3)
                    elif selected_components[l] == 'U3':
                        selected_columns.append(4)
                    else:
                        pass
                depth = defo_matrix.iloc[:, -3] / 1e3
                lat = defo_matrix.iloc[:, -2]
                lon = defo_matrix.iloc[:, -1]
                depth_condition = (depth > min_depth) & (depth < max_depth)
                lat_condition = (lat > min_lat) & (lat < max_lat)
                lon_condition = (lon > min_lon) & (lon < max_lon)
                condition_matrix = defo_matrix[depth_condition & lat_condition & lon_condition]
                for column in range(len(selected_columns)):
                    variable = condition_matrix.iloc[:, selected_columns[column]]
                    variable_extremes = [np.min(variable), np.max(variable)]
                    defo_extremes_matrix[k, 2*column:2*column+2] = variable_extremes
            minima = np.min(defo_extremes_matrix, axis=0)
            maxima = np.max(defo_extremes_matrix, axis=0)
            scale_limits = np.hstack((minima[0::2], maxima[1::2]))
        elif viscosity_input == 1:
            e_path = os.path.join(run_folder, 'e.dat')
            e = pd.read_csv(e_path, delimiter=",", skiprows=1)
            elements = e.iloc[:, 0]
            alin = e.iloc[:, 1]
            a = e.iloc[:, 2]
            an = 3.5
            strain_rate = np.zeros((len(e), 1))
            exponent = an
            visco_extremes_matrix = np.zeros((len(full_stress_files), 2))
            for l in range(len(full_stress_files)):
                opened_viscosity_matrix = pd.read_csv(full_stress_files[l], delimiter=",", skiprows=1)
                mises = opened_viscosity_matrix.iloc[:, 1]
                depth = opened_viscosity_matrix[:, -3] / 1e3
                lat = opened_viscosity_matrix[:, -2]
                lon = opened_viscosity_matrix[:, -1]
                depth_condition = depth > min_depth & depth < max_depth
                lat_condition = lat > min_lat & lat < max_lat
                lon_condition = lon > min_lon & lon < max_lon
                if b_input == 0:
                    for j in range(len(strain_rate)):
                        strain_rate[j, 1] = alin[j] * mises[j] + a[j] * mises[j] ^ exponent
                    viscosity = mises / (3 * strain_rate)
                    viscosity_matrix = np.zeros((len(opened_viscosity_matrix), 4))
                    viscosity_matrix[:, 0] = elements
                    viscosity_matrix[:, -3] = viscosity
                    viscosity_matrix[:, -2] = strain_rate
                    viscosity_matrix[:, -1] = mises
                    complete_matrix_with_viscosity = np.hstack((opened_viscosity_matrix, viscosity_matrix))
                    condition_matrix = complete_matrix_with_viscosity[depth_condition & lat_condition & lon_condition]
                    variable = condition_matrix[:, -3]
                    variable_extremes = [np.min(variable), np.max(variable)]
                    visco_extremes_matrix[l, 1: 2] = variable_extremes
                    scale_limits = [np.min(visco_extremes_matrix[:, 1]), np.max(visco_extremes_matrix[:, 2])]
                else:
                    b_matrix = np.zeros((len(full_stress_files), 2))
                    for m in range(len(full_stress_files)):
                        complete_b_matrix = np.hstack((opened_viscosity_matrix, alin, a))
                        data_points_indices = complete_b_matrix[depth_condition & lat_condition & lon_condition]
                        selected_columns = np.hstack((alin[data_points_indices, :], a[data_points_indices, :]))
                        for n in range(len(selected_columns[0])):
                            variable = selected_columns[:, n]
                            variable_extremes = [np.min(variable), np.max(variable)]
                            b_matrix[n, 2 * m - 1: 2 * n] = variable_extremes
                        scale_limits = [np.min(b_matrix[:, 1: 2:len(b_matrix[0]) - 1]),
                                        np.max(b_matrix[:, 2: 2:len(b_matrix[0])])]
        else:
            pass
        extremes_matrix[run_index, :] = scale_limits
        run_minima = np.min(extremes_matrix, axis=0)
        run_maxima = np.max(extremes_matrix, axis=0)
        matrix_half = int(len(extremes_matrix[0])/2)
        new_scale_limits = np.hstack((run_minima[0:matrix_half], run_maxima[matrix_half:]))
    return new_scale_limits
