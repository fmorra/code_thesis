import os
import numpy as np
import pdb
import re
import pandas as pd


def element_tracker(base_path, area_params):
    files_to_check = []
    headers = ['Label', 'Mises', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23', 'Viscosity']
    cycle_dir_counter = 0
    min_lat = area_params[0]
    max_lat = area_params[1]
    min_lon = area_params[2]
    max_lon = area_params[3]
    depth_range = area_params[4] * 1000
    min_depth_name = str(depth_range[0]/1000)
    max_depth_name = str(depth_range[1]/1000)
    # table_bool = input('Enter 1 to generate the element tracking matrix for te selected latitude, longitude and depth '
    #                    'ranges, 0 to skip this action: ')
    table_bool = 1
    if table_bool == 1:
        if not os.path.exists(os.path.join(base_path, 'Element_tracking_table' + min_depth_name + '_' + max_depth_name +
                                                      '.csv')):
            for iter_dir in os.listdir(base_path):
                iter_path = os.path.join(base_path, iter_dir)
                if os.path.isdir(iter_path) and 'difference' not in iter_path:
                    for step_dir in os.listdir(iter_path):
                        step_path = os.path.join(iter_path, step_dir)
                        for cycle_dir in os.listdir(step_path):
                            cycle_dir_counter = cycle_dir_counter + 1
                            stress_path = os.path.join(step_path, cycle_dir, 'stress_processing_results',
                                                       'Complete_Files')
                            for f in os.listdir(stress_path):
                                if f == 'Complete_file_EARTH.csv':
                                    files_to_check.append(os.path.join(stress_path, f))
                                elif f == 'Complete_file_coupled_EARTH.csv':
                                    files_to_check.append(os.path.join(stress_path, f))
                                else:
                                    pass
                                # if os.path.isdir(os.path.join(step_path, cycle_dir, 'stress_matrices', f)):
                                #     coupled_path = os.path.join(step_path, cycle_dir, 'stress_matrices', f)
                                #     for cf in os.listdir(coupled_path):
                                #         if cf == 'EARTH.csv':
                                #             files_to_check.append(os.path.join(coupled_path, cf)
            element_table = np.zeros((len(files_to_check), 9))
            iter_numbers = np.zeros((len(files_to_check), 3))
            count = 0
            for name in files_to_check:
                e_path = os.path.join(base_path, 'e.dat')
                with open(name, 'r') as no_arr_matrix:
                    large_matrix = no_arr_matrix.readlines()
                    large_matrix = [line.strip() for line in large_matrix[1:]]
                    large_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in large_matrix]
                    large_matrix = np.array(large_matrix)
                    # Add here the part where only the points in a certain depth and latlon range are considered
                    # and saved in a matrix called matrix
                    depth = large_matrix[:, -3]
                    lat = large_matrix[:, -2]
                    lon = large_matrix[:, -1]
                    min_depth = depth_range[0]
                    max_depth = depth_range[1]
                    depth_condition = (depth > min_depth) & (depth < max_depth)
                    lat_condition = (lat > min_lat) & (lat < max_lat)
                    lon_condition = (lon > min_lon) & (lon < max_lon)
                    matrix = large_matrix[depth_condition * lat_condition * lon_condition]
                    mises = matrix[:, 1]
                    max_stress = max(mises)
                    max_index = np.where(mises == max_stress)
                    mises = matrix[max_index, 1]
                    with open(e_path, 'r') as no_arr_e:
                        e = no_arr_e.readlines()
                        e = [line.strip() for line in e]
                        e = [np.array([float(i) for i in line.split(" ")[:]]) for line in e]
                        e = np.array(e)
                        b_diff = e[:, 1]
                        b_disl = e[:, 2]
                    b_index = int(matrix[max_index, 0]) - 1
                    viscosity = mises/(3*(b_diff[b_index] * mises + b_disl[b_index] * mises ** 3.5))

                    stresses = matrix[max_index]
                    stresses = stresses[:, 0:8]
                    value_vector = np.hstack((stresses, viscosity))
                    element_table[count, :] = value_vector
                data = [float(s) for s in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", name)]
                data = np.array([float(x) for x in data])
                data = data[1:]
                iter_numbers[count, :] = data
                count = count + 1
            if len(element_table) % 2 != 0:
                for i in range(0, len(element_table)-1, 2):
                    element_table[[i, i+1]] = element_table[[i+1, i]]
            else:
                for i in range(0, len(element_table), 2):
                    element_table[[i, i+1]] = element_table[[i+1, i]]
            indices = []
            for i in range(len(iter_numbers)):
                if i % 2 != 0:
                    indices.append(str(int(iter_numbers[i, 0])) + ' ' + str(int(iter_numbers[i, 1])) + ' ' + str(
                        int(iter_numbers[i, 2])) + ' t.m.')
                else:
                    indices.append(str(int(iter_numbers[i, 0])) + ' ' + str(int(iter_numbers[i, 1])) + ' ' + str(
                        int(iter_numbers[i, 2])) + ' f.m.')
            table = pd.DataFrame(data=element_table, index=indices, columns=headers)
            # table.round(decimals=1)
            table.to_csv(os.path.join(base_path, 'Element_tracking_table_' + str(min_lat) + '_' + str(max_lat) +
                                      str(min_lon) + '_' + str(max_lon) + '_' +
                                      min_depth_name + '_' + max_depth_name +
                                      '.csv'), float_format='%.1f')


