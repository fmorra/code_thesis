from natsort import natsorted
from python_discretizer import *
import os
import numpy as np
import re
import matplotlib.pyplot as plt
import csv
import pdb
import cartopy.crs as ccrs


def depth_classifier(sd_input, deflection_processing_path, individual_path, components_to_plot, headers_on,
                     file_extension, cartesian_classified_depth_path, geographical_classified_depth_path,
                     files_to_classify_path, files_to_classify, histogram_path):

    # Here the data is discretized in a certain number of bins, and stored in different matrices for each of those bins.
    # Distinguish two cases, one for the stresses and one for the deflections, and define relevant variables
    if sd_input == 0:
        upper_path = individual_path
        save_name = 'Stresses_until_depth_'
        histogram_name_part = 'centroids'
        matrix_to_save_part = 'Stresses_layer_'
        earth_name_part = 'Stresses_Earth_'
        discretized_headers = ['Label', 'Mises', 'S_11', 'S_22', 'S_33', 'S_12', 'S_13', 'S_23', 'X', 'Y', 'Z', 'R',
                               'Depth', 'Lat', 'Lon']
    else:
        upper_path = deflection_processing_path
        save_name = 'Deflections_until_depth_'
        histogram_name_part = 'points'
        matrix_to_save_part = 'Deflections_layer_'
        earth_name_part = 'Deflections_Earth_'
        discretized_headers = ['Label', 'U_Magn', 'U_1', 'U_2', 'U_3', 'X', 'Y', 'Z', 'R', 'Depth', 'Lat', 'Lon']

    # Extract filenames and create folders based on the components we are working with
    if components_to_plot == 'geographical':
        classified_path = geographical_classified_depth_path
        if not os.path.exists(classified_path):
            os.mkdir(classified_path)
        if not os.path.exists(histogram_path):
            os.mkdir(histogram_path)
    elif components_to_plot == 'cartesian':
        classified_path = cartesian_classified_depth_path
        if not os.path.exists(classified_path):
            os.mkdir(classified_path)
        if not os.path.exists(histogram_path):
            os.mkdir(histogram_path)
    else:
        print('Incorrect input, select either geographical or cartesian.')
    # Only leave EARTH_POINT and EARTH_POINT_LOW as files to bin, because those are the ones that interest us. Do not
    # consider the region files as they have been all saved in a single Earth part file.

    for file_to_evaluate in range(len(files_to_classify)):
        if 'I0' in files_to_classify[file_to_evaluate]:
            files_to_classify[file_to_evaluate] = 'to_delete'
        if 'Large' in files_to_classify[file_to_evaluate]:
            files_to_classify[file_to_evaluate] = 'to_delete'
        if 'Region' in files_to_classify[file_to_evaluate]:
            files_to_classify[file_to_evaluate] = 'to_delete'
        if 'LOW' in files_to_classify[file_to_evaluate]:
            files_to_classify[file_to_evaluate] = 'to_delete'
        if 'coupled' in files_to_classify[file_to_evaluate]:
            files_to_classify[file_to_evaluate] = 'to_delete'
    files_to_classify = [a for a in files_to_classify if a != 'to_delete']
    for file_to_evaluate in range(len(files_to_classify)):
        with open(os.path.join(files_to_classify_path, files_to_classify[file_to_evaluate])) as matrix:
            matrix_to_evaluate = matrix.readlines()
            matrix_to_evaluate = [line.strip() for line in matrix_to_evaluate[1:]]
            matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in matrix_to_evaluate]
            matrix_to_evaluate = np.array(matrix_to_evaluate)
            if not matrix_to_evaluate[1:, 1:].any():
                files_to_classify[file_to_evaluate] = 'to_delete'
    files_to_classify = [a for a in files_to_classify if a != 'to_delete']
    maximum_depth = 0
    earth_data = []
    for i in range(len(files_to_classify)):
        with open(os.path.join(files_to_classify_path, files_to_classify[i])) as matrix:
            data_matrix = matrix.readlines()
            data_matrix = [line.strip() for line in data_matrix[1:]]
            data_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in data_matrix]
            depth_matrix = np.array(data_matrix)
            maximum_depth = np.amax(depth_matrix[:, -3])
            earth_data.append(depth_matrix)
    earth_data = np.concatenate(earth_data, axis=0)
    if maximum_depth > 2.886e6:
        layer_depth_km = [410, 660, 2886, maximum_depth / 1000]
    else:
        layer_depth_km = [410, 660, maximum_depth / 1000]
    layer_depth_km = np.array(layer_depth_km)
    layer_depth = layer_depth_km * 1000
    check_1 = 0
    check_2 = 1
    run_input = 1
    if os.path.exists(os.path.join(histogram_path, 'Complete_Earth_distribution.png')):
        while check_1 == 1:
            run_input = eval(input('Do you want to rerun the depth discretization algorithm? 1 if yes, 0 if no: \n'))
            if run_input == 0:
                check_1 = 0
            elif run_input == 1:
                while check_2 == 1:
                    bin_input = eval(input('Do you want to change the bin number? 1 if yes, 0 if no: \n'))
                    if bin_input == 0:
                        check_2 = 0
                        check_1 = 0
                    elif bin_input == 1:
                        n_bins = eval(input('Enter number of layer bins:  \n'))
                        earth_bins = eval(input('Enter number of model bins: \n'))
                        check_2 = 0
                        check_1 = 0
                    else:
                        print('Incorrect input, select either 1 or 0. \n')
            else:
                print('incorrect input, select either 1 or 0. \n')
    fig_counter = 1

    if run_input == 0:
        print('Skipping this stage to generate MATLAB handles file.')
    else:
        print('Starting discretization of data into different depths')
        large_depth_matrix = []
        for i in range(len(files_to_classify)):
            with open(os.path.join(files_to_classify_path, files_to_classify[i])) as matrix:
                data_matrix = matrix.readlines()
                data_matrix = [line.strip() for line in data_matrix[1:]]
                data_matrix = [np.array([eval(j) for j in line.split(",")[:]]) for line in data_matrix]
                complete_matrix = np.array(data_matrix)
                large_depth_matrix.append(complete_matrix)
        large_depth_matrix = np.concatenate(large_depth_matrix, axis=0)
        for depth_element in range(len(layer_depth_km)):
            layer_depth_km[depth_element] = float(layer_depth_km[depth_element])

        for i in range(len(layer_depth)):
            layer_group_indices = [j for j, v in enumerate(large_depth_matrix[:, -3]) if v < layer_depth[i]]
            layer_matrix = large_depth_matrix[layer_group_indices, :]
            layer_depth_km_str = str(int(round(layer_depth_km[i])))
            if headers_on == 1:
                with open(os.path.join(classified_path, save_name + layer_depth_km_str + '_km.csv'), 'wb') as \
                        f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(discretized_headers)
                    writer.writerows(layer_matrix)
            else:
                with open(os.path.join(classified_path, save_name + layer_depth_km_str + '_km.csv'), 'wb') as \
                        f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(layer_matrix)
            large_depth_matrix[layer_group_indices, :] = 0
            large_depth_matrix = large_depth_matrix[~np.all(large_depth_matrix == 0, axis=1)]
        n_bins = 30
        layer_values = os.listdir(classified_path)
        layer_values_names = [filename for filename in layer_values if filename.endswith(file_extension)]
        layer_values_names = natsorted(layer_values_names)
        fig_counter = 1
        for i in range(len(layer_values_names)):
            subclassified_path = os.path.join(classified_path, 'depth_' +
                                              str(re.search(r'\d+', layer_values_names[i]).group(0)) + '_km')
            if not os.path.exists(subclassified_path):
                os.mkdir(subclassified_path)
            with open(os.path.join(classified_path, layer_values_names[i])) as layer_matrix_to_discretize:
                layer_data_matrix = layer_matrix_to_discretize.readlines()
                layer_data_matrix = [line.strip() for line in layer_data_matrix[1:]]
                layer_data_matrix = [np.array([eval(j) for j in line.split(",")[:]]) for line in layer_data_matrix]
                layer_data_matrix = np.array(layer_data_matrix)
                depth_data = layer_data_matrix[:, -3]/1000
                plt.figure(fig_counter)
                plt.hist(depth_data, n_bins, edgecolor='k', linewidth=1)
                plt.xlabel('Depth [km]')
                plt.ylabel('Number of points')
                plt.grid(axis='both', alpha=0.75)
                if i == 0:
                    plt.title('Distribution of ' + str(histogram_name_part) +
                              ' from the ABAQUS model at different depths for a depth range of 0 to ' +
                              list(filter(str.isdigit, layer_values_names[i])) + ' km')
                    plt.savefig(os.path.join(histogram_path, 'Depth_range_0_' +
                                             list(filter(str.isdigit, layer_values_names[i]))
                                             + '_km_distribution.png'), bbox_inches='tight')
                    fig_counter += 1
                else:
                    plt.title('Distribution of ' + str(histogram_name_part) +
                              ' from the ABAQUS model at different depths for a depth range of ' +
                              list(filter(str.isdigit, layer_values_names[i - 1])) + ' to ' +
                              list(filter(str.isdigit, layer_values_names[i])) + ' km')
                    plt.savefig(os.path.join(histogram_path, 'Depth_range_' +
                                             list(filter(str.isdigit, layer_values_names[i - 1])) + '-' +
                                             list(filter(str.isdigit, layer_values_names[i])) +
                                             '_km_distribution.png'), bbox_inches='tight')
                    fig_counter += 1
                indices, bin_edges = python_discretizer(depth_data, i, layer_depth, n_bins)
                for j in range(n_bins):
                    subdivision_matrix = layer_data_matrix[indices == j+1, :]
                    smaller_edge = str(int(round(bin_edges[j] / 1000)))
                    if j == n_bins-1:
                        larger_edge = str(int(round(bin_edges[-1] / 1000)))
                    else:
                        larger_edge = str(int(round(bin_edges[j+1] / 1000)))
                    if headers_on == 1:
                        with open(os.path.join(subclassified_path, matrix_to_save_part +
                                        str(re.search(r'\d+', layer_values_names[i]).group(0)) + '_km_bin_' + str(j+1) +
                                        '_edges_' + smaller_edge + '_' + larger_edge + '.csv'), 'wb') as f_write:
                            writer = csv.writer(f_write)
                            writer.writerow(discretized_headers)
                            writer.writerows(subdivision_matrix)
                    else:
                        with open(os.path.join(subclassified_path, matrix_to_save_part +
                                        str(re.search(r'\d+', layer_values_names[i]).group(0)) + '_km_bin_' + str(j+1) +
                                        '_edges_' + smaller_edge + '_' + larger_edge + '.csv'), 'wb') as f_write:
                            writer = csv.writer(f_write)
                            writer.writerows(subdivision_matrix)

    #################################################################################################################

        earth_bins = 50
        data = earth_data[:, -3]
        bins = np.arange(data.min(), data.max()+(data.max()-data.min())/(2*earth_bins), (data.max() - data.min())
                         / earth_bins)
        earth_indices = np.digitize(earth_data[:, -3], bins)
        earth_data_path = os.path.join(classified_path, 'No_layers_subdivision')
        if not os.path.exists(earth_data_path):
            os.mkdir(earth_data_path)

        for j in range(earth_bins):
            subdivision_matrix = earth_data[earth_indices == j+1, :]
            smaller_edge = str(int(round(bins[j] / 1000)))
            if j == n_bins - 1:
                larger_edge = str(int(round(bins[-1] / 1000)))
            else:
                larger_edge = str(int(round(bins[j + 1] / 1000)))
            if headers_on == 1:
                with open(os.path.join(earth_data_path, earth_name_part + 'bin_' + str(j+1) + '_edges_' + smaller_edge +
                                                        '_' + larger_edge + '.csv'), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(discretized_headers)
                    writer.writerows(subdivision_matrix)
            else:
                with open(os.path.join(earth_data_path, earth_name_part + 'bin_' + str(j+1) + '_edges_' + smaller_edge +
                                                        '_' + larger_edge + '.csv'), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(subdivision_matrix)

        if headers_on == 1:
            with open(os.path.join(earth_data_path, earth_name_part + '.csv'), 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerow(discretized_headers)
                writer.writerows(earth_data)
        else:
            with open(os.path.join(earth_data_path, earth_name_part + '.csv'), 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerows(earth_data)

        plt.figure(fig_counter)
        plt.hist(earth_data[:, -3], earth_bins, edgecolor='k', linewidth=1)
        plt.xlabel('Depth [m]')
        plt.ylabel('Number of points')
        plt.title('Distribution of points from the entire ABAQUS model at different depths, n_bins=' + str(earth_bins))
        plt.grid(axis='both', alpha=0.75)
        plt.savefig(os.path.join(histogram_path, 'Complete_Earth_distribution.png'), bbox_inches='tight')
        fig_counter += 1

