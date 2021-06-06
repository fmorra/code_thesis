from region_counter import *
from rotate_tensor import *
from cart2geo import *
import pandas as pd


def associate_stress_coord(individual_element_paths, stress_part_values, large_node_matrix_path, headers_on,
                           stress_matrices_path, coord_file, coupled_stress_folder):
    # This function associates each stress value vector to the corresponding centroid coordinates, both of them in
    # Cartesian and geographical reference systems, then saves them in separate .csv files for each part.
    import os
    import numpy as np
    import csv
    import pdb

    # Define the extension we have to read data from and gather the csv names containing the element nodes and stress
    # values to associate
    file_extension = '.csv'
    individual_stress_paths = stress_matrices_path
    dir_csv_elem_files = os.listdir(individual_element_paths)
    dir_csv_stress_files = os.listdir(individual_stress_paths)
    csv_elem_files = [filename for filename in dir_csv_elem_files if filename.endswith(file_extension)]
    csv_stress_files = [filename for filename in dir_csv_stress_files if filename.endswith(file_extension)]

    # Define the relevant folders to save centroid coordinates and complete tables and create them if they do not
    # exist yet
    centroid_files_path = os.path.join(stress_part_values, 'Centroids')
    if not os.path.exists(centroid_files_path):
        os.mkdir(centroid_files_path)
    complete_files_path = os.path.join(stress_part_values, 'Complete_files')
    if not os.path.exists(complete_files_path):
        os.mkdir(complete_files_path)
    geographical_complete_files_path = os.path.join(stress_part_values, 'Geographical_complete_files')
    if not os.path.exists(geographical_complete_files_path):
        os.mkdir(geographical_complete_files_path)

    # Define headers for the complete files and the centroid files
    complete_headers = ['Label', 'S_Mises', 'S_11', 'S_22', 'S_33', 'S_12', 'S_13', 'S_23', 'X_centroid', 'Y_centroid',
                        'Z_centroid', 'R', 'Depth', 'Lat', 'Lon']
    centroid_headers = ['Centr_label', 'X_centroid', 'Y_centroid', 'Z_centroid']

    # Initialize lists to fill
    large_centroids = []
    large_complete_file = []

    # For all the next operations we will have to refer to the large matrix containing all the part nodes, because
    # for each part and for each element we will have to extract the coordinates of each node

    # Decide whether to run the main algorithm or not based on the presence of the last file to be generated
    if os.path.isfile(os.path.join(complete_files_path, 'Stress_association_completion_certificate.txt')):
        print('The files containing centroid stresses associated to the relative coordinates already exist, moving on '
              'to classification of stress values based on depth.')
    else:
        print('Associating stress components to corresponding centroids...')
        with open(large_node_matrix_path) as file_to_read:
            # Read the matrix containing the list of nodes for all parts
            all_nodes = pd.read_csv(file_to_read, delimiter=",")

            # Filter out unwanted parts and if some matrices are completely empty, do not process them
            for file_to_evaluate in range(len(csv_stress_files)):
                if 'Region' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'LOW' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'I0' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'I1' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
            csv_stress_files = [a for a in csv_stress_files if a != 'to_delete']
            for file_to_evaluate in range(len(csv_stress_files)):
                with open(os.path.join(individual_stress_paths, csv_stress_files[file_to_evaluate])) as matrix:
                    matrix_to_evaluate = pd.read_csv(matrix, delimiter=",")
                    if len(matrix_to_evaluate) == 0:
                        csv_stress_files[file_to_evaluate] = 'to_delete'

        csv_stress_files = [a for a in csv_stress_files if a != 'to_delete']
        # Iterate over the files containing the stress values for each part or region
        for i in range(len(csv_stress_files)):
            # Eliminate the empty elements from the list at each iteration, because after we process each file they
            # will be taken away from the list
            csv_stress_files = [_f for _f in csv_stress_files if _f]
            csv_elem_files = [_f for _f in csv_elem_files if _f]
            # Select the stress file to process
            part_file = csv_stress_files[0]
            print('Processing the following stress file: ' + part_file)

            file_to_read_logical = [i.split('Elements_Part_')[1] == part_file for i in csv_elem_files]
            file_to_read_index = np.array([file_to_read_logical.index(i) for i in file_to_read_logical
                                           if i == True])
            # Open the element file and transform both the element and stress files into arrays for easier handling
            file_to_read_index = file_to_read_index[0]
            all_elems_to_open = open(os.path.join(individual_element_paths, csv_elem_files[file_to_read_index]))
            # all_elems = all_elems_to_open.readlines()
            # all_elems = [line.strip() for line in all_elems[1:]]
            # all_elems = [np.array([eval(i) for i in line.split(",")[:]]) for line in all_elems]
            # all_elems = np.array(all_elems)
            all_elems = pd.read_csv(all_elems_to_open, delimiter=",")
            part_matrix_to_open = open(os.path.join(individual_stress_paths, part_file))
            part_matrix = pd.read_csv(part_matrix_to_open, delimiter=",")
            centroid_coord = np.zeros((len(part_matrix), 4))

            # For each element, fill a vector containing all of its nodes; extract the coordinates for each of them,
            # then calculate the mean across these nodes for all three coordinates
            for j in range(len(part_matrix)):
                elem_nodes_index = [np.where(all_elems.iloc[:, 0] == part_matrix.iloc[j, 0])]
                elem_nodes = all_elems.iloc[elem_nodes_index, 1:]
                elem_nodes = elem_nodes.iloc[np.nonzero(elem_nodes)]
                elem_coord_matrix = np.zeros((len(elem_nodes), 3))
                for k in range(0, len(elem_nodes)):
                    index = [np.where(all_nodes.iloc[:, 0] == elem_nodes.iloc[k])]
                    nodes_coord = all_nodes.iloc[index, 1:]
                    if not len(nodes_coord) == 0:
                        elem_coord_matrix[k, :] = nodes_coord
                for m in range(0, 3):
                    centroid_coord[j, m+1] = np.mean(elem_coord_matrix[:, m])
                centroid_coord[j, 0] = part_matrix[j, 0]
                if np.all(centroid_coord[j, 1:]) == 0:
                    centroid_coord = np.delete(centroid_coord, j, axis=0)

            # Because of ABAQUS' reference system, swap the Y and Z coordinate first and then the X and Y
            # coordinates

            centroid_coord[:, 1], centroid_coord[:, 3] = centroid_coord[:, 3].copy(), centroid_coord[:, 1].copy()
            centroid_coord[:, 2], centroid_coord[:, 3] = centroid_coord[:, 3].copy(), centroid_coord[:, 2].copy()
            part_matrix[:, 2], part_matrix[:, 4] = part_matrix[:, 4].copy(), part_matrix[:, 2].copy()
            part_matrix[:, 3], part_matrix[:, 4] = part_matrix[:, 4].copy(), part_matrix[:, 3].copy()
            part_matrix[:, 5], part_matrix[:, 7] = part_matrix[:, 7].copy(), part_matrix[:, 5].copy()
            part_matrix[:, 5], part_matrix[:, 6] = part_matrix[:, 6].copy(), part_matrix[:, 5].copy()

            # Create the complete file for each part and region containing the stresses and centroid coordinates,
            # appending each region file of the Earth part to the Earth large matrix so that we can also have the
            # values for the part and not just for its sub-regions
            complete_file = np.hstack((part_matrix, centroid_coord[:, 1:]))
            earth_radius = 6371000
            depths = np.zeros((len(complete_file), 1))
            radial_centroid_distance = np.zeros((len(complete_file), 1))
            for k in range(len(complete_file)):
                radial_centroid_distance[k] = np.sqrt(complete_file[k, -3] ** 2 + complete_file[k, -2] ** 2 +
                                                      complete_file[k, -1] ** 2)
                depths[k] = earth_radius - radial_centroid_distance[k]
            complete_file = np.column_stack((complete_file, radial_centroid_distance, depths))
            cartesian_coordinates = complete_file[:, 8:11]
            # Calculated the geographical coordinates and append them to the previous matrix
            [lat, lon] = cart2geo(cartesian_coordinates)
            complete_file = np.column_stack((complete_file, lat, lon))

            # Define name used to save the files
            centroid_individual_path = os.path.join(centroid_files_path, 'Centroid_file_' + part_file)
            complete_individual_path = os.path.join(complete_files_path, 'Complete_file_' + part_file)
            # Write the centroid coordinates and the complete files in csvs with headers or not
            if headers_on == 1:
                with open(centroid_individual_path, 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(centroid_headers)
                    writer.writerows(centroid_coord)
                with open(complete_individual_path, 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(complete_headers)
                    writer.writerows(complete_file)
            else:
                with open(centroid_individual_path, 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(centroid_coord)
                with open(complete_individual_path, 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(complete_file)

            # Calculate and save the geographical stress components
            rotate_tensor(part_matrix, geographical_complete_files_path, headers_on, complete_headers,
                          complete_individual_path, part_file)

        # Only save
        for csv_file in range(len(csv_stress_files)):
            if csv_stress_files[csv_file] == 'EARTH.csv':
                if not os.path.exists(os.path.join(complete_files_path, 'Geographical_complete_file_coupled_' +
                                                                        csv_stress_files[csv_file])):
                    with open(os.path.join(complete_files_path, 'Complete_file_' + csv_stress_files[csv_file]),
                              'r') as uncoupled_coordinate_read_obj:
                        uncoupled_coordinate_matrix = uncoupled_coordinate_read_obj.readlines()
                        uncoupled_coordinate_matrix = [line.strip() for line in uncoupled_coordinate_matrix[1:]]
                        uncoupled_coordinate_matrix = [np.array([eval(i) for i in line.split(",")[:]])
                                                       for line in uncoupled_coordinate_matrix]
                        uncoupled_coordinate_matrix = np.array(uncoupled_coordinate_matrix)
                    with open(os.path.join(coupled_stress_folder, csv_stress_files[csv_file]), 'r') as coordinate_read_obj:
                        individual_coordinate_matrix = coordinate_read_obj.readlines()
                        individual_coordinate_matrix = [line.strip() for line in individual_coordinate_matrix[1:]]
                        individual_coordinate_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                                        individual_coordinate_matrix]
                        individual_coordinate_matrix = np.array(individual_coordinate_matrix)
                    final_coupled_matrix = np.hstack((individual_coordinate_matrix, uncoupled_coordinate_matrix[:, 8:]))
                    complete_individual_path_coupled = os.path.join(complete_files_path, 'Complete_file_coupled_' +
                                                                    csv_stress_files[csv_file])
                    with open(complete_individual_path_coupled, 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(complete_headers)
                        writer.writerows(final_coupled_matrix)
                    rotate_tensor(final_coupled_matrix, geographical_complete_files_path, headers_on, complete_headers,
                                  complete_individual_path_coupled, 'coupled_' + csv_stress_files[csv_file])

        with open(os.path.join(complete_files_path, 'Stress_association_completion_certificate.txt'), 'wb') as f_write:
            f_write.write('Stress association completed.')

    return centroid_files_path, complete_files_path, geographical_complete_files_path


