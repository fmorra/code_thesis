from region_counter import *
from rotate_tensor import *
from cart2geo import *


def associate_stress_coord(individual_element_paths, stress_part_values, large_node_matrix_path, headers_on,
                           stress_matrices_path, coupled_stress_folder):
    # This function calculates element centroids and associates each of them to the corresponding stress components.
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

    # Define the folders where to save centroid coordinates and the complete files with stresses and centroid
    # coordinates
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

    # For all the next operations it is necessary to refer to the large matrix containing all the part nodes, because
    # for each part and for each element we will have to extract the coordinates of each node

    # Only run the script if the last files to be generated are not there already
    if os.path.isfile(os.path.join(complete_files_path, 'Stress_association_completion_certificate.txt')):
        print 'The files containing centroid stresses associated to the relative coordinates already exist, ' \
              'moving on to classification of stress values based on depth.'
    else:
        print 'Associating stress components to corresponding centroids...'
        with open(large_node_matrix_path) as file_to_read:
            # Transform the list into a traditional matrix
            all_nodes = file_to_read.readlines()
            all_nodes = [line.strip() for line in all_nodes[1:]]
            all_nodes = [np.array([eval(i) for i in line.split(",")[:]]) for line in all_nodes]
            all_nodes = np.array(all_nodes)

            # If some matrices are completely empty, do not process them; this check can be done after transforming
            # each matrix into an array
            for file_to_evaluate in range(len(csv_stress_files)):
                if 'Region' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'LOW' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'I0' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'I1' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
            csv_stress_files = filter(lambda a: a != 'to_delete', csv_stress_files)
            for file_to_evaluate in range(len(csv_stress_files)):
                with open(os.path.join(individual_stress_paths, csv_stress_files[file_to_evaluate])) as matrix:
                    matrix_to_evaluate = matrix.readlines()
                    matrix_to_evaluate = [line.strip() for line in matrix_to_evaluate[1:]]
                    matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                          matrix_to_evaluate]
                    matrix_to_evaluate = np.array(matrix_to_evaluate)
                    if len(matrix_to_evaluate) == 0:
                        csv_stress_files[file_to_evaluate] = 'to_delete'

        csv_stress_files = filter(lambda a: a != 'to_delete', csv_stress_files)
        # Iterate over the files containing the stress values for each part or region
        for i in range(len(csv_stress_files)):
            # Eliminate the empty elements from the list at each iteration, because after we process each file they
            # will be taken away from the list
            csv_stress_files = filter(None, csv_stress_files)
            csv_elem_files = filter(None, csv_elem_files)
            # Select the stress file to process
            part_file = csv_stress_files[0]
            print 'Processing the following stress file: ', part_file
            part_matrix_to_open = open(os.path.join(individual_stress_paths, part_file))
            part_matrix = part_matrix_to_open.readlines()

            file_to_read_logical = [i.split('Elements_Part_')[1] == part_file for i in csv_elem_files]
            file_to_read_index = np.array([file_to_read_logical.index(i) for i in file_to_read_logical
                                           if i == True])
            # Open the element file and transform both the element and stress files into arrays for easier handling
            file_to_read_index = file_to_read_index[0]
            all_elems_to_open = open(os.path.join(individual_element_paths, csv_elem_files[file_to_read_index]))
            all_elems = all_elems_to_open.readlines()
            all_elems = [line.strip() for line in all_elems[1:]]
            all_elems = [np.array([eval(i) for i in line.split(",")[:]]) for line in all_elems]
            all_elems = np.array(all_elems)
            part_matrix = [line.strip() for line in part_matrix[1:]]
            part_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in part_matrix]
            part_matrix = np.array(part_matrix)
            centroid_coord = np.zeros((len(part_matrix), 4))

            # For each element, fill a vector containing all of its nodes; extract the coordinates for each of them,
            # then calculate the mean across these nodes for all three coordinates
            for j in range(len(part_matrix)):
                elem_nodes_index = [np.where(all_elems[:, 0] == part_matrix[j, 0])]
                elem_nodes = all_elems[elem_nodes_index, 1:]
                elem_nodes = elem_nodes[np.nonzero(elem_nodes)]
                elem_coord_matrix = np.zeros((len(elem_nodes), 3))
                for k in range(0, len(elem_nodes)):
                    index = [np.where(all_nodes[:, 0] == elem_nodes[k])]
                    nodes_coord = all_nodes[index, 1:]
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
            # then calculate radial distance, depth, latitude and longitude for each element and also add them to the
            # matrix
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

            # Calculate and save the geographical stress components in a complete matrix
            rotate_tensor(part_matrix, geographical_complete_files_path, headers_on, complete_headers,
                          complete_individual_path, part_file)

        # Save the complete coupled stresses matrix
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

        # Write file attesting compeltion of the program
        with open(os.path.join(complete_files_path, 'Stress_association_completion_certificate.txt'), 'wb') as f_write:
            f_write.write('Stress association completed.')

    return centroid_files_path, complete_files_path, geographical_complete_files_path


