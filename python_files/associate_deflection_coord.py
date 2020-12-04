import os
import csv
from rotate_deflections import *
import numpy as np
from cart2geo import *
import pdb


def associate_deflection_coord(deflection_path, individual_node_path, headers_on, deflection_processing_path):
    # This function associates deflections at each node to the relative node coordinates.

    # Define the extension to read and extract all filenames needed
    extension_to_read = 'csv'
    dir_deflection_csv = os.listdir(deflection_path)
    csv_deflection_files = [filename for filename in dir_deflection_csv if filename.endswith(extension_to_read)]
    complete_deflection_path = os.path.join(deflection_processing_path, 'Complete_deflection_files')
    # Create relevant folders if they are not already there
    if not os.path.exists(complete_deflection_path):
        os.mkdir(complete_deflection_path)
    geographical_complete_deflection_path = os.path.join(deflection_processing_path,
                                                         'Geographical_complete_deflection_files')
    if not os.path.exists(geographical_complete_deflection_path):
        os.mkdir(geographical_complete_deflection_path)
    # Define the headers and initialize the matrices containing values for all parts
    complete_headers = ['Label', 'U_Magn', 'U_X', 'U_Y', 'U_Z', 'X', 'Y', 'Z', 'R', 'Depth', 'Lat', 'Lon']

    # Filter for empty csv files after transforming them into matrices to analyze their content
    for file_to_evaluate in range(len(csv_deflection_files)):
        if 'LOW' in csv_deflection_files[file_to_evaluate]:
            csv_deflection_files[file_to_evaluate] = 'to_delete'
        if 'I0' in csv_deflection_files[file_to_evaluate]:
            csv_deflection_files[file_to_evaluate] = 'to_delete'
    csv_deflection_files = filter(lambda a: a != 'to_delete', csv_deflection_files)
    for file_to_evaluate in range(len(csv_deflection_files)):
        with open(os.path.join(deflection_path, csv_deflection_files[file_to_evaluate])) as matrix:
            matrix_to_evaluate = matrix.readlines()
            matrix_to_evaluate = [line.strip() for line in matrix_to_evaluate[1:]]
            matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in matrix_to_evaluate]
            matrix_to_evaluate = np.array(matrix_to_evaluate)
            if len(matrix_to_evaluate) == 0:
                csv_deflection_files[file_to_evaluate] = 'to_delete'
    csv_deflection_files = filter(lambda a: a != 'to_delete', csv_deflection_files)
    # Decide whether to run the algorithm or not based on whether the last file to be created is already there
    if os.path.isfile(os.path.join(complete_deflection_path, 'Deflection_association_completion_certificate.txt')):
        print 'The files containing deflections associated to the relative coordinates' \
              ' already exist, moving on to classification of stress values based on depth.'
    else:
        print 'Associating deflection components to corresponding nodes...'
        # Iterate over the files containing the deflections for each part
        for csv_file in range(len(csv_deflection_files)):
            print 'Processing the following deflection file: ', csv_deflection_files[csv_file]
            # Open the file containing the deflections, transform it into an array for easy handling and swap the y and
            # Z columns first and then the X and Y columns to take into account the different reference system that
            # ABAQUS uses.
            with open(os.path.join(deflection_path, csv_deflection_files[csv_file]), 'r') as deflection_read_obj:
                individual_deflection_matrix = deflection_read_obj.readlines()
                individual_deflection_matrix = [line.strip() for line in individual_deflection_matrix[1:]]
                individual_deflection_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                                individual_deflection_matrix]
                individual_deflection_matrix = np.array(individual_deflection_matrix)
                individual_deflection_matrix[:, 2], individual_deflection_matrix[:, 4] = \
                    individual_deflection_matrix[:, 4], individual_deflection_matrix[:, 2].copy()
                individual_deflection_matrix[:, 3], individual_deflection_matrix[:, 4] = \
                    individual_deflection_matrix[:, 4], individual_deflection_matrix[:, 3].copy()
            # Open the file containing the coordinates, transform it into an array for easy handling and swap the Y and
            # Z columns first and then the X and Y columns to take into account the different reference system that
            # ABAQUS uses. We can perform the swap here for both the deformation tensor and coordinate vector as, unlike
            # the stresses, they are both vectors. Moreover, we don't need a complex file identification system because
            # here there is a direct connection between the number of deformation and coordinate files, since the
            # deflection report file does not differentiate into regions like the stress report file.
            with open(os.path.join(individual_node_path, 'Nodes', 'Nodes_Part_' + csv_deflection_files[csv_file]), 'r') as \
                    coordinate_read_obj:
                individual_coordinate_matrix = coordinate_read_obj.readlines()
                individual_coordinate_matrix = [line.strip() for line in individual_coordinate_matrix[1:]]
                individual_coordinate_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                                individual_coordinate_matrix]
                # Swap coordinates
                individual_coordinate_matrix = np.array(individual_coordinate_matrix)
                individual_coordinate_matrix[:, 1], individual_coordinate_matrix[:, 3] = \
                    individual_coordinate_matrix[:, 3], individual_coordinate_matrix[:, 1].copy()
                individual_coordinate_matrix[:, 2], individual_coordinate_matrix[:, 3] = \
                    individual_coordinate_matrix[:, 3], individual_coordinate_matrix[:, 2].copy()

            # Create the complete deflection matrix, rotate the deflection components, append cartesian and geographical
            # complete files to the larger matrices for all parts
            complete_deflection_matrix = np.column_stack((individual_deflection_matrix,
                                                          individual_coordinate_matrix[:, 1:]))
            geographical_complete_matrix = rotate_deflections(complete_deflection_matrix)

            earth_radius = 6371000
            depths = np.zeros((len(complete_deflection_matrix), 1))
            radial_centroid_distance = np.zeros((len(complete_deflection_matrix), 1))
            for k in range(len(complete_deflection_matrix)):
                radial_centroid_distance[k] = np.sqrt(complete_deflection_matrix[k, -3] ** 2 + complete_deflection_matrix[k, -2] ** 2 +
                                                      complete_deflection_matrix[k, -1] ** 2)
                depths[k] = earth_radius - radial_centroid_distance[k]
            complete_deflection_matrix = np.column_stack((complete_deflection_matrix, radial_centroid_distance, depths))
            geographical_complete_matrix = np.column_stack((geographical_complete_matrix, radial_centroid_distance, depths))
            cartesian_coordinates = complete_deflection_matrix[:, 5:8]
            [lat, lon] = cart2geo(cartesian_coordinates)
            complete_deflection_matrix = np.column_stack((complete_deflection_matrix, lat, lon))
            geographical_complete_matrix = np.column_stack((geographical_complete_matrix, lat, lon))
            # Save the complete matrices for each part with or without headers
            if headers_on == 1:
                with open(os.path.join(complete_deflection_path, 'Complete_file_' + csv_deflection_files[csv_file]),
                          'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(complete_headers)
                    writer.writerows(complete_deflection_matrix)
                with open(os.path.join(geographical_complete_deflection_path, 'Geographical_Complete_file_' +
                                       csv_deflection_files[csv_file]), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(complete_headers)
                    writer.writerows(geographical_complete_matrix)
            else:
                with open(os.path.join(complete_deflection_path, 'Complete_file_' + csv_deflection_files[csv_file]),
                          'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(complete_deflection_matrix)
                with open(os.path.join(geographical_complete_deflection_path, 'Geographical_Complete_file_',
                                       csv_deflection_files[csv_file]), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(geographical_complete_matrix)

        with open(os.path.join(complete_deflection_path, 'Deflection_association_completion_certificate.txt'), 'wb') \
                as f_write:
            f_write.write('Deflection association completed.')
    return complete_deflection_path, geographical_complete_deflection_path
