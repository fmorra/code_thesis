def rotate_tensor(part_matrix, geographical_centroids_files_path, geographical_components_files_path,
                  geographical_complete_files_path, headers_on, complete_headers, centroid_headers,
                  complete_individual_path, part_file):

    # This function converts the stress components from catesian to spherica, leaving the coordinates as cartesian.
    # Conversion to lat and lon is carried out in the depth_classifier function.
    import numpy as np
    import os
    import csv
    import pdb

    # Allocate matrix for the geographical stresses and define headers
    geographical_components = np.zeros((len(part_matrix), 7))
    stress_headers = ['S_Mises', 'S_S11', 'S_S22', 'S_S33', 'S_S12', 'S_S13', 'S_S23']

    # Open the complete file whose stress components have to be changed
    with open(complete_individual_path, "r") as file_to_read:
        # Transform the opened file into a matrix for easy handling
        complete_file = file_to_read.readlines()
        complete_file = [line.strip() for line in complete_file[1:]]
        complete_file = [np.array([eval(i) for i in line.split(",")[:]]) for line in complete_file]
        complete_file = np.array(complete_file)

        # R, Depth, lat and lon are the same, so read the values to concatenate them to the transformed tensor
        # components
        geographical_coords = complete_file[:, 8:]
        # Define the stress tensor components with the necessary component swaps to take into account the change in
        # reference system from ABAQUS to a normal cartesian one (ZXY) -> (XYZ)
        for line in range(len(complete_file)):
            # Allocate the stress tensor
            S = np.zeros((3, 3))
            S[0, 0] = complete_file[line, 2]
            S[1, 0] = complete_file[line, 5]
            S[2, 0] = complete_file[line, 6]
            S[0, 1] = S[1, 0]
            S[1, 1] = complete_file[line, 3]
            S[2, 1] = complete_file[line, 7]
            S[0, 2] = S[2, 0]
            S[1, 2] = S[2, 1]
            S[2, 2] = complete_file[line, 4]

            # Define transformation cosines and then create the transformation matrix
            centroid_coords = complete_file[line, 8:]
            r_3d = np.sqrt(centroid_coords[0] ** 2 + centroid_coords[1] ** 2 + centroid_coords[2] ** 2)
            r_azimuth = np.sqrt(centroid_coords[0] ** 2 + centroid_coords[1] ** 2)
            cos_theta = centroid_coords[0] / r_azimuth
            sin_theta = centroid_coords[1] / r_azimuth
            cos_phi = r_azimuth / r_3d
            sin_phi = centroid_coords[2] / r_3d
        
            # T = np.zeros((3, 3))
            # T[0, 0] = cos_lon * cos_lat
            # T[1, 0] = cos_lon * sin_lat
            # T[2, 0] = -sin_lon
            # T[0, 1] = sin_lon * cos_lat
            # T[1, 1] = sin_lon * sin_lat
            # T[2, 1] = cos_lon
            # T[0, 2] = sin_lat
            # T[1, 2] = -cos_lat
            # T[2, 2] = 0.0
            transformation_tensor = np.zeros((3, 3))
            unit_x = np.array([-sin_theta, -cos_theta * sin_phi, cos_theta * cos_phi])
            unit_y = np.array([cos_theta, -sin_theta * sin_phi, sin_theta * cos_phi])
            unit_z = np.array([0, cos_phi, sin_phi])
            unit_theta = np.array([-sin_theta, cos_theta, 0])
            unit_phi = np.array([-cos_theta * sin_phi, -sin_theta * sin_phi, cos_phi])
            unit_r = np.array([cos_theta * cos_phi, sin_theta * cos_phi, sin_phi])
            transformation_tensor[0, 0] = np.dot(np.transpose(unit_theta), unit_x)
            transformation_tensor[1, 0] = np.dot(np.transpose(unit_phi), unit_x)
            transformation_tensor[2, 0] = np.dot(np.transpose(unit_r), unit_x)
            transformation_tensor[0, 1] = np.dot(np.transpose(unit_theta), unit_y)
            transformation_tensor[1, 1] = np.dot(np.transpose(unit_phi), unit_y)
            transformation_tensor[2, 1] = np.dot(np.transpose(unit_r), unit_y)
            transformation_tensor[0, 2] = np.dot(np.transpose(unit_theta), unit_z)
            transformation_tensor[1, 2] = np.dot(np.transpose(unit_phi), unit_z)
            transformation_tensor[2, 2] = np.dot(np.transpose(unit_r), unit_z)

            # Save the geographical components for each element in a line of the file allocated previously
            # geographic_tensor = np.dot(transformation_tensor, np.dot(S, np.transpose(transformation_tensor)))
            geographic_tensor = transformation_tensor.dot(S).dot(np.transpose(transformation_tensor))
            geographical_components[line, 0] = complete_file[line, 1]
            geographical_components[line, 1] = geographic_tensor[0, 0]
            geographical_components[line, 2] = geographic_tensor[1, 1]
            geographical_components[line, 3] = geographic_tensor[2, 2]
            geographical_components[line, 4] = geographic_tensor[0, 1]
            geographical_components[line, 5] = geographic_tensor[0, 2]
            geographical_components[line, 6] = geographic_tensor[1, 2]

    # Define the complete files and their labels
    labels_complete_file = complete_file[:, 0]
    labels_part_matrix = part_matrix[:, 0]
    labels_complete_file = labels_complete_file[np.newaxis]
    labels_part_matrix = labels_part_matrix[np.newaxis]
    labels_complete_file = np.transpose(labels_complete_file)
    labels_part_matrix = np.transpose(labels_part_matrix)
    complete_file_geographical = np.column_stack((labels_complete_file, geographical_components, geographical_coords))
    geographical_coords = np.column_stack((labels_part_matrix, geographical_coords[:]))

    # Save the complete files based on whether we need headers or not
    if headers_on == 1:
        # with open(os.path.join(geographical_centroids_files_path, 'Geographical_coordinates_' + part_file),
        #           'wb') as f_write:
        #     writer = csv.writer(f_write)
        #     writer.writerow(centroid_headers)
        #     writer.writerows(geographical_coords)
        # with open(os.path.join(geographical_components_files_path, 'Geographical_components_' + part_file),
        #           'wb') as f_write:
        #     writer = csv.writer(f_write)
        #     writer.writerow(complete_headers)
        #     writer.writerows(geographical_components)
        with open(os.path.join(geographical_complete_files_path, 'Geographical_complete_file_' + part_file), 'wb') as \
                f_write:
            writer = csv.writer(f_write)
            writer.writerow(complete_headers)
            writer.writerows(complete_file_geographical)
    else:
        # with open(os.path.join(geographical_centroids_files_path, 'Geographical_coordinates_' + part_file),
        #           'wb') as f_write:
        #     writer = csv.writer(f_write)
        #     writer.writerows(geographical_coords)
        with open(os.path.join(geographical_complete_files_path, 'Geographical_complete_file_' + part_file), 'wb') as \
                f_write:
            writer = csv.writer(f_write)
            writer.writerows(complete_file_geographical)

    return complete_file_geographical, geographical_coords

