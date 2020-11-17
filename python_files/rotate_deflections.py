import numpy as np
import pdb


def rotate_deflections(complete_deflection_matrix):

    # This function rotates the deflections and then saves them, the node coordinates are cartesian and are converted
    # to latitude and longitude in th depth_classifier function.

    # Extract deflections and node coordinates for each part, then allocate the necessary vectors
    displacements = complete_deflection_matrix[:, 2:5]
    node_coordinates = complete_deflection_matrix[:, 5:8]
    geographical_displacements = np.zeros((len(displacements), len(displacements[0])))
    U_magn = np.zeros((len(displacements), 1))
    # Iterate over the files containing all displacements for that part
    for i in range(len(displacements)):
        # Define the transformation tensor components and fill it with the necessary values
        R_3D = np.sqrt(node_coordinates[i, 0] ** 2 + node_coordinates[i, 1] ** 2 + node_coordinates[i, 2] ** 2)
        R_azimuth = np.sqrt(node_coordinates[i, 0] ** 2 + node_coordinates[i, 1] ** 2)
        cos_theta = node_coordinates[i, 0] / R_azimuth
        sin_theta = node_coordinates[i, 1] / R_azimuth
        cos_phi = R_azimuth / R_3D
        sin_phi = node_coordinates[i, 2] / R_3D

        transformation_tensor = np.zeros((3, 3))
        transformation_tensor[0, 0] = cos_theta * cos_phi
        transformation_tensor[1, 0] = cos_theta * sin_phi
        transformation_tensor[2, 0] = -sin_theta
        transformation_tensor[0, 1] = sin_theta * cos_phi
        transformation_tensor[1, 1] = sin_theta * sin_phi
        transformation_tensor[2, 1] = cos_theta
        transformation_tensor[0, 2] = sin_phi
        transformation_tensor[1, 2] = -cos_phi
        transformation_tensor[2, 2] = 0.0
        # Calculate the geographical displacements and then use them to calculate the deflection magnitude
        geographical_displacements[i, :] = np.transpose(transformation_tensor.dot(np.transpose(displacements[i, :])))
        U_magn[i] = np.sqrt(geographical_displacements[i, 0] ** 2 + geographical_displacements[i, 1] ** 2 +
                            geographical_displacements[i, 2] ** 2)
    # Create and return the complete matrix for this part containing the geographical displacement components
    complete_geographical_displacements = np.column_stack((U_magn, geographical_displacements))
    geographical_complete_matrix = np.column_stack((complete_deflection_matrix[:, 0], complete_geographical_displacements,
                                                 node_coordinates))

    return geographical_complete_matrix
