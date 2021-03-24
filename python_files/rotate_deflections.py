import numpy as np
import pdb


def rotate_deflections(complete_deflection_matrix):
    # This function rotates the deflections and creates a new complete matrix to be saved.

    # Extract deflections and node coordinates for each part, then allocate the necessary vectors
    displacements = complete_deflection_matrix[:, 2:5]
    node_coordinates = complete_deflection_matrix[:, 5:8]
    geographical_displacements = np.zeros((len(displacements), len(displacements[0])))
    u_magnitude = np.zeros((len(displacements), 1))
    # Iterate over the files containing all displacements for that part
    for i in range(len(displacements)):
        # Define the transformation tensor components and fill it with the necessary values
        r_3d = np.sqrt(node_coordinates[i, 0] ** 2 + node_coordinates[i, 1] ** 2 + node_coordinates[i, 2] ** 2)
        r_azimuth = np.sqrt(node_coordinates[i, 0] ** 2 + node_coordinates[i, 1] ** 2)
        cos_theta = node_coordinates[i, 0] / r_azimuth
        sin_theta = node_coordinates[i, 1] / r_azimuth
        cos_phi = r_azimuth / r_3d
        sin_phi = node_coordinates[i, 2] / r_3d

        # Create the transformation matrix for every node by defining the unit vectors for each of them
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

        # Calculate the new deformation magnitude
        new_displacements = displacements[i, :][:, np.newaxis]
        geographical_displacements[i, :] = np.transpose(np.dot(transformation_tensor, new_displacements))
        u_magnitude[i] = np.sqrt(geographical_displacements[i, 0] ** 2 + geographical_displacements[i, 1] ** 2 +
                                 geographical_displacements[i, 2] ** 2)

    # Create and return the complete matrix for this part containing the geographical displacement components
    complete_geographical_displacements = np.column_stack((u_magnitude, geographical_displacements))
    geographical_complete_matrix = np.column_stack((complete_deflection_matrix[:, 0],
                                                    complete_geographical_displacements,  node_coordinates))

    return geographical_complete_matrix
