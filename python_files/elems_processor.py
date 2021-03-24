import os
import numpy as np
import re
import csv
import pdb


def elems_processor(node_lines, elem_lines, nset_lines, stress_processing_path,
                    dat_path, file_identifiers, headers_on):

    # Define the lines which denote the parts of the file where we have to search for the element nodes; however,
    # towards the end of the part of the file we are interested in there is a series of lines containing the keyword
    # which is supposed to tell the algorithm to stop looking for element nodes. Therefore we have to filter out all
    # of these lines and then we transform the lists containing the lines into arrays.
    large_elem_matrix_unsorted = []
    large_elem_matrix_unsorted_reset = []
    new_elem_lines = elem_lines[0:len(node_lines) - 1]
    new_nset_lines = np.zeros((len(node_lines) - 1))
    for i in range(0, len(new_elem_lines)):
        for j in range(i, len(nset_lines)):
            if nset_lines[j] > new_elem_lines[i]:
                new_nset_lines[i] = nset_lines[j]
                break
    new_elem_lines = np.array(new_elem_lines)
    new_elem_lines = new_elem_lines.astype(int)
    new_nset_lines = new_nset_lines.astype(int)
    # Define all relevant paths and create the relevant folders if necessary
    large_elem_matrix_path = os.path.join(stress_processing_path, 'Large_Element_Matrix' + '.csv')
    individual_element_paths = os.path.join(stress_processing_path, 'Elements')
    if not os.path.exists(individual_element_paths):
        os.makedirs(individual_element_paths)
    # Define the headers and decide whether to run the code or not based on the presence of the last file to be created
    headers = ['Label', 'Node_1', 'Node_2', 'Node_3', 'Node_4', 'Node_5', 'Node_6', 'Node_7', 'Node_8']
    if os.path.isfile(large_elem_matrix_path):
        print('The files containing the element nodes already exist, moving on to centroid calculation and stress '
              'association.')
    else:
        print('Processing nodes to extract coordinates...')
        # Open the dat file, this time to save all element nodes 
        with open(dat_path, "r") as read_dat:
            dat_file = read_dat.readlines()
            for k in range(len(new_elem_lines)):
                # Only focus on the areas of the file where information about the elements is present
                elem_vector = np.zeros((new_nset_lines[k] - new_elem_lines[k], 9))
                # Iterate over each of the lines in the interval found and scan for float values
                for line in range(new_elem_lines[k], new_nset_lines[k]):
                    s = dat_file[line].strip()
                    data = [float(s) for s in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)]
                    data = np.array([float(x) for x in data])
                    # 4 cases have to be taken into account: two because elements can be made up of 8 or 6 elements, and
                    # 2 more because to this we have to occasionally take into account and filter out the counter every
                    # 5 lines that ABAQUS leaves in the dat file.
                    if len(data) > 6:
                        if len(data) == 10:
                            data = data[1:]
                        elif len(data) == 8:
                            data = data[1:]
                        if len(data) == 9:
                            elem_vector[line - new_elem_lines[k] + 1, :] = data
                        else:
                            elem_vector[line - new_elem_lines[k] + 1, 0:7] = data
                        # Filter out possible numerical errors
                        if elem_vector[line - new_elem_lines[k] + 1, 1] == elem_vector[line - new_elem_lines[k] + 1, 2]:
                            elem_vector[line - new_elem_lines[k] + 1, :] = 0

                elem_vector = elem_vector[~np.all(elem_vector == 0, axis=1)]
                individual_elem_matrix = elem_vector

                # Save inndividual element matrices in csv format and append them to the large element lists
                if headers_on == 1:
                    with open(os.path.join(individual_element_paths, 'Elements_Part_' + file_identifiers[k] + ".csv"),
                              'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(individual_elem_matrix)
                else:
                    with open(os.path.join(individual_element_paths, 'Elements_Part_' + file_identifiers[k] + ".csv"),
                              'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(individual_elem_matrix)

                large_elem_matrix_unsorted_reset.append(individual_elem_matrix)
                large_elem_matrix_unsorted.append(elem_vector)

        # Sort the matrix containing all elements
        element_dictionary = {}
        for elem_matrix in large_elem_matrix_unsorted:
            for j in range(len(elem_matrix)):
                dictionary_matrix = elem_matrix[j]
                element_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
        element_ids = element_dictionary.keys()
        sorted_matrix = np.zeros((len(element_ids), 9))
        sorted_element_ids = list(sorted(element_ids))
        for element in range(len(sorted_element_ids)):
            sorted_matrix[element, 0] = sorted_element_ids[element]
            sorted_matrix[element, 1:] = element_dictionary[sorted_element_ids[element]]
        # Save the larger matrices with or without headers
        if headers_on == 1:
            with open(large_elem_matrix_path, 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerow(headers)
                writer.writerows(sorted_matrix)
        else:
            with open(large_elem_matrix_path, 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerows(sorted_matrix)

    return individual_element_paths
