def nodes_processor(common_files_path, node_lines, elem_lines, file_identifiers, dat_path, headers_on):
    # Read all the nodes for each part and save them in separate .csv files, then save a sorted list for all parts
    # using a dictionary
    import os
    import numpy as np
    import re
    import csv
    import pdb

    # here define the .csv path and the file headers
    large_node_matrix_path = os.path.join(common_files_path, 'Large_Node_Matrix.csv')
    individual_path = os.path.join(common_files_path, 'Individual_part_values')
    if not os.path.exists(individual_path):
        os.makedirs(individual_path)
    headers = ['Label', 'X', 'Y', 'Z']
    large_node_matrix = []

    # Skip the operations if the last file to be created is already there
    if os.path.isfile(large_node_matrix_path):
        print ('The files containing the node coordinates already exist, moving on'
               ' to element nodes extraction for each part.')
    else:
        print('Processing nodes to extract coordinates...')
        # Open the dat file as a list of lines
        with open(dat_path, "r") as read_dat:
            dat_file = read_dat.readlines()
            # Iterate over the lines where, in the dat file, the definition of all the node coordinates for each part
            # begins
            for i in range(0, len(node_lines) - 1):
                # Allocate the matrix where all node values for that part will be stored
                individual_matrix = []
                # Iterate only over the lines of the dat file where the node coordinates are defined
                for line in range(node_lines[i], elem_lines[i]):
                    s = dat_file[line].strip()
                    # Search for floats at each line
                    data = [float(s) for s in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)]
                    data = np.array([float(x) for x in data])
                    # Delete the first element of the data vector if the vectr is longer than 4, because every 5 lines
                    # ABAQUS leaves a line counter which is not needed, then append the data to the individual matrix
                    if len(data) > 3:
                        if data[0] is not 0 and data[0] == data[1] and data[1] == data[2] and data[2] == data[3]:
                            pass
                        else:
                            if len(data) == 5:
                                data = data[1:]
                            individual_matrix.append(data)
                # Append each individual matrix to the matrix containing the nodes of all the parts
                large_node_matrix.append(individual_matrix)

                # Save the file for every single part differentiating between headers or not
                if headers_on == 1:
                    with open(os.path.join(individual_path, 'Nodes_Part_' + file_identifiers[i] + ".csv"), 'wb') \
                            as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(individual_matrix)
                else:
                    with open(os.path.join(individual_path, 'Nodes_Part_' + file_identifiers[i] + ".csv"), 'wb') \
                            as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(individual_matrix)

        # Use dictionaries to sort all the nodes, then save the sorted node list for all parts
        node_dictionary = {}
        for node_matrix in large_node_matrix:
            for j in range(len(node_matrix)):
                dictionary_matrix = node_matrix[j]
                node_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
        node_ids = list(node_dictionary.keys())
        sorted_matrix = np.zeros((len(node_ids), 4))
        sorted_node_ids = list(sorted(node_ids))
        for node in range(len(sorted_node_ids)):
            sorted_matrix[node, 0] = sorted_node_ids[node]
            sorted_matrix[node, 1:] = node_dictionary[sorted_node_ids[node]]
        if headers_on == 1:
            with open(os.path.join(large_node_matrix_path), 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerow(headers)
                writer.writerows(sorted_matrix)
        else:
            with open(os.path.join(large_node_matrix_path), 'wb') as f_write:
                writer = csv.writer(f_write)
                writer.writerows(sorted_matrix)

    return individual_path, large_node_matrix_path
