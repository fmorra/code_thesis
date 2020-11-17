def stress_report_reader(report_file, search_key, stress_matrices_path, headers_on):

    import os
    import csv
    import numpy as np
    import re
    import pdb

    # Try to open the report file
    opened_report = 0
    try:
        stress_report = open(report_file, "r")
        opened_report = stress_report.readlines()
    except:
        print "File not available or permission denied."
    # Initialize counters to delimit the search area for the data we are interested in in the dat file, and search for
    # the line containing the headers
    lines = []
    keyword_counter = 0
    file_identifiers = []
    headers = []
    line_counter = 0
    headers_key = 'S.Mises'

    # Loop over each line of the rpt file
    for line in opened_report:
        line_counter += 1
        # Find the keyword that is found at the start of each part of the report file that contains new stresses for
        # each part or region
        keyword_flag = line.find(search_key)
        # If the keyword is found, increase the number of counters (for verification), add the keyword line to an array
        # containing them, extract only the part of the line after : without newline and replace the dots in this
        # extracted part with an underscore. This part of the line will be used as the csv file name.
        if keyword_flag != -1:
            keyword_counter = keyword_counter + 1
            lines.append(line_counter)
            csv_name = line.split(': ')[1].rstrip().replace(".", "_")
            file_identifiers.append(csv_name)
        # If we also want to extract the headers for all the columns, find the header key in each file row and if it is
        # found extract the line after the word "Element" stripping it of the newline and then substitute the dots
        # with underscores.
        headers_index = line.find(headers_key)
        if headers_index != -1:
            headers_line = line.split('Element ')[1].rstrip().replace(".", "_")
            headers.append(headers_line)
    headers = headers[0]
    # Insert the 'Labels' string at the beginning of the list containing the headers
    headers = "Label " + headers
    headers = headers.split()
    large_stress_matrix = []
    large_stress_matrix_path = os.path.join(stress_matrices_path, 'EARTH.csv')
    print large_stress_matrix_path
    # Run the main script if the last of the csv files to be generated is not there yet
    if os.path.exists(large_stress_matrix_path):
        print('The files containing the stress components already exist, moving on'
              ' to stress coordinate extraction for each part.')
    else:
        print('Processing rpt to extract stress components...')
        # Iterate over each of the lines where a new part or region is defined
        for i in range(len(lines)):
            stress_submatrix = []
            # Determine whether we are in the last region of the file or not
            if i + 1 < len(lines):
                # Define the lines where to search for stresses
                search_lines = np.arange(lines[i], lines[i+1], 1)
                line_counter = 0
                # Iterate over these lines searching for float values
                for line in opened_report:
                    line_counter += 1
                    if line_counter in search_lines:
                        s = line.strip()
                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
                        # If the data vector is composed of 8 elements, fill the matrix containing the stress values
                        # with it
                        if len(data) > 7:
                            stress_submatrix.append(np.array([float(x) for x in data]))
                # Save the matrix containing the stress values for each part and region differentiating between whether
                # we want headers or not
                if headers_on == 1:
                    with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(stress_submatrix)
                else:
                    with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(stress_submatrix)
                if 'Region' in file_identifiers[i]:
                    large_stress_matrix.append(stress_submatrix)
            else:
                # For the last region, iterate until the end of the file, then the rest works the same as before
                search_lines = np.arange(lines[i], len(opened_report), 1)
                line_counter = 0
                for line in opened_report:
                    line_counter += 1
                    if line_counter in search_lines:
                        s = line.strip()
                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)
                        # if the data is not empty, fill the matrix containing the stress values with it
                        if len(data) > 7:
                            stress_submatrix.append(np.array([float(x) for x in data]))

                if headers_on == 1:
                    with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(stress_submatrix)
                else:
                    with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(stress_submatrix)
                if 'Region' in file_identifiers[i]:
                    large_stress_matrix.append(stress_submatrix)
            stress_dictionary = {}
            for stress_matrix in large_stress_matrix:
                for j in range(len(stress_matrix)):
                    dictionary_matrix = stress_matrix[j]
                    stress_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
            stress_ids = stress_dictionary.keys()
            sorted_matrix = np.zeros((len(stress_ids), 8))
            sorted_stress_ids = list(sorted(stress_ids))
            for stress in range(len(sorted_stress_ids)):
                sorted_matrix[stress, 0] = sorted_stress_ids[stress]
                sorted_matrix[stress, 1:] = stress_dictionary[sorted_stress_ids[stress]]
            if headers_on == 1:
                with open(os.path.join(large_stress_matrix_path), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerow(headers)
                    writer.writerows(sorted_matrix)
            else:
                with open(os.path.join(large_stress_matrix_path), 'wb') as f_write:
                    writer = csv.writer(f_write)
                    writer.writerows(sorted_matrix)
    stress_report.close()
    return file_identifiers, headers


