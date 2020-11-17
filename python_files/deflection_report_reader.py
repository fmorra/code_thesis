def deflection_report_reader(report_file, search_key, deflection_matrices_paths, headers_on):

    import os
    import csv
    import numpy as np
    import re
    import pdb

    # Try to open the report file
    opened_report = 0
    try:
        deflection_report = open(report_file, "r")
        opened_report = deflection_report.readlines()
    except:
        print "File not available or permission denied."

    # Initialize counters to delimit the search area for the data we are interested in in the dat file, and search for
    # the line containing the headers
    lines = []
    keyword_counter = 0
    file_identifiers = []
    headers = []
    line_counter = 0
    headers_key = 'U.Magnitude'

    # Loop over each line of the rpt file
    for line in opened_report:
        line_counter += 1
        # Find the keyword that is found at the start of each stress region in the stress file
        keyword_flag = line.find(search_key)
        # If the keyword is found, increase the number of counters (for verification), add the keyword line to an array
        # containing them, extract only the part of the line after : without newline and replace the dots in this
        # extracted part with an underscore.
        if keyword_flag != -1:
            keyword_counter = keyword_counter + 1
            lines.append(line_counter)
            csv_name = line.split(': ')[1].rstrip().replace(".", "_")
            file_identifiers.append(csv_name)
        # If we also want to extract the headers for all the columns, find the header key in each file row  and if it is
        # found extract the line after Element stripping it of the newline and then substitute the dots with underscores
        headers_index = line.find(headers_key)
        if headers_index != -1:
            headers_line = line.split('Node ')[1].rstrip().replace(".", "_")
            headers.append(headers_line)
    headers = headers[0]
    # Insert the 'Labels' string at the beginning of the list containing the headers
    headers = "Label " + headers
    headers = headers.split()
    # Run the main script if the last of the csv files to be generated is not there yet
    if os.path.exists(os.path.join(deflection_matrices_paths, 'I1.csv')):
        print('The files containing the stress components already exist, moving on'
              ' to nodes coordinate extraction for each part.')
    else:
        print('Processing rpt to extract deflection components...')
        # Iterate over each of the lines where a new part or region is defined
        for i in range(len(lines)):
            deflection_submatrix = []
            # Determine whether we are in the last region of the file or not
            if i + 1 < len(lines):
                # Define the lines where to search for deflections
                search_lines = np.arange(lines[i], lines[i+1], 1)
                line_counter = 0
                # Iterate over these lines searching for float values
                for line in opened_report:
                    line_counter += 1
                    if line_counter in search_lines:
                        s = line.strip()
                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
                        # If the data vector is composed of 5 elements, fill the matrix containing the stress values
                        # with it
                        if len(data) > 4:
                            deflection_submatrix.append(np.array([float(x) for x in data]))
                # Save the matrix containing the stress values for each part and region differentiating between whether
                # we want headers or not
                if headers_on == 1:
                    with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(deflection_submatrix)
                else:
                    with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(deflection_submatrix)
            else:
                # for the last region, iterate until the end of the file, then the rest works the same as before
                search_lines = np.arange(lines[i], len(opened_report), 1)
                line_counter = 0
                for line in opened_report:
                    line_counter += 1
                    if line_counter in search_lines:
                        s = line.strip()
                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)
                        # if the data is not empty, fill the matrix containing the stress values with it
                        if len(data) > 4:
                            deflection_submatrix.append(np.array([float(x) for x in data]))

                if headers_on == 1:
                    with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerow(headers)
                        writer.writerows(deflection_submatrix)
                else:
                    with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                        writer = csv.writer(f_write)
                        writer.writerows(deflection_submatrix)

    deflection_report.close()

    return deflection_matrices_paths

