def deflection_report_reader(report_file, deflection_matrices_paths, headers_on):

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
        print("File not available or permission denied.")

    lines = []
    keyword_counter = 0
    file_identifiers = []
    headers = []
    line_counter = 0
    headers_key = 'U.Magnitude'
    search_key = 'reported at nodes'

    # Loop over each line of the rpt file
    for line in opened_report:
        line_counter += 1
        keyword_flag = line.find(search_key)
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
                search_lines = np.arange(lines[i], lines[i+1], 1)
            else:
                search_lines = np.arange(lines[i], len(opened_report), 1)
            line_counter = 0
            for line in opened_report:
                line_counter += 1
                if line_counter in search_lines:
                    s = line.strip()
                    data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
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


