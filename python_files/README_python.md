

### Main file - python_reader_processing_only

This is the main file of the data processing step, and begins with the definition of a series of parameters and paths used to locate the files to analyze and to store the files that will be successively generated if they do not exist already. The base path follows a naming convention for which the run folder is called "run_n", the iteration, step and stress iteration directories are all defined based on it and therefore the user only needs to modify one line to move the location of all the files to be generated. 

```python
run = '18'
iteration = 1
step = 0
cycle = 1
if not isinstance(run, str):
    run = str(run)
stress_rpt_name = 'abaqus_run_' + run + '.rpt'
deflection_rpt_name = 'abaqus_run_' + run + '_deflections.rpt'
# e_csv_name = 'e.csv'
if run == 'point':
    dat_name = 'Earth_point.dat'
else:
    dat_name = 'Earth.dat'
e_dat_name = 'e.dat'
# Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
# windows or linux.
opersys = platform.system()
if opersys == 'Windows':
    base_path = os.path.join('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
else:
    base_path = os.path.join('/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run))

# Define the folder path where to find the report and e.dat files. The ABAQUS output to use consists in this folder that
# has to be manually downloaded using mobaXterm and the Earth.dat file, and it is associaated not only to the run number
# but also to an iteration, a step and a cycle number.

iteration_path = os.path.join(base_path, 'Iteration_' + str(iteration))
if not os.path.exists(iteration_path):
    os.mkdir(iteration_path)
step_path = os.path.join(iteration_path, 'step_' + str(step))
if not os.path.exists(step_path):
    os.mkdir(step_path)
cycle_path = os.path.join(step_path, 'cycle_' + str(cycle) + '_reports')
if not os.path.exists(cycle_path):
    os.mkdir(cycle_path)
level_dir = cycle_path

# Define the paths for the configuration .dat file and report files
dat_path = os.path.join(base_path, dat_name)
stress_report_path = os.path.join(level_dir, stress_rpt_name)
deflection_report_path = os.path.join(level_dir, deflection_rpt_name)

# Define all the other encessary directories and, if they do not exist yet, create them
stress_matrices_path = os.path.join(level_dir, 'stress_matrices')
stress_processing_path = os.path.join(level_dir, 'stress_processing_results')
deflection_matrices_paths = os.path.join(level_dir, 'deflection_matrices')
deflection_processing_path = os.path.join(level_dir, 'deflection_processing_results')
common_files_path = os.path.join(level_dir, 'common_files')
coupled_stress_folder = os.path.join(stress_matrices_path, 'coupled_stress_matrices')
base_e_path = os.path.join(level_dir, e_dat_name)

if not os.path.exists(stress_matrices_path):
    os.mkdir(stress_matrices_path)
if not os.path.exists(stress_processing_path):
    os.mkdir(stress_processing_path)
if not os.path.exists(deflection_matrices_paths):
    os.mkdir(deflection_matrices_paths)
if not os.path.exists(deflection_processing_path):
    os.mkdir(deflection_processing_path)
```

The next part is a function call to coordinate_reader, where every element or node is going to be associated to its centroid or node coordinates and depth, latitude and longitude.

```python
start_time_coordinate_reader = time.time()
individual_path, complete_files_path, geographical_complete_files_path, large_node_matrix_path = \
    coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_matrices_paths,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder)
end_time_coordinate_reader = time.time() - start_time_coordinate_reader
```
---
### Coordinate_reader and its functions 
Coordinate_reader calls a series of functions to read data and associate stresses and deformations to their respective coordinates and additional values such as the latitude and longitude of the corresponding element or node, and they will be described in order here. A timer is set and run for each of them. The first function to be run is called dat_processor. 

#### dat_processor
``` python
node_lines, elem_lines, nset_lines, file_identifiers, new_earth_dat = dat_processor(dat_path)
```

The first part of this function initializes a series of variables and counters used to delimit "regions" of the Earth.dat file (processed here) where only the nodes and elements for all parts, associated to their coordinates and composing nodes respectively, are located. First some lists are initialized where to save the lines where nodes, elements, part descriptions and node sets begin, then the counters for the amount of these lines are also initialized and finally the keywords used to identify where each "region of interest" containing the values to extract begins and ends are defined. The last counter 

```python
    with open(fname, "r+") as read_dat:
        # Read the file as a list of lines
        opened_dat = read_dat.readlines()
        print "Processing the dat file to make it more readable"
        # Initialize counters containing a series of lines in the dat file where we will perform the search for node
        # coordinates and element nodes
        node_lines = []
        elem_lines = []
        part_lines = []
        nset_lines = []
        node_keyword_counter = 0
        elem_keyword_counter = 0
        part_keyword_counter = 0
        nset_keyword_counter = 0
        file_identifiers = []
        # Define the keywords used to search for the lines delimiting the start and end of those regions
        node_search_key = '*Node'
        element_search_key = '*Element'
        part_search_key = 'PART INSTANCE'
        region_finisher_key = '*Nset'
        input_paragraph_end = 'OPTIONS BEING PROCESSED'
        line_counter = 0
        end_line = 0
        end_line_counter = 0
```
The next section of the function scans every line to find the keywords defined above, and when they are found, increases the counters for every of them and appends the line number to the node, elements etc, lists. The last elif checks whether we have reached a point in the .dat file where, after that, there are no interesting values anymore, because in that case all the lines after tat can also not be processed. 
```python
        for line in opened_dat:
            line = line.strip()
            line_counter += 1
            node_keyword_flag = line.find(node_search_key)
            elem_keyword_flag = line.find(element_search_key)
            part_keyword_flag = line.find(part_search_key)
            region_finisher_flag = line.find(region_finisher_key)
            end_keyword_flag = line.find(input_paragraph_end)
            if node_keyword_flag != -1:
                node_keyword_counter = node_keyword_counter + 1
                node_lines.append(line_counter)
            elif elem_keyword_flag != -1:
                elem_keyword_counter = elem_keyword_counter + 1
                elem_lines.append(line_counter)
            elif region_finisher_flag != -1:
                nset_keyword_counter = nset_keyword_counter + 1
                nset_lines.append(line_counter)
            elif part_keyword_flag != -1:
                part_keyword_counter = part_keyword_counter + 1
                part_lines.append(line_counter)
                file_identifiers.append(line.split(": ")[1].rstrip().replace(".", "_"))
            elif end_keyword_flag != -1:
                end_line = line
                end_line_counter = line_counter
                break
            else:
                pass
        node_lines.append(end_line_counter)
```
The last part of the file creates a new .dat file with all the necessary information but without commas and taking away the "LINE" word at the beginning of every 5 lines.

```python
    new_earth_dat = fname[0:-4] + "_new.dat"
    with open(fname, "r+") as read_dat:
        with open(new_earth_dat, "w") as new_read_dat:
            for line in read_dat:
                replacements = line.replace(',', ' ').replace('LINE', ' ')
                new_read_dat.write(replacements)

```
After displaying the time spent to run this function, the next one is nodes_processor, which stores all the nodes and coordinates for each part to process in a separate .csv file. Moreover, a large file is created with a list of nodes and their coordinates for all parts.
#### nodes_processor
The first part of the function defines, and if necessary creates, the paths where to save the files and also defines the headers to use for every .csv to file to create and initializes the large matrix of nodes for all parts, large_node_matrix.
```python
    # Define the relevant paths and if necessary create relevant directories
    large_node_matrix_path = os.path.join(common_files_path, 'Large_Node_Matrix.csv')
    individual_path = os.path.join(common_files_path, 'Individual_Part_Values')
    individual_node_path = os.path.join(individual_path, 'Nodes')
    if not os.path.exists(individual_node_path):
        os.makedirs(individual_node_path)

    # Define node labels and initialize matrix with all node data
    headers = ['Label', 'X', 'Y', 'Z']
    large_node_matrix = []
```
Next, the function checks whether the last file to be generated already exists or not, and if it does, this process is skipped to save time. However, if the .csv files have to be created, the .dat file is opened and individual_matrix is initialized as a list to store the node values. Then, the function iterates over each of the possible intervals between node_lines and elem_lines, which is where all the nodes for each part are defined, to extract all doubles in that row. THis vector of doubles is checked to throw away values that are no interest for this work such as the line counter present every 5 lines in the new Earth.dat file to process.
```python
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
```
The last part of the function defines a dictionary to sort large_node_matrix before saving it based on whether headers are desired or not.

```python
node_dictionary = {}
        for node_matrix in large_node_matrix:
            for j in range(len(node_matrix)):
                dictionary_matrix = node_matrix[j]
                node_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
        node_ids = node_dictionary.keys()
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
                writer.writerows(sorted_matrix)node_dictionary = {}
        for node_matrix in large_node_matrix:
            for j in range(len(node_matrix)):
                dictionary_matrix = node_matrix[j]
                node_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
        node_ids = node_dictionary.keys()
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
```
After displaying again the time spent to run this function, coordinate_reader checks whether the user is processing stresses or deflections based on the value of sd_input. This happens because processing stresses requires an additional function before associating stresses and centroids, while the nodes can be directly associated to their deformation values. This is the case because the centroid coordinates need to be calculated and to do this each element must be associated to all the labels of nodes making it up.
Because the functions to associates stresses to centroids and deformations to nodes are very similar, the stress one will be discussed first and then the deflection one will be described mainly through its difference with the aforementioned function.
#### elems_processor
This function starts by reprocessing the nset_lines list from dat_processor, since it contains 
```python
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
```
The necessary paths where to save the matrices containing all the nodes for each element are then defined together with the headers to use, if desired, for the .csv files. Once again the function checks whether the last file to be generated already exists or not, and if the file processing procedure has to be carried out, the .dat file is opened.
Once again, the file is scanned by intervals between the items of elem_lines and nset_lines and every line that is read in them is stripped and its doubles extracted and stored in an array. 
```python
        with open(dat_path, "r") as read_dat:
            dat_file = read_dat.readlines()
            for k in range(len(new_elem_lines)):
                # Only focus on the areas of the file where information about the elements is present
                elem_vector = np.zeros((new_nset_lines[k] - new_elem_lines[k], 9))
                elem_vector_reset = np.zeros((new_nset_lines[k] - new_elem_lines[k], 9))
                # Iterate over each of the lines in the interval found and scan for float values
                for line in range(new_elem_lines[k], new_nset_lines[k]):
                    s = dat_file[line].strip()
                    data = [float(s) for s in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)]
                    data = np.array([float(x) for x in data])
```
Additional operations are carried out on the data array to filter for unwanted values (such as, again, the counter every 5 lines of the Earth.dat file) or numerical errors in the .dat file such as a couple of rows containing all the same value. The data to process is saved in two arrays, one that resets its counter for each part and one that doesn't, called elem_vector and elem_vector_reset respectively. 
```python
                    if len(data) > 6:
                        if len(data) == 10:
                            data = data[1:]
                        elif len(data) == 8:
                            data = data[1:]
                        if len(data) == 9:
                            elem_vector[line - new_elem_lines[k] + 1, :] = data
                            elem_vector_reset[line - new_elem_lines[k] + 1, :] = data
                        else:
                            elem_vector[line - new_elem_lines[k] + 1, 0:7] = data
                            elem_vector_reset[line - new_elem_lines[k] + 1, 0:7] = data

                        # Filter out possible numerical errors
                        if elem_vector[line - new_elem_lines[k] + 1, 1] == elem_vector[line - new_elem_lines[k] + 1, 2]:
                            elem_vector[line - new_elem_lines[k] + 1, :] = 0
                        if elem_vector_reset[line - new_elem_lines[k] + 1, 1] == \
                                elem_vector_reset[line - new_elem_lines[k] + 1, 2]:
                            elem_vector_reset[line - new_elem_lines[k] + 1, :] = 0

                # Take away all zero rows from the element matrix relative to the relevant part we are processing and
                # create labels for the elements of each part so that they reset for each part in one of the large
                # matrices that we will create later.
                elem_vector = elem_vector[~np.all(elem_vector == 0, axis=1)]
                elem_vector_reset = elem_vector_reset[~np.all(elem_vector_reset == 0, axis=1)]
                reset_labels = range(1, len(elem_vector_reset) + 1)
                reset_labels = np.transpose(reset_labels)
                elem_vector_reset[:, 0] = reset_labels
                elem_vector[:, 0] = reset_labels
```
The matrix whose element labels reset for each part is saved as .csv with or without headers based on headers_input. This matrix is appended to an unsorted large matrix whose labels reset with each part, while elem_vector is appended to a large matrix whose labels do not reset with each part. 
```python
                large_elem_matrix_unsorted_reset.append(individual_elem_matrix)
                large_elem_matrix_unsorted.append(elem_vector)
```
This latter matrix is then sorted through a dictionary and stored in another variable called sorted_matrix. It has been decided to produce these three different outputs for future use, since they might be needed for other kinds of processing. After saving all the outputs, the next function starts. This is the crucial function of the process where the association between stress values and element centroids is made. 
#### associate_stress_coord
The function starts by reading all the .csv files containing stress values and element nodes for each part.
```python
    file_extension = '.csv'
    individual_stress_paths = stress_matrices_path
    dir_csv_elem_files = os.listdir(individual_element_paths)
    dir_csv_stress_files = os.listdir(individual_stress_paths)
    csv_elem_files = [filename for filename in dir_csv_elem_files if filename.endswith(file_extension)]
    csv_stress_files = [filename for filename in dir_csv_stress_files if filename.endswith(file_extension)]
```
New paths are defined where to store the necessary values and the headers for the .csv files to be generated are also defined.
```python
    centroid_files_path = os.path.join(stress_part_values, 'Centroids')
    if not os.path.exists(centroid_files_path):
        os.mkdir(centroid_files_path)
    complete_files_path = os.path.join(stress_part_values, 'Complete_Files')
    if not os.path.exists(complete_files_path):
        os.mkdir(complete_files_path)
    geographical_centroids_files_path = os.path.join(stress_part_values, 'Geographical_Centroids')
    if not os.path.exists(geographical_centroids_files_path):
        os.mkdir(geographical_centroids_files_path)
    geographical_components_files_path = os.path.join(stress_part_values, 'Geographical_Components')
    if not os.path.exists(geographical_components_files_path):
        os.mkdir(geographical_components_files_path)
    geographical_complete_files_path = os.path.join(stress_part_values, 'Geographical_Complete_Files')
    if not os.path.exists(geographical_complete_files_path):
        os.mkdir(geographical_complete_files_path)

    # Define headers for the complete files and the centroid files
    complete_headers = ['Label', 'S_Mises', 'S_11', 'S_22', 'S_33', 'S_12', 'S_13', 'S_23', 'X_centroid', 'Y_centroid',
                        'Z_centroid', 'R', 'Depth', 'Lat', 'Lon']
    centroid_headers = ['Centr_label', 'X_centroid', 'Y_centroid', 'Z_centroid']
```
Another check is performed to assess whether the operations of this function have to be carried or not based on the existence of the last file to be generated. If this file does not exist, the  association operations start being carried out. The first one is to filter out all the files that do not need to be processed (in this case, any of them that are not for part EARTH).
```python
            for file_to_evaluate in range(len(csv_stress_files)):
                if 'Region' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'LOW' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if 'I0' in csv_stress_files[file_to_evaluate]:
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
```
Then, for each of the .csv files containing stress values that are left to process, the corresponding element file is found. 
```python
            part_matrix_to_open = open(os.path.join(individual_stress_paths, part_file))
            part_matrix = part_matrix_to_open.readlines()

            file_to_read_logical = [i.split('Elements_Part_')[1] == part_file for i in csv_elem_files]
            file_to_read_index = np.array([file_to_read_logical.index(i) for i in file_to_read_logical
                                           if i == True])
            # Open the element file and transform both the element and stress files into arrays for easier handling
            file_to_read_index = file_to_read_index[0]
```
Both of the files are then opened and transformed into matrices for easier data manipulation.
```python
            all_elems_to_open = open(os.path.join(individual_element_paths, csv_elem_files[file_to_read_index]))
            all_elems = all_elems_to_open.readlines()
            all_elems = [line.strip() for line in all_elems[1:]]
            all_elems = [np.array([eval(i) for i in line.split(",")[:]]) for line in all_elems]
            all_elems = np.array(all_elems)
            part_matrix = [line.strip() for line in part_matrix[1:]]
            part_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in part_matrix]
            part_matrix = np.array(part_matrix)
            centroid_coord = np.zeros((len(part_matrix), 4))
```
The function then iterates over each line of the stress matrix and for each of them finds the index of the element matrix with the same element label as the one from the stress matrix. Then, all the nodes making up that element are read, a matrix of zeros with as many rows as the number of nodes and 3 columns called elem_coord_matrix is allocated and each of them is searched in the large node matrix opened at the beginning of the function. Their coordinates are stored in elem_coord_matrix and then the mean is performed over the rows to obtain the X, Y and Z coordinates of the centroids.
```python
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
```
The centroid coordinates need to be swapped around, together with the stress tensor components, because ABAQUS uses a Z, X, Y coordinate system and it has to be transformed to a X, Y, Z system.
```python
            centroid_coord[:, 1], centroid_coord[:, 3] = centroid_coord[:, 3].copy(), centroid_coord[:, 1].copy()
            centroid_coord[:, 2], centroid_coord[:, 3] = centroid_coord[:, 3].copy(), centroid_coord[:, 2].copy()
            part_matrix[:, 2], part_matrix[:, 3] = part_matrix[:, 3].copy(), part_matrix[:, 2].copy()
            part_matrix[:, 3], part_matrix[:, 4] = part_matrix[:, 4].copy(), part_matrix[:, 3].copy()
            part_matrix[:, 6], part_matrix[:, 7] = part_matrix[:, 7].copy(), part_matrix[:, 6].copy()
            part_matrix[:, 5], part_matrix[:, 6] = part_matrix[:, 6].copy(), part_matrix[:, 5].copy()
```
The next part computes the radial distance and depth of each centroid and then uses the centroid coordinates to calculate their latitude and longitude with the cart2geo function. These variables are all appended to the stress and centroid coordinates. 
```python
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
```
---
##### cart2geo
This function uses the cartesian to geographical angle relations to derive latitude and longitude for all centroids starting from the cartesian coordinate vectors.
```python
    X = cartesian_coordinates[:, 0]
    Y = cartesian_coordinates[:, 1]
    Z = cartesian_coordinates[:, 2]
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    lat = np.arcsin(Z/R) * 180 / np.pi
    lon = np.arctan2(Y, X) * 180 / np.pi
```
---
After this, the centroid coordinate and complete matrix for each part is saved with or without headers. The next step is the conversion of all the stress values to the geographical coordinate system, performed in rotate_tensor.

---
##### rotate_tensor
The function starts by allocating the matrix where to save the new stresses and its headers. The complete file is then opened and transformed into a matrix and an iteration starts over each of its lines. The stresses are read from the line and organized in a tensor, then the transformation tensor is also defined and filled up with the transformation cosines and the matrix operations for the reference system change is then carried out. The result is a tensor whose values are stored in the corresponding line of the matrix allocated previously.

```python
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
            R_3D = np.sqrt(centroid_coords[0] ** 2 + centroid_coords[1] ** 2 + centroid_coords[2] ** 2)
            R_azimuth = np.sqrt(centroid_coords[0] ** 2 + centroid_coords[1] ** 2)
            cos_lon = centroid_coords[0] / R_azimuth
            sin_lon = centroid_coords[1] / R_azimuth
            cos_lat = R_azimuth / R_3D
            sin_lat = centroid_coords[2] / R_3D
            T = np.zeros((3, 3))

            T[0, 0] = cos_lon * cos_lat
            T[1, 0] = cos_lon * sin_lat
            T[2, 0] = -sin_lon
            T[0, 1] = sin_lon * cos_lat
            T[1, 1] = sin_lon * sin_lat
            T[2, 1] = cos_lon
            T[0, 2] = sin_lat
            T[1, 2] = -cos_lat
            T[2, 2] = 0.0

            # Save the geographical components for each element in a line of the file allocated previously
            S_geographical = T.dot(S).dot(np.transpose(T))
            geographical_components[line, 0] = complete_file[line, 1]
            geographical_components[line, 1] = S_geographical[0, 0]
            geographical_components[line, 2] = S_geographical[1, 1]
            geographical_components[line, 3] = S_geographical[2, 2]
            geographical_components[line, 4] = S_geographical[0, 1]
            geographical_components[line, 5] = S_geographical[0, 2]
            geographical_components[line, 6] = S_geographical[1, 2]
```
These new stresses are concatenated to the corresponding labels, radial distance, depth, latitude and longitude coming from the complete file and the new matrix is saved with or without headers for the current part. The function then terminates.

---
The for cycle over each stress matrix is then exited and the last operation to carry out is the coupling of coupled stresses for the specified iteration, step and stress iteration with the corresponding depth etc., performed at the end of the script for part EARTH only. This has been necessary during the thesis work for analysis results and can be tweaked for more parts or be left out. This new matrix with all the values is then also saved in the same folder as the complete file for the other analyzed parts. The function then terminates.

#### associate_deflections_coord
This function, as said before, works in a very similar way to associate_stress_coord with the main difference that the nodes can directly be associated to the corresponding deformation values. Again, the necessary paths where to save the complete files and their headers are defined first and then the filenames are filtered to avoid processing parts of no interest and then an if cycle determines whether to run the rest of the function or not based on the existence of the last file to be generated. If the function runs, a for cycle starts iterating over the remaining files to process and opens each deformation and node coordinate file as a matrix. Their columns are swapped to again take into account ABAQUS' different cartesian reference system. 
```python
            with open(os.path.join(deflection_path, csv_deflection_files[csv_file]), 'r') as deflection_read_obj:
                individual_deflection_matrix = deflection_read_obj.readlines()
                individual_deflection_matrix = [line.strip() for line in individual_deflection_matrix[1:]]
                individual_deflection_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                                individual_deflection_matrix]
                individual_deflection_matrix = np.array(individual_deflection_matrix)
                individual_deflection_matrix[:, 2], individual_deflection_matrix[:, 3] = \
                    individual_deflection_matrix[:, 3], individual_deflection_matrix[:, 2].copy()
                individual_deflection_matrix[:, 3], individual_deflection_matrix[:, 4] = \
                    individual_deflection_matrix[:, 4], individual_deflection_matrix[:, 3].copy()
            with open(os.path.join(individual_node_path, 'Nodes', 'Nodes_Part_' + csv_deflection_files[csv_file]), 'r') as \
                    coordinate_read_obj:
                individual_coordinate_matrix = coordinate_read_obj.readlines()
                individual_coordinate_matrix = [line.strip() for line in individual_coordinate_matrix[1:]]
                individual_coordinate_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in
                                                individual_coordinate_matrix]
                # Swap coordinates
                individual_coordinate_matrix = np.array(individual_coordinate_matrix)
                individual_coordinate_matrix[:, 1], individual_coordinate_matrix[:, 2] = \
                    individual_coordinate_matrix[:, 2], individual_coordinate_matrix[:, 1].copy()
                individual_coordinate_matrix[:, 2], individual_coordinate_matrix[:, 3] = \
                    individual_coordinate_matrix[:, 3], individual_coordinate_matrix[:, 2].copy()
```
After this, the deformations are transformed into the new geographical reference system with rotate_deflections to also obtain a complete matrix in this new reference system.
```python
            complete_deflection_matrix = np.column_stack((individual_deflection_matrix,
                                                          individual_coordinate_matrix[:, 1:]))
            geographical_complete_matrix = rotate_deflections(complete_deflection_matrix)
```
--- 
##### rotate_deflections
Similarly to the rotate_tensor function, the function starts by opening the current complete file and extracting its displacements and node coordinates, then two matrices are allocated to store the displacement components and magnitude respectively in geographical coordinates.
```python
    displacements = complete_deflection_matrix[:, 2:5]
    node_coordinates = complete_deflection_matrix[:, 5:8]
    geographical_displacements = np.zeros((len(displacements), len(displacements[0])))
    U_magn = np.zeros((len(displacements), 1))
```
Again, the transformation cosines are defined and used to fill the transformation tensor, then this tensor is used for the reference transformation of the displacements which are in turn used to calculate the new displacement magnitude. The complete matrix in geographical coordinates is then created by concatenating displacement magnitude, components, node coordinates and is returned to the associate_deflection_coordinates function.

--- 
Radial distance, depth, latitude and longitude are again calculated for both Cartesian and geographical complete matrices and then concatenated to them.

```python
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
```
Both matrices are then saved with or without headers depending on the user's input. The function then terminates and coordinate_reader is also exited.

---
Back to the main function, the time to run coordinate reader is displayed and the user is then prompted to enter the reference system of the quantities to plot. Once this input has been received, another counter starts for the function call to depth_classifier. 

---
#### depth_classifier

---
Once this has been done, the time to run depth_classifier is displayed and the latitude, longitude and depth ranges to use for successive plots are defined. Then, another function called element_tracker is called. Its main importance resides in the generation of tables tracking the element with highest stress in the ranges defined above to check with more detail what is happening to the model.

---
#### element_tracker
This function creates a .csv table containing the stress and viscosity value for the element with the highest Mises stress in the latitude, longitude and depth ranges initialized above. The function begins by extracting the extremes of these ranges and the table headers, then the user is asked whether the function has to be run or not. It is recommended to do this only after processing a sufficiently large amount of stress iterations. 
```python
    fles_to_check = []
    headers = ['Label', 'Mises', 'S11', 'S22', 'S33', 'S12', 'S13', 'S23', 'Viscosity']
    cycle_dir_counter = 0
    min_lat = area_params[0]
    max_lat = area_params[1]
    min_lon = area_params[2]
    max_lon = area_params[3]
    depth_range = area_params[4] * 1000
    min_depth_name = str(depth_range[0]/1000)
    max_depth_name = str(depth_range[1]/1000)
    table_bool = input('Enter 1 to generate the element tracking matrix for te selected latitude, longitude and depth '
                       'ranges, 0 to skip this action: ')
```
If the function has to run, every stress iteration folder is scanned for the complete stress and complete coupled stress matrices generated in associate_stress_coord. All the filenames that are found are saved in a list initialized previously and called files_to_check.
```python
        if not os.path.exists(os.path.join(base_path, 'Element_tracking_table' + min_depth_name + '_' + max_depth_name +
                                                      '.csv')):
            for iter_dir in os.listdir(base_path):
                iter_path = os.path.join(base_path, iter_dir)
                if os.path.isdir(iter_path) and 'difference' not in iter_path:
                    for step_dir in os.listdir(iter_path):
                        step_path = os.path.join(iter_path, step_dir)
                        for cycle_dir in os.listdir(step_path):
                            cycle_dir_counter = cycle_dir_counter + 1
                            stress_path = os.path.join(step_path, cycle_dir, 'stress_processing_results',
                                                       'Stress_Part_Values', 'Complete_Files')
                            for f in os.listdir(stress_path):
                                if f == 'Complete_file_EARTH.csv':
                                    files_to_check.append(os.path.join(stress_path, f))
                                elif f == 'Complete_file_coupled_EARTH.csv':
                                    files_to_check.append(os.path.join(stress_path, f))
                                else:
                                    pass
```
Two variables are then initialized, one which is the table where to store the values called element_table and one which is a series of numbers used as row labels made by a trio of numbers indicating the iteration, step and stress iteration. The function loops over every file in files_to_check and opens it as a matrix, then extracts only the rows in that file which satisfy the latitude, longitude and depth range conditions, storing them in a variable called matrix.
```python
            element_table = np.zeros((len(files_to_check), 9))
            iter_numbers = np.zeros((len(files_to_check), 3))
            count = 0
            for name in files_to_check:
                e_path = os.path.join(base_path, 'e.dat')
                with open(name, 'r') as no_arr_matrix:
                    large_matrix = no_arr_matrix.readlines()
                    large_matrix = [line.strip() for line in large_matrix[1:]]
                    large_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in large_matrix]
                    large_matrix = np.array(large_matrix)
                    # Add here the part where only the points in a certain depth and latlon range are considered
                    # and saved in a matrix called matrix
                    depth = large_matrix[:, -3]
                    lat = large_matrix[:, -2]
                    lon = large_matrix[:, -1]
                    min_depth = depth_range[0]
                    max_depth = depth_range[1]
                    depth_condition = (depth > min_depth) & (depth < max_depth)
                    lat_condition = (lat > min_lat) & (lat < max_lat)
                    lon_condition = (lon > min_lon) & (lon < max_lon)
                    matrix = large_matrix[depth_condition * lat_condition * lon_condition]
```
The maximum Mises stress in this matrix is then read and its index extracted, then the e.dat file is opened as a matrix too and the index is used to find the B_diff and B_disl value of that element to calculate the corresponding viscosity. The maximum Mises index is then used to also extract the label and stresses for that element which are then concatenated with the viscosity in a vector called value_vector used to fill one row of element_table. 
```python
                    mises = matrix[:, 1]
                    max_stress = max(mises)
                    max_index = np.where(mises == max_stress)
                    mises = matrix[max_index, 1]
                    with open(e_path, 'r') as no_arr_e:
                        e = no_arr_e.readlines()
                        e = [line.strip() for line in e]
                        e = [np.array([float(i) for i in line.split(" ")[:]]) for line in e]
                        e = np.array(e)
                        b_diff = e[:, 1]
                        b_disl = e[:, 2]
                    viscosity = mises/(b_diff[max_index] * mises + b_disl[max_index] * mises ** 3.5)
                    stresses = matrix[max_index]
                    stresses = stresses[:, 0:8]
                    value_vector = np.hstack((stresses, viscosity))
                    element_table[count, :] = value_vector
```
Finally, the row labels are defined. The iteration, step and stress iteration numbers are extracted from the base path where the stress matrices are stored and then an additional string has to be added to them, f.m. or t.m, meaning from and to model respectively. This gives a good distinction between the simple stresses (coming from the model at this stress iteration) and the coupled ones used for the same stress iteration but as an input for the generation of the results of the next stress iteration. Because of the order in which files are processed with the ones coming from the model coming after the ones used as input, the rows need to be swapped 2 by 2. This happens because coupled_stresses, the folder containing the coupled stresses matrices, is processed before the simple stresses outside of it in alphabetical order. 
```python
                data = [float(s) for s in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", name)]
                data = np.array([float(x) for x in data])
                data = data[1:]
                iter_numbers[count, :] = data
                count = count + 1
            if len(element_table) % 2 != 0:
                for i in range(0, len(element_table)-1, 2):
                    element_table[[i, i+1]] = element_table[[i+1, i]]
            else:
                for i in range(0, len(element_table), 2):
                    element_table[[i, i+1]] = element_table[[i+1, i]]
            indices = []
```
After swapping the rows, f.m. and t.m. are each assigned every 2 rows.
```python
            for i in range(len(iter_numbers)):
                if i % 2 != 0:
                    indices.append(str(int(iter_numbers[i, 0])) + ' ' + str(int(iter_numbers[i, 1])) + ' ' + str(
                        int(iter_numbers[i, 2])) + ' t.m.')
                else:
                    indices.append(str(int(iter_numbers[i, 0])) + ' ' + str(int(iter_numbers[i, 1])) + ' ' + str(
                        int(iter_numbers[i, 2])) + ' f.m.')
```
After the for cycle over every file to process is completed, element_table with its headers and this column of row headers called indices are joined in a single pandas data frame which is then saved. 
```python
            table = pd.DataFrame(data=element_table, index=indices, columns=headers)
            table.to_csv(os.path.join(base_path, 'Element_tracking_table_' + min_depth_name + '_' + max_depth_name +
                                      '.csv'))
```
---
The last function to be called from this main is called matlab_variables_writer, which generates a .txt file containing all the information that MATLAB needs to find the data to plot and what kind of data it is.

---
#### matlab_variables_writer
This last function to run writes all the handles that MATLAB need s to generate plots into a .txt file that is afterwards read by the MATLAB main function. This operation is very quick as every variable to write has already been defined somewhere in the main, so all the variables are simply concatenated in a long list whose elements are written one for each line after being converted to strings if necessary. This file is then closed and the function terminates. 

```python
    if sd_input == 0:
        quantity_to_plot = 'stresses'
    else:
        quantity_to_plot = 'deflections'
    second_list = [quantity_to_plot, iteration_path, classified_path, components_to_plot, complete_files_path]
    depth = area_params[-1]
    min_depth = depth[0]
    max_depth = depth[1]
    variable_list = files_to_classify + second_list + area_params[:-1] + [min_depth, max_depth]
    for i in range(len(variable_list)):
        if isinstance(variable_list[i], int):
            variable_list[i] = str(variable_list[i])

    print iteration_path
    with open(os.path.join(iteration_path, 'matlab_' + components_to_plot + '_' + quantity_to_plot + '_' + 'variables.txt'), 'wb') as f_write:
        f_write.write('\n'.join(variable_list))
```
---
The total time to run the main is then calculated and displayed and the Python process is finished. The plotting operations can start in MATLAB.