# Data processing 
## Main file - python_reader_processing_only

This is the main file of the data processing step and begins with the definition of a series of parameters and paths used to locate the result files to analyze and to store the files that will be successively generated if they do not exist already. This script receives as input the base result folder coming from ABAQUS, which is also on GitHub, and gives a series of files where each element and node is associated to its stress or deformation components respectively. This is only done for part Earth, but there is the possibility of expanding the code to also process the others if necessary. 

The first function to be called is called coordinate_reader and includes different functions calls inside of it.

---
## Coordinate_reader and its functions 
Coordinate_reader calls a series of functions to read the data stored in the sample folder and associate stresses and deformations to their respective elements or nodes and their coordinates, both Cartesian and geographical. All of the functions called in this one are timed in order to check the system performance as some operations, especially the calculation of element centroids, can take a very long time. Its inputs are:
- sd_input, bool deciding whether to plot stresses or deflections;
- dat_path, path to the Earth.dat file;
- stress_matrices_path, path for the .csv files containing stress values;
- stress_processing_path, path where to save the stress data after stress values processing;
- deflection_path, path for the .csv files containing deformation values;
- deflection_processing_path, path where to save the deformation data after deformation values processing;
- headers_on, bool to decide whether to generate .csv files with or without headers;
- dat_name, name of the .dat file with model information;
- common_files_path, path where to save files used for both stress and deformation processing;
- coupled_stress_folder, path where to find the coupled stress values (in ABAQUS used as model input).
Its outputs are:
- complete_files_path, path where to find the file containing stress and deformations with the respective node or centroid coordinates in Cartesian coordinates;
- geographical_complete_files_path, same as before but in geographical stress or deformation components;
- large_node_matrix_path, path for the node list for the entire model.
The first function called here is dat_processor. 

### dat_processor

This function uses as input the Earth.dat file where the information for all the node coordinates and element for all parts is stored. Its inputs are:
- fname, path of the .dat file to open.

Its outputs are:
- node_lines, line indices in the .dat file where the definition of node coordinates for each part begins;
- elem_lines, line indices in the .dat file where the definition of elements with the nodes making them up begins and the node definition ends;
- nset_lines, line indices where the element definition for each part ends;
- file_identifiers, list of part names extracted from the .dat file;
- new_earth_dat, new .dat file path which is easier to process for the next scripts.

### nodes_processor

This function extracts node coordinates and labels for each part from the new Earth.dat file and saves them in a large node list for the entire model and a separate node list for every part. Its inputs are:
- common_files_path;
- node_lines;
- elem_lines;
- file_identifiers;
- dat_path;
- headers_on.

The outputs are instead:
- individual_path, path for the files to be read by both stress and deformations divided into each part.
- large_node_matrix_path.

Because the functions to associates stresses to centroids and deformations to nodes are very similar, the stress one will be discussed first and then the deflection one will be described mainly through its difference with the aforementioned function.

### elems_processor
This function associates each element to the nodes making it up, an additional step compared to node processing because element centroids need to be calculated as they are not given in an output file like for the node coordinates. 
The function inputs are: 
- node_lines;
- elem_lines;
- nset_lines;
- stress_processing_path;
- dat_path;
- file_identifiers;
- headers_on.

The function outputs are:
individual_element_paths, path where to save the elements associated to all the nodes making them up.

### associate_stress_coord
This function associates the stress components to the corresponding element centroids. It also calls cart2geo and rotate_tensor to convert the centroid Cartesian coordinates to geographical and rotate_tensor to convert the stress tensor components to a geographical reference system. Its inputs are:
- individual_element_paths; 
- stress_part_values; 
- large_node_matrix_path; 
- headers_on; 
- stress_matrices_path;
- coupled_stress_folder.

Its outputs are:
- centroid_files_path, path where to save every element label associated to its centroid coordinates;
- complete_files_path;
- geographical_complete_files_path.

#### cart2geo
This function uses the Cartesian to geographical angle relations to derive latitude and longitude for all centroids starting from the Cartesian coordinate vectors.  Its inputs are:
- cartesian_coordinates, the Cartesian coordinate vectors for the selected part.
Its outputs are:
- lat, lon, the latitude and longitude vectors corresponding to the Cartesian coordinates. 

#### rotate_tensor
This function performs the coordinate system transformation from Cartesian to geographical on the stress tensor components. Its inputs are:
- part_matrix, the matrix for the model part being analyzed containing stresses or deformations associated to centroids or nodes;
- geographical_complete_files_path, the file path where to save the matrix with transformed stress components;
- headers_on;
- complete_headers, headers for the matrix with stresses associated to centroids;
- complete_individual_path, path to the matrix with components to convert;
- part_file, the part name.
Its outputs are:
- complete_file_geographical, the complete matrix with geographical components;
- geo_coordinates, the radial distance, depth, latitude and longitude for each element.

### associate_deflection_coord

This function, as said before, works in a very similar way to associate_stress_coord with the main difference that the nodes can directly be associated to the corresponding deformation values. As it can be expected, this function associates every node label to its deformation values and both Cartesian and geographical coordinates. 
Similarly to the rotate_tensor function, this one performs the coordinate system transformation on the deformation components. Its inputs are: 
- deflection_path;
- individual_node_path, path to the node matrix for the current part;
- headers_on;
- deflection_processing_path.

Its outputs are:
- complete_deflection_path, path where to save the Cartesian complete deflection matrices with deformations associated to their nodes;
- geographical_complete_deflection_path, path where to save the previous matrix but in a geographical reference system.

#### rotate_deflections
Similarly to rotate_tensor, this function converts deformations to the geographical coordinate system. Its inputs are:
- complete_deflection_matrix, the matrix with stresses and deformations in Cartesian coordinates;
Its outputs are:
- geographical_complete_matrix, the matrix with stresses and deformations in geographical coordinates.

---
### depth_classifier
This function is not a necessary part of the data processing workflow anymore but it is useful as it can give an idea of the depth ranges which actually contain points for plotting. The function bins the stress and deformation values in a series of depth ranges which are determined based on a number of bins entered by the user. Its inputs are: 
- sd_input;
- deflection_processing_path;
- individual_path;
- components_to_plot, 
- headers_on;
- file_extension, defining the .csv file to save and search for;
- Cartesian_classified_depth_path, path where to save the complete data binned by depth in a Cartesian reference frame;
- geographical_classified_depth_path, same as the former variable but in a geographical reference frame;
- files_to_classify_path, path for the complete data files;
- files_to_classify, list of files in the path defined by the previous variable;
- histogram_path, path where to save the histograms for the point distribution by depth.

Its outputs are histograms with point distributions by depth and .csv files containing the stresses or deformations binned by depth.

---
### element_tracker
This function creates a .csv table containing the stress and viscosity value for the element with the highest Mises stress in the latitude, longitude and depth ranges initialized above. This is done for every stress iteration and for both the model inputs and outputs. The function inputs are: 
- components_to_plot, the reference system of the data to process;
- base_path, path for all the data for the simulation run;
- area_params, AOI bounds in lat, lon and depth ranges.

The function output is a table with the evolution of the element with highest Mises stress with time.


---
### matlab_variables_writer
This last function to run writes all the handles that MATLAB needs to generate plots  into a .txt file that is afterwards read by the MATLAB main function. This operation is very quick as every variable to write has already been defined somewhere in the main, so all the variables are simply concatenated in a long list whose elements are written one for each line after being converted to strings if necessary. This file is then closed and the function terminates. The function inputs are:
- sd_input; 
- complete_files_path;
- files_to_classify; 
- iteration_path, the folder relative to the sea level iteration; 
- components_to_plot.
The function output is a handle file to pass to MATLAB to locate all necessary files for plotting.

This function has been removed in the Python 3 workflow as plotting is also being ported to Python.
The total time to run the main is then calculated and displayed and the Python process is finished. The plotting operations can start in MATLAB.

## Additional scripts
### mises_verifier
This scripts calculates the Mises stress with the formula given in the thesis report in order to verify that it can be used as an approximation of how ABAQUS calculates the Mises stress. 