from dat_processor import *
from nodes_processor import *
from elems_processor import *
from associate_stress_coord import *
from associate_deflection_coord import *
import pdb
import time


def coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_path,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder):

    # Process the .dat files containing node coordinates and save a new copy, then calculate the time spent to run the
    # code. The time counting operation is done for all the other function calls too.
    start_time_dat = time.time()
    node_lines, elem_lines, nset_lines, file_identifiers, new_earth_dat = dat_processor(dat_path)
    end_time_dat = time.time() - start_time_dat
    # Time count
    if end_time_dat < 1:
        print('The time spent to reprocess the dat file to make it more readable is less than 1 second')
    else:
        print('The time spent to reprocess the dat file to make it more readable is ', str(int(round(end_time_dat))) +
              ' seconds')

    # Create files with nodes for each part and then a large node list for all parts, saving them in separate .csv
    start_time_nodes = time.time()
    individual_node_path, large_node_matrix_path = nodes_processor(common_files_path, node_lines, elem_lines,
                                                                   file_identifiers, dat_path, headers_on)
    end_time_nodes = time.time() - start_time_nodes
    # Time count
    if end_time_nodes < 1:
        print('The time spent to extract the node coordinates for each part is less than 1 second')
    else:
        print('The time spent to extract the node coordinates for each part is ', str(int(round(end_time_nodes))) +
              ' seconds')

    # Create a separate file for each part containing the elements associated to the corresponding nodes and then
    # associate them to the centroid coordinates that need to be calculated
    stress_part_values = stress_processing_path
    if sd_input == 0:
        start_time_elements = time.time()
        individual_element_paths = elems_processor(node_lines, elem_lines, nset_lines, stress_processing_path,
                                                   stress_part_values, dat_path, file_identifiers, headers_on)
        end_time_elements = time.time() - start_time_elements
        # Time count
        if end_time_elements < 1:
            print('The time spent to extract the element nodes for each part is less than 1 second')
        else:
            print('The time spent to extract the element nodes for each part is ', str(int(round(end_time_elements))) +
                  ' seconds')
        start_time_complete_stresses = time.time()
        centroid_files_path, complete_files_path, geographical_complete_files_path = \
            associate_stress_coord(individual_element_paths, stress_part_values, large_node_matrix_path, headers_on,
                                   stress_matrices_path, dat_name, coupled_stress_folder)
        end_time_complete_stresses = time.time() - start_time_complete_stresses
        # Time count
        if end_time_complete_stresses < 1:
            print('The time spent to associate the stresses to each corresponding element for each part is less than ' \
                  '1 second')
        else:
            print('The time spent to extract the element nodes for each part is ' +
                  str(int(round(end_time_complete_stresses))),  ' seconds')

    # For the deformations, directly associate them to the node coordinates as nothing else needs to be recalculated
    # since the nodes are directly given by the model and the centroids aren't
    else:
        start_time_complete_deflections = time.time()
        complete_files_path, geographical_complete_files_path = associate_deflection_coord(deflection_path,
                                                                                           individual_node_path,
                                                                                           headers_on,
                                                                                           deflection_processing_path)
        end_time_complete_deflections = time.time() - start_time_complete_deflections
        # Time count
        if end_time_complete_deflections < 1:
            print('The time spent to create the complete deflection files for each part and region is less than ' +
                  '1 second')
        else:
            print('The time spent to create the complete deflection files for each part and region is ' +
                  str(end_time_complete_deflections) + ' seconds')

    return stress_part_values, complete_files_path, geographical_complete_files_path, large_node_matrix_path
