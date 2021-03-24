from dat_processor import *
from nodes_processor import *
from elems_processor import *
from associate_stress_coord import *
from associate_deflection_coord import *
import pdb
import time


def coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_path,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder):

    start_time_dat = time.time()
    node_lines, elem_lines, nset_lines, file_identifiers, new_earth_dat = dat_processor(dat_path)
    end_time_dat = time.time() - start_time_dat
    if end_time_dat < 1:
        print 'The time spent to reprocess the dat file to make it more readable is less than 1 second'
    else:
        print 'The time spent to reprocess the dat file to make it more readable is ', str(int(round(end_time_dat))), \
            ' seconds'
    start_time_nodes = time.time()
    individual_node_path, large_node_matrix_path = nodes_processor(common_files_path, node_lines, elem_lines,
                                                                   file_identifiers, dat_path, headers_on)
    end_time_nodes = time.time() - start_time_nodes
    if end_time_nodes < 1:
        print 'The time spent to extract the node coordinates for each part is less than 1 second'
    else:
        print 'The time spent to extract the node coordinates for each part is ', str(int(round(end_time_nodes))), \
            ' seconds'
    if sd_input == 0:
        start_time_elements = time.time()
        individual_element_paths = elems_processor(node_lines, elem_lines, nset_lines, stress_processing_path,
                                                   dat_path, file_identifiers, headers_on)
        end_time_elements = time.time() - start_time_elements
        if end_time_elements < 1:
            print 'The time spent to extract the element nodes for each part is less than 1 second'
        else:
            print 'The time spent to extract the element nodes for each part is ', str(int(round(end_time_elements))),\
                ' seconds'
        start_time_complete_stresses = time.time()
        centroid_files_path, complete_files_path, geographical_complete_files_path = \
            associate_stress_coord(individual_element_paths, stress_processing_path, large_node_matrix_path, headers_on,
                                   stress_matrices_path, coupled_stress_folder)
        end_time_complete_stresses = time.time() - start_time_complete_stresses
        if end_time_complete_stresses < 1:
            print 'The time spent to associate the stresses to each corresponding element for each part is less than ' \
                  '1 second'
        else:
            print 'The time spent to extract the element nodes for each part is ', \
                str(int(round(end_time_complete_stresses))),  ' seconds'

    else:
        start_time_complete_deflections = time.time()
        complete_files_path, geographical_complete_files_path = associate_deflection_coord(deflection_path,
                                                                                           individual_node_path,
                                                                                           headers_on,
                                                                                           deflection_processing_path)
        end_time_complete_deflections = time.time() - start_time_complete_deflections
        if end_time_complete_deflections < 1:
            print 'The time spent to create the complete deflection files for each part and region is less than ' \
                  '1 second'
        else:
            print 'The time spent to create the complete deflection files for each part and region is ', \
                end_time_complete_deflections, ' seconds'

    return stress_part_values, complete_files_path, geographical_complete_files_path, large_node_matrix_path
