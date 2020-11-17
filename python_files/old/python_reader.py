import os
import platform
import time
import pdb
import shutil
import numpy as np
from stress_report_reader import *
from deflection_report_reader import *
from first_iter_actions import *
from coordinate_reader import *
from depth_classifier import *
from stress_coupling import *
from matlab_variables_writer import *

# General parameters: based on how the files have been produced in ABAQUS, create the name of the dat and rpt files by
# modifying the run number
program_start_time = time.time()
run = '6'
if not isinstance(run, str):
    run = str(run)
stress_rpt_name = 'abaqus_run_' + run + '.rpt'
deflection_rpt_name = 'abaqus_run_' + run + '_deflections.rpt'
# e_csv_name = 'e.csv'
if run == 'point':
    dat_name = 'Earth_point.dat'
else:
    dat_name = 'Earth.dat'
iteration = 1
e_dat_name = 'e.dat'
# Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
# windows or linux. This is the only path that has to be modified
opersys = platform.system()
if opersys == 'Windows':
    base_path = os.path.join('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
else:
    base_path = os.path.join('/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run))
# Define the path of the dat and rpt files to read
stress_report_path = os.path.join(base_path, stress_rpt_name)
deflection_report_path = os.path.join(base_path, deflection_rpt_name)
dat_path = os.path.join(base_path, dat_name)
base_e_path = os.path.join(base_path, e_dat_name)

iteration_path = os.path.join(base_path, 'Iteration_' + str(iteration))
if not os.path.exists(iteration_path):
    os.mkdir(iteration_path)
# Define the paths of the stress and coordinate files and if they are not found create these folders
stress_matrices_path = os.path.join(iteration_path, 'stress_matrices')
stress_processing_path = os.path.join(iteration_path, 'stress_processing_results')
deflection_matrices_paths = os.path.join(iteration_path, 'deflection_matrices')
deflection_processing_path = os.path.join(iteration_path, 'deflection_processing_results')
common_files_path = os.path.join(iteration_path, 'common_files')
coupled_stress_folder = os.path.join(stress_matrices_path, 'coupled_stress_matrices')

if not os.path.exists(stress_matrices_path):
    os.mkdir(stress_matrices_path)
if not os.path.exists(stress_processing_path):
    os.mkdir(stress_processing_path)
if not os.path.exists(deflection_matrices_paths):
    os.mkdir(deflection_matrices_paths)
if not os.path.exists(deflection_processing_path):
    os.mkdir(deflection_processing_path)

# Call the functions used to read data from the deformation and stress report files and count their running time
headers_on = 1
stress_search_key = 'reported at element'
deflection_search_key = 'reported at nodes'
start_time_reports = time.time()
file_identifiers, headers = stress_report_reader(stress_report_path, stress_search_key, stress_matrices_path,
                                                 headers_on)
deflection_path = deflection_report_reader(deflection_report_path, deflection_search_key, deflection_matrices_paths,
                                           headers_on)
end_time_reports = time.time() - start_time_reports
print 'The time spent to read the stress and deflection reports and extract data from them is: ', end_time_reports, \
    ' seconds'
new_e_iter_path = os.path.join(iteration_path, e_dat_name)
old_e_iter_path = os.path.join(os.path.join(base_path, 'Iteration_' + str(iteration - 1)), e_dat_name)

new_stress_test = 500000
if not os.path.exists(coupled_stress_folder):
    os.mkdir(coupled_stress_folder)
large_coupled_matrix = stress_coupling(stress_matrices_path, coupled_stress_folder, headers_on, new_stress_test,
                                       headers, iteration, base_path, base_e_path, old_e_iter_path, new_e_iter_path)
stress_matrices_path = coupled_stress_folder

# Decide whether to save the stresses or deflections
# sd_input = 1
check_0 = 1
processing_input = 1
# while check_0 == 1:
#     processing_input = input('Do you wish to process data for plotting? Enter 1 for yes or 0 for no: \n')
#     if processing_input == 1 or processing_input == 0:
#         check_0 = 0
#     else:
#         print('incorrect input, select either yes or no. \n')

if processing_input == 0:
    print 'No files have been generated. The Mises stress has been modified and saved. Going to the next iteration.'
else:
    check_1 = 1
    sd_input = 1
    # while check_1 == 1:
    #     sd_input = input('Enter 0 to discretize stress components based on depth, 1 for the deflections: \n')
    #     if sd_input == 0 or sd_input == 1:
    #         check_1 = 0
    #     else:
    #         print('Incorrect input, select either 1 or 0. \n')

    # Call the function to create the complete stresses or deflection files and count running time
    start_time_coordinate_reader = time.time()
    individual_path, complete_files_path, spherical_complete_files_path, large_node_matrix_path = \
        coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_path,
                          deflection_processing_path, headers_on, dat_name, common_files_path)
    end_time_coordinate_reader = time.time() - start_time_coordinate_reader

    if sd_input == 0:
        print 'The time spent to create the complete files containing the stress components and centroid coordinates ' \
              'is: ', end_time_coordinate_reader, ' seconds'
    else:
        print 'The time spent to create the complete files containing the deflection components and node coordinates ' \
              'is: ', end_time_coordinate_reader, ' seconds'

    # Discretize the data into different depths based on whether we want to work with spherical or cartesian components
    check_2 = 1
    components_to_plot = 'spherical'
    # while check_2 == 1:
    #     components_to_plot = input('Enter "cartesian" to work with cartesian components, "spherical" to work with '
    #                                'spherical ones: \n')
    #     if components_to_plot == 'cartesian' or components_to_plot == 'spherical':
    #         check_2 = 0
    #     else:
    #         print('incorrect input, select either "cartesian" or "spherical". \n')

    classified_path, files_to_classify, maximum_depth, fig_counter = \
        depth_classifier(sd_input, deflection_processing_path, individual_path, components_to_plot,
                         complete_files_path, spherical_complete_files_path, headers_on)
    # Count the time spent by the entire program to run
    program_end_time = time.time() - program_start_time

    matlab_variables_writer(sd_input, files_to_classify, iteration_path, classified_path,
                            components_to_plot)
    #
    # python_plotter(sd_input, files_to_classify, base_path, iteration_path, classified_path, components_to_plot, run,
    #                coupled_stress_folder, e_dat_name, fig_counter)

    print 'The time spent to run the entire program is: ', end_time_coordinate_reader, ' seconds'
    print 'Data processing files generated. Going to the next iteration.'

