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
dat_path = os.path.join(base_path, dat_name)
base_e_path = os.path.join(base_path, e_dat_name)

iteration_path = os.path.join(base_path, 'Iteration_' + str(iteration))
if not os.path.exists(iteration_path):
    os.mkdir(iteration_path)
# Define the paths of the stress and coordinate files and if they are not found create these folders
stress_report_path = os.path.join(iteration_path, stress_rpt_name)
deflection_report_path = os.path.join(iteration_path, deflection_rpt_name)
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
# file_identifiers, headers = stress_report_reader(stress_report_path, stress_search_key, stress_matrices_path,
#                                                  headers_on)
deflection_path = deflection_report_reader(deflection_report_path, deflection_search_key, deflection_matrices_paths,
                                           headers_on)
end_time_reports = time.time() - start_time_reports
print 'The time spent to read the stress and deflection reports and extract data from them is: ', end_time_reports, \
    ' seconds'
new_e_iter_path = os.path.join(iteration_path, e_dat_name)
old_e_iter_path = os.path.join(os.path.join(base_path, 'Iteration_' + str(iteration - 1)), e_dat_name)

# new_stress_test = 1000000
# if not os.path.exists(coupled_stress_folder):
#     os.mkdir(coupled_stress_folder)
# large_coupled_matrix = stress_coupling(stress_matrices_path, coupled_stress_folder, headers_on, new_stress_test,
#                                        headers, iteration, base_path, base_e_path, old_e_iter_path, new_e_iter_path)
# stress_matrices_path = coupled_stress_folder
