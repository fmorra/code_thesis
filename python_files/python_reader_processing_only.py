# This is a new version of the main Python file, which needs to be used on the already read and processed report files
# to plot the results. All the data for each iteration is stored in an iteration folder in ABAQUS which has to be
# downloaded on the PC and then this script can be run on it.

# The setup for this file is still the same, with a folder for the run and folders from the different iterations. The
# main difference this time will be that the report and e files will be in the iteration folder, only the dat file will
# be out of it. We can do this because the large dat file is only used by us to find the node coordinates and element
# labels, which should not change based on the iteration: only the stress and deformation values in that file are going
# to vary (ask Bas or Wouter to make sure of this). This is convenient also because downloading the earth.dat file can
# take a really long time, so doing that if we wanted to plot each iteration would be really time-consuming. We could
# avoid that and use this script on ABAQUS if I knew how to call and modify functions in the ABAQUS PDE.

# Import relevant libraries
import os
import platform
import time
import pdb
from stress_report_reader import *
from deflection_report_reader import *
from coordinate_reader import *
from depth_classifier import *
from matlab_variables_writer import *
from element_tracker import *

# Define general parameters which are necessary to find the source file to extract data from to organize it. We need the
# name of the e.dat file and the Earth.dat file containing the node and element information for the entire model, and
# the run, iteration, step, cycle numbers together with the report filenames.
program_start_time = time.time()
run = '26'
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

# Choose whether to save files with headers or not and  choose what quantity to work with
headers_on = 1
check_1 = 1
sd_input = 1
# while check_1 == 1:
#     sd_input = input('Enter 0 to work with stress components, 1 for the deflections: \n')
#     if sd_input == 0 or sd_input == 1:
#         check_1 = 0
#     else:
#         print('Incorrect input, select either 1 or 0. \n')

# Call the function to create stress and deflection files containing the components together with the centroid or node
# coordinates, their depth and latitude and longitude
start_time_coordinate_reader = time.time()
individual_path, complete_files_path, geographical_complete_files_path, large_node_matrix_path = \
    coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_matrices_paths,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder)
end_time_coordinate_reader = time.time() - start_time_coordinate_reader
if sd_input == 0:
    print 'The time spent to create the complete files containing the stress components and centroid coordinates ' \
          'is ', str(int(round(end_time_coordinate_reader))), ' seconds'
else:
    print 'The time spent to create the complete files containing the deflection components and node coordinates ' \
          'is ', str(int(round(end_time_coordinate_reader))), ' seconds'
# Bin the data into different depth bins based on whether we want to work with geographical or cartesian components
check_2 = 1
components_to_plot = 1
while check_2 == 1:
    # components_to_plot = input('Enter 0 to work with cartesian components, 1 to work with '
    #                            'geographical ones: \n')
    if components_to_plot == 0 or components_to_plot == 1:
        if components_to_plot == 0:
            components_to_plot = 'cartesian'
        else:
            components_to_plot = 'geographical'
        check_2 = 0
    else:
        print('incorrect input, select either "cartesian" or "geographical". \n')

depth_start_time = time.time()
classified_path, files_to_classify = depth_classifier(sd_input, deflection_processing_path, individual_path,
                                                      components_to_plot, complete_files_path,
                                                      geographical_complete_files_path, headers_on)
depth_end_time = time.time() - depth_start_time
print 'The time spent to run depth discretization algorithm is: ', str(int(round(depth_end_time))), ' seconds'

# Moved from MATLAB: depths for plotting and latlon range, because the element with the highest stress has to be
# followed in this range only
min_lat = -90
max_lat = -65
min_lon = -150
max_lon = -60
# 145, 150
# 543, 552
# 775, 780
# 979, 995
depths_to_plot = np.array([545, 550])
area_params = [min_lat, max_lat, min_lon, max_lon, depths_to_plot]
element_tracker(base_path, area_params)
# Generate text files to be read by MATLAB containing all the relevant information for plotting purposes
if components_to_plot == 'cartesian':
    matlab_variables_writer(sd_input, complete_files_path, files_to_classify, level_dir, classified_path,
                            components_to_plot, area_params)
else:
    matlab_variables_writer(sd_input, geographical_complete_files_path, files_to_classify, level_dir, classified_path,
                            components_to_plot, area_params)
program_end_time = time.time() - program_start_time

# Count the time spent by the entire program to run

print 'The time spent to run the entire program is: ', str(int(round(program_end_time))), ' seconds'
print 'Data processing files generated. Going to the next iteration.'

