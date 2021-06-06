# This is the main file without the functions used to save stresses and deformations in .csv files for each part. More
# details are in the markdown file.

import os
import time
import platform
import pdb
from stress_report_reader import *
from deflection_report_reader import *
from coordinate_reader import *
from depth_classifier import *
from matlab_variables_writer import *
from element_tracker import *
from main_plot import *

# Define the run number () and time instant defined by the model sea-level iteration, step and viscosity iteration.
program_start_time = time.time()
run = '21'
iteration = 1
step = 0
cycle = 1
if not isinstance(run, str):
    run = str(run)
# Define the stress and deflection report files names
stress_rpt_name = 'abaqus_run_' + run + '.rpt'
deflection_rpt_name = 'abaqus_run_' + run + '_deflections.rpt'
# Define additional files for point coordinates and viscosity calculation coefficients
if run == 'point':
    dat_name = 'Earth_point.dat'
else:
    dat_name = 'Earth.dat'
e_dat_name = 'e.dat'

# Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
# windows or linux.
opersys = platform.system()
if opersys == 'Windows':
    base_path = os.path.join('C:\\Users\\fabri\\Desktop\\tu_delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
else:
    base_path = os.path.join('/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run))


iteration_path = os.path.join(base_path, 'Iteration_' + str(iteration))
if not os.path.exists(iteration_path):
    os.mkdir(iteration_path)
step_path = os.path.join(iteration_path, 'step_' + str(step))
if not os.path.exists(step_path):
    os.mkdir(step_path)
cycle_path = os.path.join(step_path, 'cycle_' + str(cycle) + '_reports')
if not os.path.exists(cycle_path):
    os.mkdir(cycle_path)
level_dir = cycle_path # Can be higher or lower level than this

# Define the paths where to read the necessary files and store the stress and deformation results
dat_path = os.path.join(base_path, dat_name)
stress_report_path = os.path.join(level_dir, stress_rpt_name)
deflection_report_path = os.path.join(level_dir, deflection_rpt_name)
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
while check_1 == 1:
    sd_input = eval(input('Enter 0 to work with stress components, 1 for the deflections: \n'))
    if sd_input == 0 or sd_input == 1:
        check_1 = 0
    else:
        print('Incorrect input, select either 1 or 0. \n')

# Process stresses and deformations here into csv files for each part here if necessary

# stress_report_reader(stress_report_path, stress_matrices_path, headers_on)
# deflection_report_reader(deflection_report_path, deflection_matrices_paths, headers_on)

# Process data and produce files containing streses and deformations associated to the corresponding centroid and node
# coordinates
start_time_coordinate_reader = time.time()
individual_path, complete_files_path, geographical_complete_files_path, large_node_matrix_path = \
    coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_matrices_paths,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder)
end_time_coordinate_reader = time.time() - start_time_coordinate_reader
if sd_input == 0:
    print(('The time spent to create the complete files containing the stress components and centroid coordinates ' \
          'is ', str(int(round(end_time_coordinate_reader))), ' seconds'))
else:
    print(('The time spent to create the complete files containing the deflection components and node coordinates ' \
          'is ', str(int(round(end_time_coordinate_reader))), ' seconds'))

# Decide whether to work with Cartesian or geographical coordinates, this is necessary for plotting
check_2 = 1
while check_2 == 1:
    components_to_plot = eval(input('Enter 0 to work with cartesian components, 1 to work with geographical ones: \n'))
    if components_to_plot == 0 or components_to_plot == 1:
        if components_to_plot == 0:
            components_to_plot = 'cartesian'
        else:
            components_to_plot = 'geographical'
        check_2 = 0
    else:
        print('incorrect input, select either "cartesian" or "geographical". \n')

#
if sd_input == 0:
    upper_path = individual_path
else:
    upper_path = deflection_processing_path
cartesian_classified_depth_path = os.path.join(upper_path, 'Depth_files')
geographical_classified_depth_path = os.path.join(upper_path, 'Geographical_depth_files')
file_extension = '.csv'
if components_to_plot == 'geographical':
    files_to_classify_path = geographical_complete_files_path
    histogram_path = os.path.join(geographical_classified_depth_path, 'histogram_plots')
    files_to_classify = [filename for filename in os.listdir(files_to_classify_path)
                         if filename.endswith(file_extension)]
else:
    files_to_classify_path = complete_files_path
    histogram_path = os.path.join(cartesian_classified_depth_path, 'histogram_plots')
    files_to_classify = [filename for filename in os.listdir(files_to_classify_path)
                         if filename.endswith(file_extension)]
if not os.path.exists(os.path.join(histogram_path, 'Complete_Earth_distribution.png')):
    depth_start_time = time.time()
    depth_classifier(sd_input, deflection_processing_path, individual_path, components_to_plot, headers_on,
                     file_extension, cartesian_classified_depth_path, geographical_classified_depth_path,
                     files_to_classify_path, files_to_classify, histogram_path)
    depth_end_time = time.time() - depth_start_time
    print('The time spent to run depth discretization algorithm is: ', str(int(round(depth_end_time))), ' seconds')

# Moved from MATLAB: depths for plotting and latlon range, because the element with the highest stress has to be
# followed in this range only
min_lat = -90
max_lat = -65
min_lon = -180
max_lon = 180
depths_to_plot = np.array([145, 160])
area_params = [min_lat, max_lat, min_lon, max_lon, depths_to_plot]

# Run function to track values over time for the element with the highest stress by producing a table containing them
# element_tracker(components_to_plot, base_path, area_params)

# Generate text files to be read by MATLAB containing all the relevant information for plotting purposes - This is
# not necessary when Python plotting is finished
matlab_variables_writer(sd_input, complete_files_path, files_to_classify, iteration_path, components_to_plot)
program_end_time = time.time() - program_start_time

# Decide whether to enter the plotting section or not
plotting_input = eval(input('Enter 1 to plot values over the AOI, 0 to terminate the program: \n'))
if plotting_input == 1:
    main_plot(sd_input, components_to_plot, iteration_path)
else:
    print('It has been decided to not generate plots, exiting the program.')

