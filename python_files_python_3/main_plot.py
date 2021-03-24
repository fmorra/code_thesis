import os
import pdb
import re
from stress_deflection_plotter import *
from viscosity_plotter import *


def main_plot(files_to_classify, report_path):
    python_base_path = 'C:\\Users\\fabri\\Desktop\\tu_delft\\Thesis\\ABAQUS\\test_run_python'
    figure_counter = 1
    run = '21'
    run_vec = []
    for run_directory in os.listdir(python_base_path):
        if re.findall(r'\d+', run_directory):
            run_vec = re.findall(r'\d+', run_directory)
    iteration = '1'
    step = '0'
    simul_time = '1 ka'
    cycle = '1'
    python_variables_path = os.path.join(python_base_path, 'run_' + run, 'Iteration_' + iteration, 'step_' + step,
                                         'cycle_' + cycle + '_reports')
    coordinate_system_input = eval(input('Enter the coordinate system to use for plotting, 0 for cartesian, 1 for '
                                         'geographical: \n'))
    sd_input = eval(input('Enter 0 to plot the stresses, 1 for the deflections: \n'))
    if sd_input == 0:
        quantity = 'stress'
    else:
        quantity = 'deflection'
    figure_folder = quantity + '_plots'
    if coordinate_system_input == 0:
        coordinate_system = 'cartesian'
        complete_matrices_name = 'Complete_files'
    else:
        coordinate_system = 'geographical'
        complete_matrices_name = 'Geographical_complete_files'
    min_lat = -90
    max_lat = -65
    min_lon = -180
    max_lon = 180
    depths_to_plot = [145, 160]
    complete_matrices_path = os.path.join(python_variables_path, quantity + '_processing_results',
                                          complete_matrices_name)
    figures_path = os.path.join(python_variables_path, figure_folder)
    if not os.path.exists(figures_path):
        os.mkdir(figures_path)
    run_folder = os.path.join(python_base_path, 'run_', run)
    if depths_to_plot[0] > depths_to_plot[1]:
        depths_to_plot[0], depths_to_plot[1] = depths_to_plot[1], depths_to_plot[0].copy()
    diff_matrix_path = os.path.join(python_base_path, 'run_' + run, 'difference_matrices_plots')
    if not os.path.exists(diff_matrix_path):
        os.mkdir(diff_matrix_path)
    while_check = 1
    while while_check == 1:
        plots_bool = eval(input('Do you want to plot the stresses or deformations? \n'
                                'Enter 1 for plots, 0 to skip this step: \n'))
        if plots_bool == 1:
            figure_counter = stress_deflection_plotter(run, sd_input, coordinate_system, complete_matrices_path,
                                                       figures_path, diff_matrix_path, iteration, step, cycle, min_lat,
                                                       max_lat, min_lon, max_lon, depths_to_plot, figure_counter,
                                                       python_base_path, run_vec, simul_time, coordinate_system_input)
            while_check = 0
        elif plots_bool == 1:
            while_check = 0
        else:
            print('Incorrect input, enter 1 for plots and 0 to skip this step.')
    pdb.set_trace()
    visco_check = 1
    while visco_check == 1:
        visco_plots_bool = eval(input('Do you want to plot the viscosity/strain rate? \n Enter 1 for plots, 0 to '
                                      'skip them: \n'))
        if visco_plots_bool == 1 or visco_plots_bool == 0:
            visco_check = 0
        else:
            print('Incorrect input, enter 1 for the viscosity or strain rate plots or 0 to skip this step: \n')
    if visco_plots_bool == 1:
        viscosity_figures_path = os.path.join(report_path, 'viscosity_plots')
        if not os.path.exists(viscosity_figures_path):
            os.mkdir(viscosity_figures_path)
        figure_counter, visco_diff_path = viscosity_plotter(sd_input, python_variables_path, coordinate_system,
                                                            complete_matrices_path, report_path, min_lat, max_lat,
                                                            min_lon, max_lon, depths_to_plot, iteration, step, cycle,
                                                            diff_matrix_path, run_folder, run, figure_counter,
                                                            python_base_path, run_vec, simul_time)
    diff_check = 1
    while diff_check == 1:
        diff_input = input(
            'Enter 1 to generate plots of differences between different cycles and steps, 0 to skip this: \n')
        if diff_input == 1 or diff_input == 0:
            diff_check = 0
        else:
            print('Incorrect input, enter 1 to generate difference plots or 0 to skip: \n')
    b_diff_plots_folder = os.path.join(diff_matrix_path, 'plots', 'viscosity')
    if not os.path.exists(b_diff_plots_folder):
        os.mkdir(b_diff_plots_folder)
    b_diff_calc(b_diff_plots_folder, python_base_path, min_lat, max_lat, min_lon, max_lon, depths_to_plot[0],
                depths_to_plot[1], figure_counter)
    if diff_input == 1:
        diff_calculator(diff_matrix_path, min_lat, max_lat, min_lon, max_lon,run, depths_to_plot, iteration, step,
                        figure_counter)
