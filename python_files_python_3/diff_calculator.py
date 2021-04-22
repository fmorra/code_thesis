import os
import pdb
import numpy as np
import cartopy.crs as chart
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pandas as pd
import re
from matplotlib import cm
from plot_scale_extremes import *
from geo2cart import *
import math


def diff_calculator(diff_matrix_path_incomplete, min_lat, max_lat, min_lon, max_lon, run, depths_to_plot, iteration,
                    step, figure_counter):
    # [1 0 1 145 150]
    # [1 0 2 145 150]
    resolution = 0.25
    latlim = [min_lat, max_lat]
    lonlim = [min_lon, max_lon]
    min_depth = depths_to_plot[0]
    max_depth = depths_to_plot[1]
    earth_radius = 6371000
    lat_lin = np.arange(max_lat, min_lat, -resolution)
    lon_lin = np.arange(min_lon, max_lon, resolution)
    [gridded_lon, gridded_lat] = np.meshgrid(lon_lin, lat_lin)
    r_out = np.array([])
    lon_plot = np.array([])
    lat_plot = np.array([])
    depthrange = np.arange(min_depth, max_depth, 1)
    for dd in depthrange:
        lon_plot = np.vstack((lon_plot, gridded_lon)) if lon_plot.size else gridded_lon
        lat_plot = np.vstack((lat_plot, gridded_lat)) if lat_plot.size else gridded_lat
        temp = (earth_radius - dd * 1e3) * np.ones((len(gridded_lon), len(gridded_lon[0])))
        r_out = np.vstack([r_out, temp]) if r_out.size else temp

    visco_choice = 0
    if visco_choice == 0:
        stress_deflection_input = 1
        ref_system_input = 1
        if stress_deflection_input == 0:
            quantity = 'stresses'
        else:
            quantity = 'deflections'
        if ref_system_input == 0:
            ref_system = 'cartesian'
        else:
            ref_system = 'geographical'
        diff_matrix_path = os.path.join(diff_matrix_path_incomplete, quantity, ref_system)
        diff_plots_folder = os.path.join(diff_matrix_path_incomplete, 'plots', quantity, ref_system)
    else:
        visco_strain_input = 0
        if visco_strain_input == 0:
            quantity = 'viscosity'
        else:
            quantity = 'strain_rate'
        diff_matrix_path = os.path.join(diff_matrix_path_incomplete, quantity)
        diff_plots_folder = os.path.join(diff_matrix_path_incomplete, 'plots', quantity)

    dir_diff_files = os.path.join(diff_matrix_path, '*.csv')
    diff_files = []
    if not os.path.exists(diff_plots_folder):
        os.mkdir(diff_plots_folder)
    for i in range(len(dir_diff_files)):
        file_name = dir_diff_files[i].split('\\')
        diff_files.append(file_name[-1])
    matrix_1_flag = 0
    found_values_1 = []
    quintuplet_1 = np.arr([])
    while matrix_1_flag == 0:
        quintuplet_1 = input('Enter the iteration, step and cycle number, and the minimum and maximum depth of the first '
                                'matrix used to plot the difference with respect to another moment in time as a vector of '
                                'square brackets of five values: \n')
        for i in range(len(diff_files)):
            found_values_1 = re.findall(r"\((\d+(?:,\d+)*)\)", quintuplet_1)
            if (quintuplet_1 == found_values_1).all():
                matrix_1_flag = 1
                file_to_read_1 = diff_files[i]
                iteration_1 = quintuplet_1[0]
                step_1 = quintuplet_1[1]
                cycle_1 = quintuplet_1[2]
                break
        if matrix_1_flag == 1:
            break
        else:
            print('No match for the selected values, the matrix does not exist. Exiting search.')
            break

    found_values_2 = []
    matrix_2_flag = 0
    quintuplet_2 = np.arr([])
    while matrix_2_flag == 0:
        quintuplet_2 = input('Enter the iteration, step and cycle number, and the minimum and maximum depth of the first '
                                'matrix used to plot the difference with respect to another moment in time as a vector of '
                                'square brackets of five values: \n')
        for i in range(len(diff_files)):
            found_values_2 = re.findall(r"\((\d+(?:,\d+)*)\)", quintuplet_2)
            if (quintuplet_2 == found_values_2).all():
                matrix_1_flag = 1
                file_to_read_2 = diff_files[i]
                iteration_2 = quintuplet_2[0]
                step_2 = quintuplet_2[1]
                cycle_2 = quintuplet_2[2]
                break
        if matrix_2_flag == 1:
            break
        else:
            print('No match for the selected values, the matrix does not exist. Exiting search.')
            break

    if matrix_1_flag == 1 and matrix_2_flag == 1:
        min_depth = quintuplet_2[3]
        max_depth = quintuplet_2[4]
        headers_1 = pd.read_csv(os.path.join(diff_matrix_path, file_to_read_1), nrows = 1)
        headers_2 = pd.read_csv(os.path.join(diff_matrix_path, file_to_read_2), nrows = 1)
        values_1 = pd.read_csv(os.path.join(diff_matrix_path, file_to_read_1), skiprows = 1)
        values_2 = pd.read_csv(os.path.join(diff_matrix_path, file_to_read_2), skiprows = 1)
        possible_variables_for_plotting = headers_1.intersection(headers_2, 'stable')
        possible_variables_for_plotting[-3] = []
        possible_variables_for_plotting[-2] = []
        possible_variables_for_plotting[-1] = []
        print('The possible variables to generate the difference plots between these two cycles are: \n')
        print(possible_variables_for_plotting)
        diff_variables_to_plot = possible_variables_for_plotting
        filtered_lat = values_1[:, -3]
        filtered_lon = values_1[:, -2]
        depth_out = values_1[:, -1]
        filtered_r = earth_radius - 1e3 * depth_out
        for i in range(len(diff_variables_to_plot)):
            diff_variable_1 = values_1[:, (headers_1 == diff_variables_to_plot[i]).all()]
            diff_variable_2 = values_2[:, (headers_2 == diff_variables_to_plot[i]).all()]
            plot_variable = diff_variable_2 - diff_variable_1
            [x_in, y_in, z_in] = geo2cart(filtered_lat, filtered_lon, filtered_r)
            [x_out, y_out, z_out] = geo2cart(np.deg2rad(lon_plot), np.deg2rad(lat_plot), r_out)
            plot_variable_out = interp.griddata(x_in, y_in, z_in, plot_variable, x_out, y_out, z_out, 'nearest')

            # Open image
            diff_figure = plt.figure()
            # Geographic map
            ax = plt.axes(projection=chart.Orthographic(0, -90))
            ax.coastlines(resolution='50m', color='orange', linewidth=1)
            # Add surface
            ax.plot_surface(lon_plot, lat_plot, plot_variable_out, cmap=cm.summer, linewidth=0, antialiased=False)
            # Colorbar settings
            scale = plt.colorbar(diff_figure)
            diff_figure.clim('auto')
            scale = plt.colorbar(diff_figure)
            if visco_choice == 1:
                if visco_strain_input == 0:
                    scale.set_label(diff_variables_to_plot[i] + ' [N \cdot s/m^2]')
                else:
                    scale.set_label(diff_variables_to_plot[i])
            else:
                if stress_deflection_input == 0:
                    scale.set_label(diff_variables_to_plot[i] + ' [Pa]')
                else:
                    scale.set_label(diff_variables_to_plot[i] + ' [m]')
            if visco_choice == 0:
                if diff_variables_to_plot[i] == 'Mises':
                    plt.title('Map of value differences for the Mises stresses')
                elif diff_variables_to_plot[i] == 'Magnitude':
                    plt.title('Map of value differences for the deformation magnitude')
                else:
                    plt.title('Map of value differences for ' +  quantity + ' component ' + diff_variables_to_plot[i])
            else:
                plt.title('Map of differences for the ' + quantity)
            if visco_choice == 1:
                plt.savefig(os.path.join(diff_plots_folder, 'diff_plot_' + quantity + '_' + '_depth_range_[' + str(min_depth)
                + '-' + str(max_depth) + ']_km_' + '_' + run + '_[' + str(iteration_1) + '_' + str(iteration_2) + ']_['
                + str(step_1) + '_' + str(step_2) + ']_[' + str(cycle_1) + '_' + str(cycle_2) + '].png'))
            else:
                plt.savefig(os.path.join(diff_plots_folder, 'diff_plot_' + quantity + '_' + '_depth_range_[' + str(min_depth)
                + '-' + str(max_depth) + ']_km_' + '_' + run + '_[' + str(iteration_1) + '_' + str(iteration_2) + ']_['
                + str(step_1) + '_' + str(step_2) + ']_[' + str(cycle_1) + '_' + str(cycle_2) + '].png'))

        figure_counter = figure_counter + 1