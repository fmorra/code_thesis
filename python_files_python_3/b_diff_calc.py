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


def b_diff_calc(b_diff_plots_folder, python_base_path, min_lat, max_lat, min_lon, max_lon, min_depth,
                    max_depth, run_vec):

    b_runs = [run_vec[0], run_vec[1]]
    diff_files = []
    for i in range(len(b_runs)):
        diff_matrix_path = os.path.join(python_base_path, 'run_' + str(b_runs[i]),
                                  'difference_matrices_plots\\viscosity\\B_coeffs')
        dir_diff_files = os.path.join(diff_matrix_path, '*.csv')
        for j in range(len(dir_diff_files)):
            file_name = dir_diff_files[i].split('\\')
            diff_files.append(file_name[-1])

    resolution = 0.25
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
            if i == 0:
                scale.set_label(diff_variables_to_plot[i])
            else:
                scale.set_label(diff_variables_to_plot[i])
            plt.savefig(os.path.join(b_diff_plots_folder, 'diff_plot_' + possible_variables_for_plotting[i] + '_' +
                                     '_depth_range_[' + str(min_depth) + '-' + str(max_depth) + ']_km.png'))

