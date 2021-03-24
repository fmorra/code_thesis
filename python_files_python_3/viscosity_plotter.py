import os
import pdb
import numpy as np
import cartopy.crs as chart
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from plot_scale_extremes import *
from geo2cart import *
from plot_scale_extremes import *
import math


def viscosity_plotter(sd_input, python_variables_base_path, coordinate_system, complete_matrices_path, report_path,
                      min_lat, max_lat, min_lon, max_lon, depths_to_plot, iteration, step, cycle, diff_matrix_path,
                      run_folder, run, figure_counter, python_base_path, run_vec, simul_time):

    parts_to_plot = ['EARTH']
    visco_strain_input = 0
    if visco_strain_input == 0:
        quantity = 'viscosity'
    else:
        quantity = 'strain_rate'
    viscosity_input = 1
    visco_diff_path = os.path.join(diff_matrix_path, quantity)
    if not os.path.exists(visco_diff_path):
        os.mkdir(visco_diff_path)
    e_file_name = 'e.dat'
    viscosity_figures_path = os.path.join(report_path, 'viscosity_plots')
    if not os.path.exists(viscosity_figures_path):
        os.mkdir(viscosity_figures_path)
    b_input = 0
    latlim = [min_lat, max_lat]
    lonlim = [min_lon, max_lon]
    min_depth = depths_to_plot[0]
    max_depth = depths_to_plot[1]
    resolution = 0.25
    earth_radius = 6371000
    lat_lin = np.arange(max_lat, min_lat, resolution)
    lon_lin = np.arange(min_lon, max_lon, resolution)
    [gridded_lon, gridded_lat] = np.meshgrid(lon_lin, lat_lin)
    r_out = []
    lon_plot_2 = []
    lat_plot_2 = []
    depthrange = np.arange(min_depth, max_depth, 1)
    for dd in depthrange:
        lon_plot_2 = np.vstack((lon_plot_2, gridded_lon[:]))
        lat_plot_2 = np.vstack((lat_plot_2, gridded_lat[:]))
        temp = (earth_radius - dd * 1e3) * np.ones(len(gridded_lon))
        r_out = np.vstack((r_out, temp[:]))

    for i in range(len(parts_to_plot)):
        iter_path = python_variables_base_path
        e_path = os.path.join(iter_path, e_file_name)
        if coordinate_system == 'cartesian':
            matrix_to_open_path = os.path.join(complete_matrices_path, 'Complete_file_' + parts_to_plot[i] + '.csv')
        else:
            matrix_to_open_path = os.path.join(complete_matrices_path, 'geographical_complete_file_' + parts_to_plot[i]
                                               + '.csv')
        colorbarlimits = plot_scale_extremes(sd_input, min_lat, max_lat, min_lon, max_lon, depths_to_plot,
                                             quantity, viscosity_input, b_input, python_base_path,
                                             run_vec, coordinate_system)
        e = open(e_path)
        e = e.readlines()
        e = [line.strip() for line in e[1:]]
        e = [np.array([eval(i) for i in line.split(",")[:]]) for line in e]
        e = np.array(e)
        elements = e[:, 1]
        alin = e[:, 2]
        a = e[:, 3]
        an = 3.5
        strain_rate = np.zeros((len(e), 1))
        opened_visco_matrix = open(matrix_to_open_path)
        opened_visco_matrix = opened_visco_matrix.readlines()
        opened_visco_matrix = [line.strip() for line in opened_visco_matrix[1:]]
        opened_visco_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in opened_visco_matrix]
        opened_visco_matrix = np.array(opened_visco_matrix)
        mises = opened_visco_matrix[:, 2]
        power = (an)
        for j in range(len(strain_rate)):
            strain_rate[j, 0] = alin[j] * mises[j] + a[j] * mises[j] ^ power
        viscosity = mises / (3 * strain_rate)
        viscosity_matrix = np.zeros((len(opened_visco_matrix), 4))
        viscosity_matrix[:, 0] = elements
        viscosity_matrix[:, -3] = viscosity
        viscosity_matrix[:, -2] = strain_rate
        viscosity_matrix[:, -1] = mises
        complete_matrix_with_viscosity = np.hstack((opened_visco_matrix, viscosity_matrix))
        matrix_to_read = complete_matrix_with_viscosity
        depth = matrix_to_read[:, -7] / 1e3
        lat = matrix_to_read[:, -6]
        lon = matrix_to_read[:, -5]
        depth_condtion = depth > min_depth & depth < max_depth
        lat_condition = lat > min_lat & lat < max_lat
        lon_condition = lon > min_lon & lon < max_lon
        data_points_indices = matrix_to_read[depth_condtion * lat_condition * lon_condition]
        matrix_for_difference = np.zeros((len(data_points_indices), 4))
        matrix_for_difference_headers = [quantity, 'Latitude', 'Longitude', 'Depth']
        if visco_strain_input == 0:
            plot_variable = matrix_to_read[data_points_indices, -3]
            plot_variable = math.log10(plot_variable)
        else:
            plot_variable = matrix_to_read[data_points_indices, -2]

        filtered_lat = lat[data_points_indices]
        filtered_lon = lon[data_points_indices]
        depth_out = depth[data_points_indices]
        filtered_r = earth_radius - 1e3 * depth_out
        matrix_for_difference[:, -4] = plot_variable
        matrix_for_difference[:, -3] = filtered_lat
        matrix_for_difference[:, -2] = filtered_lon
        matrix_for_difference[:, -1] = depth_out
        [x_in, y_in, z_in] = geo2cart(filtered_lat, filtered_lon, filtered_r)
        [x_out, y_out, z_out] = geo2cart(np.deg2rad(lon_plot_2), np.deg2rad(lat_plot_2), r_out)
        plot_variable_out = interp.griddata(x_in, y_in, z_in, plot_variable, x_out, y_out, z_out, 'nearest')

        # Open image
        stress_defo_figure = plt.figure()
        # Geographic map
        ax = plt.axes(projection=chart.Orthographic(0, -90))
        ax.coastlines(resolution='50m', color='orange', linewidth=1)
        # Add surface
        ax.plot_surface(lon_plot_2, lat_plot_2, plot_variable_out, cmap=cm.summer, linewidth=0, antialiased=False)
        # Colorbar settings
        scale = plt.colorbar(stress_defo_figure)
        stress_defo_figure.clim(math.log10(colorbarlimits[0]), math.log10(colorbarlimits[1]))  # same as caxis
        if sd_input == 0:
            scale.set_label(quantity[j] + ' (Pa)', rotation=270)
        else:
            scale.set_label(quantity[j] + ' (m)', rotation=270)
        # Title settings
        if run % 2 == 0:
            rheology = ', wet rheology'
        else:
            rheology = ', dry rheology'
        if sd_input == 0:
            plt.title('Time ' + simul_time + ', stress iteration ' + str(cycle) + rheology)
            # plt.title({['Time ' simul_time], [' ']});
        else:
            plt.title('Time ' + simul_time + ', stress iteration ' + str(cycle) + rheology)
            # plt.title({['Time ' simul_time], [' ']});

        if visco_strain_input == 0:
            plt.savefig(os.path.join(viscosity_figures_path, 'Viscosity_' + parts_to_plot[i] + '[' + str(min_depth)
                                     + '-' + str(max_depth) + ']_km_' + quantity + '_' + run + '_' + iteration + '_'
                                     + step + '_' + cycle + '.png'))
        else:
            plt.savefig(os.path.join(viscosity_figures_path, 'Strain_rate_' + parts_to_plot[i] + '[' + str(min_depth)
                                     + '-' + str(max_depth) + ']_km_' + quantity + '_' + run + '_' + iteration + '_'
                                     + step + '_' + cycle + '.png'))
        diff_matrix = pd.DataFrame(data=matrix_for_difference, columns=matrix_for_difference_headers[1:])
        diff_matrix.to_csv(os.path.join(visco_diff_path, 'Iteration_' + iteration + '_step_' +
                                        step + '_cycle_' + cycle + '_' + str(min_depth) + '_' +
                                        str(max_depth) + '_km.csv'), index=False)

    figure_counter = figure_counter+3

    return figure_counter, visco_diff_path
