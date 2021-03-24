import os
import pdb
import re
import numpy as np
import cartopy.crs as chart
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pandas as pd
import cartopy
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from plot_scale_extremes import *
from geo2cart import *


def stress_deflection_plotter(run, sd_input, coordinate_system, complete_matrices_path, figures_path,
                              diff_matrix_path, iteration, step, cycle, min_lat, max_lat, min_lon, max_lon,
                              depths_to_plot, figure_counter, python_base_path, run_vec, simul_time,
                              coordinate_system_input):
    if sd_input == 0:
        components_to_plot = 'stresses'
    else:
        components_to_plot = 'deflections'
    viscosity_input = 0
    b_input = 0
    complete_diff_path = os.path.join(diff_matrix_path, components_to_plot, coordinate_system)
    if not os.path.exists(complete_diff_path):
        os.makedirs(complete_diff_path)
    parts_to_plot = ['EARTH']
    resolution = 0.25
    latlim = [min_lat, max_lat]
    lonlim = [min_lon, max_lon]
    min_depth = depths_to_plot[0]
    max_depth = depths_to_plot[1]
    earth_radius = 6371000
    lat_lin = np.arange(max_lat, min_lat-resolution, -resolution)
    lon_lin = np.arange(min_lon, max_lon+resolution, resolution)
    [gridded_lon, gridded_lat] = np.meshgrid(lon_lin, lat_lin)
    r_out = np.array([])
    lon_plot_2 = np.array([])
    lat_plot_2 = np.array([])
    depthrange = np.arange(min_depth, max_depth, 1)
    for dd in depthrange:
        lon_plot_2 = np.vstack((lon_plot_2, gridded_lon.reshape(-1, 1))) if lon_plot_2.size else \
            gridded_lon.reshape(-1, 1)
        lat_plot_2 = np.vstack((lat_plot_2, gridded_lat.reshape(-1, 1))) if lat_plot_2.size else \
            gridded_lat.reshape(-1, 1)
        temp = (earth_radius - dd * 1e3) * np.ones((len(gridded_lon), len(gridded_lon[0])))
        r_out = np.vstack([r_out, temp.reshape(-1, 1)]) if r_out.size else temp.reshape(-1, 1)
    if sd_input == 0:
        selected_component = eval(input('Enter the desired stress component(s) to plot, possible values are Mises, '
                                        'S11, S22, S33, S12, S13, S23: \n'))
        selected_components = selected_component.split(" ")
        colorbarlimits = plot_scale_extremes(sd_input, min_lat, max_lat, min_lon, max_lon, depths_to_plot,
                                             selected_components, viscosity_input, b_input, python_base_path,
                                             run_vec, coordinate_system_input)
        selected_columns = []
        for k in range(len(selected_components)):
            if selected_components[k] == 'Mises':
                selected_columns.append(1)
            elif selected_components[k] == 'S11':
                selected_columns.append(2)
            elif selected_components[k] == 'S22':
                selected_columns.append(3)
            elif selected_components[k] == 'S33':
                selected_columns.append(4)
            elif selected_components[k] == 'S12':
                selected_columns.append(5)
            elif selected_components[k] == 'S13':
                selected_columns.append(6)
            else:
                selected_columns.append(7)
    else:
        selected_component = eval(input('Enter the desired deflection components(s) to plot, possible values are '
                                        'Magnitude, U1, U2, U3: \n'))
        selected_components = selected_component.split(" ")
        colorbarlimits = plot_scale_extremes(sd_input, min_lat, max_lat, min_lon, max_lon, depths_to_plot,
                                             selected_components, viscosity_input, b_input, python_base_path,
                                             run_vec, coordinate_system)
        selected_columns = []
        for k in range(len(selected_components)):
            if selected_components[k] == 'Magnitude':
                selected_columns.append(1)
            elif selected_components[k] == 'U1':
                selected_columns.append(2)
            elif selected_components[k] == 'U2':
                selected_columns.append(3)
            else:
                selected_columns.append(4)

    for i in range(len(parts_to_plot)):
        min_depth = depths_to_plot[0]
        max_depth = depths_to_plot[1]
        if coordinate_system == 'cartesian':
            matrix_to_open = open(os.path.join(complete_matrices_path, 'Complete_file_' + parts_to_plot[i] + '.csv'))
        else:
            matrix_to_open = open(os.path.join(complete_matrices_path, 'Geographical_complete_file_' +
                                               parts_to_plot[i] + '.csv'))
        # data_matrix = matrix_to_open.readlines()
        # data_matrix = [line.strip() for line in data_matrix[1:]]
        # data_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in data_matrix]
        # matrix_to_read = np.array(data_matrix)
        matrix_to_read = pd.read_csv(matrix_to_open, delimiter=",")
        depth = matrix_to_read.iloc[:, -3] / 1e3
        lat = matrix_to_read.iloc[:, -2]
        lon = matrix_to_read.iloc[:, -1]
        depth_condtion = (depth >= min_depth) & (depth <= max_depth)
        lat_condition = (lat >= min_lat) & (lat <= max_lat)
        lon_condition = (lon >= min_lon) & (lon <= max_lon)
        filtered_matrix = matrix_to_read[depth_condtion & lat_condition & lon_condition]
        matrix_for_difference = np.zeros((len(filtered_matrix), len(selected_columns) + 3))
        matrix_for_difference_headers = [selected_components, 'Latitude', 'Longitude', 'Depth']
        depth_out = filtered_matrix.iloc[:, -3]
        filtered_lat = filtered_matrix.iloc[:, -2]
        filtered_lon = filtered_matrix.iloc[:, -1]
        filtered_r = earth_radius - depth_out
        matrix_for_difference[:, -3] = filtered_lat
        matrix_for_difference[:, -2] = filtered_lon
        matrix_for_difference[:, -1] = depth_out
        for j in range(len(selected_columns)):
            plot_variable = filtered_matrix.iloc[:, selected_columns[j]]
            matrix_for_difference[:, j] = plot_variable
            [x_in, y_in, z_in] = geo2cart(filtered_lon, filtered_lat, filtered_r)
            [x_out, y_out, z_out] = geo2cart(np.deg2rad(lon_plot_2), np.deg2rad(lat_plot_2), r_out)
            visual_check = plt.figure()
            ax = Axes3D(visual_check)
            ax.scatter(x_in, y_in, z_in, alpha=0.8)
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel('Z [m]')
            ax.set_title('Distribution of points for the [' + str(min_depth) + '-' + str(max_depth) + '] km range')
            plt.savefig(os.path.join(figures_path, 'visual_mesh_check_depth_[' + str(min_depth) + '-' + str(max_depth) +
                                     ']_km.png'), bbox_inches='tight')
            plt.close(visual_check)
            plot_variable_out = interp.griddata((x_in, y_in, z_in), plot_variable, (x_out, y_out, z_out), 'nearest')
            plot_variable_out = plot_variable_out.reshape(len(lon_plot_2), len(lat_plot_2))
            # Open image
            stress_defo_figure = plt.figure()
            # Geographic map
            ax = plt.axes(projection=chart.Orthographic(0, -90))
            ax.add_feature(cartopy.feature.LAND, edgecolor='orange')
            # ax = plt.axes(projection=chart.SouthPolarStereo())
            # ax.set_extent([-180, 180, -90, -65], crs=chart.PlateCarree())
            # ax.add_feature(cartopy.feature.LAND, edgecolor='orange')
            # Add surface
            # ax.plot_surface(lon_plot_2, lat_plot_2, plot_variable_out, cmap=cm.summer, linewidth=0, antialiased=False)
            pdb.set_trace()
            ax.contourf(lon_plot_2, lat_plot_2, plot_variable_out, transform=chart.Orthographic())
            # Colorbar settings
            scale = plt.colorbar(stress_defo_figure)
            stress_defo_figure.clim(colorbarlimits[j], colorbarlimits[j + len(colorbarlimits) / 2]) # same as caxis
            if sd_input == 0:
                scale.set_label(selected_components[j] + ' (Pa)', rotation=270)
            else:
                scale.set_label(selected_components[j] + ' (m)', rotation=270)
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
            plt.show()
            pdb.set_trace()
            plt.savefig(os.path.join(figures_path, coordinate_system + '_' + components_to_plot + '_' + parts_to_plot[i]
                                     + '[' + str(min_depth) + '-' + str(max_depth) + ']_km_' + selected_components[j] +
                                     '_' + run + '_' + iteration + '_' + step + '_' + cycle + '.png'))
            figure_counter = figure_counter + 1
            plt.close(stress_defo_figure)

        diff_matrix = pd.DataFrame(data=matrix_for_difference, columns=matrix_for_difference_headers[1:])
        diff_matrix.to_csv(os.path.join(complete_diff_path, 'Iteration_' + iteration + '_step_' +
                                                               step + '_cycle_' + cycle + '_' + str(min_depth) + '_' +
                                                               str(max_depth) + '_km.csv'), index=False)