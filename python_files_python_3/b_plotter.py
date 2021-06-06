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


def b_plotter(viscosity_figures_path, python_variables_base_path, depths_to_plot, min_lat, max_lat, min_lon,
              max_lon, run, python_base_path, coordinate_system, diff_matrix_path, quantity,
              run_vec, parts_to_plot, complete_matrices_path):

    b_input = 1
    viscosity_input = 1
    visco_strain_input = 0
    e_path = os.path.join(python_variables_base_path, 'e.dat')
    opened_e = open(e_path)
    e = pd.read_csv(opened_e, delimiter=" ", header=None)
    names = ['B_{diff}', 'B_{disl}']
    sd_input = 0
    min_depth = depths_to_plot[0]
    max_depth = depths_to_plot[1]

    resolution = 0.25
    earth_radius = 6371000
    lat_lin = np.arange(max_lat, min_lat - resolution, -resolution)
    lon_lin = np.arange(min_lon, max_lon + resolution, resolution)
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

    colorbarlimits = plot_scale_extremes(sd_input, min_lat, max_lat, min_lon, max_lon, depths_to_plot,
                                         quantity, viscosity_input, visco_strain_input, b_input, python_base_path,
                                         run_vec, coordinate_system)

    for j in range(len(parts_to_plot)):
        visco_matrix = pd.read_csv(open(os.path.join(complete_matrices_path, 'Complete_file_' + parts_to_plot[j] +
                                                     '.csv')), delimiter=",")
        matrix_to_read = np.hstack((visco_matrix, e))
        depth = matrix_to_read[:, -7] / 1e3
        lat = matrix_to_read[:, -6]
        lon = matrix_to_read[:, -5]
        depth_condtion = (depth >= min_depth) & (depth <= max_depth)
        lat_condition = (lat >= min_lat) & (lat <= max_lat)
        lon_condition = (lon >= min_lon) & (lon <= max_lon)
        filtered_matrix = matrix_to_read[depth_condtion & lat_condition & lon_condition]
        depth_out = filtered_matrix[:, -7]
        filtered_lat = filtered_matrix[:, -6]
        filtered_lon = filtered_matrix[:, -5]
        filtered_r = earth_radius - 1e3 * depth_out
        matrix_for_difference = np.zeros((len(filtered_matrix), len(names) + 3))
        matrix_for_difference[:, -3] = filtered_lat
        matrix_for_difference[:, -2] = filtered_lon
        matrix_for_difference[:, -1] = depth_out
        matrix_for_difference_headers = [names[0], names[1], 'Latitude', 'Longitude', 'Depth']
        B_path = os.path.join(diff_matrix_path, 'B_coefficients')
        if not os.path.exists(B_path):
            os.mkdir(B_path)
        variables = np.hstack((filtered_matrix[-3], filtered_matrix[-2]))

        for i in range(len(variables[0])):
            plot_variable = variables.iloc[:, i]
            matrix_for_difference[:, i] = plot_variable
            pandas_length = len(filtered_lon)
            [x_in, y_in, z_in] = geo2cart(filtered_lon, filtered_lat, filtered_r, pandas_length)
            [x_out, y_out, z_out] = geo2cart(lon_plot, lat_plot, r_out, pandas_length)
            plot_variable_out = interp.griddata((x_in, y_in, z_in), plot_variable, (x_out, y_out, z_out), 'nearest')
            viscosity_figure = plt.figure()
            # Geographic map
            ax = plt.axes(projection=chart.SouthPolarStereo())
            ax.set_extent([-180, 180, -90, -65], chart.PlateCarree())
            ax.coastlines(zorder=3)
            ax.gridlines(draw_labels=True)
            # Add surface and colorbar
            levels = np.linspace(colorbarlimits[0], colorbarlimits[-1], 10)
            surf = ax.contourf(lon_plot, lat_plot, plot_variable_out, levels=levels, cmap='summer', antialiased=False,
                               alpha=0.7, extend='min', transform=chart.PlateCarree())
            scale = viscosity_figure.colorbar(surf)
            if i == 0:
                scale.set_label('B_{diff}', labelpad=10)
            else:
                scale.set_label('B_{disl}', labelpad=10)
            # Title settings
            if int(run) % 2 == 0:
                rheology = ', wet rheology'
            else:
                rheology = ', dry rheology'
            plt.title('Map of ' + names[i] + rheology)

            plt.savefig(os.path.join(viscosity_figures_path, names[i] + '_' + parts_to_plot[j] + '_' + run + '_[' +
                                     str(min_depth) + '-' + str(max_depth) + ']_km.png'))
            diff_matrix = pd.DataFrame(data=matrix_for_difference, columns=matrix_for_difference_headers)
            diff_matrix.to_csv(os.path.join(B_path, 'B_' +  str(run) +  '_' + str(min_depth) + '_' + str(max_depth) +
                                            '_km.csv'), index=False)

