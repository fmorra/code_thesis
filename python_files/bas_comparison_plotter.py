import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib.ticker import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
import glob
import pdb

extension = '.dat'
my_files_path = 'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\bas_comparison_files'
bas_files_path = 'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\ABAQUS_3D_GIA_model\\' \
                 'Input_output_standard_run\\dat_files'
dat_names = []
for f in glob.glob(os.path.join(my_files_path, '*' + extension)):
    dat_names.append(f)
fig_counter = 0
for dat in dat_names:
    dat_name = dat.split('\\')[-1]
    my_data_path = open(os.path.join(my_files_path, dat_name))
    my_data = my_data_path.readlines()
    my_data = [line.strip() for line in my_data]
    my_data = [np.array([float(i) for i in line.split(" ")[:]]) for line in my_data]
    my_data = np.array(my_data)

    bas_data_path = open(os.path.join(bas_files_path, dat_name))
    bas_data = bas_data_path.readlines()
    bas_data = [line.strip() for line in bas_data]
    bas_data = [np.array([float(i) for i in line.split(" ")[:]]) for line in bas_data]
    bas_data = np.array(bas_data)

    diff_matrix = my_data - bas_data
    diff_matrix_percent = (diff_matrix/bas_data)*100
    for j in range(len(diff_matrix_percent)):
        for k in range(len(diff_matrix_percent[1])):
            if abs(diff_matrix_percent[j, k]) > 5:
                diff_matrix_percent[j, k] = float("NaN")
    diff_matrix_percent = np.transpose(diff_matrix_percent)
    lon_vector = np.linspace(-180, 180, len(diff_matrix_percent))
    lat_vector = np.linspace(-90, 90, len(diff_matrix_percent[1]))
    fig = plt.figure(fig_counter)
    ax = fig.gca(projection='3d')
    ax.set_xlim(lon_vector[0], lon_vector[-1])
    ax.set_ylim(lat_vector[0], lat_vector[-1])
    plt.xticks(np.arange(lon_vector[0], lon_vector[-1] + 20, 20))
    plt.yticks(np.arange(lat_vector[0], lat_vector[-1] + 10, 10))
    # ax.view_init(azim=0, elev=90)

    # lon_vector_plot, lat_vector_plot = np.mgrid[lon_vector[0]:lon_vector[-1]:1000j, lat_vector[0]:lat_vector[-1]:500j]
    # pdb.set_trace()
    interp_spline = interpolate.RectBivariateSpline(lon_vector, lat_vector, diff_matrix_percent)
    lon_vector_to_plot, lat_vector_to_plot = np.meshgrid(lon_vector, lat_vector)
    diff_matrix_percent_new = interp_spline(lon_vector, lat_vector)
    # surf = ax.plot_surface(lon_vector_to_plot, lat_vector_to_plot, diff_matrix_percent_new, cmap='summer', rstride=1, cstride=1,
    #                        alpha=None, antialiased=True)
    ax.plot_wireframe(lon_vector_to_plot, lat_vector_to_plot, diff_matrix_percent_new, color='k')
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    # ax.set_zlabel('Percentage difference (%)')
    plt.title('Percentage difference for ' + dat_name)
    plt.show()
    print diff_matrix_percent
    print diff_matrix_percent_new
    pdb.set_trace()
    fig_counter = fig_counter + 1
