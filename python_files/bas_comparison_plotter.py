import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib.ticker import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
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
    my_data_path = open(os.path.join(my_files_path, dat))
    my_data = my_data_path.readlines()
    my_data = [line.strip() for line in my_data]
    my_data = [np.array([float(i) for i in line.split(" ")[:]]) for line in my_data]
    my_data = np.array(my_data)

    bas_data_path = open(os.path.join(my_files_path, dat))
    bas_data = bas_data_path.readlines()
    bas_data = [line.strip() for line in bas_data]
    bas_data = [np.array([float(i) for i in line.split(" ")[:]]) for line in bas_data]
    bas_data = np.array(bas_data)

    lat_vector = np.linspace(-90, 90, len(my_data))
    lon_vector = np.linspace(0, 360, len(my_data[1]))
    lon_vector_plot = lon_vector - 180
    lat_vector_plot = lat_vector
    diff_matrix = my_data - bas_data
    diff_matrix_percent = (diff_matrix/bas_data)*100
    for j in range(len(diff_matrix_percent)):
        for k in range(len(diff_matrix_percent[1])):
            if abs(diff_matrix_percent[j, k]) > 5:
                diff_matrix_percent[j, k] = float("NaN")

    fig = plt.figure(fig_counter)
    ax = fig.gca(projection='3d')
    ax.set_xlim(lon_vector_plot[0], lon_vector_plot[-1])
    ax.set_ylim(lat_vector_plot[0], lat_vector_plot[-1])
    plt.xticks(np.arange(lon_vector_plot[0], lon_vector_plot[-1] + 20, 20))
    plt.yticks(np.arange(lat_vector_plot[0], lat_vector_plot[-1] + 10, 10))
    ax.view_init(azim=90, elev=90)
    plt.show()
    pdb.set_trace()
    lon_vector_plot, lat_vector_plot = np.meshgrid(lon_vector_plot, lat_vector_plot)
    surf = ax.plot_surface(lon_vector_plot, lat_vector_plot, diff_matrix_percent, edgecolor='k', linewidth=1)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # ax.set_zlabel('Percentage difference (%)')
    plt.title('Percentage difference for ' + dat_name)
    # pdb.set_trace()

    # plt.grid(axis='both', alpha=0.75)
    # bar = fig.colorbar(surf)
    # bar.set_label('Percentage difference (%)')

    if fig_counter == 1:
        print diff_matrix
        pdb.set_trace()
        plt.show()
        pdb.set_trace()
    fig_counter = fig_counter + 1
