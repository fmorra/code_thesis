from matplotlib import *


def geographic_plots(files_to_classify,report_path,classified_path):

    import numpy as np
    import math
    import glob
    import os
    from natsort import natsorted

    figure_folder = 'stress_plots'
    figures_path = os.path.join(report_path, figure_folder)
    if ~os.path.exist(figures_path):
        os.mkdir(figures_path)
    dfolders = [f.path for f in os.scandir(classified_path) if f.is_dir()]
    dnames = natsorted(dfolders)
    possible_depths = np.zeros(len(dnames), 1)
    regex = re.match(r'\d+')

    for i in range(len(dnames)):
        if ~math.isnan(re.match(regex, dnames[i])):
            possible_depths[i] = re.match(regex, dnames[i])
        else:
            possible_depths[i] = 0
            possible_depths[possible_depths == 0] = []

    selected_parts = possible_depths(1)

    for i in range(len(selected_parts)):
        classified_files_path = glob.glob(os.path.join(classified_path, 'Subclassification_depth_', str(selected_parts[i]), '_km', '*.csv'))
        classified_files = {classified_files_path.name}
        all_depths = np.zeros(len(classified_files), 1)
        for j in range(len(classified_files)):
            # all_depths(j) = str2double(extractAfter(classified_files{j}, 'bin_'));
            all_depths[j] = re.match(regex, classified_files[j])

        print('The layer bins for the selcted depth, ' + str(selected_parts[i]) + ' km, are:')
        print(all_depths)
        selected_depths = input('Enter the desired layer depth(s) for stress plotting, separated by commas:\n').split(',')
        if len(np.intersect1d(selected_depths, all_depths)) < len(selected_depths):
            print 'One or more values are outside the possible bins, program will be terminated.'
            break
        selected_component = input('Enter the desired stress component(s) separated by commas for stress plotting, possible values are S11, S22, S33, S12, S13, S23:\n')
        selected_components = selected_component.split(',')
        selected_columns = np.arr(np.zeros(len(selected_components), 1))
        for k in range(len(selected_components)):
            if selected_components[i] == 'S11':
                selected_columns[k] = 3
            elif selected_components[i] == '22':
                selected_columns[k] = 4
            elif selected_components[i] == '33':
                selected_columns[k] = 5
            elif selected_components[i] == 'S12':
                selected_columns[k] = 6
            elif selected_components[i] == 'S13':
                selected_columns[k] = 7
            else:
                selected_columns[k] = 8
        for k in range(len(selected_depths)):
            matrix_to_read = open(os.path.join(classified_path, 'subclassification_depth_', str(selected_parts[i]), '_km', 'Stresses_layer_',
                                               str(selected_parts[i]), '_km_bin_', str(selected_depths[k]), '.csv'))
        for l in range(len(selected_columns)):
            stress_to_plot = matrix_to_read[:, selected_columns[l]]
            lon_to_plot = matrix_to_read[:, -1]  # latitude[deg]
            lat_to_plot = matrix_to_read[:, -2]  # longitude[deg]
            # Figures here, continue when you know how to plot maps in python