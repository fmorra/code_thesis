def python_plotter(sd_input, files_to_classify, base_path, iteration_path, classified_path, components_to_plot, run,
                   coupled_stress_folder, e_dat_name, fig_counter):
    import os
    from _scandir import scandir
    from natsort import natsorted
    import numpy as np
    import re
    import matplotlib.pyplot as plt
    from mpl_toolkits import basemap
    import pdb

    # file_extension = '.csv'
    # while_check = 1
    # while while_check == 1:
    #     plots_bool = input('Do you want to plot the stress results or go to deflection processing? \n'
    #                        'Enter 1 for plots, 0 to skip this step: \n')
    #     if plots_bool == 1 or plots_bool == 0:
    #         while_check = 0
    #     else:
    #         print 'Incorrect input, enter 1 for the stress plots and 0 to go to deflections'
    #
    # if plots_bool == 1:
    #     if sd_input == 0:
    #         matrix_to_read_part = 'Stresses_layer_'
    #         figure_folder = 'stress_plots'
    #     else:
    #         matrix_to_read_part = 'Deflections_layer_'
    #         figure_folder = 'deflection_plots'
    #
    #     figures_path = os.path.join(iteration_path, figure_folder)
    #     if not os.path.exists(figures_path):
    #         os.mkdir(figures_path)
    #     depth_folders = [f.name for f in scandir(classified_path) if f.is_dir()]
    #     depth_names = natsorted(depth_folders)
    #     possible_depths = np.zeros((len(depth_names), 1))
    #     pdb.set_trace()
    #     for i in range(len(depth_names)):
    #         print i
    #         if re.compile(r'\d+').findall(depth_names[i]):
    #         #     if len(re.compile(r'\d+').findall(depth_names[i])) > 1:
    #         #         if float(re.compile(r'\d+').findall(depth_names[i])[1][0]) > 5 or float\
    #         #                     (re.compile(r'\d+').findall(depth_names[i])[1][0]) == 5:
    #         #             possible_depths[i] = float(re.compile(r'\d+').findall(depth_names[i])[0])
    #         #         else:
    #         #             possible_depths[i] = float(re.compile(r'\d+').findall(depth_names[i])[0])
    #         #     else:
    #             possible_depths[i] = re.compile(r'\d+').findall(depth_names[i])[0]
    #             print possible_depths[i]
    #         else:
    #             possible_depths[i] = 0
    #     possible_depths = filter(lambda b: b != 0, possible_depths)
    #
    #     print 'Select model layers for stress plotting, possible values are (km):'
    #     print possible_depths
    #     selected_parts = input('\nInput them a vector in square brackets: \n')
    #     print 0
    #     for i in range(len(selected_parts)):
    #         print 1
    #         classified_files_path = os.path.join(classified_path, "Subclassification_depth_['" + str(selected_parts[i])
    #                                              + "']_km")
    #         print 2
    #         classified_files = [filename for filename in os.listdir(classified_files_path) if
    #                             filename.endswith(file_extension)]
    #         print 3
    #         classified_files = natsorted(classified_files)
    #         print 4
    #         all_depths = np.zeros(len(classified_files), 1)
    #         print 5
    #         for j in range(len(classified_files)):
    #             all_depths[j] = re.compile(r'\d+').findall(classified_files[j].split('bin')[1])
    #         print 6
    #
    #         print('The layer bins for the selected depth(s), ', str(selected_parts[i]), ' km, are:')
    #         print all_depths
    #         selected_depths = input('\nEnter the desired bin(s) for stress plotting, as a vector in square brackets:\n')
    #         if len(np.intersect1d(selected_depths, all_depths)) < len(selected_depths):
    #             print 'One or more values are outside the possible bins, program will be terminated.'
    #             break
    #     pdb.set_trace()
    #     if sd_input == 0:
    #         selected_component = input('Enter the desired stress component(s) to plot, possible values are Mises, S11, '
    #                                    'S22, S33, S12, S13, S23:\n')
    #         selected_components = selected_component.split(" ")
    #         selected_columns = np.zeros(len(selected_components), 1)
    #         for k in range(len(selected_components)):
    #             if selected_components[k] == 'Mises':
    #                 selected_columns[k] = 1
    #             elif selected_components[k] == 'S11':
    #                 selected_columns[k] = 2
    #             elif selected_components[k] == 'S22':
    #                 selected_columns[k] = 3
    #             elif selected_components[k] == 'S33':
    #                 selected_columns[k] = 4
    #             elif selected_components[k] == 'S12':
    #                 selected_columns[k] = 5
    #             elif selected_components[k] == 'S13':
    #                 selected_columns[k] = 6
    #             else:
    #                 selected_columns[k] = 7
    #     else:
    #         selected_component = input('Enter the desired deflection components(s) to plot, possible values are '
    #                                     'Magnitude, U1, U2, U3:\n')
    #         selected_components = selected_component.split(" ")
    #         selected_columns = np.zeros(len(selected_components), 1)
    #         for k in range(len(selected_components)):
    #             if selected_components[k] == 'Magnitude':
    #                 selected_columns[k] = 2
    #             elif selected_components[k] == 'U1':
    #                 selected_columns[k] = 3
    #             elif selected_components[k] == 'U2':
    #                 selected_columns[k] = 4
    #             else:
    #                 selected_columns[k] = 5
    #
    #     ### Plot definition - stresses and deflections
    #


        ### Plot definition - viscosity

        # All the values we need are contained in the e.dat file. I should compare the results from run 6 before
        # combining, run 6 after combining, and run 7. What will change is the Mises stress in each of the files.

    plot_type = 'run_6_not_coupled'
    viscosity_folder = os.path.join(iteration_path, 'viscosity_plots')
    if not os.path.exists(viscosity_folder):
        os.mkdir(viscosity_folder)

    with open(os.path.join('C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_6\\Iteration_1\\stress_matrices\\Large_stress_matrix.csv')) as opened_mises_stress:
        mises_stress = opened_mises_stress.readlines()
        mises_stress = [line.strip() for line in mises_stress[1:]]
        mises_stress = [np.array([eval(j) for j in line.split(",")[:]]) for line in mises_stress]
        mises_stress = np.array(mises_stress)
    with open(os.path.join(base_path, 'no_iter_e.dat')) as opened_e:
        e = opened_e.readlines()
        e = [line.strip() for line in e]
        e = [np.array([eval(j) for j in line.split(",")[:]]) for line in e]
        e = np.array(e)
    elements = e[:, 0]
    alin = e[:, 1]
    a = e[:, 2]
    an = 3.5
    qtild = mises_stress[:, 1]
    strain_rate = np.dot(alin, qtild) + np.dot(a, (qtild**an))
    viscosity = np.divide(qtild, strain_rate)

    plt.figure(fig_counter)
    plt.plot(elements, viscosity, 'r-')
    plt.xlabel('Element N.')
    plt.ylabel('Viscosity [Pa s]')
    plt.grid(axis='both', alpha=0.75)
    plt.savefig(os.path.join(viscosity_folder, 'Viscosity_' + plot_type), bbox_inches='tight')













