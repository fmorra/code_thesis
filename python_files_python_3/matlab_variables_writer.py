def matlab_variables_writer(sd_input, complete_files_path, files_to_classify, iteration_path, components_to_plot):

    import csv
    import os
    import pdb

    if sd_input == 0:
        quantity_to_plot = 'stresses'
    else:
        quantity_to_plot = 'deflections'
    second_list = [quantity_to_plot, iteration_path, components_to_plot, complete_files_path]
    # depth = area_params[-1]
    # min_depth = depth[0]
    # max_depth = depth[1]
    variable_list = files_to_classify + second_list
    for i in range(len(variable_list)):
        if isinstance(variable_list[i], int):
            variable_list[i] = str(variable_list[i])

    # with open(os.path.join(iteration_path, 'matlab_variable_names' + '.csv'), 'wb') as f_write:
    #     writer = csv.writer(f_write)
    #     writer.writerow(variable_list)
    print(iteration_path)
    with open(os.path.join(iteration_path, 'matlab_' + components_to_plot + '_' + quantity_to_plot + '_' +
                                           'variables.txt'), 'wb') as f_write:
        f_write.write('\n'.join(variable_list))

    print('File with variable names necessary for MATLAB plotting has been written.')
