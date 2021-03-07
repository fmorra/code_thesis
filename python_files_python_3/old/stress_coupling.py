from first_iter_actions import *


def stress_coupling(stress_matrices_path, coupled_stress_folder, headers_on, new_stress_test, stress_headers,
                    iteration, base_path, base_e_path, old_e_iter_path, new_e_iter_path):
    import os
    import numpy as np
    import csv
    import pdb

    # Copy what happens in MATLAB, and add the possibility of modifying the csv or dat file with Bdiff and Bdisl.
    # if iteration == 1:
    #     no_iter_path = first_iter_actions(base_path, base_e_path, new_e_iter_path)
    no_iter_path = first_iter_actions(base_path, base_e_path, new_e_iter_path)
    file_extension = 'csv'
    dir_csv_stress_files = os.listdir(stress_matrices_path)
    csv_stress_files = [filename for filename in dir_csv_stress_files if filename.endswith(file_extension)]
    large_new_matrix_path = os.path.join(stress_matrices_path, 'Large_Stress_Matrix.csv')
    large_new_matrix = []
    if os.path.exists(new_e_iter_path):
        pass
    else:
        selected_component = input('Enter the desired stress component(s) to modify, possible values are Mises, '
                                   'S11, S22, S33, S12, S13, S23:\n')
        selected_components = selected_component.split()
        selected_columns = np.zeros((len(selected_components), 1))
        for k in range(len(selected_components)):
            if selected_components[k] == 'Mises':
                selected_columns[k] = 1
            elif selected_components[k] == 'S11':
                selected_columns[k] = 2
            elif selected_components[k] == 'S22':
                selected_columns[k] = 3
            elif selected_components[k] == 'S33':
                selected_columns[k] = 4
            elif selected_components[k] == 'S12':
                selected_columns[k] = 5
            elif selected_components[k] == 'S13':
                selected_columns[k] = 6
            else:
                selected_columns[k] = 7

        for file_to_evaluate in range(len(csv_stress_files)):
            with open(os.path.join(stress_matrices_path, csv_stress_files[file_to_evaluate])) as matrix:
                matrix_to_evaluate = matrix.readlines()
                matrix_to_evaluate = [line.strip() for line in matrix_to_evaluate[1:]]
                matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in matrix_to_evaluate]
                matrix_to_evaluate = np.array(matrix_to_evaluate)
                if 'Region' in csv_stress_files[file_to_evaluate]:
                    csv_stress_files[file_to_evaluate] = 'to_delete'
                if not matrix_to_evaluate[1:, 1:].any():
                    csv_stress_files[file_to_evaluate] = 'to_delete'
        csv_stress_files = filter(lambda a: a != 'to_delete', csv_stress_files)

        for i in range(len(csv_stress_files)):
            with open(os.path.join(stress_matrices_path, csv_stress_files[i])) as matrix:
                stresses_to_modify = matrix.readlines()
                stresses_to_modify = [line.strip() for line in stresses_to_modify[1:]]
                stresses_to_modify = [np.array([eval(elem) for elem in line.split(",")[:]])
                                      for line in stresses_to_modify]
                stresses_to_modify = np.array(stresses_to_modify)
                new_stress_length = np.ones((len(stresses_to_modify), 1))
                new_column = np.zeros((len(stresses_to_modify), ))
                for j in range(len(selected_columns)):
                    column_to_add = new_stress_length * new_stress_test
                    for k in range(len(stresses_to_modify[:, int(selected_columns[j])])):
                        new_column[k] = stresses_to_modify[k, int(selected_columns[j])] + int(column_to_add[k])
                    stresses_to_modify[:, int(selected_columns[j])] = new_column
                print matrix
                new_matrix = stresses_to_modify
                mises_stress = (1 / np.sqrt(2)) * np.sqrt((new_matrix[:, 2] - new_matrix[:, 3]) ** 2 +
                                                          (new_matrix[:, 3] - new_matrix[:, 4]) ** 2 +
                                                          (new_matrix[:, 4] - new_matrix[:, 2])
                                                          ** 2 + 6 * (new_matrix[:, 5] ** 2 + new_matrix[:, 7] +
                                                                      new_matrix[:, 6] ** 2))
                new_matrix[:, 1] = mises_stress
                large_new_matrix.append(new_matrix)
                print i
                print csv_stress_files[i]
                with open(os.path.join(coupled_stress_folder, csv_stress_files[i]), 'wb') as f_write:
                    if headers_on == 1:
                        writer = csv.writer(f_write)
                        writer.writerow(stress_headers)
                        writer.writerows(new_matrix)
                    else:
                        writer = csv.writer(f_write)
                        writer.writerows(new_matrix)

        # new_dictionary = {}
        # for new_matrix in large_new_matrix:
        #     for j in range(len(new_matrix)):
        #         dictionary_matrix = new_matrix[j]
        #         new_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
        # new_ids = new_dictionary.keys()
        # sorted_matrix = np.zeros((len(new_ids), 8))
        # sorted_new_ids = list(sorted(new_ids))
        # for new in range(len(sorted_new_ids)):
        #     sorted_matrix[new, 0] = sorted_new_ids[new]
        #     sorted_matrix[new, 1:] = new_dictionary[sorted_new_ids[new]]
        # if headers_on == 1:
        #     with open(os.path.join(coupled_stress_folder, 'Large_stress_matrix.csv'), 'wb') as f_write:
        #         writer = csv.writer(f_write)
        #         writer.writerow(stress_headers)
        #         writer.writerows(sorted_matrix)
        # else:
        #     with open(os.path.join(coupled_stress_folder, 'Large_stress_matrix.csv'), 'wb') as f_write:
        #         writer = csv.writer(f_write)
        #         writer.writerows(sorted_matrix)

        # new_e_column = sorted_matrix[:, 1]
        # if iteration == 1:
        #     old_e_iter_path = no_iter_path
        # with open(old_e_iter_path, "r") as f:
        #     data = f.readlines()
        #     data = [line.strip() for line in data]
        #     data = [line.strip().replace(" ", ",") for line in data]
        #     data = [np.array([eval(i) for i in line.split(",")[:]]) for line in data]
        #     data = np.array(data)
        # data[:, -1] = new_e_column
        # data = np.matrix(data)
        # with open(new_e_iter_path, 'wb') as f_save:
        #     np.savetxt(f_save, data, fmt="%s")
        # with open(base_e_path, 'wb') as f_save:
        #     np.savetxt(f_save, data, fmt="%s")

        with open(os.path.join(coupled_stress_folder, 'EARTH.csv')) as matrix:
            e_stresses = matrix.readlines()
            e_stresses = [line.strip() for line in e_stresses[1:]]
            e_stresses = [np.array([eval(elem) for elem in line.split(",")[:]])
                          for line in e_stresses]
            e_stresses = np.array(e_stresses)

        new_e_column = e_stresses[:, 1]
        if iteration == 1:
            old_e_iter_path = no_iter_path
        with open(old_e_iter_path, "r") as f:
            data = f.readlines()
            data = [line.strip() for line in data]
            data = [line.strip().replace(" ", ",") for line in data]
            data = [np.array([eval(i) for i in line.split(",")[:]]) for line in data]
            data = np.array(data)
        data[:, -1] = new_e_column
        data = np.matrix(data)
        with open(new_e_iter_path, 'wb') as f_save:
            np.savetxt(f_save, data, fmt="%s")
        with open(base_e_path, 'wb') as f_save:
            np.savetxt(f_save, data, fmt="%s")

    # with open(new_e_iter_path, 'wb') as f_write:
    #     writer = csv.writer(f_write)
    #     writer.writerows(e_file)



            

