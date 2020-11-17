import os
import numpy as np

data_matrix = open('C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_11\\'
                   'force_data_matrix_v2').readlines()
mises_stress = (1/np.sqrt(2))*np.sqrt((data_matrix[: ,2] - data_matrix[: ,3]).^2+(data_matrix(:, 4)-data_matrix(:, 5)).^2+(data_matrix(:, 5)-data_matrix(:, 3)).^2+...
    6*(new_matrix(:, 6).^2 + new_matrix(:, 8).^2+new_matrix(:, 7).^2))
diff = mises_stress - data_matrix[:, 2]