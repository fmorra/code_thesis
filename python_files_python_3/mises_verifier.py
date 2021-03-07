import os
import numpy as np

data_matrix = open('C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_21\\Iteration_1\\step_1\\'
                   'cycle_1_reports\\stress_matrices\\EARTH.csv').readlines()
data_matrix = [line.strip() for line in data_matrix[1:]]
data_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in data_matrix]
data_matrix = np.array(data_matrix)
mises_stress = (1 / np.sqrt(2)) * np.sqrt((data_matrix[:, 2] - data_matrix[:, 3]) ** 2 +
                                          (data_matrix[:, 3] - data_matrix[:, 4]) ** 2 +
                                          (data_matrix[:, 4] - data_matrix[:, 2]) ** 2 +
                                          6 * (data_matrix[:, 5] ** 2 + data_matrix[:, 7] ** 2 +
                                               data_matrix[:, 6] ** 2))
diff = abs(mises_stress - data_matrix[:, 1])
rel_diff = diff / data_matrix[:, 1] * 100
print(max(rel_diff))
print(max(diff))
print(min(diff))
