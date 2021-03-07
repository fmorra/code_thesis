# e_0 =
import os
import numpy as np
import pandas
import pdb
import re
import csv

e_1 = open('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\New folder\\e_inp_1.dat').readlines()
e_2 = open('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\New folder\\e_inp_2.dat').readlines()
e_1_matrix = np.zeros((len(e_1), 4))
e_2_matrix = np.zeros((len(e_2), 4))
for line in range(len(e_1)):
    s_1 = e_1[line].strip()
    s_2 = e_2[line].strip()
    # Search for floats at each line
    e_1_data = [float(s_1) for s_1 in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s_1)]
    e_1_matrix[line, :] = np.array([float(x) for x in e_1_data])
    e_2_data = [float(s_2) for s_2 in re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s_2)]
    e_2_matrix[line, :] = np.array([float(x) for x in e_2_data])
inp_0 = np.zeros((len(e_1), 1))
inp_0_larger = np.ones((len(e_1), 1))*1e4
inp_1 = np.transpose(e_1_matrix[:, -1][np.newaxis, :])
inp_2 = np.transpose(e_2_matrix[:, -1][np.newaxis, :])
diff_1_0 = inp_1 - inp_0
diff_1_0_larger = inp_1 - inp_0_larger
diff_2_1 = inp_2 - inp_1
data_matrix = np.column_stack((diff_1_0, diff_1_0_larger, diff_2_1))
with open('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\New folder\\diff_file.csv', 'wb') as f_write:
    writer = csv.writer(f_write)
    writer.writerows(data_matrix)
