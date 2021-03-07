import os
import numpy as np
import pdb
import math


e_path = 'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_15\\Iteration_1\\step_0' \
         '\\cycle_1_reports\\e.dat'
matrix_to_open_path = ''
with open(e_path, "r+") as e:
    opened_dat = e.readlines()
    for line in opened_dat:
        line = line.replace('e', 'E')
pdb.set_trace()
opened_dat = [line.strip() for line in opened_dat]
opened_dat = [np.array([eval(i) for i in line.split(",")[:]]) for line in opened_dat]
opened_dat = np.array(opened_dat)
elements = opened_dat[:, 0]
alin = opened_dat[:, 1]
a = opened_dat[:, 2]
mises = opened_dat[:, 3]
threshold = 1E-35
for i in alin:
    if i < threshold:
        i = i * (math.floor(math.log(i, 10))/threshold)
pdb.set_trace()
a_coeff = 1E40
an = 3.5
strain_rate = np.zeros((len(opened_dat), 1))
n = an
for j in range(len(strain_rate)):
    strain_rate[j, 1] = alin[j] * mises[j] + a[j] * mises[j] ** pow
viscosity = [1.0 / [3 * sr] for sr in strain_rate]