import os
from _scandir import scandir
from natsort import natsorted
import numpy as np
import re
import matplotlib.pyplot as plt
import pdb

# iteration_path = 'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_6\\Iteration_1'
base_path = 'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_6'
plot_type = 'run_6_not_coupled'
viscosity_folder = os.path.join(base_path, 'viscosity_plots')
if not os.path.exists(viscosity_folder):
    os.mkdir(viscosity_folder)

with open(os.path.join(
        'C:\\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_6\\Iteration_1\\stress_matrices\\Large_stress_matrix.csv')) as opened_mises_stress:
    mises_stress = opened_mises_stress.readlines()
    mises_stress = [line.strip() for line in mises_stress[1:]]
    mises_stress = [np.array([eval(j) for j in line.split(",")[:]]) for line in mises_stress]
    mises_stress = np.array(mises_stress)
new_e = 'new_e.csv'

with open(os.path.join(base_path, 'no_iter_e.dat'), 'r+') as opened_e:
    e = opened_e.readlines()
    # e = [line.strip() for line in e]
    # e = [np.array([eval(j) for j in line.split(",")[:]]) for line in e]
    # e = np.array(e)
    # e = map(float, opened_e.read().split('\n'))
    with open(os.path.join(base_path, new_e), "w") as new_read_dat:
        for e_line in opened_e:
            replacements = e_line.replace(',', ' ').strip()
            print replacements
            new_read_dat.write(replacements)

pdb.set_trace()
elements = e[:, 0]
alin = e[:, 1]
a = e[:, 2]
an = 3.5
qtild = mises_stress[:, 1]
strain_rate = np.dot(alin, qtild) + np.dot(a, (qtild ** an))
viscosity = np.divide(qtild, strain_rate)

plt.figure()
plt.plot(elements, viscosity, 'r-')
plt.xlabel('Element N.')
plt.ylabel('Viscosity [Pa s]')
plt.grid(axis='both', alpha=0.75)
plt.savefig(os.path.join(viscosity_folder, 'Viscosity_' + plot_type), bbox_inches='tight')