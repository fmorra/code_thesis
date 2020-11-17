import os
import platform
import time
import pdb
import numpy as np
import math as m
from stress_report_reader import *
from deflection_report_reader import *
from coordinate_reader import *
from depth_classifier import *
from matlab_variables_writer import *
from element_tracker import *

# General parameters: based on how the files have been produced in ABAQUS, create the name of the dat and rpt files by
# modifying the run number
program_start_time = time.time()
run = '11'
iteration = 1
step = 0
cycle = 2
if not isinstance(run, str):
    run = str(run)
stress_rpt_name = 'abaqus_run_' + run + '.rpt'
deflection_rpt_name = 'abaqus_run_' + run + '_deflections.rpt'
dat_name = 'Earth.dat'
e_dat_name = 'e.dat'
# Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
# windows or linux. This is the only path that has to be modified
# base_path = '/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run)
base_path = os.path.join('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
# Define the path of the dat and rpt files to read
dat_path = os.path.join(base_path, dat_name)

# Add additional levels because we need to work on while cycles on the timesteps contained in each iteration. The report
# files and all the relative files will be saved, from now on, in that folder. level_dir is the variable that determines
# at what path level the files will be stored.

iteration_path = os.path.join(base_path, 'Iteration_' + str(iteration))
if not os.path.exists(iteration_path):
    os.mkdir(iteration_path)
step_path = os.path.join(iteration_path, 'step_' + str(step))
if not os.path.exists(step_path):
    os.mkdir(step_path)
cycle_path = os.path.join(step_path, 'cycle_' + str(cycle) + '_reports')
if not os.path.exists(cycle_path):
    os.mkdir(cycle_path)
level_dir = cycle_path
# Define the paths of the stress and coordinate files and if they are not found create these folders

stress_matrices_path = os.path.join(level_dir, 'stress_matrices')
stress_report_path = os.path.join(level_dir, stress_rpt_name)
stress_matrix = os.path.join(stress_matrices_path, 'EARTH.csv')
base_e_path = os.path.join(level_dir, e_dat_name)
if not os.path.exists(stress_matrices_path):
    os.mkdir(stress_matrices_path)
headers_on = 1


#######################################################################################################################
# Here begins the part to eb transferred to ABAQUS, the previous code is already there
stress_search_key = 'reported at element'
deflection_search_key = 'reported at nodes'
file_identifier = 'EARTH'
node_vector = []

with open(dat_path, "r+") as read_dat:
    # Read the file as a list of lines
    opened_dat = read_dat.readlines()
    elem_vector = np.zeros((len(opened_dat), 9))
    print "Processing the dat file to make it more readable"
    # Initialize counters containing a series of lines in the dat file where we will perform the search for node
    # coordinates and element nodes
    node_lines = []
    elem_lines = []
    part_lines = []
    nset_lines = []
    earth_lines = []
    node_keyword_counter = 0
    elem_keyword_counter = 0
    part_keyword_counter = 0
    nset_keyword_counter = 0
    file_identifiers = []
    # Define the keywords used to search for the lines delimiting the start and end of those regions
    node_search_key = '*Node'
    element_search_key = '*Element'
    part_search_key = 'PART INSTANCE'
    region_finisher_key = '*Nset'
    input_paragraph_end = 'OPTIONS BEING PROCESSED'
    part_earth_search_key = '** PART INSTANCE: EARTH'
    line_counter = 0
    end_line = 0
    end_line_counter = 0
    for line in opened_dat:
        line = line.strip()
        line_counter += 1
        node_keyword_flag = line.find(node_search_key)
        elem_keyword_flag = line.find(element_search_key)
        part_keyword_flag = line.find(part_search_key)
        region_finisher_flag = line.find(region_finisher_key)
        part_earth_search_flag = line == part_earth_search_key
        if part_earth_search_flag is True:
            earth_lines.append(line_counter)
        end_keyword_flag = line.find(input_paragraph_end)
        if node_keyword_flag != -1:
            node_keyword_counter = node_keyword_counter + 1
            node_lines.append(line_counter)
        elif elem_keyword_flag != -1:
            elem_keyword_counter = elem_keyword_counter + 1
            elem_lines.append(line_counter)
        elif region_finisher_flag != -1:
            nset_keyword_counter = nset_keyword_counter + 1
            nset_lines.append(line_counter)
        elif part_keyword_flag != -1:
            part_keyword_counter = part_keyword_counter + 1
            part_lines.append(line_counter)
        elif end_keyword_flag != -1:
            end_line = line
            end_line_counter = line_counter
            break
        else:
            pass
    node_lines.append(end_line_counter)
    for i in range(0, len(node_lines) - 1):
        individual_matrix = []
        for line_save in range(node_lines[i], elem_lines[i]):
            sentence = opened_dat[line_save].strip()
            a = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
            a = np.array([float(x) for x in a])
            if len(a) > 3:
                if a[0] is not 0 and a[0] == a[1] and a[1] == a[2] and a[2] == a[3]:
                    pass
                else:
                    if len(a) == 5:
                        a = a[1:]
                    individual_matrix.append(a)
        node_vector.append(individual_matrix)
    node_dictionary = {}
    for node_matrix in node_vector:
        for j in range(len(node_matrix)):
            dictionary_matrix = node_matrix[j]
            node_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
    node_ids = node_dictionary.keys()
    sorted_matrix = np.zeros((len(node_ids), 4))
    sorted_node_ids = list(sorted(node_ids))
    for node in range(len(sorted_node_ids)):
        sorted_matrix[node, 0] = sorted_node_ids[node]
        sorted_matrix[node, 1:] = node_dictionary[sorted_node_ids[node]]
    for i in range(0, len(node_lines) - 1):
        individual_matrix = []
        for line_save in range(node_lines[i], elem_lines[i]):
            sentence = opened_dat[line_save].strip()
            a = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
            a = np.array([float(x) for x in a])
            if len(a) > 3:
                if a[0] is not 0 and a[0] == a[1] and a[1] == a[2] and a[2] == a[3]:
                    pass
                else:
                    if len(a) == 5:
                        a = a[1:]
                    individual_matrix.append(a)
    new_elem_lines = elem_lines[0:len(node_lines) - 1]
    new_nset_lines = np.zeros((len(node_lines) - 1))
    for i in range(0, len(new_elem_lines)):
        for j in range(i, len(nset_lines)):
            if nset_lines[j] > new_elem_lines[i]:
                new_nset_lines[i] = nset_lines[j]
                break
    new_elem_lines = np.array(new_elem_lines)
    new_elem_lines = new_elem_lines.astype(int)
    new_nset_lines = new_nset_lines.astype(int)

    for line_save_2 in range(new_elem_lines[-1], new_nset_lines[-1]):
        sentence = opened_dat[line_save_2].strip()
        a = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
        a = np.array([float(x) for x in a])
        if len(a) > 6:
            if len(a) == 10:
                a = a[1:]
            elif len(a) == 8:
                a = a[1:]
            if len(a) == 9:
                elem_vector[line_save_2 - new_elem_lines[-1], :] = a
            else:
                elem_vector[line_save_2 - new_elem_lines[-1], 0:7] = a

node_matrix = sorted_matrix[~np.all(sorted_matrix == 0, axis=1)]
elem_matrix = elem_vector[~np.all(elem_vector == 0, axis=1)]
part_matrix = open(stress_matrix).readlines()
part_matrix = [line.strip() for line in part_matrix[1:]]
part_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in part_matrix]
part_matrix = np.array(part_matrix)
centroid_coord = np.zeros((len(part_matrix), 4))
pdb.set_trace()
if os.path.exists(os.path.join(base_path, 'force_data_matrix.csv')):
    print('Centroid coordinate file already created.')
else:
    for j in range(len(elem_matrix)):
        elem_nodes_index = [np.where(elem_matrix[:, 0] == part_matrix[j, 0])]
        elem_nodes = elem_matrix[elem_nodes_index, 1:]
        elem_nodes = elem_nodes[np.nonzero(elem_nodes)]
        elem_coord_matrix = np.zeros((len(elem_nodes), 3))
        for k in range(0, len(elem_nodes)):
            index = [np.where(node_matrix[:, 0] == elem_nodes[k])]
            nodes_coord = node_matrix[index, 1:]
            if not len(nodes_coord) == 0:
                elem_coord_matrix[k, :] = nodes_coord
        for m in range(0, 3):
            centroid_coord[j, m + 1] = np.mean(elem_coord_matrix[:, m])
        centroid_coord[j, 0] = part_matrix[j, 0]
        if np.all(centroid_coord[j, 1:]) == 0:
            centroid_coord = np.delete(centroid_coord, j, axis=0)
    data_matrix = np.hstack((part_matrix, centroid_coord[:, 1:]))
    with open(os.path.join(base_path, 'force_data_matrix.csv'), 'wb') as f_write:
        writer = csv.writer(f_write)
        writer.writerows(data_matrix)
# Calculate R and depth
pdb.set_trace()
data_matrix = open(os.path.join(base_path, 'force_data_matrix.csv')).readlines()
data_matrix = [line.strip() for line in data_matrix]
data_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in data_matrix]
data_matrix = np.array(data_matrix)
earth_radius = 6371000
depths = np.zeros((len(data_matrix), 1))
radial_centroid_distance = np.zeros((len(data_matrix), 1))
X = data_matrix[:, -3]
Y = data_matrix[:, -2]
Z = data_matrix[:, -1]
for k in range(len(data_matrix)):
    radial_centroid_distance[k] = np.sqrt(X[k] ** 2 + Y[k] ** 2 + Z[k] ** 2)
    depths[k] = earth_radius - radial_centroid_distance[k]
data_matrix = np.column_stack((data_matrix, radial_centroid_distance, depths))

# Discretize based on a certain number of bins
mantle_discontinuity = 660000
N_bins = 19
depth = data_matrix[:, -1]
bins = np.arange(0, depth.max()+(depth.max()-depth.min())/(2*N_bins), (depth.max()-depth.min())/(2*N_bins))
indices = np.digitize(depth, bins)
viscosity_edge_change = np.argmax(bins >= mantle_discontinuity) - 1  # Bin up until viscosity has the upper mantle value
data_matrix = np.column_stack((data_matrix, indices))

# Definition and discretization of the scaling factor

u_0_cm = 2.5  # cm/yr
u_0 = u_0_cm*(1/(86400*365))*(1/100)  # m/s
plate_thickness = 30000  # As an average from the figure using the Moho depicted as black dashes
a_r = 8  # Aspect ratio
D = 2.786E6
mu_m = [3.5e20, 2e21]  # viscosity, two different values for upper and lower mantle
u_vec = np.linspace(-u_0, u_0, len(bins))
tau_h_scaling_factor = np.zeros((len(data_matrix), 1))
tau_n_scaling_factor = np.zeros((len(data_matrix), 1))
# pdb.set_trace()
for i in range(len(data_matrix)):
    bin_index = int(data_matrix[i, -1])
    if bin_index < viscosity_edge_change:
        mu = mu_m[0]
    else:
        mu = mu_m[1]
    tau_h_scaling_factor[i] = (mu*2*u_vec[bin_index])/D
    tau_n_scaling_factor[i] = (tau_h_scaling_factor[i]*a_r)/plate_thickness
pdb.set_trace()
data_matrix = np.column_stack((data_matrix, tau_h_scaling_factor, tau_n_scaling_factor))
# After calculating the scaling factors, define traction in geographical coordinates
T_module = 4E6  # MPa, taken from Tutu
alpha = 60  # Average angle, as appears from Tutu, between a meridian and the average traction direction over West
# Antarctica
alpha_cos = m.cos(m.radians(60))
alpha_sin = m.sin(m.radians(60))
# T = np.zeros((1, 3))
# T[0, 0] = T_module*alpha_cos
# T[0, 1] = T_module*alpha_sin
T = np.transpose(([[T_module*alpha_cos, T_module*alpha_sin, 0]]))
# This data is already in the ABAQUS reference system
sphere_norm = np.sqrt((2*X)**2 + (2*Y)**2 + (2*Z)**2)
unit_normal_matrix = np.column_stack((2*X/sphere_norm, 2*Y/sphere_norm, 2*Z/sphere_norm))

# Convert traction to ABAQUS coordinate system and obtain stresses, then extract components

S_matrix = np.zeros((len(data_matrix), 7))
for element in range(len(data_matrix[:, 0])):
    R_3D = np.sqrt(X[element] ** 2 + Y[element] ** 2 + Z[element] ** 2)
    R_azimuth = np.sqrt(X[element] ** 2 + Y[element] ** 2)
    cos_theta = X[element] / R_azimuth
    sin_theta = Y[element] / R_azimuth
    cos_phi = R_azimuth / R_3D
    sin_phi = Z[element] / R_3D

    transformation_tensor = np.zeros((3, 3))
    transformation_tensor[0, 0] = cos_theta * cos_phi
    transformation_tensor[1, 0] = cos_theta * sin_phi
    transformation_tensor[2, 0] = -sin_theta
    transformation_tensor[0, 1] = sin_theta * cos_phi
    transformation_tensor[1, 1] = sin_theta * sin_phi
    transformation_tensor[2, 1] = cos_theta
    transformation_tensor[0, 2] = sin_phi
    transformation_tensor[1, 2] = -cos_phi
    transformation_tensor[2, 2] = 0.0

    # T_abaqus = np.zeros((3, 1))
    T_abaqus = np.array(np.matmul(np.transpose(transformation_tensor), T))
    unit_normal = np.transpose(unit_normal_matrix[element][np.newaxis, :])
    S = np.matmul(unit_normal, np.transpose(T_abaqus))
    S_matrix[element, 1] = S[0, 0] * tau_n_scaling_factor[element]
    S_matrix[element, 2] = S[1, 1] * tau_n_scaling_factor[element]
    S_matrix[element, 3] = S[2, 2] * tau_n_scaling_factor[element]
    S_matrix[element, 4] = S[0, 1] * tau_h_scaling_factor[element]
    S_matrix[element, 5] = S[0, 2] * tau_h_scaling_factor[element]
    S_matrix[element, 6] = S[1, 2] * tau_h_scaling_factor[element]
    # data_matrix[element, 2:8] = data_matrix[element, 2:8] + S_matrix
    S_matrix[element, 0] = (1 / np.sqrt(2)) * np.sqrt((S_matrix[element, 1] - S_matrix[element, 2]) ** 2 +
                                                      (S_matrix[element, 2] - S_matrix[element, 3]) ** 2 +
                                                      (S_matrix[element, 3] - S_matrix[element, 1]) ** 2 + 6 *
                                                      (S_matrix[element, 4] ** 2 + S_matrix[element, 6] ** 2 +
                                                       S_matrix[element, 5] ** 2))

data_matrix = np.column_stack((data_matrix, S_matrix))
pdb.set_trace()
# Finally, save the indices and corresponding scaling factor next to the stress values
with open(os.path.join(base_path, 'complete_force_data_matrix.csv'), 'wb') as f_write:
    writer = csv.writer(f_write)
    writer.writerows(data_matrix)

