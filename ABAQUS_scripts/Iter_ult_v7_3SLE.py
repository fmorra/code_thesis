## A multi-layer model iterator for Abaqus-CAE, apply after the the model generated
## Author: Bas Blank & Haiyang Hu
## Build : 19-05-2014
## Latest Build: 7-1-2019


import numpy as np
import math
import os
import platform
import sys
import pdb
import shutil
import time
import csv
import glob
import re
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import*
print 'abaqus modules successfully imported before iteration starts'
Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)

import Model_data2_top
from Model_data2_top import *
import Model_data2
from Model_data2 import *
os.chdir(complete_dir)
print os.getcwd()

# Number of Iteration
# Total number of iterations, it has been proven that convergence with a number already between 3 and 5
Res=3
## Iteration
# Starting iteration, incremented at each step, needed for the while condition
Num=1
CPUs = 16
mdb=openMdb(Model_name)

# If else cycle is only for testing phase, remove it afterwards

while Num<=Res:   # We do not give a convergence condition because we know convergence is reached in 4 or 5 steps

    ####################################################################################
    # Check if the job for each iteration has already been completed or not to save time
    if os.path.exists('Iteration_' + str(Num) + '_job_completion_certificate.txt'):
        print 'Job already submitted and completed, skipping this action.'
    else:
        print 'Start FEM,round '+str(Num)
        with open(os.path.join(complete_dir, 'iteration_counter.txt'), 'w') as writefile:
            writefile.write(str(Num) + '\n')
        if Num == 1:
            fortran_path = os.path.join(Dir2, 'user_2.f')
        else:
            fortran_path = os.path.join(Dir2, 'user_2_mises.f')
        mdb=openMdb(Model_name)
        ## Remove lock-file if present
        if os.path.isfile(Model_name+'.lck'):
            os.remove(Model_name+'.lck') 
        mdb.Job(name=Model_name, model=Model_name, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=ON, 
        modelPrint=ON, contactPrint=ON, historyPrint=ON, userSubroutine=fortran_path, 
        scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        ## Define Output
        mdb.models[Model_name].fieldOutputRequests['F-Output-1'].setValues(variables=('U','S','MISES','COORD'))
        ## Save model database (this was commented at first)
        mdb.saveAs(pathName=os.path.join(complete_dir,Model_name+'.cae'))
        print 'Model job has been defined.'
        mdb.jobs[Model_name].submit(consistencyChecking=OFF)
        mdb.jobs[Model_name].waitForCompletion()
        with open('Iteration_' + str(Num) + '_job_completion_certificate.txt', 'wb') as file:
            file.write('Job completed for iteration ' + str(Num))
##        if Num == 1:
##            os.mkdir(complete_dir, 'base_files')
##            if not os.path.exists(os.path.join(complete_dir, 'base_files', Model_name + '.odb' )):
##                os.mkdir(os.path.join(complete_dir, 'base_files')
##                shutil.copyfile(os.path.join(complete_dir, Model_name + '.odb' ), os.path.join(complete_dir, 'base_files', Model_name + '.odb' ))
##            if os.path.exists(os.path.join(complete_dir, 'base_files')):
##                os.mkdir(os.path.join(complete_dir, 'base_files'))
        print 'Finish FEM\n'
    ## This part is usualy skipped
    if initialize_interpol_lists==1 and Num==1:
        print 'Start initializing interpolation lists\n'
        # The output database containing the results is read
        odb=openOdb(Model_name+ '.odb')

        # We extract all the node sets for all layers and assign them to a new vector based on their name, which is in
        # turn based on the presence or not of coarse layers.
        Node_set=[0]*N_layer
        for i in range(N_layer):
            if i > Coarse_Layers:
                Node_set[i]=odb.rootAssembly.instances[Model_name.upper()].nodeSets['IF'+str(i)]
            else:
                Node_set[i]=odb.rootAssembly.instances[Model_name.upper()+'_LOW'].nodeSets['IF'+str(i)]
                
        ## horizontal speeds are defined for an y-north axis system!!
        nodes_length=np.zeros(N_layer)
        fileout = open('Node_list_new.dat','w')

        # For each layer, we check whether we are at the surface or not. If not, we extract the nodes of the layer
        # and save the coordinates of each node in a file called data_out. For the surface, the data_out remains full of
        # zeros. We then close the file and save it in a dat file.
        for j in range(N_layer):
            if Density[j+1]-Density[j]<0:
                L1=len(Node_set[j].nodes)
                nodes_length[j]=L1
                Data_out=np.zeros([L1,3])
                for k in range(L1):
                    Data_out[k,0:3]=Node_set[j].nodes[k].coordinates                
            else:
                L1=len(Node_set[j].nodes)
                Data_out=np.zeros([L1,3])
                nodes_length[j]=L1
            np.savetxt(fileout,Data_out)
        fileout.close() 
        np.savetxt(os.path.join(complete_dir, 'Nodes_sets.dat'), nodes_length)
        print 'Start matlab process\n'
 #        os.system('cd C:\Program Files\MATLAB\R2014b\bin')
        CHECK=os.system("matlab.exe -nodisplay -nosplash -nodesktop -r run('ConvertNodeLocationList_v2.m');") ##1
        print 'Finish matlab process\n'
        print 'Start initializing interpolation lists. Check if MATLAB is finished before the sub-python process has started\n'
              
    ##############################################################################################
    # Run preliminary operations to calculate initial deflections at each step if not done already
    if os.path.exists('Iteration_' + str(Num) + '_preliminary_operations_completion_certificate.txt'):
        print ('The preliminary operations on the odb for iteration number ' + str(Num) + ' have already been carried out, moving to iterative procedure.')
    else:
        print ('Start odb process\n')
        
    ##        print ('Start matlab process\n')
    ##        os.system('cd C:\Program Files\MATLAB\R2014b\bin')
    ##        CHECK=os.system("matlab.exe -nodisplay -nosplash -nodesktop -r run('ConvertNodeLocationList_v2.m');") ##1
    ##        print ('Finish matlab process\n')
    ##        print ('Start initializing interpolation lists. Check if MATLAB is finished before the sub-python process has started\n')
    
        odb=openOdb(path=Model_name+ '.odb')   
        Node_set=[0]*N_layer
        for i in range(N_layer):
            if i > Coarse_Layers:
                Node_set[i]=odb.rootAssembly.instances[Model_name.upper()].nodeSets['IF'+str(i)]
            else:
                Node_set[i]=odb.rootAssembly.instances[Model_name.upper()+'_LOW'].nodeSets['IF'+str(i)]
        #### con2
        nodes_length=np.zeros(N_layer)
        fileout=open('Data_output.dat','w')
        S_0 = np.loadtxt('S0.dat')
        fileout2 = open('Data_output_hor.dat','w')
        
        for i in range(N_step):        
            ### Compute initial deflection
            # All of this is done more than once at every iteration,  because we iterate for every larger iteration over the
            # ice history.
            U=odb.steps['Step'+str(i+1)].frames[-1].fieldOutputs['U']
            ## horizontal speeds are defined for an y-north axis system!!
            # Then for each layer we check if it is the surface or not. If not, we extract the deflection on the region U1
            # which is the subset of nodes pertaining to that specified layer, then we create two save files. Data_out
            # now has 4 columns and te first three contain the nodes coordinates while the fourth one the sum of the product
            # between the coordinates of the nodes and vector a, divided by the radius vector of that node. Data_output_hor
            # has 6 columns, where the first 3 are the node coordinates and the last 3 are the deflections at each node.
            # For the surface, as before, they are full of zeros. We then save these outputs
            # These results are stored for each node.
            U1=[0]*N_layer
            for j in range(N_layer):
                if Density[j+1]-Density[j]<0:
                    U1=U.getSubset(region=Node_set[j])
                    L1=len(Node_set[j].nodes)
                    nodes_length[j]=L1
                    Data_out=np.zeros([L1,4])
                    Data_out_Hor=np.zeros([L1,6])
                    for k in range(L1):
                        Data_out[k,0:3]=Node_set[j].nodes[k].coordinates
                        Data_out_Hor[k,0:3]=Node_set[j].nodes[k].coordinates
                        X=Node_set[j].nodes[k].coordinates
                        a=U1.values[k].data
                        Data_out_Hor[k,3:6]=a
                        radius = np.linalg.norm(X)
                        Data_out[k,3]=(X*a).sum()/radius                    
                else:
                    L1=len(Node_set[j].nodes) # amount of nodes for each layer
                    Data_out=np.zeros([L1,4])
                    Data_out_Hor=np.zeros([L1,6])
                    nodes_length[j]=L1
                np.savetxt(fileout,Data_out)
                np.savetxt(fileout2,Data_out_Hor)
        fileout.close()
        fileout2.close()                                                                                                                   
        np.savetxt(os.path.join(complete_dir, 'Nodes_config.dat'), nodes_length)
        session.odbs[Model_name+'.odb'].close()
        with open('Iteration_' + str(Num) + '_preliminary_operations_completion_certificate.txt','wb') as file:
            file.write('Preliminary operations on the odb for iteration number ' + str(Num) + ' carried out.')
    
    #################################################################################################################
    ####################################New part added: stress pairing at each timestep##############################
    # stop
    stress_rpt_name = 'abaqus_run_'+ str(run_number) + '.rpt'
    deflection_rpt_name = 'abaqus_run_'+ str(run_number) + '_deflections.rpt'
    iter_dir = os.path.join(complete_dir, 'Iteration_'+ str(Num))

    if not os.path.exists(iter_dir):
        os.mkdir(iter_dir)
    for step_counter in range(N_step): # for step_counter in range(1, N_step + 1):
        step_dir = os.path.join(iter_dir, 'step_' + str(step_counter))
        if not os.path.exists(step_dir):
            os.mkdir(step_dir)
        diff = 1e6
        cycle_counter = 1
        while diff > 100 and cycle_counter < 4: 
        # This section is the main file on the local machine, and defines many parameters used in the next sections, such as file paths and file names.
        # General parameters: based on how the files have been produced in ABAQUS, create the name of the dat and rpt files by modifying the run number
            headers_on = 1
            stress_search_key = 'reported at element'
            deflection_search_key = 'reported at nodes'
            start_time_reports = time.time()
            # Define the subfolder relative to the current cycle
            cycle_dir = os.path.join(step_dir, 'cycle_' + str(cycle_counter) + '_reports')
            if not os.path.exists(cycle_dir):
                os.mkdir(cycle_dir)
            level_dir = cycle_dir
            run = run_number
            if not isinstance(run, str):
                run = str(run)
            if run == 'point':
                dat_name = 'Earth_point.dat'
            else:
                dat_name = 'Earth.dat'
            iteration = Num
            e_dat_name = 'e.dat'
            # Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
            # windows or linux. This is the only path that has to be modified
            opersys = platform.system()
            if opersys == 'Windows':
                base_path = os.path.join('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
            else:
        ##        base_path = os.path.join('/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run))
                base_path = complete_dir # was complete_dir - go down 1 level
            # Specify where to save and read the e.dat file for the FORTRAN routine and specify which folder to use as a base for this iteration, step and cycle
            dat_path = os.path.join(base_path, dat_name)
            base_e_path = os.path.join(base_path, e_dat_name)
            iteration_path = level_dir 
            
            # Define the paths inside the folder relative to the current iteration, step and cycle where to save the stresses and deflection matrices 
            stress_report_path = os.path.join(iteration_path, stress_rpt_name)
            deflection_report_path = os.path.join(iteration_path, deflection_rpt_name)
            stress_matrices_path = os.path.join(iteration_path, 'stress_matrices')
            stress_processing_path = os.path.join(iteration_path, 'stress_processing_results')
            deflection_matrices_paths = os.path.join(iteration_path, 'deflection_matrices')
            deflection_processing_path = os.path.join(iteration_path, 'deflection_processing_results')
            common_files_path = os.path.join(iteration_path, 'common_files')
            coupled_stress_folder = os.path.join(stress_matrices_path, 'coupled_stress_matrices')
            
            if not os.path.exists(stress_matrices_path):
                os.mkdir(stress_matrices_path)
            if not os.path.exists(stress_processing_path):
                os.mkdir(stress_processing_path)
            if not os.path.exists(deflection_matrices_paths):
                os.mkdir(deflection_matrices_paths)
            if not os.path.exists(deflection_processing_path):
                os.mkdir(deflection_processing_path)
            
            # If they do not exist, generate the report files for stresses and deflections.
            if os.path.exists(os.path.join(level_dir, deflection_rpt_name)):
                print 'Skipping report generation, they are already there.'
            else:
                os.chdir(level_dir)
                base_model = 'Earth'
                # Define output parameters
                step_to_access = step_counter + 1
                Model_name = base_model
                if os.path.isfile(Model_name+'.lck'):
                    os.remove(Model_name+'.lck')
                # In this case we can only generate the field output reports, as the job definition and all the forces are already defined in the iter file. 
                odb = openOdb(path=os.path.join(complete_dir, Model_name+'.odb'))
                session.viewports['Viewport: 1'].setValues(displayedObject = odb)
                selected_step = odb.steps['Step'+str(step_to_access)]
                selected_frame = selected_step.frames[-1]
                print('Generating report files')
                # Generate stress report
                session.writeFieldReport(fileName = stress_rpt_name, append = OFF, sortItem = "Element Label", odb = odb, step = step_to_access, frame = selected_frame, outputPosition = ELEMENT_CENTROID, 
                                         variable = (('S',INTEGRATION_POINT,((INVARIANT,'Mises'),(COMPONENT,'S11'),(COMPONENT,'S22'),(COMPONENT,'S33'),(COMPONENT,'S12'),(COMPONENT,'S13'),(COMPONENT,'S23'),)),)) 
                print("Stress report generated.")
                # Generate deflection report
                session.writeFieldReport(fileName = deflection_rpt_name, append = OFF, sortItem = "Node Label", odb = odb, step = step_to_access, frame = selected_frame, outputPosition = NODAL, 
                                         variable = (('U',NODAL,),))
                print("Deflection report generated.")
                os.chdir(complete_dir)
            # stop
            if cycle_counter == 2 and step_counter == 1:
                stop_step_2
    ############################# Combine the two stresses together and, if wanted, process this data to plot the results of this coupling ####################################
    ############ Read stress reports #######################################################################################################
        # Here we read the stress report file and create csvs containing the values for each part, then we combine the csvs for part EARTH 
        # into a larger file (because that part is saved as divided in regions)
        
        ##    Here we do not need to change anything, the paths are updated in the first section of the modifications.
            
            report_file = stress_report_path
            lines = []
            keyword_counter = 0
            file_identifiers = []
            headers = []
            line_counter = 0
            headers_key = 'S.Mises'
            # Loop over each line of the rpt file to save the filenames and csv file headers
            
            with open(report_file) as opened_report: 
                for line in opened_report:
                    line_counter += 1
                    # Find the keyword that is found at the start of each part of the report file that contains new stresses for
                    # each part or region
                    keyword_flag = line.find(stress_search_key)
                    # If the keyword is found, increase the number of counters (for verification), add the keyword line to an array
                    # containing them, extract only the part of the line after : without newline and replace the dots in this
                    # extracted part with an underscore. This part of the line will be used as the csv file name.
                    if keyword_flag != -1:
                        keyword_counter = keyword_counter + 1
                        lines.append(line_counter)
                        csv_name = line.split(': ')[1].rstrip().replace(".", "_")
                        file_identifiers.append(csv_name)
                    # If we also want to extract the headers for all the columns, find the header key in each file row and if it is
                    # found extract the line after the word "Element" stripping it of the newline and then substitute the dots
                    # with underscores.
                    headers_index = line.find(headers_key)
                    if headers_index != -1:
                        headers_line = line.split('Element ')[1].rstrip().replace(".", "_")
                        headers.append(headers_line)
        
            headers = headers[0]
            # Insert the 'Labels' string at the beginning of the list containing the headers
            headers = "Label " + headers
            headers = headers.split()
            stress_headers = headers
            if os.path.exists(os.path.join(stress_matrices_path, 'I1.csv')):
                pass
            else:
                print('Stress matrices generation checkpoint')
                large_stress_matrix = []
                large_stress_matrix_path = os.path.join(stress_matrices_path, 'EARTH.csv')
                print large_stress_matrix_path
                # Run the main script if the last of the csv files to be generated is not there yet
                if os.path.exists(large_stress_matrix_path):
                    print('The files containing the stress components already exist, moving on'
                          ' to stress coordinate extraction for each part.')
                else:
                    print('Processing rpt to extract stress components...')
                    # Iterate over each of the lines where a new part or region is defined
                    opening_report_start_time = time.time() 
                    for i in range(len(lines)):
                        stress_submatrix = []
                        # Determine whether we are in the last region of the file or not
                        if i + 1 < len(lines):
                            # Define the lines where to search for stresses
                            search_lines = np.arange(lines[i], lines[i+1], 1)
                            line_counter = 0
                            # Iterate over these lines searching for float values
                            with open(report_file) as opened_report: 
                                for line in opened_report:
                                    line_counter += 1
                                    if line_counter in search_lines:
                                        s = line.strip()
                                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
                                        # If the data vector is composed of 8 elements, fill the matrix containing the stress values
                                        # with it
                                        if len(data) > 7:
                                            stress_submatrix.append(np.array([float(x) for x in data]))
                            # Save the matrix containing the stress values for each part and region differentiating between whether
                            # we want headers or not
                            print 'Saving csv file for part: ' + file_identifiers[i]
                            if headers_on == 1:
                                with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerow(headers)
                                    writer.writerows(stress_submatrix)
                            else:
                                with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerows(stress_submatrix)
                            if 'Region' in file_identifiers[i]:
                                large_stress_matrix.append(stress_submatrix)
                        else:
                            # For the last region, iterate until the end of the file, then the rest works the same as before
                            with open(report_file) as opened_report: 
                                line_counter = 0
                                for line in opened_report:
                                    line_counter += 1
                                last_line = line_counter
                                line_counter = 0
                                search_lines = np.arange(lines[i], last_line, 1)
                                for line in opened_report:
                                    line_counter += 1
                                    if line_counter in search_lines:
                                        s = line.strip()
                                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)
                                        # if the data is not empty, fill the matrix containing the stress values with it
                                        if len(data) > 7:
                                            stress_submatrix.append(np.array([float(x) for x in data]))
                            print 'Saving csv file for part: ' + file_identifiers[i]
                            if headers_on == 1:
                                with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerow(headers)
                                    writer.writerows(stress_submatrix)
                            else:
                                with open(os.path.join(stress_matrices_path, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerows(stress_submatrix)
                            if 'Region' in file_identifiers[i]:
                                large_stress_matrix.append(stress_submatrix)
                    opening_report_end_time = time.time() - opening_report_start_time
                    print 'Time elapsed to open the reports is ' + str(opening_report_end_time)      
                    
                    stress_dictionary = {}
                    for stress_matrix in large_stress_matrix:
                        for j in range(len(stress_matrix)):
                            dictionary_matrix = stress_matrix[j]
                            stress_dictionary[dictionary_matrix[0]] = dictionary_matrix[1:]
                    stress_ids = stress_dictionary.keys()
                    sorted_matrix = np.zeros((len(stress_ids), 8))
                    sorted_stress_ids = list(sorted(stress_ids))
                    for stress in range(len(sorted_stress_ids)):
                        sorted_matrix[stress, 0] = sorted_stress_ids[stress]
                        sorted_matrix[stress, 1:] = stress_dictionary[sorted_stress_ids[stress]]
                    if headers_on == 1:
                        with open(os.path.join(large_stress_matrix_path), 'wb') as f_write:
                            writer = csv.writer(f_write)
                            writer.writerow(headers)
                            writer.writerows(sorted_matrix)
                    else:
                        with open(os.path.join(large_stress_matrix_path), 'wb') as f_write:
                            writer = csv.writer(f_write)
                            writer.writerows(sorted_matrix)
            
        ## Here, recalculate the difference vale for the validity of the while cycle. The difference is calculated based on the
        ## values in the stress report at two different iterations; this has to be skipped for the first while cycle, because the report is generated 
        ## before the stresses are coupled. EARTH is the part used to calculate the sum of the Mises stresses.
            
            
            if not cycle_counter == 1:
                diff_matrix_previous_iteration = os.path.join(os.path.join(step_dir, 'cycle_' + str(cycle_counter-1) + '_reports/stress_matrices'), 'EARTH.csv')
                with open(diff_matrix_previous_iteration) as matrix:
                    previous_matrix_to_evaluate = matrix.readlines()
                    previous_matrix_to_evaluate = [line.strip() for line in previous_matrix_to_evaluate[1:]]
                    previous_matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in previous_matrix_to_evaluate]
                    previous_matrix_to_evaluate = np.array(previous_matrix_to_evaluate)
                    mises_previous_iteration = previous_matrix_to_evaluate[:,1]
                diff_matrix_current_iteration = os.path.join(stress_matrices_path, 'EARTH.csv')
                with open(diff_matrix_current_iteration) as matrix:
                    current_matrix_to_evaluate = matrix.readlines()
                    current_matrix_to_evaluate = [line.strip() for line in current_matrix_to_evaluate[1:]]
                    current_matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in current_matrix_to_evaluate]
                    current_matrix_to_evaluate = np.array(current_matrix_to_evaluate)
                    mises_current_iteration = current_matrix_to_evaluate[:,1]
                diff = abs(np.sum(mises_current_iteration) - np.sum(mises_previous_iteration))
            elif (cycle_counter == 1 and step_counter == 0 and Num == 1):
                pass
            else:
                # Extract the highest cycle number in the previous step
                previous_step = step_counter-1
                step_sentence = 'step_' + str(previous_step)
                all_files = os.listdir(complete_dir)
                txt_files = filter(lambda x: x[-4:] == '.txt', all_files)
                step_txt_files = []
                for file in txt_files:
                    if step_sentence in file:
                        step_txt_files.append(file)
                previous_max_cycle = 0
                for filename in step_txt_files:
                    triplet = re.compile(r'\d+').findall(filename)
                    cycle = triplet[-1]
                    if cycle > previous_max_cycle:
                        previous_max_cycle = cycle
                diff_matrix_previous_iteration = os.path.join(os.path.join(iter_dir, 'step_' + str(previous_step), 'cycle_' + str(previous_max_cycle)
                                                                           + '_reports/stress_matrices'), 'EARTH.csv')
                #diff_matrix_previous_iteration = os.path.join(os.path.join(iter_dir, 'step_' + str(step_counter-1), 'cycle_3_reports/stress_matrices'), 'EARTH.csv')
                # Proceed as in the first condition
                with open(diff_matrix_previous_iteration) as matrix:
                    previous_matrix_to_evaluate = matrix.readlines()
                    previous_matrix_to_evaluate = [line.strip() for line in previous_matrix_to_evaluate[1:]]
                    previous_matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in previous_matrix_to_evaluate]
                    previous_matrix_to_evaluate = np.array(previous_matrix_to_evaluate)
                    mises_previous_iteration = previous_matrix_to_evaluate[:,1]
                diff_matrix_current_iteration = os.path.join(stress_matrices_path, 'EARTH.csv')
                with open(diff_matrix_current_iteration) as matrix:
                    current_matrix_to_evaluate = matrix.readlines()
                    current_matrix_to_evaluate = [line.strip() for line in current_matrix_to_evaluate[1:]]
                    current_matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in current_matrix_to_evaluate]
                    current_matrix_to_evaluate = np.array(current_matrix_to_evaluate)
                    mises_current_iteration = current_matrix_to_evaluate[:,1]
                diff = abs(np.sum(mises_current_iteration) - np.sum(mises_previous_iteration))
            
            # Read deflection reports ###################################################################################################

            ## Here, too, we do not have to change anything, because the reports will always be in the same format and the paths are 
            ## updated in the first section of the modifications.
            if os.path.exists(os.path.join(deflection_matrices_paths, 'I1.csv')):
                pass
            else:
                print('Deflection matrices generation checkpoint')
                
                report_file_deflections = deflection_report_path
                lines = []
                keyword_counter = 0
                file_identifiers = []
                headers = []
                line_counter = 0
                headers_key = 'U.Magnitude'
            
                # Loop over each line of the rpt file
                with open(report_file_deflections) as opened_report:
                    for line in opened_report:
                        line_counter += 1
                        # Find the keyword that is found at the start of each stress region in the stress file
                        keyword_flag = line.find(deflection_search_key)
                        # If the keyword is found, increase the number of counters (for verification), add the keyword line to an array
                        # containing them, extract only the part of the line after : without newline and replace the dots in this
                        # extracted part with an underscore.
                        if keyword_flag != -1:
                            keyword_counter = keyword_counter + 1
                            lines.append(line_counter)
                            csv_name = line.split(': ')[1].rstrip().replace(".", "_")
                            file_identifiers.append(csv_name)
                        # If we also want to extract the headers for all the columns, find the header key in each file row  and if it is
                        # found extract the line after Element stripping it of the newline and then substitute the dots with underscores
                        headers_index = line.find(headers_key)
                        if headers_index != -1:
                            headers_line = line.split('Node ')[1].rstrip().replace(".", "_")
                            headers.append(headers_line)
                headers = headers[0]
                # Insert the 'Labels' string at the beginning of the list containing the headers
                headers = "Label " + headers
                headers = headers.split()
                # Run the main script if the last of the csv files to be generated is not there yet
                if os.path.exists(os.path.join(deflection_matrices_paths, 'I1.csv')):
                    print('The files containing the stress components already exist, moving on'
                          ' to nodes coordinate extraction for each part.')
                else:
                    print('Processing rpt to extract deflection components...')
                    # Iterate over each of the lines where a new part or region is defined
                    for i in range(len(lines)):
                        deflection_submatrix = []
                        # Determine whether we are in the last region of the file or not
                        if i + 1 < len(lines):
                            # Define the lines where to search for deflections
                            search_lines = np.arange(lines[i], lines[i+1], 1)
                            line_counter = 0
                            # Iterate over these lines searching for float values
                            with open(report_file_deflections) as opened_report:
                                print opened_report
                                for line in opened_report:
                                    line_counter += 1
                                    if line_counter in search_lines:
                                        s = line.strip()
                                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", s)
                                        # If the data vector is composed of 5 elements, fill the matrix containing the stress values
                                        # with it
                                        if len(data) > 4:
                                            deflection_submatrix.append(np.array([float(x) for x in data]))
                            # Save the matrix containing the stress values for each part and region differentiating between whether
                            # we want headers or not
                            if headers_on == 1:
                                with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerow(headers)
                                    writer.writerows(deflection_submatrix)
                            else:
                                with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerows(deflection_submatrix)
                        else:
                            # for the last region, iterate until the end of the file, then the rest works the same as before
                            line_counter = 0
                            with open(report_file_deflections) as opened_report:
                                for line in opened_report:
                                    line_counter += 1
                                last_line = line_counter
                                line_counter = 0
                                search_lines = np.arange(lines[i], last_line, 1)
                                for line in opened_report:
                                    line_counter += 1
                                    if line_counter in search_lines:
                                        s = line.strip()
                                        data = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[ee][+-]?\d+)?", s)
                                        # if the data is not empty, fill the matrix containing the stress values with it
                                        if len(data) > 4:
                                            deflection_submatrix.append(np.array([float(x) for x in data]))
            
                            if headers_on == 1:
                                with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerow(headers)
                                    writer.writerows(deflection_submatrix)
                            else:
                                with open(os.path.join(deflection_matrices_paths, file_identifiers[i] + ".csv"), 'wb') as f_write:
                                    writer = csv.writer(f_write)
                                    writer.writerows(deflection_submatrix)
            
            end_time_reports = time.time() - start_time_reports
            print 'The time spent to read the stress and deflection reports and extract data from them is: ', end_time_reports, \
                ' seconds'
##            if cycle_counter == 3:
##                stop
            ######################################################################################################################
            #################################### Define the stresses to couple ###################################################
##            stress_search_key = 'reported at element'
##            deflection_search_key = 'reported at nodes'
##            file_identifier = 'EARTH'
##            node_vector = []
##            
##            with open(dat_path, "r+") as read_dat:
##                # Read the file as a list of lines
##                opened_dat = read_dat.readlines()
##                elem_vector = np.zeros((len(opened_dat), 9))
##                node_vector = np.zeros((len(opened_dat), 4))
##                print "Processing the dat file to make it more readable"
##                # Initialize counters containing a series of lines in the dat file where we will perform the search for node
##                # coordinates and element nodes
##                elem_lines = []
##                nset_lines = []
##                # Define the keywords used to search for the lines delimiting the start and end of those regions
##                element_search_key = '*Element'
##                region_finisher_key = '*Nset'
##                input_paragraph_end = 'OPTIONS BEING PROCESSED'
##                part_earth_search_key = '** PART INSTANCE: EARTH'
##                line_counter = 0
##                end_line = 0
##                end_line_counter = 0
##                for line in opened_dat:
##                    line = line.strip()
##                    line_counter += 1
##                    elem_keyword_flag = line.find(element_search_key)
##                    region_finisher_flag = line.find(region_finisher_key)
##                    part_earth_search_flag = line == part_earth_search_key
##                    end_keyword_flag = line.find(input_paragraph_end)
##                    if part_earth_search_flag is True:
##                        starting_earth_line = line_counter
##                    elif elem_keyword_flag != -1:
##                        elem_lines.append(line_counter)
##                    elif region_finisher_flag != -1:
##                        nset_lines.append(line_counter)
##                    elif end_keyword_flag != -1:
##                        end_line = line
##                        end_line_counter = line_counter
##                        break
##                    else:
##                        pass
##                for e in range(len(elem_lines)):
##                    if elem_lines[e] < starting_earth_line:
##                        elem_lines[e] = 'a'
##                elem_lines = filter(lambda b: b != 'a', elem_lines)
##                for line_save in range(starting_earth_line, elem_lines[0]):
##                    sentence = opened_dat[line_save].strip()
##                    a = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
##                    a = np.array([float(x) for x in a])
##                    if len(a) > 3:
##                        if a[0] is not 0 and a[0] == a[1] and a[1] == a[2] and a[2] == a[3]:
##                            pass
##                        else:
##                            if len(a) == 5:
##                                a = a[1:]
##                            node_vector[line_save - starting_earth_line, :] = a
##                for e in range(len(nset_lines)):
##                    if nset_lines[e] < elem_lines[0]:
##                        nset_lines[e] = 'a'
##                nset_lines = filter(lambda b: b != 'a', nset_lines)
##                for line_save_2 in range(elem_lines[0], nset_lines[0]):
##                    sentence = opened_dat[line_save_2].strip()
##                    a = [float(s) for s in re.findall(r'-?\d+\.?\d*', sentence)]
##                    a = np.array([float(x) for x in a])
##                    if len(a) > 6:
##                        if len(a) == 10:
##                            a = a[1:]
##                        elif len(a) == 8:
##                            a = a[1:]
##                        if len(a) == 9:
##                            elem_vector[line_save_2 - elem_lines[0], :] = a
##                        else:
##                            elem_vector[line_save_2 - elem_lines[0], 0:7] = a
##            node_matrix = node_vector[~np.all(node_vector == 0, axis=1)]
##            elem_matrix = elem_vector[~np.all(elem_vector == 0, axis=1)]
##            part_matrix = open(stress_matrix).readlines()
##            part_matrix = [line.strip() for line in part_matrix[1:]]
##            part_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in part_matrix]
##            part_matrix = np.array(part_matrix)
##            centroid_coord = np.zeros((len(part_matrix), 4))
##            if os.path.exists(os.path.join(level_dir, 'force_data_matrix.csv')):
##                print('Centroid coordinate file already created.')
##            else:
##                for j in range(len(elem_matrix)):
##                    elem_nodes_index = [np.where(elem_matrix[:, 0] == part_matrix[j, 0])]
##                    elem_nodes = elem_matrix[elem_nodes_index, 1:]
##                    elem_nodes = elem_nodes[np.nonzero(elem_nodes)]
##                    elem_coord_matrix = np.zeros((len(elem_nodes), 3))
##                    for k in range(0, len(elem_nodes)):
##                        index = [np.where(node_matrix[:, 0] == elem_nodes[k])]
##                        nodes_coord = node_matrix[index, 1:]
##                        if not len(nodes_coord) == 0:
##                            elem_coord_matrix[k, :] = nodes_coord
##                    for m in range(0, 3):
##                        centroid_coord[j, m + 1] = np.mean(elem_coord_matrix[:, m])
##                    centroid_coord[j, 0] = part_matrix[j, 0]
##                    if np.all(centroid_coord[j, 1:]) == 0:
##                        centroid_coord = np.delete(centroid_coord, j, axis=0)
##                data_matrix = np.hstack((part_matrix, centroid_coord[:, 1:]))
##                with open(os.path.join(level_dir, 'force_data_matrix.csv'), 'wb') as f_write:
##                    writer = csv.writer(f_write)
##                    writer.writerows(data_matrix)
##            # Calculate R and depth
##            data_matrix = open(os.path.join(level_dir, 'force_data_matrix.csv')).readlines()
##            data_matrix = [line.strip() for line in data_matrix]
##            data_matrix = [np.array([eval(i) for i in line.split(",")[:]]) for line in data_matrix]
##            data_matrix = np.array(data_matrix)
##            earth_radius = 6371000
##            depths = np.zeros((len(data_matrix), 1))
##            radial_centroid_distance = np.zeros((len(data_matrix), 1))
##            X = data_matrix[:, -3]
##            Y = data_matrix[:, -2]
##            Z = data_matrix[:, -1]
##            for k in range(len(data_matrix)):
##                radial_centroid_distance[k] = np.sqrt(X[k] ** 2 + Y[k] ** 2 + Z[k] ** 2)
##                depths[k] = earth_radius - radial_centroid_distance[k]
##            data_matrix = np.column_stack((data_matrix, radial_centroid_distance, depths))
##            
##            # Discretize based on a certain number of bins
##            mantle_discontinuity = 660000
##            N_bins = 50
##            depth = data_matrix[:, -1]
##            bins = np.arange(0, depth.max()+(depth.max()-depth.min())/(2*N_bins), (depth.max()-depth.min())/(2*N_bins))
##            indices = np.digitize(depth, bins)
##            viscosity_edge_change = np.argmax(bins >= mantle_discontinuity)  # Bin up until viscosity has the upper mantle value
##            data_matrix = np.column_stack((data_matrix, indices))
##            # Definition and discretization of the scaling factor
##            
##            u_0_cm = 2.5  # cm/yr
##            u_0 = u_0_cm*(1.0/(86400*365))*(1.0/100)  # m/s
##            plate_thickness = 30000  # As an average from the figure using the Moho depicted as black dashes
##            a_r = 8  # Aspect ratio
##            D = 2.786E6
##            mu_m = [3.5e20, 2e21]  # viscosity, two different values for upper and lower mantle
##            u_vec = np.linspace(-u_0, u_0, len(bins))
##            tau_h = np.zeros((len(data_matrix), 1))
##            tau_n = np.zeros((len(data_matrix), 1))
##            # pdb.set_trace()
##            for i in range(len(data_matrix)):
##                bin_index = int(data_matrix[i, -1])
##                if bin_index < viscosity_edge_change:
##                    mu = mu_m[0]
##                else:
##                    mu = mu_m[1]
##                tau_h[i] = (mu*2*u_vec[bin_index])/D
##                tau_n[i] = (tau_h[i]*a_r)/plate_thickness
##            data_matrix = np.column_stack((data_matrix, tau_h, tau_n))
##            # After calculating the scaling factors, define traction in geographical coordinates
##            alpha = 60  # Average angle, as appears from Tutu, between a meridian and the average velocity direction over WA
##            theta = 90 - alpha
##            theta_cos = m.cos(m.radians(theta))
##            theta_sin = m.sin(m.radians(theta))
##            S_matrix = np.zeros((len(data_matrix), 7))
##            
##            u0_reference_taus = np.zeros((2, 2))
##            u0_T = np.zeros((2, 2))
##            u0_T[0, 0] = theta_cos
##            u0_T[1, 1] = theta_cos
##            u0_T[0, 1] = theta_sin
##            u0_T[1, 0] = -theta_sin
##            for element in range(len(data_matrix[:, 0])):
##                print element
##                u0_reference_taus = np.zeros((2, 2))
##                S_geo = np.zeros((3, 3))
##                u0_reference_taus[0, 0] = tau_n[element]
##                u0_reference_taus[1, 0] = tau_h[element]
##                u0_reference_taus[0, 1] = tau_h[element]
##                u0_reference_taus[1, 1] = tau_n[element]
##            
##                geographical_reference_taus = np.dot(u0_T, np.dot(u0_reference_taus, np.transpose(u0_T)))
##                R_3D = np.sqrt(X[element] ** 2 + Y[element] ** 2 + Z[element] ** 2)
##                R_azimuth = np.sqrt(X[element] ** 2 + Y[element] ** 2)
##                cos_theta = X[element] / R_azimuth
##                sin_theta = Y[element] / R_azimuth
##                cos_phi = R_azimuth / R_3D
##                sin_phi = Z[element] / R_3D
##            
##                transformation_tensor = np.zeros((3, 3))
##                transformation_tensor[0, 0] = cos_theta * cos_phi
##                transformation_tensor[1, 0] = cos_theta * sin_phi
##                transformation_tensor[2, 0] = -sin_theta
##                transformation_tensor[0, 1] = sin_theta * cos_phi
##                transformation_tensor[1, 1] = sin_theta * sin_phi
##                transformation_tensor[2, 1] = cos_theta
##                transformation_tensor[0, 2] = sin_phi
##                transformation_tensor[1, 2] = -cos_phi
##                transformation_tensor[2, 2] = 0.0
##            
##                S_geo[0, 0] = u0_reference_taus[0, 0]
##                S_geo[1, 0] = u0_reference_taus[1, 0]
##                # S_geo[2, 0] =
##                S_geo[0, 1] = u0_reference_taus[0, 1]
##                S_geo[1, 1] = u0_reference_taus[1, 1]
##                # S_geo[2, 1] =
##                # S_geo[0, 2] =
##                # S_geo[1, 2] =
##                S_geo[2, 2] = 1
##                # To cartesian reference frame
##                S = np.dot(np.transpose(transformation_tensor), np.dot(S_geo, transformation_tensor))
##                S_matrix[element, 1] = S[0, 0]
##                S_matrix[element, 2] = S[1, 1]
##                S_matrix[element, 3] = S[2, 2]
##                S_matrix[element, 4] = S[0, 1]
##                S_matrix[element, 5] = S[0, 2]
##                S_matrix[element, 6] = S[1, 2]
##                # To ABAQUS reference frame
##                S_matrix[element, 2], S_matrix[element, 3] = S_matrix[element, 3].copy(), S_matrix[element, 2].copy()
##                S_matrix[element, 1], S_matrix[element, 2] = S_matrix[element, 2].copy(), S_matrix[element, 1].copy()
##                S_matrix[element, 4], S_matrix[element, 5] = S_matrix[element, 5].copy(), S_matrix[element, 4].copy()
##                S_matrix[element, 5], S_matrix[element, 6] = S_matrix[element, 6].copy(), S_matrix[element, 5].copy()
##                S_matrix[element, 0] = (1 / np.sqrt(2)) * np.sqrt((S_matrix[element, 1] - S_matrix[element, 2]) ** 2 +
##                                                                  (S_matrix[element, 2] - S_matrix[element, 3]) ** 2 +
##                                                                  (S_matrix[element, 3] - S_matrix[element, 1]) ** 2 + 6 *
##                                                                  (S_matrix[element, 4] ** 2 + S_matrix[element, 6] ** 2 +
##                                                                   S_matrix[element, 5] ** 2))
##                
##            
##            data_matrix = np.column_stack((data_matrix, S_matrix))
##            # Finally, save the indices and corresponding scaling factor next to the stress values
##            with open(os.path.join(level_dir, 'complete_force_data_matrix.csv'), 'wb') as f_write:
##                writer = csv.writer(f_write)
##                writer.writerows(data_matrix)
            
            ######################################################################################################################
            
            step_counter_name = str(int(float(step_counter)))
            no_iter_path = os.path.join(base_path, 'no_iter_e.dat')
            if (Num == 1 and step_counter == 0 and cycle_counter == 1) and not os.path.exists('Iteration_' + str(Num) + '_step_' + str(step_counter_name) + '_cycle_' + str(cycle_counter) + '_e_creation_completion_certificate.txt'):
                shutil.copyfile(os.path.join(Dir2, e_dat_name), no_iter_path)
                with open(no_iter_path, "r") as f:
                    no_iter_data = f.readlines()
                    no_iter_data = [line.strip() for line in no_iter_data]
                    no_iter_data = [line.strip().replace("  ", " ").replace(" ", ",") for line in no_iter_data]
                    no_iter_data = [np.array([eval(i) for i in line.split(",")[:]]) for line in no_iter_data]
                    no_iter_data = np.array(no_iter_data)
                    new_column = np.zeros((len(no_iter_data), 1))
                    no_iter_data = np.column_stack((no_iter_data, new_column))
                no_iter_data = np.matrix(no_iter_data)
                # with open(new_e_iter_path, 'wb') as f_save:
                #     np.savetxt(f_save, data, fmt="%s")
                with open(base_e_path, 'wb') as f_save:
                    np.savetxt(f_save, no_iter_data, fmt="%s")
                with open('Iteration_' + str(Num) + '_step_' + str(step_counter_name) + '_cycle_' + str(cycle_counter) + '_e_creation_completion_certificate.txt', 'wb') as file:
                    file.write('Created e file.' + '\n')
            else:
                print('The base e file does not need to be genearted again')
            new_e_iter_path = os.path.join(cycle_dir, e_dat_name)
            
            # Define the new stress to couple 
            
            new_stress_test = 10000 #0.01MPa
            if not os.path.exists(coupled_stress_folder):
                os.mkdir(coupled_stress_folder)
     
            ##########################################################################################################################
            # Select the components for stress coupling
            
            if os.path.exists('Iteration_' + str(Num) + '_step_' + str(step_counter_name) + '_cycle_' + str(cycle_counter) + '_coupling_completion_certificate.txt'):
                print ('Coupling process already done for cycle ' + str(cycle_counter) + ', skipping this step')
            else:
                print('Stress coupling checkpoint')          
                file_extension = 'csv'
                dir_csv_stress_files = os.listdir(stress_matrices_path)
                csv_stress_files = [filename for filename in dir_csv_stress_files if filename.endswith(file_extension)]
                selected_component = 'S11 S12' # Only in one shear direction
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
                
            
                # Filter to find the files where to pair stresses: in this case, only EARTH 
                for file_to_evaluate in range(len(csv_stress_files)):
                    with open(os.path.join(stress_matrices_path, csv_stress_files[file_to_evaluate])) as matrix:
                        matrix_to_evaluate = matrix.readlines()
                        matrix_to_evaluate = [line.strip() for line in matrix_to_evaluate[1:]]
                        matrix_to_evaluate = [np.array([eval(i) for i in line.split(",")[:]]) for line in matrix_to_evaluate]
                        matrix_to_evaluate = np.array(matrix_to_evaluate)
                        if not matrix_to_evaluate.any():
                            csv_stress_files[file_to_evaluate] = 'to_delete'
                        if 'Region' in csv_stress_files[file_to_evaluate]:
                            csv_stress_files[file_to_evaluate] = 'to_delete'
                        if 'LOW' in csv_stress_files[file_to_evaluate]:
                            csv_stress_files[file_to_evaluate] = 'to_delete'
                        if 'I0' in csv_stress_files[file_to_evaluate]:
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
                        
##                        for j in range(len(S_matrix[:, 1:][0])):
##                            stresses_to_modify[:, j+1] = stresses_to_modify[:, j+1] + S_matrix[:, j]
                        new_matrix = stresses_to_modify
                        mises_stress = (1 / np.sqrt(2)) * np.sqrt((new_matrix[:, 2] - new_matrix[:, 3]) ** 2 +
                                                                  (new_matrix[:, 3] - new_matrix[:, 4]) ** 2 +
                                                                  (new_matrix[:, 4] - new_matrix[:, 2]) ** 2 + 
                                                                  6 * (new_matrix[:, 5] ** 2 + new_matrix[:, 7] ** 2 +
                                                                              new_matrix[:, 6] ** 2))
                        new_matrix[:, 1] = mises_stress
                        with open(os.path.join(coupled_stress_folder, csv_stress_files[i]), 'wb') as f_write:
                            if headers_on == 1:
                                writer = csv.writer(f_write)
                                writer.writerow(stress_headers)
                                writer.writerows(new_matrix)
                            else:
                                writer = csv.writer(f_write)
                                writer.writerows(new_matrix)
            
                with open(os.path.join(coupled_stress_folder, 'EARTH.csv')) as matrix:
                    e_stresses = matrix.readlines()
                    e_stresses = [line.strip() for line in e_stresses[1:]]
                    e_stresses = [np.array([eval(elem) for elem in line.split(",")[:]]) for line in e_stresses]
                    e_stresses = np.array(e_stresses)
            
                new_e_column = e_stresses[:, 1]
                with open(base_e_path, "r") as f:
                    data = f.readlines()
                    data = [line.strip() for line in data]
                    data = [line.strip().replace(" ", ",") for line in data]
                    data = [np.array([eval(i) for i in line.split(",")[:]]) for line in data]
                    data = np.array(data)
                data[:, -1] = new_e_column
                data = np.matrix(data)
                with open(new_e_iter_path, 'wb') as f_save:
                    np.savetxt(f_save, data, fmt="%s")
                os.remove(base_e_path)
                with open(base_e_path, 'wb') as f_save:
                    np.savetxt(f_save, data, fmt="%s")
                with open('Iteration_' + str(Num) + '_step_' + step_counter_name + '_cycle_' + str(cycle_counter) + '_coupling_completion_certificate.txt', 'wb') as file:
                    file.write('Coupling completed for iteration ' + str(Num) + ', step ' + step_counter_name + ' and cycle ' + str(cycle_counter) + '\n')
                    
            # stop
##            if cycle_counter == 3:
##                stop
            
            if os.path.exists('Iteration_' + str(Num) + '_step_' + str(step_counter_name) + '_cycle_' + str(cycle_counter) + '_job_completion_certificate.txt'):
                print 'Job already submitted and completed, skipping this action.'
            else:
                job_start_time = time.time()
                mdb.Job(name=Model_name, model=Model_name, description='', 
                type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=ON, 
                modelPrint=ON, contactPrint=ON, historyPrint=ON, userSubroutine= os.path.join(Dir2, 'user_2_mises.f'), 
                scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)        
                if os.path.isfile(Model_name+'.lck'):
                    os.remove(Model_name+'.lck') 
                mdb.jobs[Model_name].submit(consistencyChecking=OFF)
                mdb.jobs[Model_name].waitForCompletion()
                job_end_time = time.time() - job_start_time
                with open('Iteration_' + str(Num) + '_step_' + str(step_counter_name) + '_cycle_' + str(cycle_counter) + '_job_completion_certificate.txt', 'wb') as file:
                    file.write('Job completed for iteration ' + str(Num) + ', step ' + step_counter_name + ' and cycle ' + str(cycle_counter) + '\n')
                    if not (cycle_counter == 1 and step_counter == 0 and Num == 1):
                        file.write('Mises difference between this and previous cycle is: ' + str(diff) + '\n')
                    file.write('The time to complete the job was: ' + str(int(float(job_end_time))) + ' seconds' + '\n')
                    for line in range(0,20):
                        np.savetxt(file, data[line, :], fmt="%s")
                   
            cycle_counter = cycle_counter + 1
            
    ###################################################################################################################
    ###################################################################################################################
    
    # process the output
    if Sphere_int==1:
        print 'Start sub-Python process\n'
        CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2/sph_tools_SphereInt_v1_0_SLE.py') ##1
        print 'Finish sub-Python process\n'
	
    else:
        print 'Start sub-Python process\n'
        CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2/sph_tools_TPW7_4SLE.py') ##1
        print 'Finish sub-Python process\n'
    if CHECK!=0:
        sys.exit(0)
    
    Coor_sph=mdb.models[Model_name].rootAssembly.datums.values()[0]
    ## Regenerate new load
    for i in range(N_layer):
        exec("Data"+str(i)+"=np.loadtxt('Data_in'+str(i)+'.dat')")
        if i == (N_layer-1):
            DataIce=np.loadtxt(Dir+'DataIce_in.dat')                    
          
    exec("grid1=Data"+str(j)+".shape[1]/2") 
	
    theta=np.linspace(0,180,grid1+1)
    theta=np.hstack((0,theta))
    theta=theta.reshape(theta.size,1)
            
    lam=np.linspace(0,360,grid1*2+1)
    lam=np.hstack((lam[lam<=180],np.array([J-360 for J in lam if J>180])))   #in Abaqus, Th start from 0~pi then -pi~0.	

    for i in range(N_step):
        ### Sea level equation
        if SLE_true == 1:
            # Extract different loads based on whether we compensate for the spherical harmonics of degree 1 and 0 or not
            if DisableSH_0_1 ==1:
                print 'Start iterate sub-Python process to remove SHdegree 0 and 1\n'
                CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/sph_tools_initial_load_v2.py') ##1
                print 'Finish iterate sub-Python process\n'

                DataLoad=np.loadtxt('DataLoad_in.dat')
                
            else:  
                exec("SeaLevel = np.loadtxt('Data_in_Slvl.dat')")
                Sea_lvl_allTimes_ice_eq = SeaLevel*Density_water/Density_Ice

                DataLoad =  DataIce  +  Sea_lvl_allTimes_ice_eq

            # As seed in modelgen, we need to distinguish if the iteration is larger than 0, extracting loads in a
            # different part of the matrix and with a different force field for the 0th iteration. The force fields are
            # then alwys saved in a tuple.
            if i>0:
                if filter_elastic==1 and i==(N_step-1):
                    New_Force_field=np.vstack((lam,DataLoad[(grid1+1)*(i-1):(grid1+1)*i,:]))
                else:
                    New_Force_field=np.vstack((lam,DataLoad[(grid1+1)*i:(grid1+1)*(i+1),:]))
                New_Force_field=np.hstack((theta,New_Force_field))
    
                
                Old_Force_field = New_Force_field
                Force_field= tuple([tuple(J) for J in New_Force_field]) 
            else:
                Force_field=np.vstack((lam,DataLoad[(grid1+1)*i:(grid1+1)*(i+1),:]))
                Force_field=np.hstack((theta,Force_field))
                Old_Force_field = Force_field
                Force_field=tuple([tuple(J) for J in Force_field])

            # Again create a mapped field with the total load
            mdb.models[Model_name].MappedField(name='TotalLoad' +str(i)+'_field', description='',
               regionType=POINT, partLevelData=False, localCsys=Coor_sph,
               pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32,
               gridPointData={str(Radius[N_layer-1]):Force_field})
                 
                 
            ## Create iteration-step ice+sea load
            a = mdb.models[Model_name].rootAssembly
            region = a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]
            # Here we only work with the surface and with slightly different commands based on the timestep
            if i==0:
                mdb.models[Model_name].Pressure(name='TotalLoad_'+str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='TotalLoad'+str(i)+'_field', magnitude=Density_Ice*Gacc[N_layer-1],amplitude=UNSET)
        
            else:
                mdb.models[Model_name].TabularAmplitude(name='Amp-'+str(i), timeSpan=TOTAL, smooth=SOLVER_DEFAULT, data=((
                   Time[i-1], 1.0), (Time[i], 0.0)))
                    
                mdb.models[Model_name].Pressure(name='TotalLoad_'+str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='TotalLoad'+str(i)+'_field', magnitude=Density_Ice*Gacc[N_layer-1],amplitude=UNSET)
        
            if i>0:
                mdb.models[Model_name].loads['TotalLoad_'+str(i-1)].deactivate('Step'+str(i+1))

            

        for j in range(N_layer):
            if i==0:
                Gap=Time[0]
            else:
                Gap=Time[i]-Time[i-1]
            ## delete extra computation statement
            exec("phi1_matrix = np.loadtxt('Data_in"+str(j)+".dat')")
            phi1 = phi1_matrix[(grid1+1)*i:(grid1+1)*(i+1),:]
                
            if Num == Res:
                # geoid is calculated and saved
                fileout=open('Geoid_output.dat','w')
                geoid = phi1/Gacc[j]
                np.savetxt(fileout,geoid)
                fileout.close()
			
            exec("Data1_field=np.vstack((lam,Data"+str(j)+"[(grid1+1)*i:(grid1+1)*(i+1),:]))")
            
            if Compensate_cog ==1:
                # Again as in modelgen we define matrices for theta and phi, then we take and calculate the norm of the
                # center of gravity vector
                theta_mat           = np.tile(np.linspace(0,2*np.pi,(2*grid1+1)) ,[grid1+1, 1])
                phi_mat             = np.tile(np.transpose([np.linspace(np.pi/2.0,-np.pi/2.0,(grid1+1))]) ,[1, 2*grid1+1])
                r_cog_vec_steps     = np.loadtxt('CoG_vector.dat')
                #r_cog_vec = [r_cog_lm[0], r_cog_lm[1], r_cog_lm[2]]
                r_cog_vec   = r_cog_vec_steps[i,:]
                r_cog_norm  = np.linalg.norm(r_cog_vec)
            
                theta_cog   = np.arctan2(r_cog_vec[0],r_cog_vec[2])
                Delta_theta = theta_cog - theta_mat    
                # Correct center of gravity if the correction is different from 0 for that layer, timestep and iteration
                #cos_gamma_field[(grid1+1)*i:(grid1+1)*(i+1),:] = r_cog_vec[1]/r_cog_norm*np.sin(phi_mat)+np.sqrt(r_cog_vec[0]**2+r_cog_vec[2]**2)/r_cog_norm*np.cos(phi_mat)*np.cos(Delta_theta)
                if r_cog_norm>0:
                    cos_gamma_field_singleStep = r_cog_vec[1]/r_cog_norm*np.sin(phi_mat)+np.sqrt(r_cog_vec[0]**2+r_cog_vec[2]**2)/r_cog_norm*np.cos(phi_mat)*np.cos(Delta_theta)
                    CoG_correction = r_cog_norm*cos_gamma_field_singleStep
                else:
                    CoG_correction==0
                Data1_field[1:(grid1+2),:] = Data1_field[1:(grid1+2),:] - CoG_correction*Gacc[j]            

            # Again distinguish based on sub-iteration, if it is 0 or larger. If it is larger, we have the iterative
            # process to define the new load fields based on the sub-timestep we are at
            if i>0:
                New_Data1_field=np.hstack((theta,Data1_field))
                
                if Rampload_enabled==1:
                    Data1_field_diff = Old_Data1_field[1:len(theta)+1,1:len(lam)+1] - New_Data1_field[1:len(theta)+1,1:len(lam)+1]
                    Data1_field_diff=np.vstack((lam,Data1_field_diff))
                    Data1_field_diff=np.hstack((theta,Data1_field_diff))
                    Data1_field_diff = tuple([tuple(J) for J in Data1_field_diff])
                    
                    mdb.models[Model_name].MappedField(name='P_Diff_' +str(i)+'_'+str(j)+'_field', description='',
                       regionType=POINT, partLevelData=False, localCsys=Coor_sph,
                       pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32,
                       gridPointData={str(Radius[j]):Data1_field_diff})
                                                   
                Old_Data1_field = New_Data1_field
                Data1_field= tuple([tuple(J) for J in New_Data1_field])
                
            else:
                Data1_field=np.hstack((theta,Data1_field))
                Old_Data1_field = Data1_field
                Data1_field=tuple([tuple(J) for J in Data1_field])

            
            mdb.models[Model_name].MappedField(name='P_'+str(i)+'_'+str(j)+ '_field', description='', 
                regionType=POINT, partLevelData=False, localCsys=Coor_sph, 
                pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32, 
                gridPointData={str(Radius[j]):Data1_field})
                                               
                
            ## Create new loads
            # The loads here are created based on whether we are considering the first interface or not. If we are at
            # the 0th interface we see if we have an homogeneous or a layered Earth and then geenrate the pressure field.
            # If we are not at the 0th interface,  we check whether we are treating coarse layers or not and define the
            # pressure fields based on that.
            a = mdb.models[Model_name].rootAssembly
            
            Coor_sph=mdb.models[Model_name].rootAssembly.datums.values()[0]                ##Define load potential P20(cos(T))
            
            if j==0:
                if N_layer==1:
                    region = a.instances[Model_name].surfaces['Sout0']
                    mdb.models[Model_name].Pressure(name='Phi_'+str(i)+'_'+str(j), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='P_'+str(i)+'_'+str(j)+'_field', magnitude=Density[1]-Density[0],amplitude=UNSET)
                else:
                    if Density[j+1]-Density[j]<0:
                        region = a.instances[Model_name+'_Low'].surfaces['Sin1']
                        if i==0:
                            mdb.models[Model_name].Pressure(name='Phi_'+str(i)+'_'+str(j), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='P_'+str(i)+'_'+str(j)+'_field', magnitude=(-1)**(flag_lc+1)*(Density[0]-Density[1]),amplitude=UNSET) 
                        else:
                            mdb.models[Model_name].Pressure(name='Phi_'+str(i)+'_'+str(j), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='P_'+str(i)+'_'+str(j)+'_field', magnitude=(-1)**(flag_lc+1)*(Density[0]-Density[1]),amplitude=UNSET) 

            else:
                if Density[j+1]-Density[j]<0:
                    if j > Coarse_Layers:
                        region = a.instances[Model_name].surfaces['Sout'+str(j)]
                    else:
                        region = a.instances[Model_name+'_Low'].surfaces['Sout'+str(j)]
                        
                    if i==0:
                        
                        mdb.models[Model_name].Pressure(name='Phi_'+str(i)+'_'+str(j), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='P_'+str(i)+'_'+str(j)+'_field', magnitude=Density[j+1]-Density[j],amplitude=UNSET) 
                    else:
                        mdb.models[Model_name].Pressure(name='Phi_'+str(i)+'_'+str(j), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='P_'+str(i)+'_'+str(j)+'_field', magnitude=Density[j+1]-Density[j],amplitude=UNSET) 
            
            if i>0:
                if Density[j+1]-Density[j]<0:
                    mdb.models[Model_name].loads['Phi_'+str(i-1)+'_'+str(j)].deactivate('Step'+str(i+1))
        
        if Num == 1:
            if Compensate_cog==1:
                del mdb.models[Model_name].loads['CoG_correction_surface' + str(i)]
                del mdb.models[Model_name].loads['CoG_correction_CMB' + str(i)]
    
    print 'finished odb process\n'
    Num=Num+1
print 'Iterative process finished. Odb model generated.'