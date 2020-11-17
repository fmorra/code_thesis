# Thesis_project
This folder contains the code developed for my thesis project, both in Python (for the data processing section) and MATLAB (mainly for the data post-processing section to plot data). 
All the scripts are available on GitHub at the following links:
* [Python and FORTRAN scripts for model and results generation](https://github.com/fmorra/Thesis_project/ABAQUS_scripts)
* [Python scripts for data processing](https://github.com/fmorra/Thesis_project/pythno_files)
* [MATLAB scripts for data post-processing](https://github.com/fmorra/Thesis_project/matlab_files)


## Model and data generation

## Data processing
This section is going to describe how the data generated in ABAQUS is going to be manipulated into a series of files which are easier to read and contain all necessary information used for plotting in the next section. A sample input folder with all necessary files can be found [here](https://github.com/fmorra/Thesis_project/sample_folder). 

### Main file - python_reader_processing_only

This is the main file of the data processing step, and begins with the definition of a series of parameters and paths used to locate the files to analyze and to store the files that will be successively generated if they do not exist already. The base path follows a naming convention for which the run folder is called "run_n", the iteration, step and stress iteration directories are all defined based on it and therefore the user only needs to modify one line to move the location of all the files to be generated. 

```python
run = '18'
iteration = 1
step = 0
cycle = 1
if not isinstance(run, str):
    run = str(run)
stress_rpt_name = 'abaqus_run_' + run + '.rpt'
deflection_rpt_name = 'abaqus_run_' + run + '_deflections.rpt'
# e_csv_name = 'e.csv'
if run == 'point':
    dat_name = 'Earth_point.dat'
else:
    dat_name = 'Earth.dat'
e_dat_name = 'e.dat'
# Detect the OS which is being used and select base, rpt and dat paths accordingly. Only if/else because we either use
# windows or linux.
opersys = platform.system()
if opersys == 'Windows':
    base_path = os.path.join('C:\Users\\fabri\\Desktop\\TU Delft\\Thesis\\ABAQUS\\test_run_python\\run_' + str(run))
else:
    base_path = os.path.join('/home/fabri/Earth_model_abaqus_SLE0/results_run_' + str(run))

# Define the folder path where to find the report and e.dat files. The ABAQUS output to use consists in this folder that
# has to be manually downloaded using mobaXterm and the Earth.dat file, and it is associaated not only to the run number
# but also to an iteration, a step and a cycle number.

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

# Define the paths for the configuration .dat file and report files
dat_path = os.path.join(base_path, dat_name)
stress_report_path = os.path.join(level_dir, stress_rpt_name)
deflection_report_path = os.path.join(level_dir, deflection_rpt_name)

# Define all the other encessary directories and, if they do not exist yet, create them
stress_matrices_path = os.path.join(level_dir, 'stress_matrices')
stress_processing_path = os.path.join(level_dir, 'stress_processing_results')
deflection_matrices_paths = os.path.join(level_dir, 'deflection_matrices')
deflection_processing_path = os.path.join(level_dir, 'deflection_processing_results')
common_files_path = os.path.join(level_dir, 'common_files')
coupled_stress_folder = os.path.join(stress_matrices_path, 'coupled_stress_matrices')
base_e_path = os.path.join(level_dir, e_dat_name)

if not os.path.exists(stress_matrices_path):
    os.mkdir(stress_matrices_path)
if not os.path.exists(stress_processing_path):
    os.mkdir(stress_processing_path)
if not os.path.exists(deflection_matrices_paths):
    os.mkdir(deflection_matrices_paths)
if not os.path.exists(deflection_processing_path):
    os.mkdir(deflection_processing_path)
```

The next part is a function call to coordinate_reader, where every e
```python
start_time_coordinate_reader = time.time()
individual_path, complete_files_path, geographical_complete_files_path, large_node_matrix_path = \
    coordinate_reader(sd_input, dat_path, stress_matrices_path, stress_processing_path, deflection_matrices_paths,
                      deflection_processing_path, headers_on, dat_name, common_files_path, coupled_stress_folder)
end_time_coordinate_reader = time.time() - start_time_coordinate_reader
```
## Data post-processing (plot generation)
