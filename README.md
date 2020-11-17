# Thesis_project
This folder contains the code developed for my thesis project, both in Python (for the data processing section) and MATLAB (mainly for the data post-processing section to plot data). 
All the scripts are available on GitHub at the following links:
* [Python and FORTRAN scripts for model and results generation](https://github.com/fmorra/code_thesis/tree/main/ABAQUS_scripts)
* [Python scripts for data processing](https://github.com/fmorra/code_thesis/tree/main/python_files)
* [MATLAB scripts for data post-processing](https://github.com/fmorra/code_thesis/tree/main/matlab_files)


## Model and data generation
The scripts belonging to this section are used to generate the ABAQUS model and run simulations on them. The main modifications have taken place in the Iter_ult script and the user_2_mises.f file, which is where the iterative procedure for stress coupling and the input instructions for the viscosity calculated with these coupled stresses take place. 

## Data processing
This section is going to describe how the data generated in ABAQUS is going to be manipulated into a series of files which are easier to read and contain all necessary information used for plotting in the next section. A sample input folder with all necessary files can be found [in this GitHub  repo](https://github.com/fmorra/sample_run_folder), and a more detailed README file for that case can be found  [here](https://github.com/fmorra/python_files/README_python.md) with its corresponding .html [here](https://github.com/fmorra/sample_run_folder).

## Data post-processing (plot generation)
