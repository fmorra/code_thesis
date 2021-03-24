# Thesis project scripts overview
This folder contains the code developed for my thesis project, both in Python (for the data processing section) and MATLAB (mainly for the data post-processing section to plot data). The MATLAB code is now also being ported to Python, with the end goal of the entire workflow being in Python 3. 
All the scripts are available on GitHub at the following links:
* [Python and FORTRAN scripts for model and results generation](https://github.com/fmorra/code_thesis/tree/main/ABAQUS_scripts)
* [Python scripts for data processing](https://github.com/fmorra/code_thesis/tree/main/python_files)
* [MATLAB scripts for data post-processing](https://github.com/fmorra/code_thesis/tree/main/matlab_files)
* [Python 3 workflow](https://github.com/fmorra/code_thesis/tree/main/python_files_python_3)


## Model and data generation
The scripts belonging to this section are used to generate the ABAQUS model and run simulations on them. The main modifications have taken place in the Iter_ult script and the user_2_mises.f file, which is where the iterative procedure for stress coupling and the input instructions for the viscosity calculated with these coupled stresses take place. A more detailed README file for data generation can be found  [here](https://github.com/fmorra/python_files/README_ABAQUS.md).

## Data processing
This section is going to describe how the data generated in ABAQUS is going to be manipulated into a series of files which are easier to read and contain all necessary information used for plotting in the next section. A sample input folder with all necessary files can be found [in this separate GitHub  repository](https://github.com/fmorra/sample_run_folder), and a more detailed README file for data processing can be found  [here](https://github.com/fmorra/python_files/README_python.md).

## Data post-processing (plot generation)
This section is going to describe how the data generated in ABAQUS is going to be manipulated into a series of files which are easier to read and contain all necessary information used for plotting in the next section. As for the other process steps, a more detailed README file for data processing can be found  [here](https://github.com/fmorra/python_files/README_MATLAB.md).

## Python 3 workflow with data processing and post-processing
This section is going to document the entire workflow after generation of the results and is still under construction as the code itself is still being written. The end result for the processing section will be very similar to the Python 2 processing one, while the MATLAB section will keep its logical structure while of course everything else will change. The more detailed documentation can be found  [here](https://github.com/fmorra/python_files/README_python_3.md).