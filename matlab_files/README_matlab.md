# Data post-processing (plotting)
Plotting data has been done in MATLAB for thesis purposes, however this is being transferred to Python together with overhauling the processing scripts to Python 3. Two ways of plotting are possible: the default one is related to one very specific stress iteration, while the other one is giving run, iteration, step and stress iteration ranges to produce stress/deformation/viscosity plots in batches. This is mainly useful when a single plot change has to be implemented in all plots. The single stress iteration case will be discussed here as the functions used are all exactly the same.  
## Main file-geographic_plots
The last segment of the data workflow concerns plotting and for this treason the main function is divided into two main parts for both single and multiple runs: the first one defines all the parameters to use and the second one calls the plotting functions.

### depth_searcher
This is the first function to be called and is related to stress and deformation plotting. The function generates a grid of points, then reads the stresses and deformations and defines the boundaries for the plots to be generated. The data is then interpolated and plotted in this range and also saved in a matrix which can be used to generate plots of the value differences. The plots are generated using a scale defined by a function used to give the same scale to all the plots to be genera ted.  Function inputs:
- run, run number as a string;
- sd_input, Boolean for stresses or deformations;
- coordinate_system, string for the coordinate system to use;
- complete_matrices_path, path where to find the stress and deformation values with their coordinates;
- figures_path, path where to save the plots;
- diff_matrix_path, path where to save the value matrices for the plots of values differences;
- iteration, iteration number as a string;
- step, step number as a string;
- cycle, stress iteration number as a string;
- min_lat, minimum latitude;
- max_lat, maximum latitude;
- min_lon, minimum longitude;
- max_lon, maximum longitude;
- depths_to_plot, depth range extremes;
- run_folder, folder where to find data for the simulation run defined above;
- figure_counter, counter for the number of figures generated;
- python_base_path, path to the folder containing run folders;
- run_vec, vector of runs to used to give the same extremes to all the plots to generate using all the run folders in the base path;
- simul_time, the timestep string to use in the plot title.
Function outputs:
- figure_counter, index for the figure to generate used to generate new figures from the next functions if necessary without overwriting the old ones.

#### caxisextremes

This function is used to give the same colorbar scale to all plots. It reads all the possible files for the selected quantity and reference system and then reads the extremes for each of them, and after that calculates the minimum and maximum over all the matrices. This is done for all the selected components to plot. This function also works with viscosity, strain rate and B coefficients. Its inputs are:
- sd_input;
- min_lat;
- max_lat;
- min_lon;
- max_lon,;
- depths_to_plot;
- selected_components;
- viscosity_input, Boolean to decide whether to plot the viscosity or not;
- b_input, Boolean to decide whether to plot the B coefficients or not;
- python_base_path;
- run_vec;
- coordinate_system;
Its outputs are:
- new_colorbarlimits, the colorbar limits based on the extremes of all the other simulation runs considered.

### visco_plotter 
This function is similar to the one used to plot stresses and deformations in structure and functionality, with the main difference of generating the viscosity or strain rate plots. The order in which these functionalities are implemented is the same as the one for stresses and deformations, so they do not need to be discussed again. Its inputs are:
- sd_input; 
- python_variables_base_path, path for the current iteration;
-coordinate_system;
- complete_matrices_path; 
- report_path; 
- min_lat;
- max_lat;
- min_lon;
- max_lon;
- depths_to_plot;
- iteration;
- step;
- cycle;
- diff_matrix_path;
- run_folder;
- run; 
- figure_counter;
- python_base_path; 
- run_vec;
- simul_time.

Its outputs are:
- figure_counter;
- visco_diff_path, path where to save the matrices to calculate difference plots for viscosity or strain rate.

#### B_plots
This function is called inside visco_plotter . Its inputs are:
- viscosity_figures_path, path to the viscosity plots;
- alin, the B_diff coefficient in the creep law;
- a, the B_disl coefficient in the creep law;
- data_points_indices, indices of the nodes or centroids in the latitude, longitude and depth ranges;
- part, name of the part being processed;
- min_depth;
- max_depth; 
- min_lat;
- max_lat; 
- min_lon; 
- max_lon; 
- lat, latitude of the points belonging to the selected plotting ranges;
- lon, longitude of the points belonging to the selected plotting ranges;
- run;
- run_folder;
- viscosity_input;
- python_base_path;
- depth, depth of the points belonging to the selected plotting ranges;
- run_vec;
- coordinate_system; 
- visco_diff_path.
    
The function has no outputs aside from the plots to generate.

### diff_calculator
This function plots differences between the common quantities existing for two different stress iterations in the latitude, longitude and depth ranges selected. A similar function to this one is used to plot differences in the B coefficients, and it is separate because they are defined for every run in this case (actually only for even and odd-numbered runs) so their file location and handling is less detailed than what is required, for example, for stresses and deformations. Its inputs are:
- diff_matrix_path_incomplete, the base path for all the difference matrices, for all quantities and reference systems;
- min_lat;
- max_lat;
- min_lon;
- max_lon;
- run;
- depths_to_plot;
- iteration;
- step;
- figure_counter.

Its outputs are:
- figure_counter.

### b_diff_calc
This function also calculates value differences and plots them between two stress iterations, but only does so for B coefficients. Its inputs are:
- b_diff_plots_folder, path where to save the B coefficients difference plots;
- python_base_path;
- min_lat;
- max_lat;
- min_lon;
- max_lon;
- min_depth;
- max_depth;
- figure_counter.
    
The function has no outputs aside from the B coefficients difference plots.

