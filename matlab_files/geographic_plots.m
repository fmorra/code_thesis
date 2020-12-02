% function [] = geographic_plots(sd_input,files_to_classify,report_path,classified_path,components_to_plot)

% This function differentiates between the two cases concerning deflections
% and stresses t plot the data that has been generated in the previous
% functions. In the case of stresses, the function reads the centroid
% coordinates for each element of and the stress values contained in each
% of those elements, while in the case of deflections it reads the node
% coordinates and relative deflections and afterwards plots them.

% First of all, check if we actually want to plot the deflection or stress
% maps.

% New part added here, because it is not possible for now to generate plots
% with Python. This script wil be used for geographic and viscosity plots.
% This first part of the code is used to read all the relevant variables
% used to access the data processed with Python.
close all; clear all; clc; 
%set(0,'DefaultFigureVisible','off')
%% Stress and deflection
% 'Mises S11 S12'
% 'Magnitude U1 U2 U3'
python_base_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python';
run = '21';
iteration = '1';
step = '0';
cycle = '1';
coordinate_system = 'geographical'; % cartesian or geographical
quantity = 'deflections'; %stresses or deflections
% python_variables_path = [python_variables_base_path '\run_' run '\iteration_' iteration];
python_variables_path = [python_base_path '\run_' run '\Iteration_' iteration...
    '\step_' step '\cycle_' cycle '_reports'];
directory = dir(python_variables_path);
directory_filenames = {directory.name};
for i=1:length(directory_filenames)
    if ~contains(directory_filenames{i},'matlab')==1
        directory_filenames{i} = {};
    end
end
directory_filenames = directory_filenames(~cellfun('isempty',directory_filenames));
for i=1:length(directory_filenames)
    if contains(directory_filenames{i},coordinate_system)==1 && ...
            contains(directory_filenames{i},quantity)==1
        file_to_open = directory_filenames{i};
    end
end

file_to_open_path = [python_variables_path '\' file_to_open];
opened_file = fopen(file_to_open_path);

i = 1;
j = 1;
files_to_classify = {};
plotting_params = zeros(1,6);
while ~feof(opened_file)
    % Read a line without including the newline
    s=fgetl(opened_file);
    if contains(s, '\') ~= 1 && contains(s, '.csv') == 1
        files_to_classify{i} = strtrim(s);
    end
%     if isnumeric(str2double(s)) == 1 && ~isnan(str2double(s))
%         sd_input = str2double(s);
%     end
    if strcmp(s, 'stresses') == 1 || strcmp(s, 'deflections') == 1
        if strcmp(s, 'stresses') == 1
            sd_input = 0;
        else
            sd_input = 1;
        end
    end
    if isempty(extractAfter(s, '_reports')) && contains(s, 'Iteration_')
        report_path  = s;
    end
    if contains(s, 'Depth') == 1
        classified_path = s;
    end
    if strcmp(s, 'cartesian')== 1 || strcmp(s, 'geographical')== 1
        components_to_plot = s;
    end
    % if contains(s, 'Complete') == 1 || contains(s, 'geographical_complete')
    if contains(s, 'omplete') == 1 
        complete_matrices_path = s;
    end
    if abs(str2double(s)) > 1 && ~isnan(str2double(s))
        plotting_params(j) = str2double(s);
        j = j + 1;
    end
    i = i + 1;
end

if sd_input == 0
    matrix_to_read_part = 'Stresses_layer_';
    figure_folder = 'stress_plots';
else
    matrix_to_read_part = 'Deflections_layer_';
    figure_folder = 'deflection_plots';
end
figures_path = [report_path '\' figure_folder];
if ~exist(figures_path, 'dir')
    mkdir(figures_path)
end

while_check = 1;
while while_check == 1
    plots_bool = input(['Do you want to plot the stresses or deformations? \n'...
        'Enter 1 for plots, 0 to skip this step: \n']);
    if plots_bool == 1 || plots_bool == 0
        while_check = 0;
    else
        disp('Incorrect input, enter 1 for plots and 0 to skip this step.');
    end
end

diff_matrix_path = [python_base_path '\run_' run '\difference_matrices_plots\matrices'];
if ~exist(diff_matrix_path, 'dir')
    mkdir(diff_matrix_path)
end
min_lat = plotting_params(1);
max_lat = plotting_params(2);
min_lon = plotting_params(3);
max_lon = plotting_params(4);
depths_to_plot = [plotting_params(5) plotting_params(6)]; 
if strcmp(quantity, 'deflections') == 1
    min_lat = -90;
    max_lat = -65;
    min_lon = -180;
    max_lon = 180;
    depths_to_plot = [100 150];
end
if depths_to_plot(1) > depths_to_plot (2)
    depths_to_plot([1 2]) = depths_to_plot([2 1]);
end
if plots_bool == 1
    depth_searcher(sd_input,coordinate_system,complete_matrices_path,...
        files_to_classify,figures_path,diff_matrix_path,iteration,step,cycle,...
        min_lat,max_lat,min_lon,max_lon,depths_to_plot);
end

%% Viscosity 
    
visco_check = 1;
while visco_check == 1
    visco_plots_bool = input(['Do you want to plot the viscosity/strain rate? \n'...
        'Enter 1 for plots, 0 to skip them: \n']);
    if visco_plots_bool == 1 || visco_plots_bool == 0
        visco_check = 0;
    else
        disp('Incorrect input, enter 1 for the viscosity or strain rate plots or 0 to skip this step:\n');
    end
end

% Plots - viscosity
if visco_plots_bool == 1
    viscosity_figures_path = [report_path '\viscosity_plots'];
    if ~exist(viscosity_figures_path, 'dir')
        mkdir(viscosity_figures_path)
    end
    visco_diff_path = visco_plotter(python_variables_path,coordinate_system, complete_matrices_path,...
        report_path,min_lat,max_lat,min_lon,max_lon,depths_to_plot,iteration,step,cycle,...
        diff_matrix_path);
end

%% Time difference plots 
diff_check = 1;
while diff_check == 1
    diff_input = input('Enter 1 to generate plots of differences between different cycles and steps, 0 to skip this:\n');
    if diff_input == 1 || diff_input == 0
        diff_check = 0;
    else
        disp('Incorrect input, enter 1 to generate difference plots or 0 to skip:\n');
    end
end

if diff_input == 1
    diff_calculator(diff_matrix_path,min_lat,max_lat,min_lon,max_lon);
end
