# Data post-processing (plotting)
## Main file-geographic_plots
The last segment of the data workflow concerns it plotting. the iteration, step, stress iteration, quantity and its reference system of the data to plot are defined and then the corresponding handle .txt file produced by Python is searched for. 
```MATLAB
python_base_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python';
run = '17';
iteration = '1';
step = '0';
cycle = '1';
coordinate_system = 'geographical'; % cartesian or geographical
quantity = 'stresses'; %stresses or deflections
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
```

This file is opened and all of its variables are extracted and stored into the corresponding variables which are going to be used next for the plot generation procedure. 

```MATLAB
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
    if contains(s, 'Complete') == 1 || contains(s, 'geographical_complete')
        complete_matrices_path = s;
    end
    if abs(str2double(s)) > 1 && ~isnan(str2double(s))
        plotting_params(j) = str2double(s);
        j = j + 1;
    end
    i = i + 1;
end
```
A few more variables used for the plot names and to save the figures are also defined, then the script asks for the user whether these plots need to be generated or not, reads the latitude, longitude and depth ranges extremes and finally the first function for the stresses or deformation plot generation is called.
```MATLAB
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
if depths_to_plot(1) > depths_to_plot (2)
    depths_to_plot([1 2]) = depths_to_plot([2 1]);
end
if plots_bool == 1
    depth_searcher(sd_input,coordinate_system,complete_matrices_path,...
        files_to_classify,figures_path,diff_matrix_path,iteration,step,cycle,...
        min_lat,max_lat,min_lon,max_lon,depths_to_plot);
end
```
---
### depth_searcher
This function starts by determining the quantity to plot and by initializing a path where to save the matrix of data for the current stress iteration, which can be used to generate difference plots. The, the function selects the stresses or deformations data columns to read based on the user's input and then starts iterating over the parts to plot, opening the corresponding matrix and then filtering its data based on the latitude, longitude and depth ranges.
```MATLAB
if sd_input == 0
    components_to_plot = 'stresses';
else
    components_to_plot = 'deflections';
end
complete_diff_path = [diff_matrix_path '\' components_to_plot  '\' coordinate_system];
if ~exist(complete_diff_path, 'dir')
    mkdir(complete_diff_path)
end

parts_to_plot = {'EARTH'};
if sd_input == 0
    selected_component = input(['Enter the desired stress component(s) to plot,' ...
        ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
    selected_components = split(selected_component);
    selected_columns = zeros(length(selected_components),1);
    for k = 1:length(selected_components)
        if strcmp(selected_components{k}, 'Mises') == 1
            selected_columns(k) = 2;
        elseif strcmp(selected_components{k}, 'S11') == 1
            selected_columns(k) = 3;
        elseif strcmp(selected_components{k}, 'S22') == 1
            selected_columns(k) = 4;
        elseif strcmp(selected_components{k}, 'S33') == 1
            selected_columns(k) = 5;
        elseif strcmp(selected_components{k}, 'S12') == 1
            selected_columns(k) = 6;
        elseif strcmp(selected_components{k}, 'S13') == 1
            selected_columns(k) = 7;
        else
            selected_columns(k) = 8;
        end
    end
else
    selected_component = input(['Enter the desired deflection components(s) to plot,' ...
        ' possible values are Magnitude, U1, U2, U3:\n']);
    selected_components = split(selected_component);
    selected_columns = zeros(length(selected_components),1);
    for k = 1:length(selected_components)
        if strcmp(selected_components{k}, 'Magnitude') == 1
            selected_columns(k) = 2;
        elseif strcmp(selected_components{k}, 'U1') == 1
            selected_columns(k) = 3;
        elseif strcmp(selected_components{k}, 'U2') == 1
            selected_columns(k) = 4;
        else
            selected_columns(k) = 5;
        end
    end
end
for i=1:length(parts_to_plot)
    min_depth = depths_to_plot(1);
    max_depth = depths_to_plot(2);

    if strcmp(coordinate_system, 'cartesian') == 1
        matrix_to_read = readmatrix([complete_matrices_path '\Complete_file_'...
            parts_to_plot{i} '.csv']);
    else
        matrix_to_read = readmatrix([complete_matrices_path '\geographical_complete_file_'...
            parts_to_plot{i} '.csv']);
    end
    depth = matrix_to_read(:,end-2)/1e3;
    lat = matrix_to_read(:,end-1);
    lon = matrix_to_read(:,end);
    depth_condtion = depth>min_depth & depth<max_depth;
    lat_condition = lat>min_lat & lat<max_lat;
    lon_condition = lon>min_lon & lon<max_lon;
    data_points_indices = matrix_to_read(depth_condtion...
        & lat_condition & lon_condition); % Gives the indices
```
The length of indices vector of the filtered data is used to allocate the matrix called matrix_for_difference which could potentially be used for the difference plots generation, with a number of columns equal to the variables to plot plus two for latitude and longitude.
```MATLAB
    matrix_for_difference = zeros(length(data_points_indices),length(selected_columns)+2);
    matrix_for_difference_headers = [selected_components','Latitude','Longitude'];
```
A for cycle then iterates over the data columns to plot each selected component, matrix_for_difference is filled with each data column and latitude and longitude for that depth range, then data is gridded using the geoloc2grid function and then saved using a variety of additional options to have a more detailed plot title and filename to save. 
```matlab
for j = 1:length(selected_columns)
        plot_variable = matrix_to_read(data_points_indices,selected_columns(j));
        matrix_for_difference(:,j) = plot_variable;
        matrix_for_difference(:,end-1) = lat(data_points_indices);
        matrix_for_difference(:,end) = lon(data_points_indices);
        [Z, refvec] = geoloc2grid(lat(data_points_indices),wrapTo360(lon(data_points_indices)),...
            plot_variable,0.5);
        load coastlines;
        cmap = colormap('jet');
        alpha 0.7;
        colormap(cmap);
        caxis('auto'); % was [0 1e6]
        h = colorbar('v'); %set(h, 'ylim', [0 1e6]);
        if sd_input == 0
            set(get(h,'ylabel'),'string',[selected_components{j} ' (Pa)'])
        else
            set(get(h,'ylabel'),'string',[selected_components{j} ' (m)'])
        end
        latlim = [min_lat max_lat];lonlim = [min_lon max_lon];
        ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
        set(ax,'Visible','off');
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        geoshow(Z, refvec, 'DisplayType', 'texture');
        plotm(coastlat,coastlon);
        if sd_input == 0
            title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j}]});
        else 
            title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j}]});
        end
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        grid on;
        % Save the figures based on whether we are working with
        % stresses or deflections
        saveas(gcf,[figures_path '\' coordinate_system '_' components_to_plot '_' parts_to_plot{i}...
            '_depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km_Component_'...
            selected_components{j} '.png']);
    end
```
Once this iteration is finished, matrix_for difference is also saved since it is now completely filled with all the variables to plot.
```matlab
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[complete_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
```
---
The next step is asking the user whether they want to generate viscosity plots and, if so, the function visco_plotter is called.
```MATLAB
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
```

---
### visco_plotter 
The function starts by deciding what quantity to plot and and defining and creating the paths where to save the .csv fiels for difference plots and then the viscosity or strain rate plots. 
```matlab
parts_to_plot = {'EARTH'};
visco_strain_input = input('Enter 0 to plot the viscosity, 1 to plot the strain rate:\n');
if visco_strain_input == 0
    quantity = 'viscosity';
else
    quantity = 'strain_rate';
end
visco_diff_path = [diff_matrix_path '\' quantity];
if ~exist(visco_diff_path, 'dir')
    mkdir(visco_diff_path)
end
e_file_name = 'e.dat';
viscosity_figures_path = [report_path '\viscosity_plots'];
if ~exist(viscosity_figures_path, 'dir')
    mkdir (viscosity_figures_path)
end
```
A for cycle starts over the parts to plot (in this case EARTH) and for each of them the e.dat file for the current stress iteration is opened followed by the complete matrix containin stresses, coordinates, depth, latitude and longitude. 
```matlab
for i=1:length(parts_to_plot)
    % for i = 1:length(iteration_subfolders)
    iter_path = python_variables_base_path;
    e_path = [iter_path '\' e_file_name];
    if strcmp(coordinate_system, 'cartesian') == 1
        matrix_to_open_path = [complete_matrices_path '\Complete_file_'...
            parts_to_plot{i} '.csv'];
    else
        matrix_to_open_path = [complete_matrices_path '\geographical_complete_file_'...
            parts_to_plot{i} '.csv'];
    end
```
The viscosity and strain rates are then calculated for all the matrix elements and then saved in a separate matrix called viscosity_matrix which is concatenated to the stress matrix. 
```matlab
    e = readmatrix(e_path);
    elements = e(:,1);
    alin = e(:,2);
    a = e(:,3);
    an = 3.5;
    strain_rate = zeros(length(e),1);
    opened_visco_matrix = readmatrix(matrix_to_open_path);
    mises = opened_visco_matrix(:,2);
    pow = (an);
    for j = 1:length(strain_rate)
        strain_rate(j,1) = alin(j)*mises(j) + a(j)*mises(j)^pow;
    end
    viscosity = mises./(3*strain_rate);
    viscosity_matrix = zeros(length(opened_visco_matrix),4);
    viscosity_matrix(:,1) = elements;
    viscosity_matrix(:,end - 2) = viscosity;
    viscosity_matrix(:,end - 1) = strain_rate;
    viscosity_matrix(:,end) = mises;
    complete_matrix_with_viscosity = [opened_visco_matrix, viscosity_matrix];
```
Once again, This matrix is filtered to only leave the points contained in the latitude, longitude and depth range passed from Python and the remaining data is plotted and saved after being gridded with geoloc2grid. The matrix used to calculate stress iteration differences for this part is also saved and the function terminates. This part is very similar to the generation of stress and deformation plots and that is why the code is not reported here.
---
Finally,the user is asked whether they want to generate plots of the value differences between two stress iterations for the same depth range, and if so, the diff_calculator function is called.
```MATLAB
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
```
---
### diff_calculator
This function plots differences between the common quantities existing for two different stress iterations in the latitude, longitude and depth ranges read from Python. It starts by choosing which quantity and in what reference system (if applicable) to plot and defines and creates where necessary all the paths used to read data and store plots.
```matlab
visco_choice = 0;
if visco_choice == 0
    stress_deflection_input = 0;
    ref_system_input = 1;
    if stress_deflection_input == 0
        quantity = 'stresses';
    else
        quantity = 'deflections';
    end
    if ref_system_input == 0
        ref_system = 'cartesian';
    else
        ref_system = 'geographical';
    end
    diff_matrix_path = [diff_matrix_path_incomplete '\' quantity '\' ref_system];
    diff_plots_folder = [diff_matrix_path_incomplete '\plots\' quantity '\' ref_system];
else
    visco_strain_input = 0;
    if visco_strain_input == 0
        quantity = 'viscosity';
    else
        quantity = 'strain_rate';
    end
    diff_matrix_path = [diff_matrix_path_incomplete '\' quantity];
    diff_plots_folder = [diff_matrix_path_incomplete '\plots\' quantity];
end
dir_diff_files = dir([diff_matrix_path '\*.csv']);

if ~exist(diff_plots_folder, 'dir')
    mkdir(diff_plots_folder)
end
diff_files = {dir_diff_files.name};
```
The next part analyzes the two matrices for which the differences have to be plotted. The matrices are found by entering two different quintuplets of values, representing the iteration, step, stress iteration and minimum and maximum depth used for previously generated plots. The depth range has to be the same for both of them and is the one given by the Python scripts, but the user has to enter it. Using these values, the script searches for the matrix containing these numbers in its name in the folder containing the matrices for difference plotting, and returns an error message if nothing is found. 
```matlab
matrix_1_flag = 0;
while matrix_1_flag == 0
    quintuplet_1 = input(['Enter the iteration, step and cycle number, and the minimum and maximum depth '...
    'of the first matrix used to plot the difference with respect to another moment in time '...
    'as a vector of square brackets of five values:\n']);
    for i=1:length(diff_files)
        found_values_1 = str2double(regexp(diff_files{i}, '\d+', 'match'));
        if isequal(quintuplet_1,found_values_1)
            matrix_1_flag = 1;
            file_to_read_1 = diff_files{i};
            iteration_1 = quintuplet_1(1);
            step_1 = quintuplet_1(2);
            cycle_1 = quintuplet_1(3);
            break;
        end
    end
    if matrix_1_flag == 1
        break;
    else
        disp('No match for the selected values, the matrix does not exist. Exiting search.');
        break;
    end
end

matrix_2_flag = 0;
while matrix_2_flag == 0
    quintuplet_2 = input(['Enter the run, step and cycle number, and the minimum and maximum depth '...
    'of the second matrix  used to plot the difference with respect to another moment in time '...
    'as a vector of square brackets of five values:\n']);
    for i=1:length(diff_files)
        found_values_2 = str2double(regexp(diff_files{i}, '\d+', 'match'));
        if isequal(quintuplet_2,found_values_2)
            matrix_2_flag = 1;
            file_to_read_2 = diff_files{i};
            iteration_2 = quintuplet_2(1);
            step_2 = quintuplet_2(2);
            cycle_2 = quintuplet_2(3);
            break;
        end
    end
    if matrix_2_flag == 1
        break;
    else
        disp('No match for the selected values, the matrix does not exist. Exiting search.');
        break;
    end
end
min_depth = quintuplet_2(4);
max_depth = quintuplet_2(5);
```
If two matrices at two different stress iterations are found for the selected depth range, they are opened and the variables common for each of them are read together with the latitude or longitude of one of them (it does not matter, since they are the same).
```matlab
if matrix_1_flag == 1 && matrix_2_flag == 1
    % Extract the headers of both matrices, intersect and display the
    % possible common variables
    [values_1,headers_1,~] = xlsread([diff_matrix_path '\' file_to_read_1]);
    [values_2,headers_2,~] = xlsread([diff_matrix_path '\' file_to_read_2]);
    possible_variables_for_plotting = intersect(headers_1,headers_2, 'stable');
    possible_variables_for_plotting{end-1} = {};
    possible_variables_for_plotting{end} = {};
    possible_variables_for_plotting = possible_variables_for_plotting(~cellfun('isempty',possible_variables_for_plotting));
    disp(['The possible variables to generate the difference plots between these'...
        ' two cycles are:\n']);
    disp(possible_variables_for_plotting);
    diff_variables_to_plot = possible_variables_for_plotting;
    lat = values_1(:,end-1);
    lon = values_1(:,end);
```
Then, in a similar fashion to the stress and deformation and viscosity plot generation functions, a for cycle iterates over each common variable found, calculates the difference between these values and plots them after using geoloc2grid. This time, the data does not need to be filtered because the function works on already manipulated data. The difference plots are then saved in the corresponding folder.

---
the main script then terminates. 