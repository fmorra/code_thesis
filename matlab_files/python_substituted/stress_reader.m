% Script used to read the stress vectors and plot their values over a world
% map based on how ABAQUS associates each element to a coordinate on the
% sphere. 
% Author: Fabrizio Morra
% Date: 18/2/2018
% clear all; close all; clc;

% Run numbering begins with 5 because this is the first run for which there
% has been a successful stress output. Generate all the necessary paths
% which will be used in the next calculations. 

close all; clc;
total_time_tic = tic;
report_processing_tic = tic;
run = 7;
e_name = 'e.csv';
if isstring(run) == 0
    run = num2str(run);
end
stress_report_name = ['abaqus_run_' run '.rpt'];
deflection_report_name = ['abaqus_run_' run '_deflections.rpt'];
if (strcmp(run,'point') == 1)
    coord_file = 'Earth_point.dat';
else
    coord_file = 'Earth.dat';
end
iteration = 1;
report_path = ['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\stress_output_files\run_' run];
% For testing, put the dat and rpt files in a folder and copy its path here
% iteration_path = ['C:\Users\fabri\Desktop\TU Delft\algorithm_trial_folder'];
stress_file = [report_path '\' stress_report_name];
deflection_file = [report_path '\' deflection_report_name];
file_coord = [report_path '\' coord_file];
e_path = [report_path '\' e_name];
iteration_path = [report_path '\Iteration_' num2str(iteration)];
stress_matrices_path = [iteration_path '\stress_matrices'];
stress_processing_path = [iteration_path '\stress_data_processing'];
deflection_matrices_path = [iteration_path '\deflection_matrices'];
deflection_processing_path = [iteration_path '\deflection_data_processing'];
common_files_path = [iteration_path '\common_values'];
coupled_stress_path = [stress_matrices_path '\coupled_stress_matrices'];
if ~exist(stress_matrices_path, 'dir')
    mkdir(stress_matrices_path)
end
if ~exist(stress_processing_path, 'dir')
    mkdir(stress_processing_path)
end
if ~exist(deflection_matrices_path, 'dir')
    mkdir(deflection_matrices_path)
end
if ~exist(deflection_processing_path, 'dir')
    mkdir(deflection_processing_path)
end
if ~exist(iteration_path, 'dir')
    mkdir(iteration_path)
end
% Decide if we just need to read the file or read parts of it based on a
% word key. Note: using the headers for the first time will take longer for
% the program to run, around 5 minutes.

headers_on = 1;
search_key_stress = 'reported at element';
search_key_deflections = 'reported at nodes';

% Read the report files for both stresses and deflections and store the
% values for each part, and, if applicable, region, in a CSV

headers = read_stress_report(stress_file,search_key_stress,stress_matrices_path,headers_on);
deflection_paths = read_deflections_report(deflection_file,search_key_deflections, ...
    deflection_matrices_path, headers_on);
report_processing_toc = toc(report_processing_tic);
disp(['The time spent to process the report files is: ' ...
    num2str(report_processing_toc) ' seconds.']);

% Add additional stresses (test)

if iteration == 1
    % First iteration: simply calculate the stresses and classify them
    copyfile(e_path, [iteration_path '\' e_name]);
else
    % Second iteration onwards: add a new stress and recalculate the Mises
    % stress
    new_stress_test = 1000000; %Pa - 1MPa
    copyfile [report_path '\Iteration_' num2str(iteration - 1)] ...
        [iteration_path '\e.csv'];
    if ~exist(coupled_stress_path, 'dir')
        mkdir(coupled_stress_path)
    end
    large_new_matrix = stress_coupling(stress_matrices_path, coupled_stress_path, headers_on, ...
        new_stress_test, headers);
    stress_matrices_path = coupled_stress_path; % To work with coupled matrices
end

% Decide whether we want to work on the stresses or the deflections.
check_0 = 1;
while check_0 == 1
    processing_input = input('Do you wish to process data for plotting? Enter yes or no: \n');
    if strcmp(processing_input, 'yes') == 1 || strcmp(processing_input, 'no') == 1
        check_0 = 0;
    else
        disp('Incorrect input, select either yes or no. \n')
    end
end

if strcmp(processing_input, 'no') == 1
    disp('No files have been generated. The Mises stress has been modified and saved. Going to the next iteration.')
else
    check_1 = 1;
    % sd_input = 1;
    while check_1 == 1
        sd_input = input(['Enter 0 to discretize stress components based on depth,'...
            ' 1 for the deflections: \n']);
        if sd_input == 0 || sd_input == 1
            check_1 = 0;
        else
            disp('incorrect input, select either 1 or 0. \n');
        end
    end
    
    % Divide the nodes and elements in csv files for each part and, if
    % applicable, region, after modifying the dat file to perform the
    % operation. After this, associate all the stresses and deflections
    % components to the relative centroids or nodes.
    
    [individual_path,complete_files_path,spherical_complete_files_path,figure_counter,...
        large_node_matrix_path] = read_coordinates(sd_input,file_coord,stress_processing_path,...
        deflection_paths,headers_on,stress_matrices_path,deflection_processing_path,coord_file,common_files_path);
    
    % Classify the results for each part and region based on depth as radial
    % distance from the center, we pass the folders with both unspherical and
    % unspherical stress or deflection components. Unspherical means cartesian,
    % spherical spherical.
    
    check_2 = 1;
    % components_to_plot = 'spherical';
    while check_2 == 1
        components_to_plot = input(['Enter "cartesian" to work with cartesian components, ' ...
            '"spherical" to work with spherical ones: \n']);
        if strcmp(components_to_plot, 'cartesian') ==1 || strcmp(components_to_plot, 'spherical') == 1
            check_2 = 0;
        else
            disp('Incorrect input, select either "cartesian" or "spherical".');
        end
    end
    
    % Here we use a function to classify the different stress values based on
    % their depth, because we can see that the centroids are calculated at
    % different "layers" and each centroid value has a stress. Clssifying and
    % saving them is better for plotting purposes, because it means we do not
    % have to search for the points at the same depths every time we want to
    % produce plots. Ths function also accepts the deflection files, based on
    % the values of sd_input.
    
    depth_classification_tic = tic;
    [classified_path,files_to_classify,maximum_depth] = depth_classifier(sd_input,...
        deflection_processing_path,individual_path,components_to_plot,complete_files_path,...
        spherical_complete_files_path,headers_on);
    depth_classification_toc = toc(depth_classification_tic);
    disp(['The time spent to discretize the files based on depth is: ' ...
        num2str(depth_classification_toc) ' seconds.']);
    
    %--------------------------------------------------------------------------
    % Test section to understand what causes the the shift in lon
    % stress_scatter_test = readmatrix(['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\stress_output_files\run_6\stress_matrices_composition\Individual_Part_Values\Depth_classified_files\subclassification_depth_410_km\Stresses_layer_410_km_bin_1.csv']);
    % figure(2)
    % earth_example
    % scatter3(stress_scatter_test(:,9),stress_scatter_test(:,10),stress_scatter_test(:,11),2,stress_scatter_test(:,2));
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    %
    % matrix_to_read = readmatrix('C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\stress_output_files\run_6\stress_matrices_composition\Individual_Part_Values\Depth_classified_files\subclassification_depth_410_km\Stresses_layer_410_km_bin_1.csv');
    % variable_to_plot =  matrix_to_read(:,2);
    % lon_to_plot = matrix_to_read(:,end); %longitude [deg]
    % lat_to_plot = matrix_to_read(:,end-1); %latitude [deg]
    % % Create figure
    % figure()
    % % Load world map with coastlines
    % %worldmap world;
    % worldmap ([-90 -70],[-180 -45])
    % load coastlines;
    % plotm(coastlat, coastlon);
    % % Grid the data before plotting, and generate a surface
    % [LatGrid, LonGrid] = meshgrid(linspace(min(lat_to_plot), max(lat_to_plot),1000), ...
    %     linspace(min(lon_to_plot), max(lon_to_plot),1000));
    % stressgrid = griddata(lat_to_plot, lon_to_plot, variable_to_plot, LatGrid, LonGrid);
    % surfm(LatGrid, LonGrid, stressgrid);
    % alpha 0.5;
    
    % matrix_to_read = readmatrix('C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\stress_output_files\run_6\stress_matrices_composition\Individual_Part_Values\Depth_classified_files\subclassification_depth_410_km\Stresses_layer_410_km_bin_1.csv');
    % variable_to_plot =  matrix_to_read(:,3);
    % lon_to_plot = matrix_to_read(:,end); %longitude [deg]
    % lat_to_plot = matrix_to_read(:,end-1); %latitude [deg]
    % figure()
    % % worldmap world;
    % worldmap ([-90 -70],[-180 -45])
    % load coastlines;
    % plotm(coastlat, coastlon);
    % % Grid the data before plotting, and generate a surface
    % [LatGrid, LonGrid] = meshgrid(linspace(min(lat_to_plot), max(lat_to_plot),1000), ...
    %     linspace(min(lon_to_plot), max(lon_to_plot),1000));
    % stressgrid = griddata(lat_to_plot, lon_to_plot, variable_to_plot, LatGrid, LonGrid);
    % surfm(LatGrid, LonGrid, stressgrid);
    % % Figure propeties definition
    % if sd_input == 0
    %     title('Map of the stresses');
    % else
    %     title(['Map of the deflections']);
    % end
    % alpha 0.5;
    % c = colorbar;
    % if sd_input == 0
    %     c.Label.String = 'Stress [Pa]';
    % else
    %     c.Label.String = 'Deflections  [m]';
    % end
    % c.Location = 'southoutside';
    % grid on;
    
    %--------------------------------------------------------------------------
    % Plot the selected results in 2-D lat/lon graphs
    
    geographic_plots_tic = tic;
%     geographic_plots(sd_input,files_to_classify,iteration_path,classified_path,components_to_plot);
    geographic_plots_toc = toc(geographic_plots_tic);
    disp(['The time spent to plot the results is: ' ...
        num2str(geographic_plots_toc) ' seconds.']);
    total_time_toc = toc(total_time_tic);
    disp(['The total elapsed time is: ' ...
        num2str(total_time_toc) ' seconds.']);
end

