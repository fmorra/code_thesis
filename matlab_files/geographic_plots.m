% This is the main script for the generation of stress, deformation,
% viscosity and their differences in time plots.
set(groot,'DefaultFigureVisible','off')

python_base_path = 'C:\Users\fabri\Desktop\tu_delft\Thesis\ABAQUS\test_run_python';
figure_counter = 1;

%% Singe run case
% Define the run and then the lit of runs. This is necessary to give the
% same colorbar extremes to all plots
run = '21';
runs_directory = dir(python_base_path);
run_folders_flags = [runs_directory.isdir];
run_folders = runs_directory(run_folders_flags);
run_vec = [];
run_counter = 0;
for i = 1:length(run_folders)
    run_folder_name = run_folders(i).name;
    if isempty(regexp(run_folder_name, '\d+', 'match')) == 0
        run_counter = run_counter + 1;
        run_vec(run_counter) = str2double(regexp(run_folder_name, '\d+', 'match'));
    end
end

% Define iteration, step and corresponding ice history time, stress
% iteration, plotting parameters such as what to plot and in what reference
% system, the lat, lon and depth ranges 
iteration = '1';
step = '0';
simul_time = '1 ka'; % '31 ka'
cycle = '2';
coordinate_system = 'geographical'; % cartesian or geographical
quantity = 'stresses'; %stresses or deflections
min_lat = -90;
max_lat = -65;
min_lon = -180;
max_lon = 180;
depths_to_plot = [145 160];
python_variables_path = [python_base_path '\run_' run '\Iteration_' iteration...
    '\step_' step '\cycle_' cycle '_reports'];
run_folder = [python_base_path '\run_' run];
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

% Read the necessary handles from the MATLAB text file produced in Python
file_to_open_path = [python_variables_path '\' file_to_open];
opened_file = fopen(file_to_open_path);
i = 1;
files_to_classify = {};
while ~feof(opened_file)
    % Read a line without including the newline
    s=fgetl(opened_file);
    % Read the parts, as of now only EARTH
    if contains(s, '\') ~= 1 && contains(s, '.csv') == 1
        files_to_classify{i} = strtrim(s);
    end
    % Read whether to plot stresses or deflections, this will agree with
    % the previous user input as that one determines the MATLAB file to
    % read
    if strcmp(s, 'stresses') == 1 || strcmp(s, 'deflections') == 1
        if strcmp(s, 'stresses') == 1
            sd_input = 0;
        else
            sd_input = 1;
        end
    end
    % Read the path where the stress iteration files are
    if isempty(extractAfter(s, '_reports')) && contains(s, 'Iteration_')
        report_path  = s;
    end
    % Read the path for the complete matrices
    if contains(s, 'omplete') == 1 
        complete_matrices_path = s;
    end
    i = i + 1;
end

% Set up the depth range, quantities to plot and figure path
if depths_to_plot(1) > depths_to_plot (2)
    depths_to_plot([1 2]) = depths_to_plot([2 1]);
end
if sd_input == 0
    figure_folder = 'stress_plots';
else
    figure_folder = 'deflection_plots';
end
figures_path = [report_path '\' figure_folder];
if ~exist(figures_path, 'dir')
    mkdir(figures_path)
end

% Choose whether to plot stresses and deflection and, if so, set up the
% path for the matrices used to store data on which to bace the difference
% plots

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

diff_matrix_path = [python_base_path '\run_' run '\difference_matrices_plots'];
if ~exist(diff_matrix_path, 'dir')
    mkdir(diff_matrix_path)
end

% Call the function for stress and deformation plotting
if plots_bool == 1
    figure_counter = depth_searcher(run,sd_input,coordinate_system,complete_matrices_path,...
        files_to_classify,figures_path,diff_matrix_path,iteration,step,cycle,min_lat,max_lat,...
        min_lon,max_lon,depths_to_plot,run_folder,figure_counter,python_base_path,run_vec,simul_time);
end

%% Viscosity 
% Again, decide whether to plot the viscosity or not and if so run the
% corresponding function after allocating the figures path
if strcmp(quantity,'stresses') == 1
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
        [figure_counter, visco_diff_path] = visco_plotter(sd_input,python_variables_path,...
            coordinate_system, complete_matrices_path, report_path,min_lat,max_lat,min_lon,...
            max_lon,depths_to_plot,iteration,step,cycle,diff_matrix_path,run_folder,run,...
            figure_counter,python_base_path,run_vec,simul_time);
    end
end
%% Time difference plots 
% Decide whether to generate the difference plots or not and if so geenarte
% the B coefficients and other quantities difference plots
diff_check = 1;
while diff_check == 1
    diff_input = input('Enter 1 to generate plots of differences between different cycles and steps, 0 to skip this:\n');
    if diff_input == 1 || diff_input == 0
        diff_check = 0;
    else
        disp('Incorrect input, enter 1 to generate difference plots or 0 to skip:\n');
    end
end
b_diff_plots_folder = [diff_matrix_path '\plots\viscosity\'];
if ~exist(b_diff_plots_folder, 'dir')
        mkdir(b_diff_plots_folder)
end
if diff_input == 1
    b_diff_calc(b_diff_plots_folder,python_base_path,min_lat,max_lat,...
    min_lon,max_lon,depths_to_plot(1),depths_to_plot(2),figure_counter)
    figure_counter = diff_calculator(diff_matrix_path,min_lat,max_lat,min_lon,max_lon,...
                run,depths_to_plot,iteration,step,figure_counter);
end


%% Multiple runs case
% The same descriptions as for the single run case apply. The main
% differences are iteration cycles to loop over, for example, runs and
% steps.

% run_vec = [25];
% iteration = '1';
% step_vec = [0 1];
% simul_time_vec = {'1 ka','31 ka'}; % Update for new timesteps
% for i = 1:length(run_vec)
%     run = num2str(run_vec(i));
%     for j=1:length(step_vec)
%         simul_time = simul_time_vec{j};
%         step = num2str(step_vec(j));
%         step_folder = [python_base_path '\run_' run '\Iteration_' iteration...
%             '\step_' step];
%         files = dir(step_folder);
%         dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
%         subFolders = files(dirFlags);
%         subfolder_names = extractfield(subFolders,'name');
%         cycle_vec = regexp(subfolder_names,'[\d\.]+','match');
%         for k=1:length(cycle_vec)
%             
%             cycle = char(cycle_vec{k});
%             disp([run iteration step cycle]);
%             coordinate_system = 'geographical'; % cartesian or geographical
%             quantity = 'stresses'; %stresses or deflections
%             % python_variables_path = [python_variables_base_path '\run_' run '\iteration_' iteration];
%             python_variables_path = [python_base_path '\run_' run '\Iteration_' iteration...
%                 '\step_' step '\cycle_' cycle '_reports'];
%             run_folder = [python_base_path '\run_' run];
%             
%             directory = dir(python_variables_path);
%             directory_filenames = {directory.name};
%             for a=1:length(directory_filenames)
%                 if ~contains(directory_filenames{a},'matlab')==1
%                     directory_filenames{a} = {};
%                 end
%             end
%             directory_filenames = directory_filenames(~cellfun('isempty',directory_filenames));
%             for a=1:length(directory_filenames)
%                 if contains(directory_filenames{a},coordinate_system)==1 && ...
%                         contains(directory_filenames{a},quantity)==1
%                     file_to_open = directory_filenames{a};
%                 end
%             end
%             
%             file_to_open_path = [python_variables_path '\' file_to_open];
%             opened_file = fopen(file_to_open_path);
%             
%             b = 1;
%             c = 1;
%             files_to_classify = {};
%             plotting_params = zeros(1,6);
%             while ~feof(opened_file)
%                 % Read a line without including the newline
%                 s=fgetl(opened_file);
%                 if contains(s, '\') ~= 1 && contains(s, '.csv') == 1
%                     files_to_classify{b} = strtrim(s);
%                 end
%                 %     if isnumeric(str2double(s)) == 1 && ~isnan(str2double(s))
%                 %         sd_input = str2double(s);
%                 %     end
%                 if strcmp(s, 'stresses') == 1 || strcmp(s, 'deflections') == 1
%                     if strcmp(s, 'stresses') == 1
%                         sd_input = 0;
%                     else
%                         sd_input = 1;
%                     end
%                 end
%                 if isempty(extractAfter(s, '_reports')) && contains(s, 'Iteration_')
%                     report_path  = s;
%                 end
%                 if contains(s, 'Depth') == 1
%                     classified_path = s;
%                 end
%                 if strcmp(s, 'cartesian')== 1 || strcmp(s, 'geographical')== 1
%                     components_to_plot = s;
%                 end
%                 % if contains(s, 'Complete') == 1 || contains(s, 'geographical_complete')
%                 if contains(s, 'omplete') == 1
%                     complete_matrices_path = s;
%                 end
%                 if abs(str2double(s)) > 1 && ~isnan(str2double(s))
%                     plotting_params(c) = str2double(s);
%                     c = c + 1;
%                 end
%                 b = b + 1;
%             end
%             
%             if sd_input == 0
%                 figure_folder = 'stress_plots';
%             else
%                 figure_folder = 'deflection_plots';
%             end
%             figures_path = [report_path '\' figure_folder];
%             if ~exist(figures_path, 'dir')
%                 mkdir(figures_path)
%             end
%             
% %             while_check = 1;
% %             while while_check == 1
% %                 plots_bool = input(['Do you want to plot the stresses or deformations? \n'...
% %                     'Enter 1 for plots, 0 to skip this step: \n']);
% %                 if plots_bool == 1 || plots_bool == 0
% %                     while_check = 0;
% %                 else
% %                     disp('Incorrect input, enter 1 for plots and 0 to skip this step.');
% %                 end
% %             end
%             plots_bool = 1;
%             diff_matrix_path = [python_base_path '\run_' run '\difference_matrices_plots\matrices'];
%             if ~exist(diff_matrix_path, 'dir')
%                 mkdir(diff_matrix_path)
%             end
%             % min_lat = plotting_params(1);
%             % max_lat = plotting_params(2);
%             % min_lon = plotting_params(3);
%             % max_lon = plotting_params(4);
%             %depths_to_plot = [plotting_params(5) plotting_params(6)];
%             min_lat = -90;
%             max_lat = -65;
%             min_lon = -180;
%             max_lon = 180;
%             depths_to_plot = [145 160];
%             if depths_to_plot(1) > depths_to_plot (2)
%                 depths_to_plot([1 2]) = depths_to_plot([2 1]);
%             end
% %             if plots_bool == 1
% %                 figure_counter = depth_searcher(run,sd_input,coordinate_system,complete_matrices_path,...
% %                     files_to_classify,figures_path,diff_matrix_path,iteration,step,cycle,...
% %                     min_lat,max_lat,min_lon,max_lon,depths_to_plot,run_folder,...
% %                     figure_counter,python_base_path,run_vec,simul_time);
% %             end
%             viscosity_figures_path = [report_path '\viscosity_plots'];
%             if ~exist(viscosity_figures_path, 'dir')
%                 mkdir(viscosity_figures_path)
%             end
% %             if sd_input == 0
% %                 [figure_counter, visco_diff_path] = visco_plotter(sd_input,python_variables_path,coordinate_system, complete_matrices_path,...
% %                     report_path,min_lat,max_lat,min_lon,max_lon,depths_to_plot,iteration,step,cycle,...
% %                     diff_matrix_path,run_folder,run,figure_counter,python_base_path,run_vec,simul_time);
% %             end
%             figure_counter = diff_calculator(diff_matrix_path,min_lat,max_lat,min_lon,max_lon,...
%                 run,depths_to_plot,iteration,step,figure_counter);
%         end
%     end
% end
