 function [new_colorbarlimits] = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,...
     depths_to_plot,selected_components,~,viscosity_input,b_input,python_base_path,run_vec,...
     coordinate_system)
    % This function calculates the colorbar extremes by reading all the
    % simulation run folder in the directory where they are saved in order
    % to give the same range to all plots.
    
    % read all the runs and allocate a matrix for all the extremes, meaning
    % for every component of every run, then a vector for the extremes of
    % every component after calculating min and max on every column in the
    % previous matrix
    runs = run_vec;
    extremes_matrix = zeros(length(runs),2*length(selected_components));
    new_colorbarlimits = zeros(1,size(extremes_matrix,2));
    % Iterate over every simulation
    for run=1:length(runs)
        % Extract all the possible files for every iteration and step in
        % every run folder, for stresses or deflections
        run_folder = [python_base_path '\run_' num2str(runs(run))];
        if strcmp(coordinate_system,'cartesian') == 1
            list = dir([run_folder '\**\Complete_file_EARTH.csv']);
        else
            list = dir([run_folder '\**\Geographical_complete_file_EARTH.csv']);
        end
        names = extractfield(list,'name');
        names_paths = extractfield(list,'folder');
        full_stress_files = cell(length(names), 1);
        full_defo_files = cell(length(names), 1);
        for i=1:length(names)
            full_file = [names_paths{i} '\' names{i}];
            if contains(full_file,'stress') == 1
                full_stress_files{i} = full_file;
            else
                full_defo_files{i} = full_file;
            end
        end
        full_stress_files = full_stress_files(~cellfun('isempty',full_stress_files));
        full_defo_files = full_defo_files(~cellfun('isempty',full_defo_files));
        
        % Stress case
        if sd_input == 0 && viscosity_input == 0
            % Read the right data based on selected columns
            stress_extremes_matrix = zeros(length(full_stress_files),length(selected_components)*2);
            for j = 1:length(full_stress_files)
                stress_matrix = readmatrix(full_stress_files{j});
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
                    elseif strcmp(selected_components{k}, 'S23') == 1
                        selected_columns(k) = 8;
                    else
                        
                    end
                end
                % Only read the points of the selected variable which are
                % in the desired ranges
                depth = stress_matrix(:,end-2)/1e3;
                lat = stress_matrix(:,end-1);
                lon = stress_matrix(:,end);
                min_depth = depths_to_plot(1);
                max_depth = depths_to_plot(2);
                depth_condtion = depth>min_depth & depth<max_depth;
                lat_condition = lat>min_lat & lat<max_lat;
                lon_condition = lon>min_lon & lon<max_lon;
                data_points_indices = stress_matrix(depth_condtion...
                    & lat_condition & lon_condition);
                condition_matrix = stress_matrix(data_points_indices,:);
                % Calculate min and max of the current data columna nd
                % store them in conscutive cells in the current stress
                % iteration row
                for l = 1:length(selected_columns)
                    variable = condition_matrix(:,selected_columns(l));
                    variable_extremes = [min(variable),max(variable)];
                    stress_extremes_matrix(j,2*l-1:2*l) = variable_extremes;
                end
            end
            % Calculate the colorbar limits with min and max over every 2
            % matrix cells in order to do this for every component
            colorbarlimits = [min(stress_extremes_matrix(:,1:2:size(stress_extremes_matrix,2)-1))...
                max(stress_extremes_matrix(:,2:2:size(stress_extremes_matrix,2)))];
        % Deformations
        elseif sd_input == 1 && viscosity_input == 0
            % The procedure is exactly the same as for the stresses.
            defo_extremes_matrix = zeros(length(full_defo_files),length(selected_components)*2);
            for k = 1:length(full_defo_files)
                defo_matrix = readmatrix(full_defo_files{k});
                selected_columns = zeros(length(selected_components),1);
                for l = 1:length(selected_components)
                    if strcmp(selected_components{l}, 'Magnitude') == 1
                        selected_columns(l) = 2;
                    elseif strcmp(selected_components{l}, 'U1') == 1
                        selected_columns(l) = 3;
                    elseif strcmp(selected_components{l}, 'U2') == 1
                        selected_columns(l) = 4;
                    elseif strcmp(selected_components{l}, 'U3') == 1
                        selected_columns(l) = 5;
                    else
                        
                    end
                end
                depth = defo_matrix(:,end-2)/1e3;
                lat = defo_matrix(:,end-1);
                lon = defo_matrix(:,end);
                min_depth = depths_to_plot(1);
                max_depth = depths_to_plot(2);
                depth_condtion = depth>min_depth & depth<max_depth;
                lat_condition = lat>min_lat & lat<max_lat;
                lon_condition = lon>min_lon & lon<max_lon;
                data_points_indices = defo_matrix(depth_condtion...
                    & lat_condition & lon_condition);
                condition_matrix = defo_matrix(data_points_indices,:);
                % 
                for l = 1:length(selected_columns)
                    variable = condition_matrix(:,selected_columns(l));
                    variable_extremes = [min(variable),max(variable)];
                    defo_extremes_matrix(k,2*l-1:2*l) = variable_extremes;
                end
                colorbarlimits = [min(defo_extremes_matrix(:,1:2:size(defo_extremes_matrix,2)-1))...
                    max(defo_extremes_matrix(:,2:2:size(defo_extremes_matrix,2)))];
            end
        % Viscosity
        elseif viscosity_input == 1 && b_input == 0
            visco_extremes_matrix = zeros(length(full_stress_files),2);
            for l = 1:length(full_stress_files)
                % Read all the values necessary to calculate the viscosity
                % from the Mises stress
                e_path = [run_folder '\e.dat'];
                e = readmatrix(e_path);
                elements = e(:,1);
                alin = e(:,2);
                a = e(:,3);
                an = 3.5;
                strain_rate = zeros(length(e),1);
                opened_visco_matrix = readmatrix(full_stress_files{l});
                mises = opened_visco_matrix(:,2);
                pow = (an);
                % Calculate the viscosity and strain rate and filter them
                % for the ranges we are interested in
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
                                min_depth = depths_to_plot(1);
                max_depth = depths_to_plot(2);
                matrix_to_read = complete_matrix_with_viscosity;
                depth = matrix_to_read(:,end-6)/1e3;
                lat = matrix_to_read(:,end-5);
                lon = matrix_to_read(:,end-4);
                depth_condtion = depth>min_depth & depth<max_depth;
                lat_condition = lat>min_lat & lat<max_lat;
                lon_condition = lon>min_lon & lon<max_lon;
                data_points_indices = complete_matrix_with_viscosity(depth_condtion...
                    & lat_condition & lon_condition);
                condition_matrix = complete_matrix_with_viscosity(data_points_indices,:);
                % Same as before, but this time the matrix with min and max
                % only has 2 columns instead of 2*number of variables to
                % plot
                variable = condition_matrix(:,end-2);
                variable_extremes = [min(variable),max(variable)];
                visco_extremes_matrix(l,1:2) = variable_extremes;
                colorbarlimits = [min(visco_extremes_matrix(:,1))...
                    max(visco_extremes_matrix(:,2))];
            end
        % B coefficients
        elseif viscosity_input == 1 && b_input==1
            B_extremes_matrix = zeros(length(full_stress_files),2);
            for l = 1:length(full_stress_files)
                % Read the B coefficients and filter them for the desired
                % ranges
                e_path = [run_folder '\e.dat'];
                e = readmatrix(e_path);
                alin = e(:,2);
                a = e(:,3);
                opened_visco_matrix = readmatrix(full_stress_files{l});
                min_depth = depths_to_plot(1);
                max_depth = depths_to_plot(2);
                complete_matrix_B = [opened_visco_matrix, alin, a];
                matrix_to_read = complete_matrix_B;
                depth = matrix_to_read(:,end-4)/1e3;
                lat = matrix_to_read(:,end-3);
                lon = matrix_to_read(:,end-2);
                depth_condtion = depth>min_depth & depth<max_depth;
                lat_condition = lat>min_lat & lat<max_lat;
                lon_condition = lon>min_lon & lon<max_lon;
                data_points_indices = complete_matrix_B(depth_condtion...
                    & lat_condition & lon_condition);
                selected_columns = [alin(data_points_indices,:),a(data_points_indices,:)];
                % This time the matrix has 4 columns, for 2 variables
                for m = 1:size(selected_columns,2)
                    variable = selected_columns(:,m);
                    variable_extremes = [min(variable),max(variable)];
                    B_extremes_matrix(l,2*m-1:2*m) = variable_extremes;
                end
                colorbarlimits = [min(B_extremes_matrix(:,1:2:size(B_extremes_matrix,2)-1))...
                max(B_extremes_matrix(:,2:2:size(B_extremes_matrix,2)))];
            end
        else
            
        end
        extremes_matrix(run,:) = colorbarlimits;
    end
    % Calculate min and max over every column of the matrix which has a row
    % for each simulation run and the minima and maxima of the selected
    % values over every stress iteration in that run for its columns
    for minmax_counter=1:(size(extremes_matrix,2)/2)
        new_colorbarlimits(1,minmax_counter) = min(extremes_matrix(:,minmax_counter));
        new_colorbarlimits(1,minmax_counter+size(extremes_matrix,2)/2) = ...
            max(extremes_matrix(:,minmax_counter+size(extremes_matrix,2)/2));
    end
end

