function [colorbarlimits] = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder)

list = dir([run_folder '\**\Geographical_complete_file_EARTH.csv']);
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

if sd_input == 0
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
        for l = 1:length(selected_columns)
            variable = condition_matrix(:,selected_columns(l));
            variable_extremes = [min(variable),max(variable)];
            stress_extremes_matrix(j,2*l-1:2*l) = variable_extremes;
        end
    end
    colorbarlimits = [min(stress_extremes_matrix(:,1:2:size(stress_extremes_matrix,2)-1))...
        max(stress_extremes_matrix(:,2:2:size(stress_extremes_matrix,2)))];
%     colorbarlimits = stress_extremes_matrix;
else
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
        lat = stress_matrix(:,end-1);
        lon = stress_matrix(:,end);
        min_depth = depths_to_plot(1);
        max_depth = depths_to_plot(2);
        depth_condtion = depth>min_depth & depth<max_depth;
        lat_condition = lat>min_lat & lat<max_lat;
        lon_condition = lon>min_lon & lon<max_lon;
        data_points_indices = defo_matrix(depth_condtion...
            & lat_condition & lon_condition);
        condition_matrix = defo_matrix(data_points_indices,:);
        for l = 1:length(selected_columns)
            variable = condition_matrix(:,selected_columns(l));
            variable_extremes = [min(variable),max(variable)];
            defo_extremes_matrix(k,2*l-1:2*l) = variable_extremes;
        end
        colorbarlimits = [min(defo_extremes_matrix(:,1:2:size(defo_extremes_matrix,2)-1))...
        max(defo_extremes_matrix(:,2:2:size(defo_extremes_matrix,2)))];
%         colorbarlimits = defo_extremes_matrix;
    end
end
end

