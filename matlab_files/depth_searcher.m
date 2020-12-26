function [] = depth_searcher(run,sd_input,coordinate_system,complete_matrices_path,...
    ~,figures_path,diff_matrix_path,iteration,step,cycle,...
    min_lat,max_lat,min_lon,max_lon,depths_to_plot,run_folder)
% This function is used to select the data points to plot by giving an
% interval of points, as a substitute to the depth_classifier function.

if sd_input == 0
    components_to_plot = 'stresses';
else
    components_to_plot = 'deflections';
end
viscosity_input = 0;
complete_diff_path = [diff_matrix_path '\' components_to_plot  '\' coordinate_system];
if ~exist(complete_diff_path, 'dir')
    mkdir(complete_diff_path)
end

parts_to_plot = {'EARTH'};
if sd_input == 0
%     selected_component = input(['Enter the desired stress component(s) to plot,' ...
%         ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
    selected_component = 'Mises S11 S12';
    selected_components = split(selected_component);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder, viscosity_input);
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
%     selected_component = input(['Enter the desired deflection components(s) to plot,' ...
%         ' possible values are Magnitude, U1, U2, U3:\n']);
    selected_component = 'Magnitude U3';
    selected_components = split(selected_component);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder, viscosity_input);
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
        matrix_to_read = readmatrix([complete_matrices_path '\Geographical_complete_file_'...
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
    %data_points = matrix_to_read(data_points_indices,:);
    matrix_for_difference = zeros(length(data_points_indices),length(selected_columns)+2);
    matrix_for_difference_headers = [selected_components','Latitude','Longitude'];
    for j = 1:length(selected_columns)
        plot_variable = matrix_to_read(data_points_indices,selected_columns(j));
        matrix_for_difference(:,j) = plot_variable;
        matrix_for_difference(:,end-1) = lat(data_points_indices);
        matrix_for_difference(:,end) = lon(data_points_indices);
%         [Z, refvec] = geoloc2grid(lat(data_points_indices),wrapTo360(lon(data_points_indices)),...
%             plot_variable,0.5);
%         load coastlines;
%         cmap = colormap('jet');
%         alpha 0.7;
%         colormap(cmap);
%         latlim = [min_lat max_lat];lonlim = [min_lon max_lon];
%         ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
%         set(ax,'Visible','off');
%         set(findall(gca, 'type', 'text'), 'visible', 'on')
%         geoshow(Z, refvec, 'DisplayType', 'texture');
%         plotm(coastlat,coastlon);

        latlim = [-90 -65];
        lonlim = [-180 180];
        figure(j)
        worldmap(latlim,lonlim);
        load coastlines;
        [LatGrid, LonGrid] = meshgrid(linspace(min(lat(data_points_indices)), max(lat(data_points_indices)),500), ...
        linspace(min(lon(data_points_indices)), max(lon(data_points_indices)),500));
        stressgrid = griddata(lat(data_points_indices), lon(data_points_indices),...
            plot_variable, LatGrid, LonGrid,'v4');
        surfm(LatGrid, LonGrid, stressgrid);
        colormap summer;
        plotm(coastlat, coastlon, 'color', rgb('OrangeRed'));
        camroll(180)
        
        caxis([colorbarlimits(j) colorbarlimits(j+length(colorbarlimits)/2)]);
        h = colorbar('v');
        if sd_input == 0
            set(get(h,'ylabel'),'string',[selected_components{j} ' (Pa)'])
        else
            set(get(h,'ylabel'),'string',[selected_components{j} ' (m)'])
        end
        if sd_input == 0
            title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j} ...
                ', iteration ' num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)]});
        else 
            title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j} ...
                ', iteration ' num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)]});
        end
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        grid on;
        
        % Save the figures based on whether we are working with
        % stresses or deflections
        saveas(gcf,[figures_path '\' coordinate_system '_' components_to_plot '_' parts_to_plot{i}...
            '[' num2str(min_depth) '-' num2str(max_depth) ']_km_'...
            selected_components{j} '_' run '_' iteration '_' step '_' cycle '.png']);
    end
    %if strcmp(parts_to_plot,'EARTH') == 1
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[complete_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
    % end
end
end