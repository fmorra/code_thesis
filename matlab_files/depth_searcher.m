function [figure_counter] = depth_searcher(run,sd_input,coordinate_system,complete_matrices_path,...
    ~,figures_path,diff_matrix_path,iteration,step,cycle, min_lat,max_lat,...
    min_lon,max_lon,depths_to_plot,run_folder,figure_counter,...
    python_base_path,r_earth)
% This function is used to select the data points to plot by giving an
% interval of points, as a substitute to the depth_classifier function.

if sd_input == 0
    components_to_plot = 'stresses';
else
    components_to_plot = 'deflections';
end

viscosity_input = 0;
b_input = 0;
complete_diff_path = [diff_matrix_path '\' components_to_plot  '\' coordinate_system];
if ~exist(complete_diff_path, 'dir')
    mkdir(complete_diff_path)
end
parts_to_plot = {'EARTH'};

resolution = 0.25;
latlim = [min_lat max_lat];
lonlim = [min_lon max_lon];
min_depth = depths_to_plot(1);
max_depth = depths_to_plot(2);
s = referenceSphere('Earth');
lat_lin = max_lat:-resolution:min_lat;
lon_lin = min_lon:resolution:max_lon;
[gridded_lon,gridded_lat] = meshgrid(lon_lin,lat_lin);
% plot_lon = reshape(gridded_lon, [numel(gridded_lon),1]);
% plot_lat = reshape(gridded_lat, [numel(gridded_lat),1]);
r_out = []; lon_plot_2 = []; lat_plot_2 = [];
depthrange = min_depth:1:max_depth;
for dd = depthrange
    lon_plot_2 = [lon_plot_2; gridded_lon(:)];
    lat_plot_2 = [lat_plot_2; gridded_lat(:)];
    temp = (s.Radius-dd*1e3)*ones(size(gridded_lon));
    r_out = [r_out; temp(:)];
end

if sd_input == 0
%     selected_component = input(['Enter the desired stress component(s) to plot,' ...
%         ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
    selected_component = 'Mises S11 S12';
    selected_components = split(selected_component);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder, viscosity_input,b_input,python_base_path);
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
        selected_components, run_folder, viscosity_input, python_base_path);
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
    lat_condition = lat>=min_lat & lat<=max_lat;
    lon_condition = lon>=min_lon & lon<=max_lon;
    data_points_indices = matrix_to_read(depth_condtion...
        & lat_condition & lon_condition); % Gives the indices
    
%     [long,lati] = meshgrid(floor(min(lon)):resolution:ceil(max(lon)),...
%         floor(min(lat)):resolution:ceil(max(lat)));
%     lati=flipud(lati);
%     start_lat = lati(1,1);
%     start_lon = long(1,1);
%     truncate_no_lat1 = floor((start_lat-latlim(2))*(1/resolution))+1;
%     truncate_no_lat2 = floor((start_lat-latlim(1))*(1/resolution))+1;
%     truncate_no_lon1 = floor((lonlim(1)-start_lon)*(1/resolution))+1;
%     truncate_no_lon2 = floor((lonlim(2)-start_lon)*(1/resolution))+1;
%     R = r_earth-depth(data_points_indices);
    
    matrix_for_difference = zeros(length(data_points_indices),length(selected_columns)+2);
    matrix_for_difference_headers = [selected_components','Latitude','Longitude'];
    for j = 1:length(selected_columns)
        plot_variable = matrix_to_read(data_points_indices,selected_columns(j));
        filtered_lat = lat(data_points_indices);
        filtered_lon = lon(data_points_indices);
        depth_out = depth(data_points_indices);
        filtered_R = s.Radius - 1e3*depth_out;
        matrix_for_difference(:,end-2) = plot_variable;
        matrix_for_difference(:,end-1) = filtered_lat;
        matrix_for_difference(:,end) = filtered_lon;
        
%         plot_variable_grid = griddata(lon(data_points_indices),lat(data_points_indices)...
%             ,plot_variable,long,lati,'linear');
%         if latlim(1)==-90
%             truncate_no_lat2 = length(plot_variable_grid(:,1));
%         end
%         if lonlim(1)==-180
%             truncate_no_lon1 = 1;
%         end
%         if lonlim(2)==180
%             truncate_no_lon2 = length(plot_variable_grid(1,:));
%         end
%         plot_variable = reshape(plot_variable_grid(truncate_no_lat1:truncate_no_lat2,...
%             truncate_no_lon1:truncate_no_lon2),[numel(gridded_lon),1]);
%         
%         figure(figure_counter)
%         [Z, refvec] = geoloc2grid(plot_lat,plot_lon,plot_variable,resolution);
%         refvec(2) = start_lat - (truncate_no_lat1-1)*resolution;
%         refvec(3) = start_lon + (truncate_no_lon1-1)*resolution;
%         load coastlines;
%         latlim = [min_lat max_lat];lonlim = [min_lon max_lon];
%         ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
%         set(ax,'Visible','off');
%         cmap = colormap('summer');
%         alpha 0.7;
%         colormap(cmap);
%         set(findall(gca, 'type', 'text'), 'visible', 'on')
%         geoshow(Z, refvec, 'DisplayType', 'texture');
%         plotm(coastlat,coastlon, 'color', rgb('OrangeRed'));
%         R_lin = R*ones(length(plot_lon),1);
        [x_in,y_in,z_in]=sph2cart(deg2rad(filtered_lon),deg2rad(filtered_lat),filtered_R);
        [x_out,y_out,z_out]=sph2cart(deg2rad(lon_plot_2),deg2rad(lat_plot_2),r_out);

        figure(1)
        scatter3(x_in,y_in,z_in,10,z_in)
        plot_variable_out = griddata(x_in,y_in,z_in,plot_variable,x_out,y_out,z_out,'nearest');
        figure(figure_counter+1)
        colormap summer;
        load coastlines;
        [Z, refvec] = geoloc2grid(lat_plot_2,lon_plot_2,plot_variable_out,resolution);
        ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
        set(ax,'Visible','off');
        cmap = colormap('summer');
        alpha 0.7;
        colormap(cmap);
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        geoshow(Z, refvec, 'DisplayType', 'texture');
        plotm(coastlat, coastlon, 'color', rgb('OrangeRed'));
        
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
                ', iteration ' num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)],[' ']});
        else 
            title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j} ...
                ', iteration ' num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)],[' ']});
        end
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        grid on;
        
        % Save the figures based on whether we are working with
        % stresses or deflections
        saveas(gcf,[figures_path '\' coordinate_system '_' components_to_plot '_' parts_to_plot{i}...
            '[' num2str(min_depth) '-' num2str(max_depth) ']_km_'...
            selected_components{j} '_' run '_' iteration '_' step '_' cycle '.png']);
        figure_counter = figure_counter + 1;
    end
    %if strcmp(parts_to_plot,'EARTH') == 1
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[complete_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
    % end
end
end