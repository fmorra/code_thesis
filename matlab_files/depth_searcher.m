function [figure_counter] = depth_searcher(run,sd_input,coordinate_system,complete_matrices_path,...
    ~,figures_path,diff_matrix_path,iteration,step,cycle, min_lat,max_lat,...
    min_lon,max_lon,depths_to_plot,run_folder,figure_counter,python_base_path,...
    run_vec,simul_time)
% This function is very similar to the stress and deformation 

% Define the quantities to plot and where to save matrices with values used
% to genearte difference plots. Viscosity and b_input also have to be
% defined as false since caxisextremes uses them as input since it is
% called to generate many different plots.

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

% Generate the 3D plotting grid by appending equally spaced latitude and
% longitude a number of times equal to the values in the depth range spaced
% by 1 km

resolution = 0.25;
latlim = [min_lat max_lat];
lonlim = [min_lon max_lon];
min_depth = depths_to_plot(1);
max_depth = depths_to_plot(2);
s = referenceSphere('Earth');
lat_lin = max_lat:-resolution:min_lat;
lon_lin = min_lon:resolution:max_lon;
[gridded_lon,gridded_lat] = meshgrid(lon_lin,lat_lin);
r_out = []; lon_plot_2 = []; lat_plot_2 = [];
depthrange = min_depth:1:max_depth;
for dd = depthrange
    lon_plot_2 = [lon_plot_2; gridded_lon(:)];
    lat_plot_2 = [lat_plot_2; gridded_lat(:)];
    temp = (s.Radius-dd*1e3)*ones(size(gridded_lon));
    r_out = [r_out; temp(:)];
end

% Define two cases,s tresses and deflections
if sd_input == 0
    % Select the stress components to plot, find the minimum and maximum
    % across all simulation runs and find the data columns for the specific
    % case selected
    selected_component = input(['Enter the desired stress component(s) to plot,' ...
        ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
%     selected_component = 'Mises S11';
    selected_components = split(selected_component);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder, viscosity_input,b_input,python_base_path,run_vec,coordinate_system);
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
    % Process the deformations doing the same operations as for the
    % stresses
    selected_component = input(['Enter the desired deflection components(s) to plot,' ...
        ' possible values are Magnitude, U1, U2, U3:\n']);
%     selected_component = 'Magnitude';
    selected_components = split(selected_component);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components, run_folder, viscosity_input,b_input,python_base_path,run_vec,coordinate_system);
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
% Plotting starts
for i=1:length(parts_to_plot)
    % Define depth extremes and file to read
    min_depth = depths_to_plot(1);
    max_depth = depths_to_plot(2);
    if strcmp(coordinate_system, 'cartesian') == 1
        matrix_to_read = readmatrix([complete_matrices_path '\Complete_file_'...
            parts_to_plot{i} '.csv']);
    else
        matrix_to_read = readmatrix([complete_matrices_path '\Geographical_complete_file_'...
            parts_to_plot{i} '.csv']);
    end
    % Read depth, latitude and longitude and filter the file to extract
    % points in the ranges selected
    depth = matrix_to_read(:,end-2)/1e3;
    lat = matrix_to_read(:,end-1);
    lon = matrix_to_read(:,end);
    depth_condtion = depth>min_depth & depth<max_depth;
    lat_condition = lat>=min_lat & lat<=max_lat;
    lon_condition = lon>=min_lon & lon<=max_lon;
    data_points_indices = matrix_to_read(depth_condtion...
        & lat_condition & lon_condition); % Gives the indices
    % Allocate the data matrix for the difference plots
    matrix_for_difference = zeros(length(data_points_indices),length(selected_columns)+3);
    matrix_for_difference_headers = [selected_components','Latitude','Longitude','Depth'];
    % Define the depth, lat, lon for the points in the selected ranges and
    % transfer them to the data matri for plot difference
    filtered_lat = lat(data_points_indices);
    filtered_lon = lon(data_points_indices);
    depth_out = depth(data_points_indices);
    filtered_R = s.Radius - 1e3*depth_out;
    matrix_for_difference(:,end-2) = filtered_lat;
    matrix_for_difference(:,end-1) = filtered_lon;
    matrix_for_difference(:,end) = depth_out;
    % Generate one plot for each selected column
    for j = 1:length(selected_columns)
        % Select current variable to plot and
        plot_variable = matrix_to_read(data_points_indices,selected_columns(j));
        matrix_for_difference(:,j) = plot_variable;
        % Convert coordinates to Cartesian for interpolation
        [x_in,y_in,z_in]=sph2cart(deg2rad(filtered_lon),deg2rad(filtered_lat),filtered_R);
        [x_out,y_out,z_out]=sph2cart(deg2rad(lon_plot_2),deg2rad(lat_plot_2),r_out);
        % Generate a scatter plot check for the mesh
        visual_check = figure();
        scatter3(x_in,y_in,z_in,10)
        xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
        title(['Distribution of points for the [' num2str(min_depth) '-' num2str(max_depth) '] km range']);
        grid on;
        saveas(gcf,[figures_path '\visual_mesh_check_depth_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
        close(visual_check);
        % Grid data and generate plots using the colorbar extremes from the
        % function called beforehand
        plot_variable_out = griddata(x_in,y_in,z_in,plot_variable,x_out,y_out,z_out,'nearest');
        figure(figure_counter)
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
        disp([colorbarlimits(j) colorbarlimits(j+length(colorbarlimits)/2)])
        % caxis('auto');
        h = colorbar('v');
        if sd_input == 0
            set(get(h,'ylabel'),'string',[selected_components{j} ' (Pa)'])
        else
            set(get(h,'ylabel'),'string',[selected_components{j} ' (m)'])
        end
        if mod(run,2)==0
            rheology = ', wet rheology';
        else
            rheology = ', dry rheology';
        end
        if sd_input == 0
            % ', stress iteration ' num2str(cycle)
%             title({['Map of ' components_to_plot_title ' component ' selected_components{j} ','],...
%                 ['timestep ' num2str(str2double(step) + 1) ', stress iteration ' num2str(cycle)],[' ']});
            title({['Time ' simul_time ', stress iteration ' num2str(cycle) rheology],[' ']});
            % title({['Time ' simul_time],[' ']});
        else 
%             title({['Map of the ' components_to_plot ' for part ' parts_to_plot{i} ' with depth range '],...
%                 [num2str(min_depth) '-' num2str(max_depth) 'km and component ' selected_components{j} ...
%                 ', iteration ' num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)],[' ']});
%             title({['Map of ' components_to_plot_title ' component ' selected_components{j} ','],...
%                 ['timestep ' num2str(str2double(step) + 1)  ', stress iteration ' num2str(cycle)],[' ']});
            title({['Time ' simul_time ', stress iteration ' num2str(cycle) rheology],[' ']});
            % title({['Time ' simul_time],[' ']});
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
    % Save the data matrix for plot difference
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[complete_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
    % end
end
end