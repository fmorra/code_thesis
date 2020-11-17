function [visco_diff_path] = visco_plotter(python_variables_base_path,coordinate_system, ...
    complete_matrices_path,report_path,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
    iteration,step,cycle,diff_matrix_path)

disp('Viscosity can only be plotted for part EARTH, as it is the one we have the B values for. \n');
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
for i=1:length(parts_to_plot)
    e_file_name = 'e.dat';
    viscosity_figures_path = [report_path '\viscosity_plots'];
    if ~exist(viscosity_figures_path, 'dir')
        mkdir (viscosity_figures_path)
    end
    
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
    
    min_depth = depths_to_plot(1);
    max_depth = depths_to_plot(2);
    matrix_to_read = complete_matrix_with_viscosity;
    depth = matrix_to_read(:,end-6)/1e3;
    lat = matrix_to_read(:,end-5);
    lon = matrix_to_read(:,end-4);
    depth_condtion = depth>min_depth & depth<max_depth;
    lat_condition = lat>min_lat & lat<max_lat;
    lon_condition = lon>min_lon & lon<max_lon;
    data_points_indices = matrix_to_read(depth_condtion...
        & lat_condition & lon_condition);
    matrix_for_difference = zeros(length(data_points_indices),3);
    matrix_for_difference_headers = {quantity,'Latitude','Longitude'};
    if visco_strain_input == 0
        plot_variable = matrix_to_read(data_points_indices,end-2);
        plot_variable = log10(plot_variable);
    else
        plot_variable = matrix_to_read(data_points_indices,end-1);
    end
    matrix_for_difference(:,end-2) = plot_variable;
    matrix_for_difference(:,end-1) = lat(data_points_indices);
    matrix_for_difference(:,end) = lon(data_points_indices);
    load coastlines;
    latlim = [min_lat max_lat];
    lonlim = [min_lon max_lon];
    [Z, refvec] = geoloc2grid(lat(data_points_indices),lon(data_points_indices),...
        plot_variable,0.5);
    cmap = colormap('jet');
    colormap(cmap);
    caxis('auto');
    ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
    plotm(coastlat,coastlon);
    set(ax,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    geoshow(Z, refvec, 'DisplayType', 'texture');
    
%     worldmap (latlim,lonlim)
%     % Grid the data before plotting, and generate a surface
%     [LatGrid, LonGrid] = meshgrid(linspace(min(lat(data_points_indices)), max(lat(data_points_indices)),1000), ...
%         linspace(min(lon(data_points_indices)), max(lon(data_points_indices)),1000));
%     plotm(coastlat,coastlon);
%     stressgrid = griddata(lat(data_points_indices), lon(data_points_indices), plot_variable, LatGrid, LonGrid);
%     surfm(LatGrid, LonGrid, stressgrid);
    
    h = colorbar('v'); % set(h, 'ylim', [0 1e6]);
    if visco_strain_input == 0
        set(get(h,'ylabel'),'string','Viscosity [Ns/m^2]')
    else
        set(get(h,'ylabel'),'string','Strain rate')
    end
    if visco_strain_input == 0
        title({['Viscosity for part ' parts_to_plot{i} ' with depth range '], ...
               [num2str(min_depth) '-' num2str(max_depth) ' km']});
    else
        title({['Strain rate for part ' parts_to_plot{i} ' with depth range '], ...
                [num2str(min_depth) '-' num2str(max_depth) ' km']})
    end
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    grid on;
    % Save the figures based on whether we are working with
    % stresses or deflections
    if visco_strain_input == 0
        saveas(gcf,[viscosity_figures_path '\' 'Viscosity_' parts_to_plot{i}...
            '_depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
    else
        saveas(gcf,[viscosity_figures_path '\' 'Strain_rate_' parts_to_plot{i}...
            '_depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
    end
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[visco_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
end
for i=1:length(parts_to_plot)
    B_plots(viscosity_figures_path,alin,a,data_points_indices,parts_to_plot{i},...
        min_depth,max_depth,min_lat,max_lat,min_lon,max_lon,lat,lon);
end

end

