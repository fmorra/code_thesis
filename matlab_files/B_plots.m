function [] = B_plots(viscosity_figures_path,alin,a,data_points_indices,part,...
    min_depth,max_depth,min_lat,max_lat,min_lon,max_lon,lat,lon,run,run_folder,...
    viscosity_input,python_base_path,depth,run_vec,coordinate_system,visco_diff_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all; clc;
variables = [alin,a];
names = {'B_{diff}','B_{disl}'};
open_figures = findobj('type','figure');
sd_input = 0; % Only used when processing stresses
b_input = 1;
depths_to_plot = [min_depth max_depth];
selected_components = names;
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
matrix_for_difference = zeros(length(data_points_indices),5);
matrix_for_difference_headers = {'B_diff','B_disl','Latitude','Longitude','Depth'};
for dd = depthrange
    lon_plot_2 = [lon_plot_2; gridded_lon(:)];
    lat_plot_2 = [lat_plot_2; gridded_lat(:)];
    temp = (s.Radius-dd*1e3)*ones(size(gridded_lon));
    r_out = [r_out; temp(:)];
end
colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
    selected_components,run_folder,viscosity_input,b_input,python_base_path,run_vec,...
    coordinate_system);
filtered_lat = lat(data_points_indices);
filtered_lon = lon(data_points_indices);
depth_out = depth(data_points_indices);
matrix_for_difference(:,end-2) = filtered_lat;
matrix_for_difference(:,end-1) = filtered_lon;
matrix_for_difference(:,end) = depth_out;
B_path = [visco_diff_path '\B_coeffs'];
if ~exist(B_path, 'dir')
    mkdir(B_path)
end
for i=1:size(variables,2)
    plot_variable = variables(data_points_indices,i);
    matrix_for_difference(:,i) = plot_variable;
    filtered_R = s.Radius - 1e3*depth_out;
    [x_in,y_in,z_in]=sph2cart(deg2rad(filtered_lon),deg2rad(filtered_lat),filtered_R);
    [x_out,y_out,z_out]=sph2cart(deg2rad(lon_plot_2),deg2rad(lat_plot_2),r_out);
    latlim = [min_lat max_lat];
    lonlim = [min_lon max_lon];
    plot_variable_out = griddata(x_in,y_in,z_in,plot_variable,x_out,y_out,z_out,'nearest');
    
    figure(i)
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
    caxis([colorbarlimits(i) colorbarlimits(i+length(colorbarlimits)/2)]);
    disp([colorbarlimits(i) colorbarlimits(i+length(colorbarlimits)/2)])
    h = colorbar('v');
    if i == 1
        set(get(h,'ylabel'),'string','B_{diff}')
    else
        set(get(h,'ylabel'),'string','B_{disl}')
    end
    if mod(run,2)==0
        title(['Map of ' names{i} ', wet rheology']);
    else
        title(['Map of ' names{i} ', dry rheology']);
    end
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    grid on;
    % Save the figures based on whether we are working with
    % stresses or deflections
    saveas(gcf,[viscosity_figures_path '\' names{i} '_' part '_' run...
        '_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[B_path '\B_' num2str(run) '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
end
end
