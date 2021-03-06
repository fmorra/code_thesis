function [figure_counter,visco_diff_path] = visco_plotter(sd_input,python_variables_base_path,...
    coordinate_system,complete_matrices_path,report_path,min_lat,max_lat,min_lon,max_lon,...
    depths_to_plot,iteration,step,cycle,diff_matrix_path,run_folder,run,figure_counter,...
    python_base_path,run_vec,simul_time)
close all; clc;
% disp('Viscosity can only be plotted for part EARTH, as it is the one we have the B values for. \n');
parts_to_plot = {'EARTH'};
% visco_strain_input = input('Enter 0 to plot the viscosity, 1 to plot the strain rate:\n');
% if visco_strain_input == 0
%     quantity = 'viscosity';
% else
%     quantity = 'strain_rate';
% end
visco_strain_input = 0;
if visco_strain_input == 0
    quantity = 'viscosity';
else
    quantity = 'strain_rate';
end
viscosity_input = 1;
visco_diff_path = [diff_matrix_path '\' quantity];
if ~exist(visco_diff_path, 'dir')
    mkdir(visco_diff_path)
end
e_file_name = 'e.dat';
viscosity_figures_path = [report_path '\viscosity_plots'];
if ~exist(viscosity_figures_path, 'dir')
    mkdir (viscosity_figures_path)
end
b_input = 0;
latlim = [min_lat max_lat];
lonlim = [min_lon max_lon];
min_depth = depths_to_plot(1);
max_depth = depths_to_plot(2);
% mean_depth = (max_depth+min_depth)/2;
resolution = 0.25;
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
    selected_components = cellstr(quantity);
    colorbarlimits = caxisextremes(sd_input,min_lat,max_lat,min_lon,max_lon,depths_to_plot,...
        selected_components,run_folder,viscosity_input,b_input,python_base_path,run_vec,...
        coordinate_system);
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
    
    matrix_to_read = complete_matrix_with_viscosity;
    depth = matrix_to_read(:,end-6)/1e3;
    lat = matrix_to_read(:,end-5);
    lon = matrix_to_read(:,end-4);
    depth_condtion = depth>min_depth & depth<max_depth;
    lat_condition = lat>min_lat & lat<max_lat;
    lon_condition = lon>min_lon & lon<max_lon;
    data_points_indices = matrix_to_read(depth_condtion...
        & lat_condition & lon_condition);
    
    matrix_for_difference = zeros(length(data_points_indices),4);
    matrix_for_difference_headers = {quantity,'Latitude','Longitude','Depth'};
    if visco_strain_input == 0
        plot_variable = matrix_to_read(data_points_indices,end-2);
        plot_variable = log10(plot_variable);
    else
        plot_variable = matrix_to_read(data_points_indices,end-1);
    end
    filtered_lat = lat(data_points_indices);
    filtered_lon = lon(data_points_indices);
    depth_out = depth(data_points_indices);
    filtered_R = s.Radius - 1e3*depth_out;
    matrix_for_difference(:,end-3) = plot_variable;
    matrix_for_difference(:,end-2) = filtered_lat;
    matrix_for_difference(:,end-1) = filtered_lon;
    matrix_for_difference(:,end) = depth_out;
    disp(max(mises(data_points_indices,:)))
    disp(max(strain_rate(data_points_indices,:)))
    disp(max(viscosity(data_points_indices,:)))
    disp(max(plot_variable))
%     filtered_R = R*ones(length(filtered_lon),1);
%     R_lin = R*ones(length(plot_lon),1);
    [x_in,y_in,z_in]=sph2cart(deg2rad(filtered_lon),deg2rad(filtered_lat),filtered_R);
    [x_out,y_out,z_out]=sph2cart(deg2rad(lon_plot_2),deg2rad(lat_plot_2),r_out);
    
%     figure(1)
%     scatter3(x_in,y_in,z_in,10,z_in)
%     p = knnsearch([x_in,y_in,z_in],[x_out,y_out,z_out],'k',6);
%     arclentot = zeros(size(plot_lat));
%     plot_variable_out = zeros(size(plot_lat));
%     coordinates = zeros(60,3);
%     filtered_matrix = matrix_to_read(data_points_indices,:);
%     for k = 1:10
%         for l = 1:6
%             label = p(k,l);
%             % index = filtered_matrix(:,1) == label;
%             coordinates(6*(k-1)+l,:) = filtered_matrix(label,10:12);
%         end
%     end
%     figure(3)
%     scatter3(x_out(1:10),y_out(1:10),z_out(1:10),'r'); 
%     figure(4)
%     scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'b'); 
%     figure(5)
%     scatter3(x_out(1:10),y_out(1:10),z_out(1:10),'r');hold on;
%     scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'b'); hold off;
%     for ii = 1:6
%         arclentot = arclentot + distance(lat(p(:,ii)),lon(p(:,ii)),plot_lat,plot_lon,s); % same unit as S, so km
%     end
%     for ii = 1:6
%         arclen = distance(lat(p(:,ii)),lon(p(:,ii)),plot_lat,plot_lon,s);
%         plot_variable_out = plot_variable_out+viscosity(p(:,ii)).*(arclen./arclentot);
%     end
    
    
    plot_variable_out = griddata(x_in,y_in,z_in,plot_variable,x_out,y_out,z_out,'nearest');
    figure()
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
    
    %caxis('auto')
    caxis(log10([colorbarlimits(1) colorbarlimits(2)]));
    h = colorbar('v'); % set(h, 'ylim', [0 1e6]);
    if visco_strain_input == 0
        set(get(h,'ylabel'),'string','log_{10}Viscosity [Ns/m^2]')
    else
        set(get(h,'ylabel'),'string','log_{10}Strain rate')
    end
    if mod(run,2)==0
        rheology = ', wet rheology';
    else
        rheology = ', dry rheology';
    end
    if visco_strain_input == 0
%         title({['Viscosity for part ' parts_to_plot{i} ' with depth range '], ...
%             [num2str(min_depth) '-' num2str(max_depth) ' km, iteration ' ...
%             num2str(iteration) ', step ' num2str(step) ', cycle ' num2str(cycle)],[' ']});
        %  ', stress iteration ' num2str(cycle)
%         title({['Viscosity map, timestep ' num2str(str2double(step) + 1) ', stress iteration ' num2str(cycle)],[' ']});
        title({['Time ' simul_time ', stress iteration ' num2str(cycle) ],[' ']});
        %title({['Time ' simul_time],[' ']});
    else
%         title({['Strain rate map, timestep ' num2str(str2double(step) + 1) ', stress iteration ' num2str(cycle)],[' ']});
        title({['Time ' simul_time ', stress iteration ' num2str(cycle) ],[' ']});
        %title({['Time ' simul_time],[' ']});
    end
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    grid on;
    % Save the figures based on whether we are working with
    % stresses or deflections
    if visco_strain_input == 0
        saveas(gcf,[viscosity_figures_path '\' 'Viscosity_' parts_to_plot{i}...
            '[' num2str(min_depth) '-' num2str(max_depth) ']_km_'...
            quantity '_' run '_' iteration '_' step '_' cycle '.png']);
    else
        saveas(gcf,[viscosity_figures_path '\' 'Strain_rate_' parts_to_plot{i}...
            '[' num2str(min_depth) '-' num2str(max_depth) ']_km_'...
            quantity '_' run '_' iteration '_' step '_' cycle '.png']);
    end
    matrix_for_difference_wheaders = [matrix_for_difference_headers; num2cell(matrix_for_difference)];
    writecell(matrix_for_difference_wheaders,[visco_diff_path '\Iteration_' iteration '_step_' step ...
        '_cycle_' cycle '_' num2str(min_depth) '_' num2str(max_depth) '_km.csv']);
    
end
% for i=1:length(parts_to_plot)
%     B_plots(viscosity_figures_path,alin,a,data_points_indices,parts_to_plot{i},...
%         min_depth,max_depth,min_lat,max_lat,min_lon,max_lon,lat,lon,run,...
%         run_folder,viscosity_input,python_base_path,depth,run_vec,...
%         coordinate_system,visco_diff_path);
% end
figure_counter = figure_counter+3;
end

