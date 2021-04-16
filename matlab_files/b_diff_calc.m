function [] = b_diff_calc(b_diff_plots_folder,python_base_path,min_lat,max_lat,...
    min_lon,max_lon,min_depth,max_depth,figure_counter)
% This function calculates difference plots for the B coefficients
% Select runs for which they are different
runs = [25,26];

all_diff_files = {};
diff_paths = {};
count = 1;
for i=1:length(runs)
    diff_matrix_path = [python_base_path '\run_' num2str(runs(i)) ...
        '\difference_matrices_plots\viscosity\B_coeffs\'];
    diff_paths{i} = diff_matrix_path;
    dir_diff_files = dir([diff_matrix_path '\*.csv']);
    diff_files = {dir_diff_files.name};
    for j=1:length(diff_files)
        all_diff_files{count} = diff_files(j);
        count = count + 1;
    end
end

% Mesh data
resolution = 0.25;
latlim = [min_lat max_lat];
lonlim = [min_lon max_lon];
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

% Find data matrix for the first run
matrix_1_flag = 0;
while matrix_1_flag == 0
    for i=1:length(all_diff_files)
        triplet_1 = [runs(1) min_depth max_depth];
        a = regexp(all_diff_files{i}, '\d+', 'match');
        found_values_1 = str2double(a{1});
        if isequal(triplet_1,found_values_1)
            matrix_1_flag = 1;
            file_to_read_1 = all_diff_files{i};
            break;
        end
    end
    if matrix_1_flag == 1
        break;
    else
        disp('No match for the selected values, the matrix does not exist. Exiting search.');
        break;
    end
end

% Find data matrix for the second run
matrix_2_flag = 0;
while matrix_2_flag == 0
    for i=1:length(all_diff_files)
        triplet_2 = [runs(2) min_depth max_depth];
        b = regexp(all_diff_files{i}, '\d+', 'match');
        found_values_2 = str2double(b{1});
        if isequal(triplet_2,found_values_2)
            matrix_2_flag = 1;
            file_to_read_2 = all_diff_files{i};
            break;
        end
    end
    if matrix_2_flag == 1
        break;
    else
        disp('No match for the selected values, the matrix does not exist. Exiting search.');
        break;
    end
end
min_depth = triplet_2(2);
max_depth = triplet_2(3);

% Again, filter the data, interpolate and plot it. Save the figure
% afterwards.
if matrix_1_flag == 1 && matrix_2_flag == 1
    % Extract the headers of both matrices, intersect and display the
    % possible common variables
    [values_1,headers_1,~] = xlsread([char(diff_paths{1}) char(file_to_read_1)]);
    [values_2,headers_2,~] = xlsread([char(diff_paths{2}) char(file_to_read_2)]);
    possible_variables_for_plotting = intersect(headers_1,headers_2, 'stable');
    possible_variables_for_plotting{end-2} = {}; 
    possible_variables_for_plotting{end-1} = {};
    possible_variables_for_plotting{end} = {};
    possible_variables_for_plotting = possible_variables_for_plotting(~cellfun('isempty',possible_variables_for_plotting));
    disp(['The possible variables to generate the difference plots between these'...
        ' two cycles are:\n']);
    disp(possible_variables_for_plotting);
    diff_variables_to_plot = possible_variables_for_plotting;
    filtered_lat = values_1(:,end-2);
    filtered_lon = values_1(:,end-1);
    depth_out = values_1(:,end);
    filtered_R = s.Radius - 1e3*depth_out;
    for i=1:length(diff_variables_to_plot)
        diff_variable_1 = values_1(:,strcmp(headers_1,diff_variables_to_plot{i}));
        diff_variable_2 = values_2(:,strcmp(headers_2,diff_variables_to_plot{i}));
        diff = diff_variable_2-diff_variable_1;
        [x_in,y_in,z_in]=sph2cart(deg2rad(filtered_lon),deg2rad(filtered_lat),filtered_R);
        [x_out,y_out,z_out]=sph2cart(deg2rad(lon_plot_2),deg2rad(lat_plot_2),r_out);
        plot_variable_out = griddata(x_in,y_in,z_in,diff,x_out,y_out,z_out,'nearest');
        
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
        caxis('auto'); % was [0 1e6]
        h = colorbar('v'); %set(h, 'ylim', [0 1e6]);
        if i == 1
            set(get(h,'ylabel'),'string','B_{diff}')
            title('Map of differences for B_{diff}');
        else
            set(get(h,'ylabel'),'string','B_{disl}')
            title('Map of differences for B_{disl}');
        end
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        grid on;
        % Save the figures based on whether we are working with
        % stresses or deflections
        saveas(gcf,[b_diff_plots_folder 'diff_plot_' possible_variables_for_plotting{i} '_'...
            'depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
        figure_counter = figure_counter + 1;
    end
end
end


