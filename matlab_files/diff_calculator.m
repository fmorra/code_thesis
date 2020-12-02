function [] = diff_calculator(diff_matrix_path_incomplete,min_lat,max_lat,min_lon,max_lon)
% Calculate the differences between two similar plots, where possible.
% [1 0 1 145 150]
% [1 0 2 145 150]

visco_choice = 0;
if visco_choice == 0
    stress_deflection_input = 1;
    ref_system_input = 1;
    if stress_deflection_input == 0
        quantity = 'stresses';
    else
        quantity = 'deflections';
    end
    if ref_system_input == 0
        ref_system = 'cartesian';
    else
        ref_system = 'geographical';
    end
    diff_matrix_path = [diff_matrix_path_incomplete '\' quantity '\' ref_system];
    diff_plots_folder = [diff_matrix_path_incomplete '\plots\' quantity '\' ref_system];
else
    visco_strain_input = 0;
    if visco_strain_input == 0
        quantity = 'viscosity';
    else
        quantity = 'strain_rate';
    end
    diff_matrix_path = [diff_matrix_path_incomplete '\' quantity];
    diff_plots_folder = [diff_matrix_path_incomplete '\plots\' quantity];
end
dir_diff_files = dir([diff_matrix_path '\*.csv']);

if ~exist(diff_plots_folder, 'dir')
    mkdir(diff_plots_folder)
end
diff_files = {dir_diff_files.name};

matrix_1_flag = 0;
while matrix_1_flag == 0
    quintuplet_1 = input(['Enter the iteration, step and cycle number, and the minimum and maximum depth '...
    'of the first matrix used to plot the difference with respect to another moment in time '...
    'as a vector of square brackets of five values:\n']);
    for i=1:length(diff_files)
        found_values_1 = str2double(regexp(diff_files{i}, '\d+', 'match'));
        if isequal(quintuplet_1,found_values_1)
            matrix_1_flag = 1;
            file_to_read_1 = diff_files{i};
            iteration_1 = quintuplet_1(1);
            step_1 = quintuplet_1(2);
            cycle_1 = quintuplet_1(3);
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

matrix_2_flag = 0;
while matrix_2_flag == 0
    quintuplet_2 = input(['Enter the run, step and cycle number, and the minimum and maximum depth '...
    'of the second matrix  used to plot the difference with respect to another moment in time '...
    'as a vector of square brackets of five values:\n']);
    for i=1:length(diff_files)
        found_values_2 = str2double(regexp(diff_files{i}, '\d+', 'match'));
        if isequal(quintuplet_2,found_values_2)
            matrix_2_flag = 1;
            file_to_read_2 = diff_files{i};
            iteration_2 = quintuplet_2(1);
            step_2 = quintuplet_2(2);
            cycle_2 = quintuplet_2(3);
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
min_depth = quintuplet_2(4);
max_depth = quintuplet_2(5);

if matrix_1_flag == 1 && matrix_2_flag == 1
    % Extract the headers of both matrices, intersect and display the
    % possible common variables
    [values_1,headers_1,~] = xlsread([diff_matrix_path '\' file_to_read_1]);
    [values_2,headers_2,~] = xlsread([diff_matrix_path '\' file_to_read_2]);
    possible_variables_for_plotting = intersect(headers_1,headers_2, 'stable');
    possible_variables_for_plotting{end-1} = {};
    possible_variables_for_plotting{end} = {};
    possible_variables_for_plotting = possible_variables_for_plotting(~cellfun('isempty',possible_variables_for_plotting));
    disp(['The possible variables to generate the difference plots between these'...
        ' two cycles are:\n']);
    disp(possible_variables_for_plotting);
    diff_variables_to_plot = possible_variables_for_plotting;
    lat = values_1(:,end-1);
    lon = values_1(:,end);
    for i=1:length(diff_variables_to_plot)
        diff_variable_1 = values_1(:,strcmp(headers_1,diff_variables_to_plot{i}));
        diff_variable_2 = values_2(:,strcmp(headers_2,diff_variables_to_plot{i}));
        diff = diff_variable_2-diff_variable_1;
        [Z, refvec] = geoloc2grid(lat,lon,diff,0.5);
        load coastlines;
        cmap = colormap('jet');
        alpha 0.7;
        colormap(cmap);
        caxis('auto'); % was [0 1e6]
        h = colorbar('v'); %set(h, 'ylim', [0 1e6]);
        if visco_choice == 1
            if visco_strain_input == 0
                set(get(h,'ylabel'),'string',[diff_variables_to_plot{i} ' [N \cdot s/m^2]'])
            else
                set(get(h,'ylabel'),'string',[diff_variables_to_plot{i}])
            end
        else
            if stress_deflection_input == 0
                set(get(h,'ylabel'),'string',[diff_variables_to_plot{i} ' [Pa]'])
            else
                set(get(h,'ylabel'),'string',[diff_variables_to_plot{i} ' [m]'])
            end
        end
        latlim = [min_lat max_lat];lonlim = [min_lon max_lon];
        ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
        set(ax,'Visible','off');
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        geoshow(Z, refvec, 'DisplayType', 'texture');
        plotm(coastlat,coastlon);
        if visco_choice == 0
            title({['Map of the ' ref_system ' ' quantity ' difference for part with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km and component ' diff_variables_to_plot{i} ','],...
                ['iterations [' num2str(iteration_1) '-' num2str(iteration_2) '], steps ['...
                num2str(step_1) '-' num2str(step_2) '], cycles [' num2str(cycle_1)...
                '-' num2str(cycle_2) ']']});
        else
            title({['Map of the ' quantity ' difference for part with depth range '],...
                [num2str(min_depth) '-' num2str(max_depth) 'km,' ],...
                ['iterations [' num2str(iteration_1) '-' num2str(iteration_2) '], steps ['...
                num2str(step_1) '-' num2str(step_2) '], cycles [' num2str(cycle_1)...
                '-' num2str(cycle_2) ']']}); 
        end
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        grid on;
        % Save the figures based on whether we are working with
        % stresses or deflections
        if visco_choice == 1
            saveas(gcf,[diff_plots_folder '\diff_plot_' quantity '_'...
                '_depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km_'...
                '_Iterations_[' num2str(iteration_1)...
                '_' num2str(iteration_2) ']_steps_[' num2str(step_1) '_' num2str(step_2)...
                ']_cycles[' num2str(cycle_1) '_' num2str(cycle_2) '].png']);
        else
            saveas(gcf,[diff_plots_folder '\diff_plot_' ref_system '_' quantity '_'...
                '_depth_range_[' num2str(min_depth) '-' num2str(max_depth) ']_km_Component_'...
                diff_variables_to_plot{i} '_Iterations_[' num2str(iteration_1)...
                '_' num2str(iteration_2) ']_steps_[' num2str(step_1) '_' num2str(step_2)...
                ']_cycles[' num2str(cycle_1) '_' num2str(cycle_2) '].png']);
        end
    end
end

end

