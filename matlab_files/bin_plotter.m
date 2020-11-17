function [] = bin_plotter(sd_input,plots_bool,report_path,figures_path)

% Extract the possible layer depths divided into bins based on the
% subdivision carried out in the previous function, using regexp on the
% file names.

d = dir(classified_path);
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
dnames = natsortfiles({dfolders.name});
possible_depths = zeros(length(dnames),1);
for i = 1:length(dnames)
    if ~isempty(regexp(dnames{i}, '[\d\.]+', 'match','once'))
        possible_depths(i) = str2double(regexp(dnames{i}, '[\d\.]+', 'match','once'));
    else
        possible_depths(i) = 0;
    end
end
possible_depths(possible_depths==0) = [];
% Get the number of active generated figures to create new ones without
% overwriting the old plots and insert the bin number(s) relative to the
% layer we want to plot

cfn = length(findobj('type','figure'));

%% Plots - stress and deflection

if plots_bool == 1
    
    disp('Select model layers for stress plotting, possible values are (km):');
    fprintf('%4.f ', possible_depths)
    selected_parts = input('\nInput them a vector in square brackets: \n');
    
    % Here we have nested for cycles looping over the number of layers,
    % bins and components we want to plot. The outermost one is the layers
    % divided into a number of bins determined in the previous function.
    
    for i = 1:length(selected_parts)
        
        % Read all files where dat has been binned and sort them in natural
        % order
        
        classified_files_path = dir([classified_path '\Subclassification_depth_' num2str(selected_parts(i))...
            '_km\*.csv']);
        natsortfiles({classified_files_path.name});
        classified_files = {classified_files_path.name};
        all_depths = zeros(length(classified_files),1);
        
        
        % Read all the possible bins for that layer depth and ask for user
        % input to see which ones to plot
        
        for j = 1:length(classified_files)
            all_depths(j) = str2double(regexp(extractAfter(classified_files{j},'bin'), '[\d\.]+', 'match','once'));
        end
        
        disp(['The layer bins for the selected depth(s), ' num2str(selected_parts(i)) ' km, are:']);
        fprintf('%d ', sort(all_depths));
        selected_depths = input('\nEnter the desired bin(s) for stress plotting, as a vector in square brackets:\n');
        if length(intersect(selected_depths,all_depths))<length(selected_depths)
            disp('One or more values are outside the possible bins, program will be terminated.')
            break;
        end
        
        % Now ask for the component(s) to plot and selec the column
        % containing data accordingly, both for stresses and deflections
        
        if sd_input == 0
            selected_component = input(['Enter the desired stress component(s) to plot,' ...
                ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
            selected_components = split(selected_component);
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
            selected_component = input(['Enter the desired deflection components(s) to plot,' ...
                ' possible values are Magnitude, U1, U2, U3:\n']);
            selected_components = split(selected_component);
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
        
        % Plot a map for each layer, bin and component selected in the
        % nested for loops
        
        for k = 1:length(selected_depths)
            % Open data matrix
            matrix_to_read = readmatrix([classified_path '\subclassification_depth_' ...
                num2str(selected_parts(i)) '_km\' matrix_to_read_part num2str(selected_parts(i))...
                '_km_bin_' num2str(selected_depths(k)) '.csv']);
            for l = 1:length(selected_columns)
                % Read data columns to plot
                variable_to_plot =  matrix_to_read(:,selected_columns(l));
                lon_to_plot = matrix_to_read(:,end); % latitude [deg]
                lat_to_plot = matrix_to_read(:,end-1); %longitude [deg]
                % Create figure
                figure('Position', get(0, 'Screensize'))
                % Load world map with coastlines
                %worldmap world;
                worldmap ([-90 -70],[-180 -45])
                load coastlines;
                plotm(coastlat, coastlon);
                % Grid the data before plotting, and generate a surface
                [LatGrid, LonGrid] = meshgrid(linspace(min(lat_to_plot), max(lat_to_plot),1000), ...
                    linspace(min(lon_to_plot), max(lon_to_plot),1000));
                stressgrid = griddata(lat_to_plot, lon_to_plot, variable_to_plot, LatGrid, LonGrid);
                surfm(LatGrid, LonGrid, stressgrid);
                % Figure propeties definition
                if sd_input == 0
                    t = title(['Map of the ' components_to_plot ' stresses for layer ' num2str(selected_parts(i)) ...
                        ' and bin ' num2str(selected_depths(k)) ' with component ' selected_components{l}]);
                    pos = get(t, 'position');
                    set(t, 'position', [1.5e6 -7.7e6 0]);
                else
                    t = title(['Map of the ' components_to_plot ' deflections for layer ' num2str(selected_parts(i)) ...
                        ' and bin ' num2str(selected_depths(k)) ' with component ' selected_components{l}]);
                    pos =get (t, 'position');
                    set(t, 'position', [1.5e6 -7.7e6 0]);
                end
                
                alpha 0.5;
                c = colorbar;
                if sd_input == 0
                    c.Label.String = 'Stress [Pa]';
                else
                    c.Label.String = 'Deflections  [m]';
                end
                c.Location = 'southoutside';
                grid on;
                % Save the figures based on whether we are working with
                % stresses or deflections
                saveas(gcf,[figures_path '\' matrix_to_read_part num2str(selected_parts(i))...
                    '_bin_' num2str(selected_depths(k)) '_Component_' selected_components{l} ...
                    '_' components_to_plot '.png']);
            end
            
        end
        
    end
else
    
end

end

