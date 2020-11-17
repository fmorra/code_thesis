function [complete_deflection_path,spherical_complete_deflection_path] = ...
    associate_deflection_coord(deflection_paths,individual_path,headers_on, ...
    deflection_processing_path)
% Associate all deflections in all nodes to the corresponding node
% coordinates, in a similar way to what is done in the element stress
% function.

% As always, define necessary paths where to read files from and store
% files, and if the paths where to save files don't exist, create them

dir_csv_deflection_files = dir([deflection_paths '\*.csv']);
csv_deflection_files = {dir_csv_deflection_files.name};
complete_headers = {'Label','U_Magn','U_X','U_Y','U_Z','X','Y','Z'};
complete_deflection_path = [deflection_processing_path '\Complete_deflection_files'];
if ~exist(complete_deflection_path, 'dir')
    mkdir (complete_deflection_path)
end
spherical_complete_deflection_path = [deflection_processing_path '\Spherical_complete_deflection_files'];
if ~exist(spherical_complete_deflection_path, 'dir')
    mkdir (spherical_complete_deflection_path)
end
large_complete_file = [];
large_spherical_complete_file = [];

% Redefine the list of parts to work on, excluding empty files like I1
for i = 1:length(csv_deflection_files)
    size_to_evaluate = dir([deflection_paths '\' csv_deflection_files{i}]);
    if size_to_evaluate.bytes < 100
        csv_deflection_files{i} = {};
    end
end
csv_deflection_files = csv_deflection_files(~cellfun(@isempty,csv_deflection_files));

% Verify if the files are already all there based on whether the last file
% to be created is actually already in the folder

if isfile([spherical_complete_deflection_path '\Large_complete_spherical_deflection_file.csv'])
    disp(['The files containing deflections associated to the relative coordinates' ...
        ' already exist, moving on to classification of stress values based on depth.'])
else
    disp('Associating deflection components to corresponding nodes...')
    for i = 1:length(csv_deflection_files)
        % Read the individual deflection files and invert the y and Z
        % coordinate because of the difference in ABAQUS' coordiate system
        disp(['Processing the following deflection file: ', csv_deflection_files{i}]);
        individual_deflection_matrix = readmatrix([deflection_paths '\' csv_deflection_files{i}]);
        individual_deflection_matrix(:, [4 5]) = individual_deflection_matrix(:, [5 4]);
        individual_deflection_matrix(:, [3 4]) = individual_deflection_matrix(:, [4 3]);
        individual_coordinate_matrix = readmatrix([individual_path '\Nodes_Part_' csv_deflection_files{i}]);
        individual_coordinate_matrix(:, [3 4]) = individual_coordinate_matrix(:, [4 3]);
        individual_coordinate_matrix(:, [2 3]) = individual_coordinate_matrix(:, [3 2]);
        complete_deflection_matrix = [individual_deflection_matrix,individual_coordinate_matrix(:,2:4)];
        large_complete_file = [large_complete_file;complete_deflection_matrix];
        
        % Rotate deflection components in the spherical coordinate system, the
        % coordinates will be spherical in the geographic_plots function
        
        spherical_complete_matrix = rotate_deflections(complete_deflection_matrix);
        large_spherical_complete_file = [large_spherical_complete_file;spherical_complete_matrix];
        
        % Save individual matrices in different ways based on whether we
        % want headers or not
        if headers_on == 1
            complete_file = [complete_headers; num2cell(complete_deflection_matrix)];
            complete_file_spherical = [complete_headers; num2cell(spherical_complete_matrix)];
            writecell(complete_file,[complete_deflection_path '\Complete_file_' csv_deflection_files{i}]);
            writecell(complete_file_spherical,[spherical_complete_deflection_path '\Spherical_complete_file_' csv_deflection_files{i}]);
        else
            writematrix(complete_deflection_matrix,[complete_deflection_path '\Complete_file_' csv_deflection_files{i}]);
            writematrix(spherical_complete_matrix,[spherical_complete_deflection_path '\Spherical_complete_file_' csv_deflection_files{i}]);
        end
        
    end
    
    % Save the larger matrices also based on whether we want the headers or
    % not
    if headers_on == 1
        large_complete_file = [complete_headers; num2cell(large_complete_file)];
        large_spherical_complete_file = [complete_headers; num2cell(large_spherical_complete_file)];
        writecell(large_complete_file,[complete_deflection_path '\Large_complete_deflection_file.csv']);
        writecell(large_spherical_complete_file,[spherical_complete_deflection_path '\Large_complete_spherical_deflection_file.csv']);
    else
        writematrix(large_complete_file,[complete_deflection_path '\Large_complete_deflection_file.csv']);
        writematrix(large_spherical_complete_file,[spherical_complete_deflection_path '\Large_complete_spherical_deflection_file.csv']);
    end
end
end

