function [complete_files_path,spherical_complete_files_path] = ...
    associate_stress_coord(individual_element_paths,individual_path,...
    large_node_matrix_path,headers_on,stress_matrices_paths,coord_file)
% This function uses the results of both nodes_processor and
% elements_processor to calculate the centroid coordinates for each element
% and associate them to the corresponding element and stress components. 

% As before, define the relevant file paths and folders, creating them if
% they do not exist

individual_stress_paths = stress_matrices_paths;
dir_csv_elem_files = dir([individual_element_paths '\*.csv']);
csv_elem_files = {dir_csv_elem_files.name};
dir_csv_stress_files = dir([individual_stress_paths '\*.csv']);
csv_stress_files = {dir_csv_stress_files.name};
centroid_files_path = [individual_path '\Centroids'];
if ~exist(centroid_files_path, 'dir')
    mkdir (centroid_files_path)
end
complete_files_path = [individual_path '\Complete_Files'];
if ~exist(complete_files_path, 'dir')
    mkdir (complete_files_path)
end
spherical_centroids_files_path = [individual_path '\Spherical_Centroids'];
if ~exist(spherical_centroids_files_path, 'dir')
    mkdir (spherical_centroids_files_path)
end
spherical_components_files_path = [individual_path '\Spherical_Components'];
if ~exist(spherical_components_files_path, 'dir')
    mkdir (spherical_components_files_path)
end
spherical_complete_files_path = [individual_path '\Spherical_Complete_Files'];
if ~exist(spherical_complete_files_path, 'dir')
    mkdir (spherical_complete_files_path)
end

% Define headers if we want them in the main function

complete_headers = {'Label', 'S_Mises', 'S_11', 'S_22', 'S_33', 'S_12', 'S_13', 'S_23', ...
     'X_centroid', 'Y_centroid', 'Z_centroid'};
centroid_headers = {'Centr_label', 'X_centroid', 'Y_centroid', 'Z_centroid'};

% Initialize necessary variables

all_nodes = readmatrix(large_node_matrix_path);
large_centroids = [];
large_complete_file = [];
spherical_earth_point_large_complete = [];


% The first external for cycle runs over the list of csv files contaning
% the stress components for each part and region. However, some adjustments
% have to be made. Because the region stress files all refer to a single
% larger part element file in EARTH_POINT and have similar names, every
% time a file is processed its name is removed from the list of stress
% files names. Then, the first file in the list is always processed, since
% the previous ones have been deleted from the list of possible files to
% read.

for i = 1:length(csv_stress_files)
    size_to_evaluate = dir([dir_csv_stress_files(i).folder '\' csv_stress_files{i}]);
    if size_to_evaluate.bytes < 100
        csv_stress_files{i} = {};
    end
end
csv_stress_files = csv_stress_files(~cellfun(@isempty,csv_stress_files));

% Verify if the files are already all there based on whether the last file
% to be created is actually already in the folder

if isfile([spherical_complete_files_path '\Spherical_complete_file_' coord_file(1:5) '.csv'])
    disp(['The files containing centroid stresses associated to the relative coordinates' ...
        ' already exist, moving on to classification of stress values based on depth.'])
else 
    disp('Associating stress components to corresponding centroids...')
    for i = 1:length(csv_stress_files)
        % Here we create the list of stress files and extract the first file in
        % the list, which is the one to read
        
        csv_stress_files = csv_stress_files(~cellfun('isempty',csv_stress_files));
        csv_elem_files = csv_elem_files(~cellfun('isempty',csv_elem_files));
        part_file = csv_stress_files{1};
        disp(['Processing the following stress file: ', part_file]);
        part_matrix = readmatrix([individual_stress_paths '\' part_file]);
        
        % Here we find the name of the stress file to open and read
        
        region_files_number = nnz(find(contains(csv_stress_files,'Region')));
        region_flag = contains(part_file,'Region');
        if region_flag == 1
            file_to_read_logical = contains(extractAfter(csv_elem_files,'Part_'),extractBefore(part_file,'_Region')) == 1;
        else
            file_to_read_logical = strcmp(extractAfter(csv_elem_files,'Part_'),part_file) == 1;
        end
        file_to_read_index = (file_to_read_logical == 1);
        
        % Read all elements for a certain part and allocate a matrix where
        % to store centroid coordinates
        all_elems = readmatrix([individual_element_paths '\' csv_elem_files{file_to_read_index}]);
        centroid_coord = zeros(size(part_matrix,1),4);
        % complete_file = zeros(size(part_matrix,1),4);
        % This for cycle is used to calculate the centroid coordinates for each
        % element of each part. The external for cycle loops over the elements
        % of that part.
        
        for j=1:size(part_matrix,1)
            % Here we use the element label from the stress matrix of the part
            % or region to find its node labels based on its label in the
            % element matrix of the part, then we read the nodes for that
            % element
            
            elem_nodes_index = all_elems(:,1) == part_matrix(j,1);
            elem_nodes = all_elems(elem_nodes_index,2:end);
            
            % 6-node elements will have two zeros at the end that we will take
            % out, and we can do this because we process on row at a time
            
            elem_nodes = nonzeros(elem_nodes);
            
            % Allocate a matrix where to store the three coordinates of each node
            % of the considered element, which means an 8- or 6-by-3 matrix
            % depending on the element type itself, then calculate the
            % centroid coordinates as the mean of the coordinates of each
            % node
            
            elem_coord_matrix = zeros(length(elem_nodes),3);
            
            for k = 1:size(elem_nodes,1)
                index = all_nodes(:,1) == elem_nodes(k,1);
                nodes_coord = all_nodes(index,2:4);
                if ~isempty(nodes_coord)
                    elem_coord_matrix(k,1:end) = nodes_coord;
                end
            end
            for l = 1:3
                centroid_coord(j,l+1) = mean(elem_coord_matrix(1:end,l));
            end
            centroid_coord(j,1) = part_matrix(j,1);
            if centroid_coord(j,2:4) == 0
                centroid_coord(j,:) = [];
            end
        end
        
        % Swap the Y and Z coordinates because ABAQUS uses another
        % coordinate system, and for the same reason swap the stress
        % components around 
        
        centroid_coord(:, [3 4]) = centroid_coord(:, [4 3]);
        centroid_coord(:, [2 3]) = centroid_coord(:, [3 2]);
        part_matrix(:, [3 5]) = part_matrix(:, [5 3]);
        part_matrix(:, [3 4]) = part_matrix(:, [4 3]);
        part_matrix(:, [7 6]) = part_matrix(:, [6 7]);
        part_matrix(:, [6 8]) = part_matrix(:, [8 6]);
        
        % Write both the centroid coordinates and the complete file with
        % centroid coordinates and stress components for each region, then
        % delete the name of the file that has just been read from the list of
        % stress csv files.
        
        complete_file = [part_matrix,centroid_coord(:,2:4)];
        if (region_flag == 1)
            large_centroids = [large_centroids; centroid_coord];
            large_complete_file = [large_complete_file; complete_file];
        end
        
        centroid_individual_path = [centroid_files_path '\Centroid_file_' part_file];
        complete_individual_path = [complete_files_path '\Complete_file_' part_file];
        
        % As before, determine how to save the matrix based on the presence
        % of headers or not
        
        if headers_on == 1
            complete_file_headers = [complete_headers; num2cell(complete_file)];
            centroid_coord_headers = [centroid_headers; num2cell(centroid_coord)];
            writecell(complete_file_headers,complete_individual_path);
            writecell(centroid_coord_headers,centroid_individual_path);
        else
            writematrix(complete_file,complete_individual_path);
            writematrix(centroid_coord,centroid_individual_path);
        end
        
        % Here the rotation of the stress tensor into geographical coordinates
        % is performed.
        complete_file_spherical = rotate_tensor(part_matrix,spherical_centroids_files_path, ...
            spherical_components_files_path,spherical_complete_files_path,headers_on,...
            complete_headers,centroid_headers,complete_individual_path,part_file);
        
        % Create the complete matrix containing the spherical stress components
        % for part Earth too, needed for the depth classification.
        if (region_flag == 1)
            spherical_earth_point_large_complete = [spherical_earth_point_large_complete; complete_file_spherical];
        end
        
        csv_stress_files{1} = [];
        if region_flag == 0
            csv_elem_files{file_to_read_logical} = [];
        elseif region_files_number == 0
            csv_elem_files{file_to_read_logical} = [];
        end
    end
    
    if headers_on == 1
        large_centroids = [centroid_headers; num2cell(large_centroids)];
        large_complete_file = [complete_headers; num2cell(large_complete_file)];
        spherical_earth_point_large_complete = [complete_headers; num2cell(spherical_earth_point_large_complete)];
        writecell(large_centroids,[centroid_files_path '\Centroid_coordinates_' coord_file(1:5) '.csv']);
        writecell(large_complete_file,[complete_files_path '\Complete_file_' coord_file(1:5) '.csv']);
        writecell(spherical_earth_point_large_complete,[spherical_complete_files_path '\Spherical_complete_file_' coord_file(1:5) '.csv']);
    else
        writematrix(large_centroids,[centroid_files_path '\Centroid_coordinates_' coord_file(1:5) '.csv']);
        writematrix(large_complete_file,[complete_files_path '\Complete_file_' coord_file(1:5) '.csv']);
        writematrix(spherical_earth_point_large_complete,[spherical_complete_files_path '\Spherical_complete_file_' coord_file(1:5) '.csv']);
    end
    
end

end

