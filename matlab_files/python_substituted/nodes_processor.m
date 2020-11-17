function [individual_node_path,large_node_matrix_path] = ...
    nodes_processor(common_files_path,node_lines,elem_lines,g,file_identifiers,headers_on)

% This function uses the .dat file produced by ABAQUS and processed in 
% MATLAB  to extract the coordinates of all nodes and stored them in a 
% separate csv file for each model part.

% Initialize matrices and paths, and create the folder storing the nodes
% for each part if the folder is not already there.

large_node_matrix_path = [common_files_path '\' 'Large_Node_Matrix' '.csv'];

individual_node_path = [common_files_path '\Nodes'];
if ~exist(individual_node_path, 'dir')
    mkdir (individual_node_path)
end

large_node_matrix_unsorted = [];

headers = {'Label', 'X', 'Y', 'Z'};
% We loop only over the first n-1 elements of the vector containing the
% lines where the Node keyword is found, because the last line in the
% vector is the one where we stop processing as there are no interesting
% values anymore after that. 

% Verify if the files are already all there based on whether the last file
% to be created is actually already in the folder

if isfile(large_node_matrix_path)
    disp(['The files containing the node coordinates already exist, moving on' ...
        ' to element nodes extraction for each part.'])
else
    % Run the function if the files are not there
    disp('Processing nodes to extract coordinates...')
    for i = 1:length(node_lines)-1
        % Loop over the lines of the dat files between the Node and Eelement
        % keywords, because in these parts the node coordinates are stored.
        individual_matrix = zeros(elem_lines(i)-node_lines(i),4);
        disp(['Extracting coordinates for part: ' file_identifiers{i}])
        for j = node_lines(i):1:elem_lines(i)
            % Scan each line for floats
            s=char(g{1}(j));
            data= sscanf(s, '%f');
            if length(data) > 3
                if data(2) == data(3) && data(3) == data(4) && data(2)~= 0
                    
                else
                    data = data';
                    
                    % Every 5 lines, there will be an additional element which is
                    % the number right after LINE, in which we are not interested,
                    % so we must take it away.
                    
                    if length(data) == 5
                        data(1) = [];
                    end
                    
                    % Save the float values if the data vector is not empty
                    
                    % individual_matrix(j-node_lines(i)+1,:) = data;
                    individual_matrix(j,:) = data;
                end
            end
        end
        
        % Define different matrices to save> one for each region, and two
        % larger ones containing all values of every part both in sorted and
        % unsorted order. Unsorted because this way we preserve the enumeration
        % system used by ABAQUS. It is not necessary to store the nodes for
        % each region here.
        
        individual_matrix(~any(individual_matrix,2), : ) = [];
        large_node_matrix_unsorted = [large_node_matrix_unsorted; individual_matrix];
        large_node_matrix = sortrows(large_node_matrix_unsorted);
        if headers_on == 1
            individual_matrix = [headers; num2cell(individual_matrix)];
            individual_matrix(cellfun(@(x) ~x(1),individual_matrix(:,1)),:) = [];
            writecell(individual_matrix,[individual_node_path '\Nodes_Part_' file_identifiers{i} '.csv']);
        else
            writematrix(individual_matrix,[individual_node_path '\Nodes_Part_' file_identifiers{i} '.csv']);
        end
    end
    
    if headers_on == 1
        large_node_matrix = [headers; num2cell(large_node_matrix)];
        writecell(large_node_matrix,large_node_matrix_path);
    else
        writematrix(large_node_matrix,large_node_matrix_path);
    end
    
end

end

