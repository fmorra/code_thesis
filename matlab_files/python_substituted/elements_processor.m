function [individual_element_paths] = elements_processor(node_lines,...
    elem_lines,nset_lines,coord_paths,individual_path,g,file_identifiers,headers_on)
% This function processes part of the large dat file produced by ABAQUS nd
% modified by MATLAB to save the elements of each part and the nodes that
% make them up.

% Initialize the matrices containing all the elements 

large_elem_matrix_unsorted = [];
large_elem_matrix_unsorted_reset = [];

% Each element set is delimited by a line containing *Element as upper
% boundary and *Nset as lower boundary. 

new_elem_lines = elem_lines(1:length(node_lines)-1);
new_nset_lines = zeros(length(node_lines)-1,1);
for i=1:length(new_elem_lines)
    for j=i:length(nset_lines)
        if nset_lines(j) > new_elem_lines(i)
            new_nset_lines(i) = nset_lines(j);
            break;
        end
    end
end

% Initialize variables to save the element with the nodes making them up,
% then initialize the paths and folders we need

large_elem_matrix_path = [coord_paths '\' 'Large_Elem_Matrix' '.csv'];
large_elem_matrix_unsorted_path = [coord_paths '\' 'Large_Elem_Matrix_Unsorted' '.csv'];
large_elem_matrix_unsorted_path_reset = [coord_paths '\' 'Large_Elem_Matrix_Unsorted_Reset' '.csv'];
individual_element_paths = [individual_path '\Elements'];

if ~exist(individual_element_paths, 'dir')
    mkdir (individual_element_paths)
end

headers = {'Label', 'Node_1', 'Node_2', 'Node_3', 'Node_4', 'Node_5', 'Node_6', ...
    'Node_7', 'Node_8'};

% Verify if the files are already all there based on whether the last file
% to be created is actually already in the folder

if isfile(large_elem_matrix_unsorted_path_reset)
    disp(['The files containing the element nodes already exist, moving on' ...
        ' to centroid calculation and stress association.'])
else
    % Run the function if the files are not there
    disp('Processing elements to extract their nodes...')
    % Read the dat file and process the lines containing element
    % information for each part
    for k = 1:length(new_elem_lines)
        disp(['Writing file containing element nodes for part ' file_identifiers{k}])
        % Loop between the two lines containing delimiting keywords for the
        % elements and scan every line to save float values if found. This part
        % works in a similar fashion as the function scanning for the node
        % coordinates.
        %elcounter = 0;
        elem_vector = zeros(new_nset_lines(k) - new_elem_lines(k),9);
        elem_vector_reset = zeros(new_nset_lines(k) - new_elem_lines(k),9);
        
        for j = new_elem_lines(k):1:new_nset_lines(k)
            s=char(g{1}(j));
            data= sscanf(s, '%f');
            if ~isempty(data)
                %elcounter = elcounter + 1;
                data = data';
                
                % Every 5 lines, there will be an additional element which is
                % the number right after LINE, in which we are not interested,
                % so we must take it away. We need to diffeentiate based on the
                % number of nodes in an element: if we have 8 nodes, the
                % numbers of that line will be 10 because we have this first
                % number e must get rid of, then the element label and then the
                % node labels. Same reasonig if we have 6 nodes, which brings
                % the total to 8.
                
                if length(data) == 10
                    data(1) = [];
                elseif length(data) == 8
                    data(1) = [];
                end
                
                % Save two separate element vectors, one where the label count
                % does not reset for each part and one where it does. Also,
                % differentiate based on how many nodes make up the element (8
                % or 6).
                
                if length(data) == 9
                    elem_vector(j-new_elem_lines(k)+1,:) = data;
                    elem_vector_reset(j-new_elem_lines(k)+1,:) = data;
                    %elem_vector_reset(j-new_elem_lines(k)+1,1) = elcounter;
                else
                    elem_vector(j-new_elem_lines(k)+1,1:7) = data;
                    elem_vector_reset(j-new_elem_lines(k)+1,1:7) = data;
                    %elem_vector_reset(j-new_elem_lines(k)+1,1) = elcounter;
                end
                if elem_vector(j-new_elem_lines(k)+1,2) == elem_vector(j-new_elem_lines(k)+1,3)
                    elem_vector(j-new_elem_lines(k)+1,:) = 0;
                end
                if elem_vector_reset(j-new_elem_lines(k)+1,2) == elem_vector_reset(j-new_elem_lines(k)+1,3)
                    elem_vector_reset(j-new_elem_lines(k)+1,:) = 0;
                end
            end
            
        end
        
        elem_vector(~any(elem_vector,2), : ) = [];
        elem_vector_reset(~any(elem_vector_reset,2), : ) = [];
        reset_labels = (1:1:length(elem_vector_reset))';
        elem_vector_reset(:,1) = reset_labels;
        elem_vector(:,1) = reset_labels;
        % Save the element matrices in the same way it was done for the nodes:
        % we have individual matrices, one for each part, and then larger
        % matrices both in sorted and unsorted order. We also have a large
        % matrix where the element label resets for each part.
        
        % Define each individual matrix
        
        individual_elem_matrix = elem_vector_reset;
        
        % Build the larger matrices, a sorte one, an unsorted one, and an
        % unsorted one where the element count resets for each part
        
        large_elem_matrix_unsorted_reset = [large_elem_matrix_unsorted_reset; individual_elem_matrix];
        large_elem_matrix_unsorted = [large_elem_matrix_unsorted; elem_vector];
        large_elem_matrix_unsorted(~any(large_elem_matrix_unsorted,2), : ) = [];
        large_elem_matrix = sortrows(large_elem_matrix_unsorted);
        
        if headers_on == 1
            individual_elem_matrix = [headers; num2cell(individual_elem_matrix)];
            individual_elem_matrix(cellfun(@(x) ~x(1),individual_elem_matrix(:,1)),:) = [];
            writecell(individual_elem_matrix,[individual_element_paths '\Elements_Part_' file_identifiers{k} '.csv']);
        else
            writematrix(individual_elem_matrix,[individual_element_paths '\Elements_Part_' file_identifiers{k} '.csv']);
        end
        
    end
    
    % Matrices are saved here, differentiating based on whether we want
    % headers or not
    % Attach headers to the larger matrices before saving them
    
    if headers_on == 1
        large_elem_matrix = [headers; num2cell(large_elem_matrix)];
        large_elem_matrix_unsorted = [headers; num2cell(large_elem_matrix_unsorted)];
        large_elem_matrix_unsorted_reset = [headers; num2cell(large_elem_matrix_unsorted_reset)];
        writecell(large_elem_matrix,large_elem_matrix_path);
        writecell(large_elem_matrix_unsorted,large_elem_matrix_unsorted_path);
        writecell(large_elem_matrix_unsorted_reset,large_elem_matrix_unsorted_path_reset);
    else
        writematrix(large_elem_matrix,large_elem_matrix_path);
        writematrix(large_elem_matrix_unsorted,large_elem_matrix_unsorted_path);
        writematrix(large_elem_matrix_unsorted_reset,large_elem_matrix_unsorted_path_reset);
    end
    
end

end

