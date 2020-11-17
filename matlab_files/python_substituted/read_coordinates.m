function [stress_part_values,complete_files_path,spherical_complete_files_path,figure_counter,large_node_matrix_path] ...
    = read_coordinates(sd_input,fname,stress_processing_path,deflection_paths,headers_on,...
    stress_matrices_paths,deflection_processing_path,coord_file,common_files_path)
%% Opening the file and initializing the variables 
% Function to read the .dat file containing all information about nodes and
% elements. The file will always be divided into different parts for which
% the node coordinates and the element vertices coordinates are displayed.
% We are only interested in the first part of the generated .dat file.
% Author: Fabrizio Morra

if nargin < 1
    error('Function requires one input argument');
elseif ~ischar(fname)
    error('Input argument must be a string representing a filename');
end

%% Processing the dat file
dat_processing_tic = tic;
[node_lines,elem_lines,g,file_identifiers,nset_lines] = dat_processor(fname);
dat_processing_toc = toc(dat_processing_tic);
disp(['The time spent to process the dat file is: ' ...
    num2str(dat_processing_toc) ' seconds.']);
%% Processing the nodes
nodes_processing_tic = tic;
[individual_node_path,large_node_matrix_path] = nodes_processor(common_files_path,...
    node_lines,elem_lines,g,file_identifiers,headers_on);
nodes_processing_toc = toc(nodes_processing_tic);
disp(['The time spent to process the elements is: ' ...
    num2str(nodes_processing_toc) ' seconds.']);
figure_counter = 0;
stress_part_values = [stress_processing_path '\Stress_Part_Values'];
% The next part of this function is based on whether we need to work with
% stresses or deflections. 
%% Processing the elements and associating stress components to centroid coordinates
if sd_input == 0
    elements_processing_tic = tic;
    individual_element_paths = elements_processor(node_lines,elem_lines,...
        nset_lines,stress_processing_path,stress_part_values,g,file_identifiers,headers_on);
    elements_processing_toc = toc(elements_processing_tic);
    disp(['The time spent to process the elements is: ' ...
        num2str(elements_processing_toc) ' seconds.']);
    % Create files with the stresses for each centroid and their
    % coordinates:
    % Now we can open the list with all the nodes and the one with all the
    % elements. These two variables will be matrices.
    
    complete_stress_files_tic = tic;
    [complete_files_path,spherical_complete_files_path] = ...
        associate_stress_coord(individual_element_paths,stress_part_values,...
        large_node_matrix_path,headers_on,stress_matrices_paths,coord_file);
    complete_stress_files_toc = toc(complete_stress_files_tic);
    disp(['The time spent to create the complete stress files is: ' ...
        num2str(complete_stress_files_toc) ' seconds.']);
    
%% Test scatter plots
    
    stress_scatter_test_s = readmatrix([spherical_complete_files_path '\Spherical_complete_file_' coord_file(1:5) '.csv']);
    figure(1)
    scatter3(stress_scatter_test_s(:,9),stress_scatter_test_s(:,10),stress_scatter_test_s(:,11),2,stress_scatter_test_s(:,2));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Scatter plot of all centroids with corresponding Mises cartesian component')
    c = colorbar;
    c.Label.String = 'Stress [Pa]';
    c.Location = 'southoutside';
    figure_counter = figure_counter + 1;
    
    stress_scatter_test_c = readmatrix([complete_files_path '\Complete_file_' coord_file(1:5) '.csv']);
    figure(2)
    scatter3(stress_scatter_test_c(:,9),stress_scatter_test_c(:,10),stress_scatter_test_c(:,11),2,stress_scatter_test_c(:,2));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Scatter plot of all centroids with corresponding Mises spherical component')
    c = colorbar;
    c.Label.String = 'Stress [Pa]';
    c.Location = 'southoutside';
    figure_counter = figure_counter + 1;

else
%% Associating deflections to their nodes with a shorter function
    complete_deflection_files_tic = tic;
    [complete_files_path,spherical_complete_files_path] = ...
        associate_deflection_coord(deflection_paths,individual_node_path,headers_on,deflection_processing_path);
    complete_deflection_files_toc = toc(complete_deflection_files_tic);
    disp(['The time spent to create the complete deflection files is: ' ...
        num2str(complete_deflection_files_toc) ' seconds.']);
    
    deflection_scatter_test_s = readmatrix([spherical_complete_files_path '\Large_complete_spherical_deflection_file.csv']);
    figure(1)
    scatter3(deflection_scatter_test_s(:,6),deflection_scatter_test_s(:,7),deflection_scatter_test_s(:,8),2,deflection_scatter_test_s(:,2));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Scatter plot of all nodes with corresponding deflection magnitude')
    c = colorbar;
    c.Label.String = 'Deflections [m]';
    c.Location = 'southoutside';
    figure_counter = figure_counter + 1;
    deflection_scatter_test_c = readmatrix([complete_files_path '\Large_complete_deflection_file.csv']);
    figure(2)
    scatter3(deflection_scatter_test_c(:,6),deflection_scatter_test_c(:,7),deflection_scatter_test_c(:,8),2,deflection_scatter_test_c(:,2));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Scatter plot of all nodes with corresponding deflection magnitude')
    c = colorbar;
    c.Label.String = 'Deflections [m]';
    c.Location = 'southoutside';
    figure_counter = figure_counter + 1;
end

end
