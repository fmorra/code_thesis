function [complete_file_spherical] = rotate_tensor(part_matrix,spherical_centroids_files_path, ...
    spherical_components_files_path,spherical_complete_files_path,headers_on,...
    complete_headers,centroid_headers,complete_individual_path,part_file)

% Coordinate transformation of the stress tensor - geographic coordinates
% Read the line containing the 6 stress tensor elements and extract the
% corresponding 3*3 matrix.
    
spherical_components = zeros(size(part_matrix,1),7);
complete_file = readmatrix(complete_individual_path);
coords = complete_file(:,9:11);
stress_headers = {'S_Mises','S_S11','S_S22','S_S33','S_S12','S_S13','S_S23'};

for l=1:size(part_matrix,1)
    
    % Here we fill the stress tensor matrix with the components coming from
    % the values read in the report file.
    
    % We also need to swap elements of the stress tensor - this should only happen
    % after comparing to ABAQUS results because we have to keep the
    % (1,2,3) direction -> (z,x,y) that ABAQUS uses for this
    % comparison. Swapping in MATLAB is performed used linear indexing,
    % that is, the matrix indices are as follows:
    % [1 4 7]
    % [2 5 8]
    % [3 6 9]
    
    % Half of the matrix will be filled first and then the other
    % symmetrical values will be assigned after performing the element swap
    % due to the change of axes from the ABAQUS model.
    
    S(1,1) = complete_file(l,3);
    S(2,1) = complete_file(l,6);
    S(3,1) = complete_file(l,7);
    S(1,2) = S(2,1);
    S(2,2) = complete_file(l,4);
    S(3,2) = complete_file(l,8);
    S(3,3) = complete_file(l,5);
    S(1,3) = S(3,1);
    S(2,3) = S(3,2);
%     % Swap S(1,1) and S(3,3)
%     S([1 9]) = S([9 1]);
%     % Swap S(1,1) and S(2, 2)
%     S([1 5]) = S([5 1]);
%     S(3,1) = complete_file(l,7);
%     S(3,2) = complete_file(l,8);
%     % Swap S(2,1) and S(3,1)
%     S([2 3]) = S([3 2]);
%     % Swap S(2,1) and S(3,2)
%     S([2 6]) = S([6 2]);

    
    % Define the values used in the transformation matrix
    
    centroid_coords = complete_file(l,9:end);
    R_3D = sqrt(centroid_coords(1)^2+centroid_coords(2)^2+centroid_coords(3)^2);
    R_azimuth = sqrt(centroid_coords(1)^2+centroid_coords(2)^2);
    cos_theta = centroid_coords(1)/R_azimuth;
    sin_theta = centroid_coords(2)/R_azimuth;
    cos_phi = R_azimuth/R_3D;
    sin_phi = centroid_coords(3)/R_3D; 
    
    % Define the elements of the transformation matrix based on previous
    % values; these are computed based on te spherical trasnsformation of
    % the stress tensor by Allan F. Brower, after which trigonometric
    % properties are used to derive the latitude and lognitude after the
    % transformation instead of simple spherical angles. The longitude
    % stays the same, while the latitude is equal to 90 minus the secondary
    % angle of the transformation.
    
    T(1,1) = cos_theta*cos_phi;
    T(2,1) = cos_theta*sin_phi;
    T(3,1) = -sin_theta;
    T(1,2) = sin_theta*cos_phi;
    T(2,2) = sin_theta*sin_phi;
    T(3,2) = cos_theta;
    T(1,3) = sin_phi;
    T(2,3) = -cos_phi;
    T(3,3) = 0.0;
    S_spherical = T*S*T';
    
    % Save mises stress which is an invariant without any transformation
    % and then use the newly calculated spherical values to store them
    
    spherical_components(l,1) = part_matrix(l,2);
    spherical_components(l,2) = S_spherical(1,1);
    spherical_components(l,3) = S_spherical(2,2);
    spherical_components(l,4) = S_spherical(3,3);
    spherical_components(l,5) = S_spherical(1,2);
    spherical_components(l,6) = S_spherical(1,3);
    spherical_components(l,7) = S_spherical(2,3);
    
    % The coordinates are not spherical, that will happen later after
    % discretizing the stress values for differnt depths.
    
end

% Create files to save and save them

complete_file_spherical = [complete_file(:,1),spherical_components,coords];
centroid_coords = [part_matrix(:,1), coords];

if headers_on == 1
    coords = [centroid_headers; num2cell(centroid_coords)];
    spherical_components = [stress_headers; num2cell(spherical_components)];
    complete_file_spherical_headers = [complete_headers; num2cell(complete_file_spherical)];
    writecell(coords,[spherical_centroids_files_path '\Spherical_Coordinates_' part_file]);
    writecell(spherical_components,[spherical_components_files_path '\Spherical_Components_' part_file]);
    writecell(complete_file_spherical_headers,[spherical_complete_files_path '\Spherical_Complete_file_' part_file]);
else
    writematrix(coords,[spherical_centroids_files_path '\Spherical_Coordinates_' part_file]);
    writematrix(spherical_components,[spherical_components_files_path '\Spherical_Components_' part_file]);
    writematrix(complete_file_spherical,[spherical_complete_files_path '\Spherical_Complete_file_' part_file]);
end

end

