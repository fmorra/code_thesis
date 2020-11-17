function [spherical_complete_matrix] = rotate_deflections(complete_deflection_matrix)
% This function changes the reference system for the three deflection
% components from cartesian to spherical.

displacements = complete_deflection_matrix(:,3:5);
node_coords = complete_deflection_matrix(:,6:8);
spherical_displacements = zeros(size(displacements,1),size(displacements,2));
U_magn = zeros(size(displacements,1),1);
R_3D = sqrt(node_coords(:,1).^2+node_coords(:,2).^2+node_coords(:,3).^2);
R_azimuth = sqrt(node_coords(:,1).^2+node_coords(:,2).^2);
cos_theta = node_coords(:,1)./R_azimuth;
sin_theta = node_coords(:,2)./R_azimuth;
cos_phi = R_azimuth./R_3D;
sin_phi = node_coords(:,3)./R_3D;

for i = 1:length(displacements)
    T(1,1) = cos_theta(i)*cos_phi(i);
    T(2,1) = cos_theta(i)*sin_phi(i);
    T(3,1) = -sin_theta(i);
    T(1,2) = sin_theta(i)*cos_phi(i);
    T(2,2) = sin_theta(i)*sin_phi(i);
    T(3,2) = cos_theta(i);
    T(1,3) = sin_phi(i);
    T(2,3) = -cos_phi(i);
    T(3,3) = 0.0;
    spherical_displacements(i,:) = (T*displacements(i,:)')';
    U_magn(i) = sqrt(spherical_displacements(i)^2+spherical_displacements(i)^2+spherical_displacements(i)^2);
end

complete_spherical_displacements = [U_magn,spherical_displacements];
spherical_complete_matrix = [complete_deflection_matrix(:,1),complete_spherical_displacements,complete_deflection_matrix(:,end-2:end)];

end

