function [lat,lon] = cart2geo(cartesian_coordinates)
% This function uses a matrix containing X,Y,Z on the first, second and
% third columns respectively and returns two columns containing the
% corresponding latitude and longitude.

X = cartesian_coordinates(:,1);
Y = cartesian_coordinates(:,2);
Z = cartesian_coordinates(:,3);
R = sqrt(X.^2+Y.^2+Z.^2);

lat = asin(Z./R)*180/pi;
lon = atan2(Y,X)*180/pi;

end

