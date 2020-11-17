function [lon,lat] = cart2latlon(cart_coord)
% Code to convert cartesian coordinates to latitude and longitude, based on
% the work of Minjae Yoo:
% https://www.mathworks.com/matlabcentral/fileexchange/69514-ecef-x-y-z-to-longitude-and-latitude

X = cart_coord(1);
Y = cart_coord(2);
Z = cart_coord(3);

% Estimate the Earth's radius of curvature
a = 6378137;            % [m], semi-major axis
f = 1/298.257223563;    % Ellipsoid flattening
b = a*(1-f);            % Semi-minor axis

% Estimate auxiliary values
P = sqrt(X^2 + Y^2); 
e = sqrt(((a^2) - (b^2))/a^2); % eccentricity of The Earth

% Initial latitude value
initial_latitude = atan2(Z,(P*(1-e^2)));

% Iteration loop to estimate the latitude
for i = 1:10000
    N = a / sqrt(1-e^2*sin(initial_latitude).^2);        % Prime vertical
    h = (P/cos(initial_latitude)) - N;                   % Height
    initial_latitude = atan2(Z,(P*(1-e^2*(N/(N+h)))));   % Latitude
end

% Convert to degrees
lon = (atan2(Y,X))*180/pi;
lat = initial_latitude*180/pi;

end

