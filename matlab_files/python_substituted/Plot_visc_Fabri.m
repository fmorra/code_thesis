% load('Stress_Amundsen4000gr1000H20')
% element_loc = importdata('Element Centroid list_cone_model');
clear all
close all

%% viscosity parameters
d           = 4000;     % Micrometer
water       = 500.000;     %1000;        % atoms H per 1e6 atoms Si
phi_melt    = 0.00;

% latlim = [-82 -62];
% lonlim = [-170 -40];

latlim = [-80 -70];
lonlim = [-140 -80];

%% GOCE+ data
% Read Temperatures
% b = load('LM3D-Temp-geo.xyz');
% % 
b = load('ant-blocks-10g-T-geo.xyzT');
disp(['Finished reading temperatures']);
lon = b(:,1);
lat = b(:,2);
depth = b(:,3);      % in m, max depth 400 km
T_in = b(:,4) ;     % in degree Celsius - I don't have them and need to modify this

r = 6371e3+depth;
base_r = 6371e3;
fileString = ['GOCE+_' int2str(d) 'gr_' int2str(water) 'H2O_v3'];

Goce_true=1;
%%
Deg2Rad = pi/180;
pixelsize=0.25;

y_in = r.*sind(lat);
x_in = r.*cosd(lat).*sind(lon);
z_in = r.*cosd(lat).*cosd(lon);

%Set area |left ASE|TG|PIG|normal ASE?
lon_min     = 240;%240;%252;%259;%minmin245
lon_max     = 265;%252;%256;%261;%maxmax262
depth_max3  = 670e3;
depth_max2  = 400e3;
depth_max   = 200e3;
depth_min   = 50e3;
lat_min     = -76;%-76;%-76;%-75.5%-75
lat_max     = -74;%-75;%-75;%-74.5%-74

y_max = (base_r-depth_max3).*sind(lat_max);
y_min = base_r.*sind(lat_min);
x_min = base_r.*cosd(lat_max).*sind(lon_max);
x_max = (base_r-depth_max3).*cosd(lat_min).*sind(lon_min);  
z_min = base_r.*cosd(lat_max).*cosd(lon_min);
z_max = (base_r-depth_max3).*cosd(lat_min).*cosd(lon_max);

% Read element coordinates
% a = load('D:\Abaqus_Earth_version\Hipparchos files\Element Centroid list_2deg.txt');
a = load('Element Centroid list_cone_model.txt');
disp(['Finished reading coordinates']);
elno = a(:,1);
x_out = a(:,2);     % in m
y_out = a(:,3);
z_out = a(:,4);
r_out = sqrt(x_out.^2+y_out.^2+z_out.^2);
lon_out = atan2(x_out,z_out)*180/pi;
lat_out = asin(y_out./r_out)*180/pi;
depth_out = (base_r-r_out);

nodeID_long_int = find((lon_out<lon_max & lon_out>lon_min) | (lon_out<-360+lon_max & lon_out>-360+lon_min));
nodeID_lat_int = find((lat_out<lat_max & lat_out>lat_min));
nodeID_depth_int1 = find(depth_out<depth_max & depth_out>depth_min);
nodeID_depth_int2 = find(depth_out<depth_max2 & depth_out>depth_max);
nodeID_depth_int3 = find(depth_out<depth_max3 & depth_out>depth_max2);

nodeID_long = zeros(length(elno),1);
nodeID_lat = zeros(length(elno),1);
nodeID_depth1 = zeros(length(elno),1);
nodeID_depth2 = zeros(length(elno),1);
nodeID_depth3 = zeros(length(elno),1);

nodeID_long(nodeID_long_int) = 1;
nodeID_lat(nodeID_lat_int) = 1;
nodeID_depth1(nodeID_depth_int1) = 1;
nodeID_depth2(nodeID_depth_int2) = 1;
nodeID_depth3(nodeID_depth_int3) = 1;

el_amundsen_selection1 = find(nodeID_long.*nodeID_lat.*nodeID_depth1);
el_amundsen_selection2 = find(nodeID_long.*nodeID_lat.*nodeID_depth2);
el_amundsen_selection3 = find(nodeID_long.*nodeID_lat.*nodeID_depth3);

j=1;
k=1;

% Interpolate
T_out = griddata(x_in,y_in,z_in,T_in,x_out,y_out,z_out,'linear');

%% load stress
load coast

% Use your stresses and elements here (only one iteration at a time)
mises_list = importdata('W81\Data_output_stress_upper_labels.dat');
mises = importdata('G45\Data_output_stress_upper.dat');
mises2 = importdata('G45\Data_output_stress_upper_present.dat');
Elements = length(mises);
Stress_mises = zeros(Elements,1);
Stress_mises2 = zeros(Elements,1);
for i = 1:Elements
    Stress_mises(mises_list(i)) = mises(i); 
    Stress_mises2(mises_list(i)) = mises2(i); 
end
%% Computation of viscosity
        % melt fraction (default: 0)
is_linear   = 0;
pressure = 0.0333e9*depth_out/1E3;
T_out_K = T_out+273.15;               % units Kelvin for flow law

% Compute creep parameters
if water>0
    [Bdisl,Bdiff] = flowlawolivinewet(T_out_K,d,phi_melt,pressure,water);

else
    [Bdisl,Bdiff] = flowlawolivine(T_out_K,d,phi_melt,pressure);
end
	Bdiff = Bdiff*1e-6*(3./2.);% 1e-6 so stress can be in Pa, 3/2 see equation 5 van der wal et al (GJI 2013)

if is_linear
    Bdisl = 0;

else
    Bdisl = Bdisl*(1e-6^3.5)*(3./2.);% 1e-6 so stress can be in Pa, 3/2 see equation 5 van der wal et al (GJI 2013)
end
viscosity = 1./(3*(Bdiff+Bdisl.*(Stress_mises).^2.5));
viscosity_present = 1./(3*(Bdiff+Bdisl.*(Stress_mises2).^2.5));

%% Plot

kl2 = find( abs(depth_out-70E3)<10E3 & y_out<-(6371e3-70E3)*cosd(30) );
kl3 = find( abs(depth_out-150E3)<10E3 & y_out<-(6371e3-150E3)*cosd(30) );
kl4 = find( abs(depth_out-230E3)<10E3 & y_out<-(6371e3-230E3)*cosd(30) );

[Long,Lat] = meshgrid(min(lon_out(kl2)):pixelsize:max(lon_out(kl2)),min(lat_out(kl2)):pixelsize:max(lat_out(kl2)));
Lat=flipud(Lat);

Viscosity_AngularGrid70 = griddata(lon_out(kl2),lat_out(kl2),log10(viscosity(kl2)),Long,Lat,'linear');
Viscosity_AngularGrid150 = griddata(lon_out(kl3),lat_out(kl3),log10(viscosity(kl3)),Long,Lat,'linear');
Viscosity_AngularGrid230 = griddata(lon_out(kl4),lat_out(kl4),log10(viscosity(kl4)),Long,Lat,'linear');

Viscosity_AngularGrid70_present = griddata(lon_out(kl2),lat_out(kl2),log10(viscosity_present(kl2)),Long,Lat,'linear');
Viscosity_AngularGrid150_present = griddata(lon_out(kl3),lat_out(kl3),log10(viscosity_present(kl3)),Long,Lat,'linear');
Viscosity_AngularGrid230_present = griddata(lon_out(kl4),lat_out(kl4),log10(viscosity_present(kl4)),Long,Lat,'linear');

[X,Y] = meshgrid(min(x_out(kl2)):max(x_out(kl2))/50:(max(x_out(kl2))-max(x_out(kl2))*0.001),min(z_out(kl2)):max(z_out(kl2))/50:(max(z_out(kl2))-max(z_out(kl2))*0.001));
Viscosity_regularGrid120 = log10(griddata(x_out(kl2),z_out(kl2),viscosity(kl2),X,Y,'linear'));
Viscosity_regularGrid120_present = log10(griddata(x_out(kl2),z_out(kl2),viscosity_present(kl2),X,Y,'linear'));
[X,Y] = meshgrid(min(x_out(kl3)):max(x_out(kl3))/50:(max(x_out(kl3))-max(x_out(kl3))*0.001),min(z_out(kl3)):max(z_out(kl3))/50:(max(z_out(kl3))-max(z_out(kl3))*0.001));
Viscosity_regularGrid200 = log10(griddata(x_out(kl3),z_out(kl3),viscosity(kl3),X,Y,'linear'));
Viscosity_regularGrid200_present = log10(griddata(x_out(kl3),z_out(kl3),viscosity(kl3),X,Y,'linear'));

%% For checking, make plot at depth of 'showdepth' km


%lonlim = [-178 178];

Start_lat = Lat(1,1);
Start_long = Long(1,1);
truncate_no_lat1 = floor((Start_lat-latlim(2))*(1/pixelsize))+1;
truncate_no_lat2 = floor((Start_lat-latlim(1))*(1/pixelsize))+1;

if latlim(1)==-90
    truncate_no_lat2 = length(T_kelvin_angularGrid120(:,1));
end

truncate_no_lon1 = floor((lonlim(1) - Start_long)*(1/pixelsize))+1;
truncate_no_lon2 = floor((lonlim(2) - Start_long)*(1/pixelsize))+1;

if latlim(1)==-180
    truncate_no_lat1 = 1;
end

if latlim(2)==180
    truncate_no_lat2 = length(T_kelvin_angularGrid120(1,:));
end
LatVec = latlim(2):-pixelsize:latlim(1);
LonVec = lonlim(1):pixelsize:lonlim(2);

[lonMat, latMat] = meshgrid(LonVec, LatVec);
lon = reshape(lonMat, [numel(lonMat),1]);
lati = reshape(latMat, [numel(latMat),1]);

%%

Viscosity = reshape(-Viscosity_AngularGrid70(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(21);
colormap(hot);
caxis([-21.0 -18.5]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 70km depth for ASE for 1951');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
load coast;
coastlines= importdata('coastxy.dat');
[theta, rho]=cart2pol(coastlines(:,2),coastlines(:,1));
long=theta/pi*180;
lat = -90+rho/(pi*6371e3/180);
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [log_{10}(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')
%%

Viscosity = reshape(-Viscosity_AngularGrid150(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(22);
colormap(hot);
caxis([-19.0 -18.0]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 150km depth for ASE for 1951');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [10log(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')
%%

Viscosity = reshape(-Viscosity_AngularGrid230(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(23);
colormap(hot);
caxis([-19.0 -18.0]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 230km depth for ASE for 1951');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [10log(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')
%%

Viscosity = reshape(-Viscosity_AngularGrid70_present(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(24);
colormap(hot);
caxis([-21.0 -18.5]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 70km depth for ASE for 2014');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [10log(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')
%%

Viscosity = reshape(-Viscosity_AngularGrid150_present(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(25);
colormap(hot);
caxis([-19.0 -18.0]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 150km depth for ASE for 2014');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [10log(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')
%%

Viscosity = reshape(-Viscosity_AngularGrid230_present(truncate_no_lat1:truncate_no_lat2,truncate_no_lon1:truncate_no_lon2),[numel(lonMat),1]);
[Z, refvec] = geoloc2grid(lati,lon,Viscosity,pixelsize);
refvec(2) = Start_lat - (truncate_no_lat1-1)*pixelsize;
refvec(3) = Start_long + (truncate_no_lon1-1)*pixelsize;
f = figure(26);
colormap(hot);
caxis([-19.0 -18.0]);
h = colorbar('v'); %set(h, 'ylim', [18.4 18.8]);
hcb = worldmap(Z,refvec);
ttl = title('Viscosity at 230km depth for ASE for 2014');
ttl.Position = [3.816987465135753 -6648474.767631762*0.92 1.4551915228366852E-11];
geoshow(Z, refvec, 'DisplayType', 'texture');
plotm(lat,long,'b','LineWidth',2)
set(get(h,'ylabel'),'string','Viscosity [10log(Pa s)]');
setm(hcb,'mlabellocation',20)%title('Viscosity at 100km depth for the WinterC model in West-Antarctica')

