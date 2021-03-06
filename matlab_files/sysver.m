%file_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_21\Iteration_1\step_0\cycle_1_reports\stress_processing_results\Stress_Part_Values\Complete_Files\Complete_file_EARTH';
file_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_21\Iteration_1\step_0\cycle_1_reports\deflection_processing_results\Complete_deflection_files\Complete_file_EARTH';
file = readmatrix(file_path);
% Back to original
file(:,[7 8]) = file(:,[8 7]);
file(:,[6 8]) = file(:,[8 6]);
 
% file(:,[6 8]) = file(:,[8 6]);
% file(:,[6 7]) = file(:,[7 6]);
% % Current, correct transformation
file(:,[6 8]) = file(:,[8 6]);
file(:,[7 8]) = file(:,[8 7]);
x = file(:,6);
y = file(:,7);
z = file(:,8);
[az,el,r] = cart2sph(x,y,z);
geo_mat = [rad2deg(az),rad2deg(el),r];
r_condition = r<=6.371e6 & r>=6.366e6;
lonlim = [-135 -60];
latlim = [-80 -65];
lat_condition = geo_mat(:,2)>=latlim(1) & geo_mat(:,2)<=latlim(2);
lon_condition = geo_mat(:,1)>=lonlim(1) & geo_mat(:,1)<=lonlim(2);
area_condition = r_condition & lat_condition & lon_condition;
lon_plot = geo_mat(area_condition,1);
lat_plot = geo_mat(area_condition,2);
% lon_plot = geo_mat(area_indices,1);
% lat_plot = geo_mat(area_indices,2);

[lon,lat] = meshgrid(linspace(min(lon_plot),max(lon_plot),400),...
    linspace(min(lat_plot),max(lat_plot),400));
magnitude_plot = file(area_condition,2);
magnitude = griddata(lat_plot,lon_plot,magnitude_plot,lat,lon,'linear');
scattermagn = file(:,2);

figure(1)
worldmap(latlim,lonlim);
load coastlines
plotm(coastlat, coastlon);
surfm(lat,lon,magnitude);
alpha 0.6;
h = colorbar('v'); 
caxis([round(min(magnitude_plot)) round(max(magnitude_plot))])
set(get(h,'ylabel'),'string','Deformations (m)');

figure(2)
scatter3(x, y, z, 20, scattermagn);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
h = colorbar('v'); 
set(get(h,'ylabel'),'string','Deformations (m)');