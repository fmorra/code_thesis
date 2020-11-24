% Script to compare results from my run and Bas'
clear all; close all; clc;
extension = '.dat';
my_files_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\bas_comparison_files';
bas_files_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\ABAQUS_3D_GIA_model\Input_output_standard_run\dat_files';
dat_names = ...
    {strcat('Deflection0',extension),strcat('Deflection9',extension),strcat('Deflection10',extension),...
    strcat('dx_surface0',extension),strcat('dx_surface9',extension),strcat('dx_surface10',extension),...
    strcat('dy_surface0',extension),strcat('dy_surface9',extension),strcat('dy_surface10',extension),...
    strcat('dz_surface0',extension),strcat('dz_surface9',extension),strcat('dz_surface10',extension)};

for i = 1:length(dat_names)
    my_file = importdata([my_files_path '\' dat_names{i}]);
    bas_file = importdata([bas_files_path '\' dat_names{i}]);
    lat_vector = linspace(-90,90,size(my_file,1));
    lon_vector = linspace(0,360,size(my_file,2));
    lon_vector_plot = flip(lon_vector - 180);
    lat_vector_plot = flip(lat_vector);
    diff_matrix = my_file - bas_file;
    diff_matrix_percent = (diff_matrix./bas_file)*100;
    for  j = 1:size(diff_matrix_percent,1)
        for k = 1:size(diff_matrix_percent,2)
            if abs(diff_matrix_percent(j,k)) > 1
                diff_matrix_percent(j,k) = 0;
            end
        end
    end
    
    figure(i)
    surf(lon_vector_plot,lat_vector_plot,diff_matrix_percent, 'edgecolor', 'none')
    axis([lon_vector_plot(end) lon_vector_plot(1) lat_vector_plot(end) lat_vector_plot(1)])
    set(gca,'Xtick',lon_vector_plot(end):20:lon_vector_plot(1))
    set(gca,'Ytick',lat_vector_plot(end):10:lat_vector_plot(1))
    set(gca, 'Xdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
    xlabel('Longitude (deg)'); ylabel('Latitude(deg)');
    title('Percentage difference');
    view(180,90);
    h = colorbar('v');
    set(get(h,'ylabel'),'string','%')
    grid on;
    saveas(gcf,['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\bas_comparison_files\' dat_names{i}(1:end-4) '.png']);
    
end
