function [] = visco_diff_plotter(visco_diff_path)

visco_strain_input = 0;
if visco_strain_input == 0
    quantity = 'viscosity';
else
    quantity = 'strain_rate';
end
dir_diff_files = dir([visco_diff_path '\*.csv']);
diff_plots_folder = [visco_diff_path '\' quantity '_plots'];
if ~exist(diff_plots_folder, 'dir')
    mkdir(diff_plots_folder)
end
diff_files = {dir_diff_files.name};

end

