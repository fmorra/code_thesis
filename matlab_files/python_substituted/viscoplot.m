% This script will be used to plot the viscosity calculated using the same
% equation as in the fortran file. The data will be read from the dat file.
% Select an iteration to read the viscosity data
iteration = 1;
run = 6;
% e_path = ['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\stress_output_files\run_' ...
%     num2str(run) '\Iteration_' num2str(iteration)];
e_path = ['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_' ...
    num2str(run) '\Iteration_' num2str(iteration) '\e.dat'];
stress_path = ['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_' ...
    num2str(run) '\Iteration_' num2str(iteration) '\stress_matrices\coupled_stress_matrices\large_stress_matrix.csv'];
stress_path_coupled = ['C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_' ...
    num2str(run) '\Iteration_' num2str(iteration) '\stress_matrices\large_stress_matrix.csv'];
e = readmatrix(e_path);
large_stress_matrix = readmatrix(stress_path);
large_stress_matrix_coupled = readmatrix(stress_path_coupled);
elements = e(:,1);
a = e(:,2);
alin = e(:,3);
qtild = e(:,4);
an = 3.5;
strain_rate_coupled = zeros(length(e),1);
strain_rate = zeros(length(e),1);
stress = zeros(length(e),1);
stress_coupled = zeros(length(e),1);
coupled_mises = large_stress_matrix_coupled(:,2)*1e6;
mises = large_stress_matrix(:,2)*1e6;

for i = 1:length(strain_rate_coupled)
    strain_rate_coupled(i,1) = alin(i)*coupled_mises(i) + a(i)*coupled_mises(i)^an;
    strain_rate(i,1) = alin(i)*mises(i) + a(i)*mises(i)^an;
end

viscosity_coupled = coupled_mises./(3*strain_rate_coupled);
viscosity = mises./(3*strain_rate);

%% Plots
figure(1)
semilogy(elements, viscosity);
grid on;
figure(2)
semilogy(elements, viscosity_coupled);
grid on;
figure(3)
semilogy(elements, viscosity, elements, viscosity_coupled);
grid on;