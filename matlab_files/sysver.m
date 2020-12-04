%file_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_21\Iteration_1\step_0\cycle_1_reports\stress_processing_results\Stress_Part_Values\Complete_Files\Complete_file_EARTH';
file_path = 'C:\Users\fabri\Desktop\TU Delft\Thesis\ABAQUS\test_run_python\run_21\Iteration_1\step_0\cycle_1_reports\deflection_processing_results\Complete_deflection_files\Complete_file_EARTH';
file = readmatrix(file_path);
% Back to original
file(:,[7 8]) = file(:,[8 7]);
file(:,[6 8]) = file(:,[8 6]);

% file(:,[9 11]) = file(:,[11 9]);
% file(:,[10 11]) = file(:,[11 10]);
% % Current, correct transformation
file(:,[6 8]) = file(:,[8 6]);
file(:,[6 7]) = file(:,[7 6]);
x = file(:,6);
y = file(:,7);
z = file(:,8);
magnitude = file(:,2);

figure(1)
scatter3(x, y, z, 20, magnitude);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
h = colorbar('v'); 
set(get(h,'ylabel'),'string','Deformations (m)');