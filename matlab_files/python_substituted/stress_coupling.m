function [large_new_matrix] = stress_coupling(stress_matrices_paths, ...
    coupled_stress_folder, headers_on, new_stress_test, stress_headers)

dir_csv_stress_files = dir([stress_matrices_paths '\*.csv']);
csv_stress_files = {dir_csv_stress_files.name};
large_new_matrix = [];
large_new_matrix_path = [stress_matrices_path, '\Large_Stress_Matrix.csv'];

if exist(coupled_stress_folder)
    
else
    selected_component = input(['Enter the desired stress component(s) to modify,' ...
        ' possible values are Mises, S11, S22, S33, S12, S13, S23:\n']);
    selected_components = split(selected_component);
    selected_columns = zeros(length(selected_components),1);
    for k = 1:length(selected_components)
        if strcmp(selected_components{k}, 'Mises') == 1
            selected_columns(k) = 2;
        elseif strcmp(selected_components{k}, 'S11') == 1
            selected_columns(k) = 3;
        elseif strcmp(selected_components{k}, 'S22') == 1
            selected_columns(k) = 4;
        elseif strcmp(selected_components{k}, 'S33') == 1
            selected_columns(k) = 5;
        elseif strcmp(selected_components{k}, 'S12') == 1
            selected_columns(k) = 6;
        elseif strcmp(selected_components{k}, 'S13') == 1
            selected_columns(k) = 7;
        else
            selected_columns(k) = 8;
        end
    end
    
    for i = 1:length(csv_stress_files)
        size_to_evaluate = dir([dir_csv_stress_files(i).folder '\' csv_stress_files{i}]);
        if size_to_evaluate.bytes < 100
            csv_stress_files{i} = {};
        end
    end
    
    csv_stress_files = csv_stress_files(~cellfun(@isempty,csv_stress_files));
    
    for i = 1:length(csv_stress_files)
        
        %     equation_test = readmatrix([stress_matrices_paths '\' csv_stress_files{i}]);
        %     mises_stress = (1/sqrt(2))*sqrt((equation_test(:,3)-equation_test(:,4)).^2+...
        %         (equation_test(:,4)-equation_test(:,5)).^2 + (equation_test(:,5)-equation_test(:,3)).^2+...
        %         6.*(equation_test(:,6).^2+equation_test(:,8).^2+equation_test(:,7).^2))
        
        stresses_to_modify = readmatrix([stress_matrices_paths '\' csv_stress_files{i}]);
        new_stress_length = ones(length(stresses_to_modify),1);
        for j = 1:length(selected_columns)
            stresses_to_modify(:,selected_columns(j)) = stresses_to_modify(:,selected_columns(j))...
                + new_stress_length*new_stress_test;
        end
        new_matrix = stresses_to_modify;
        
        % Recalculate the Mises stress
        mises_stress = (1/sqrt(2))*sqrt((new_matrix(:,3)-new_matrix(:,4)).^2+...
            (new_matrix(:,4)-new_matrix(:,5)).^2 + (new_matrix(:,5)-new_matrix(:,3)).^2+...
            6.*(new_matrix(:,6).^2+new_matrix(:,8)+new_matrix(:,7).^2));
        new_matrix(:,2) = mises_stress;
        large_new_matrix = [large_new_matrix; new_matrix];
        if headers_on == 1
            new_matrix = [stress_headers; num2cell(new_matrix)];
            writecell(new_matrix,[coupled_stress_folder '\' csv_stress_files{i}]);
        else
            writematrix(new_matrix,[coupled_stress_folder '\' csv_stress_files{i}]);
        end
    end
    large_new_matrix = sortrows(large_new_matrix);
    if headers_on == 1
        large_new_matrix = [stress_headers; num2cell(large_new_matrix)];
        writecell(large_new_matrix, large_new_matrix_path);
    else
        writematrix(large_new_matrix, large_new_matrix_path);
    end
    new_e_column = large_new_matrix(1:end, 1);
    old_e_file_path = [base_path '\Iteration_' num2str(iteration - 1) 'e.csv'];
    e_file = readmatrix(old_e_file_path);
    e_file(1:end, end) = new_e_column;
    new_e_file_path = [base_path '\Iteration_' num2str(iteration) 'e.csv'];
    writematrix(e_file, new_e_file_path);
end
end

