function [deflection_paths] = read_deflections_report(fname,search_key,deflection_paths,...
    headers_on)
% This function reads deflections from the relative report file the same
% way as was done with the stress components.

%% Read the deflections from the report
% Check number and type of arguments

if nargin < 1
    error('Function requires one input argument');
elseif ~ischar(fname)
    error('Input argument must be a string representing a filename');
end

% Open file

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end

% Initialize the same counters as for the stress file

line=0;
tline=1;
lines = [];
keyword_counter = 0;
file_identifiers = {};

while ~feof(fid)
    % Read a line without including the newline
    s=fgets(fid);
    % Determine where the next part or region starts and store the file
    % line, also storing this line after manipulating it to use the part or
    % region name for the names of the csv files, as for the stress file.
    keyword_index = contains(s,search_key);
    line = line + 1;
    if (keyword_index == 1)
        keyword_counter = keyword_counter + 1;
        lines(keyword_counter,1) = line;
        file_identifiers{keyword_counter} = replace(strip(extractAfter(s,": ")),".","_");
    end
    % Different headers name
    headers_key = 'U.Magnitude';
    headers_index = contains(s,headers_key);
    if (headers_index == 1)
        headers = replace(strip(extractAfter(s,"Node ")),".","_");
    end
end
headers = ['Label';split(headers)];
headers = reshape(headers,[1 5]);

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end
g = textscan(fid,'%s','delimiter','\n');

% As in the stress file, check for the existence of the last file to be
% generated.
if isfile([deflection_paths '\I1.csv'])
    disp(['The files containing the deflection components already exist, moving on' ...
        ' to nodes coordinate extraction for each part.'])
else
    disp('Processing rpt to extract deflection components...')
    for i = 1:size(lines,1)
        if (i<size(lines,1))
            for j = lines(i):1:lines(i+1)
                % Access line j, scan it for float values and if they are found
                % transpose them and save them in the row of a matrix
                s=char(g{1}(j));
                [data]= sscanf(s, '%f');
                if ~isempty(data)
                    data = data';
                    deflection_submatrix(j-lines(i)+1,:) = data;
                end
            end
            % Delete all zero rows, save the complete stress component matrix
            % in a csv file and re-initialize it
            deflection_submatrix( ~any(deflection_submatrix,2), : ) = [];
            
            if (headers_on == 1)
                deflection_submatrix = [headers; num2cell(deflection_submatrix)];
                writecell(deflection_submatrix,[deflection_paths '\' file_identifiers{i} '.csv']);
                deflection_submatrix = [];
            else
                writematrix(deflection_submatrix,[deflection_paths '\' file_identifiers{i} '.csv']);
                deflection_submatrix = [];
            end
            
        else
            % Same procedure as above, but using a while cycle
            while ~feof(fid)
                s=char(g{1}(tline));
                [data, ncols, errmsg, nxtindex]= sscanf(s, '%f');
                if ~isempty(data)
                    data = data';
                    deflection_submatrix(tline,:) = data;
                    tline=tline+1;
                else
                    tline=tline+1;
                end
            end
            
            if (headers_on == 1)
                deflection_submatrix = [headers; num2cell(deflection_submatrix)];
                writecell(deflection_submatrix,[deflection_paths '\' file_identifiers{i} '.csv']);
                deflection_submatrix = [];
            else
                writematrix(deflection_submatrix,[deflection_paths '\' file_identifiers{i} '.csv']);
                deflection_submatrix = [];
            end
            
        end
    end
    % Close the file which had been opened again
    fclose(fid);
end

%% Classify the deflections based on depth, with the same convention as stresses

% coord_node_matrix = readmatrix(large_node_matrix_path);
% deflection_matrices = dir([deflection_paths '\*.csv']);
% deflection_files = {deflection_matrices.name};
% complete_deflection_path = [deflection_paths '\Complete_deflection_matrices'];
% if ~exist(complete_deflection_path, 'dir')
%     mkdir (complete_deflection_path)
% end
% 
% for i = 1:length(deflection_files)
%     if contains(deflection_files{i},'I1') == 1
%         
%     else
%         current_deflection_matrix = readmatrix([deflection_paths '\' deflection_files{i}]);
%         current_deflection_coord_matrix = zeros(length(current_deflection_matrix),3);
%         for j = 1:length(current_deflection_matrix)
%             nodes_bool = coord_node_matrix(:,1) == current_deflection_matrix(j,1);
%             current_deflection_coord_matrix(j,:) = coord_node_matrix(nodes_bool,2:4);
%         end
%         new_matrix = [current_deflection_matrix,current_deflection_coord_matrix];
%         writematrix(new_matrix,[complete_deflection_path '\Complete_deflection_files_part_' file_identifiers{i} '.csv']);
%     end
% end

end

