function [headers]=read_stress_report(fname,search_key,stress_paths,headers_on)
% This function reads a report file, extracting only the parts of data we are
% interested in based on the keywords delimiting them. These parts are then
% saved into separate .csv files. 
% Author: Fabrizio Morra
% Based on the RPTRead function from:
% Li Haixing; Email:windchaser_lhx@163.com

%% Opening the file and initializing the variables 
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

% Initialize two line counters, one for a first while cycle over the file
% to search for the keywords and the otehr when saving data from the last
% part of the file. Then, initialize the vector containing the numbers of
% the lines where the keywords are, the keyword counter (this is for
% verification) and the name of the .odb file regions which will be used to
% name the separate csv files.

line=0;
tline=1;
lines = [];
keyword_counter = 0;
file_identifiers = {};

%% Processing
% First of all, read the entire file and extract all the lines where a new
% region is defined. To do this, search for a keyword or keyphrase, in our
% case "reported at element". When a keyword is found, save increase the
% keyword counter, save its line and the name of the region defined after
% ':', stripping it of the newline character and substituting the dots in 
% its name with _. 

while ~feof(fid)
    % Read a line without including the newline
    s=fgets(fid);
    % Determine where the next part or region starts and store the file
    % line, also storing this line after manipulating it to use the part or
    % region name for the names of the csv files.
    keyword_index = contains(s,search_key);
    line = line + 1;
    if (keyword_index == 1)
        keyword_counter = keyword_counter + 1;
        lines(keyword_counter,1) = line;
        file_identifiers{keyword_counter} = replace(strip(extractAfter(s,": ")),".","_");
    end
    % Read and store the headers if we need them when we save the files
    headers_key = 'S.Mises';
    headers_index = contains(s,headers_key);
    if (headers_index == 1)
        headers = replace(strip(extractAfter(s,"Element ")),".","_");
    end
end
headers = ['Label';split(headers)];
headers = reshape(headers,[1 8]);

% After identifying all the lines containing the key phrase, we roughly
% know the intervals where all the stress values for that region are
% defined. But the we have float values and non-float elements, such as
% strings, and we are only interested in the nueric values. 

% Open the file again and scan it. This allows us to access every line as
% an element of a cell array that now the file is.

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end
g = textscan(fid,'%s','delimiter','\n');

% Loop over the number of keywords to evaluate the lines between each of
% them. Then use a for loop for the first n-1 keywords and a while loop for
% the last part of the file whichs tops once the file is finished to loop
% over the line intervals in which each region's stress components are
% defined. We can check for the existence of the last file to be generated
% and if it is already there we can skip this part.
if isfile([stress_paths '\I1.csv'])
    disp(['The files containing the stress components already exist, moving on' ...
        ' to nodes coordinate extraction for each part.'])
else
    disp('Processing rpt to extract stress components...')
    for i = 1:size(lines,1)
        if (i<size(lines,1))
            for j = lines(i):1:lines(i+1)
                % Access line j, scan it for float values and if they are found
                % transpose them and save them in the row of a matrix
                s=char(g{1}(j));
                [data, ncols, errmsg, nxtindex]= sscanf(s, '%f');
                if ~isempty(data)
                    data = data';
                    stress_submatrix(j-lines(i)+1,:) = data;
                end
            end
            % Delete all zero rows, save the complete stress component matrix
            % in a csv file and re-initialize it
            stress_submatrix( ~any(stress_submatrix,2), : ) = [];
            
            if (headers_on == 1)
                stress_submatrix = [headers; num2cell(stress_submatrix)];
                writecell(stress_submatrix,[stress_paths '\' file_identifiers{i} '.csv']);
                stress_submatrix = [];
            else
                writematrix(stress_submatrix,[stress_paths '\' file_identifiers{i} '.csv']);
                stress_submatrix = [];
            end
            
        else
            % Same procedure as above, but using a while cycle
            while ~feof(fid)
                s=char(g{1}(tline));
                [data, ncols, errmsg, nxtindex]= sscanf(s, '%f');
                if ~isempty(data)
                    data = data';
                    stress_submatrix(tline,:) = data;
                    tline=tline+1;
                else
                    tline=tline+1;
                end
            end
            
            if (headers_on == 1)
                stress_submatrix = [headers; num2cell(stress_submatrix)];
                writecell(stress_submatrix,[stress_paths '\' file_identifiers{i} '.csv']);
                stress_submatrix = [];
            else
                writematrix(stress_submatrix,[stress_paths '\' file_identifiers{i} '.csv']);
                stress_submatrix = [];
            end
            
        end
    end
    % Close the file which had been opened again
    fclose(fid);
end
end
