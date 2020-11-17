function [node_lines,elem_lines,g,file_identifiers,nset_lines] = dat_processor(fname)
% This function reads the .dat file produced by ABAQUS and reprocesses it
% so that saving the stress and coorinate values used for the following
% calculations is easier.

% Open the .dat file

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end

% Initialize different line counters: one to store the lines where '*Node'
% appears, one for '*Element', one for '*Nset', one for 'PART INSTANCE'.
% The first three are needed to narrow down the pieces of the .dat file
% from which the node coordinates and element nodes are read, the last one
% is needed to divide the data based on the different model parts. 
% The vectors where the lines containing these keywords are stored are also 
% initialized, and another vector containing the names of each part is 
% initialized too.

line=0;
node_lines = [];
elem_lines = [];
part_lines = [];
nset_lines = [];
node_keyword_counter = 0;
elem_keyword_counter = 0;
part_keyword_counter = 0;
nset_keyword_counter = 0;
file_identifiers = {};

%% Processing
% Define all the necessary key expressions to search for in the file.

node_search_key = '*Node';
element_search_key = '*Element';
part_search_key = 'PART INSTANCE';
region_finisher_key = '*Nset';
input_paragraph_end = 'OPTIONS BEING PROCESSED';

% This cycle runs until the end of the file

while ~feof(fid)
    % Read a line without including the newline
    s=fgets(fid);
    % For every new line in the iteration, check for the different keywords
    node_keyword_flag = contains(s,node_search_key);
    elem_keyword_flag = contains(s,element_search_key);
    part_keyword_flag = contains(s,part_search_key);
    region_finisher_flag = contains(s,region_finisher_key);
    end_keyword_flag = contains(s,input_paragraph_end); 
    % Count the line we are at, which is needed to identify the lines
    % containing the search keys and then used to partition the large dat
    % file and rpocess those parts
    line = line + 1;
    
    % Increase the counter for each key expression and save the line it is
    % found at
    
    if (node_keyword_flag == 1)
        node_keyword_counter = node_keyword_counter + 1;
        node_lines(node_keyword_counter,1) = line;
    end
    if (elem_keyword_flag == 1)
        elem_keyword_counter = elem_keyword_counter + 1;
        elem_lines(elem_keyword_counter,1) = line;
    end
    if (region_finisher_flag == 1)
        nset_keyword_counter = nset_keyword_counter + 1;
        nset_lines(nset_keyword_counter,1) = line;
    end
    % In the case of the part name, store it in a separate vector after
    % manipulating its string to use it as a name for the separate csv
    % files
    if (part_keyword_flag == 1)
        part_keyword_counter = part_keyword_counter + 1;
        part_lines(part_keyword_counter,1) = line;
        file_identifiers{part_keyword_counter} = replace(strip(extractAfter(s,": ")),".","_");
    end
    % We can stop processing before the end of the file; after this key
    % expression there is no data we are interested in 
    if (end_keyword_flag == 1)
        end_line = line;
        break;
    end
end

% Save the line where we stop processing as a delimiter for the last piece
% of file to process, then close the file

node_lines(end+1) = end_line;
fclose(fid);

% Open the file again and scan it. This allows us to access every line as
% an element of a cell array that now the file is.

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end

% Replace all commas with whitespaces

Data = fileread(fname);
Data = strrep(Data, ',', ' ');

% Some lines begin with "LINE", which is not needed and can be deleted
Data = strrep(Data, 'LINE', ' ');
fid = fopen(fname, 'w');
fwrite(fid, Data, 'char');

% Finally, after writing it again, reopen the file to access it as a cell
% array where each line is a cell. This allows for easy looping and element
% inspection over each line.

fid = fopen(fname);
if fid==-1
    error('File not found or permission denied.');
end
g = textscan(fid,'%s','delimiter','\n');

end

