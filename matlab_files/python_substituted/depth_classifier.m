function [classified_path,files_to_classify,maximum_depth] = depth_classifier(...
    sd_input,deflection_processing_path,individual_path,components_to_plot,...
    complete_files_path,spherical_complete_files_path,headers_on)
% This function takes as input the paths to the complete stress or
% deflection files containing the stress or deflection components and the
% corresponding centroid or node coordinates and discretizes the results
% based on depth. 

% Define parameters that change based on whether we are working with
% stresses or deflections, like plot names and file path parts. 

if sd_input == 0
    upper_path = individual_path;
    save_name = 'Stresses_until_depth_';
    histogram_name_part = 'centroids';
    matrix_to_save_part = 'Stresses_layer_';
    earth_name_part = 'Stresses_Earth_';
    discretized_headers = {'Label','Mises','S_11','S_22','S_33','S_12','S_13',...
        'S_23','X','Y','Z','R','Depth','Lat','Lon'};
else
    upper_path = deflection_processing_path;
    save_name = 'Deflections_until_depth_';
    histogram_name_part = 'nodes';
    matrix_to_save_part = 'Deflections_layer_';
    earth_name_part = 'Deflections_Earth_';
    discretized_headers = {'Label','U_Magn','U_1','U_2','U_3','X','Y','Z',...
        'R','Depth','Lat','Lon'};
end

% Define the necessary path and folders and create them if they are not
% there, differentiating between the two cases: working with cartesian or
% spherical stress components.

unspherical_classified_depth_path = [upper_path '\Depth_classified_files'];
spherical_classified_depth_path = [upper_path '\Spherical_depth_classified_files'];
    
if strcmp(components_to_plot, 'spherical') == 1
    files_to_classify_path = dir([spherical_complete_files_path '\*.csv']);
    files_to_classify = {files_to_classify_path.name};
    classified_path = spherical_classified_depth_path;
    histogram_path = [spherical_classified_depth_path '\histogram_plots'];
    if ~exist(classified_path, 'dir')
        mkdir (classified_path);
    end
    if ~exist(histogram_path, 'dir')
        mkdir (histogram_path);
    end
elseif strcmp(components_to_plot, 'cartesian') == 1
    files_to_classify_path = dir([complete_files_path '\*.csv']);
    files_to_classify = {files_to_classify_path.name};
    classified_path = unspherical_classified_depth_path;
    histogram_path = [unspherical_classified_depth_path '\histogram_plots'];
    if ~exist(classified_path, 'dir')
        mkdir (classified_path)
    end
    if ~exist(histogram_path, 'dir')
        mkdir (histogram_path)
    end
else
    disp('Incorrect input, select either spherical or unspherical.');
end
earth_data_path = [classified_path '\No_layers_subdivision'];
% Only process part files, taking away those relative to the regions.
% Region I1 will be skipped as it only contains empty values while I0 is
% not processed because it belongs to the core, which also does not
% interest us.

for i = 1:length(files_to_classify)
    size_to_evaluate = dir([files_to_classify_path(i).folder '\' files_to_classify{i}]);
    if contains(files_to_classify{i},'Region') == 1
        files_to_classify{i} = {};
    elseif size_to_evaluate.bytes < 100
        files_to_classify{i} = {};
    elseif contains(files_to_classify{i},'I0') == 1
        files_to_classify{i} = {};
    elseif contains(files_to_classify{i},'Large') == 1
        files_to_classify{i} = {};
    end
end
files_to_classify = files_to_classify(~cellfun(@isempty,files_to_classify));
    
% To classify the stress values based on layers, we must build a vector
% containing different layer depths. We know 410 and 660 from the literature, 
% but the depth of  the last one is determined by the maximum element depth 
% ABAQUS has generated, and we must find it.
    
maximum_depth = 0;
earth_radius = 6371000;
earth_data = [];
cfn = length(findobj('type','figure'));

% For every file we read, we need to calculate the radial cemtroid distance
% of every node.

for i = 1:length(files_to_classify)
    % Read the data file and allocate matrices to store depth and radial
    % distance of each centroid or node
    data_matrix = readmatrix([files_to_classify_path(i).folder '\' files_to_classify{i}]);
    layer_depth = zeros(length(data_matrix),1);
    radial_centroid_distance = zeros(length(data_matrix),1);
    % Calculate the radial distance of each centroid or node
    for j = 1:length(data_matrix)
        radial_centroid_distance(j) = sqrt(data_matrix(j,end-2)^2+data_matrix(j,end-1)^2+data_matrix(j,end)^2);
    end
    % Calculate the maximum depth, this will be the lower limit for the
    % lower mantle based on how the model is generated 
    for j = 1:length(radial_centroid_distance)
        radius_to_search = radial_centroid_distance(j);
        layer_depth(j) = earth_radius - radius_to_search;
        if layer_depth(j)>maximum_depth
            maximum_depth = layer_depth(j);
        end
    end
    depthmatrix = [data_matrix,radial_centroid_distance,layer_depth];
    earth_data = [earth_data;depthmatrix];
end

% Differentiate based on whether we study the core too or not

if maximum_depth > 2.886e6
    layer_depth_km = [410,660,2886,maximum_depth/1000];
    layer_depth = layer_depth_km*1000;
else
    layer_depth_km = [410,660,maximum_depth/1000];
    layer_depth = layer_depth_km*1000;
end

check_1 = 1;
check_2 = 1;
run_input = 1;

% Decide whether to run or not the depth classifier proper and if we want
% to change the number of bins we are using

if exist([earth_data_path '\' earth_name_part(1:end-1) '.csv'])
    while check_1 == 1
        run_input = input(['Do you want to rerun the depth discretization algorithm? ' ...
            '1 if yes, 0 if no: \n']);
        if run_input == 0 
            check_1 = 0;
        elseif run_input == 1
            while check_2 == 1
                bin_input = input(['Do you want to change the bin number?' ...
                    ' 1 if yes, 0 if no: \n']);
                if bin_input == 0
                    check_2 = 0;
                    check_1 = 0;
                elseif bin_input == 1
                    n_bins = input('Enter number of layer bins:  \n');
                    earth_bins = input('Enter number of model bins: \n');
                    check_2 = 0;
                    check_1 = 0;
                else
                    disp('incorrect input, select either 1 or 0. \n');
                end
            end
        else
            disp('incorrect input, select either 1 or 0. \n');
        end
    end
end

%--------------------------------------------------------------------------

if run_input == 0
    disp('Skipping this stage to go to geographic plots.')
else
    disp('Starting discretization of data into different depths');
    % We can now create a file for each layer depth based on the values we have
    % just found. First we iterate over the part files we still analyze, then
    % we store the values into the depths given by the vector and we save the
    % files with stresses at depth.

    large_depth_matrix = [];
    
    % Create a complete matrix using the stress or deflection matrix
    % containing stresses or deflections and centroid or nodes coordinates,
    % then appending the centroid distance and depth of each of those
    % points. Also iteratively add matrices to a larger csv file to store
    % the values for each part.
    
    for i = 1:length(files_to_classify)
        data_matrix = readmatrix([files_to_classify_path(i).folder '\' files_to_classify{i}]);
        radial_centroid_distance = zeros(length(data_matrix),1);
        depths = zeros(length(data_matrix),1);

        for k = 1:length(data_matrix)
            radial_centroid_distance(k) = sqrt(data_matrix(k,end-2)^2+data_matrix(k,end-1)^2+data_matrix(k,end)^2);
            depths(k) = earth_radius - radial_centroid_distance(k);
        end

        complete_matrix = [data_matrix, radial_centroid_distance, depths];
        large_depth_matrix = [large_depth_matrix; complete_matrix];
        
    end
    
    if sd_input == 0
        cartesian_coordinates = large_depth_matrix(:,9:11);
    else
        cartesian_coordinates = large_depth_matrix(:,6:8);
    end
    [lat,lon] = cart2geo(cartesian_coordinates);
%     [lat_lon_matrix] = ecef2lla(cartesian_coordinates);
%     lat = lat_lon_matrix(:,1);
%     lon = lat_lon_matrix(:,2);
    large_depth_matrix = [large_depth_matrix,lat,lon];
    
    % After creating the complete large matrix, divide it different layers
    % based on the centroid or nodes depth

    for i = 1:length(layer_depth)
        layer_group_indices = large_depth_matrix(:,end-2) < layer_depth(i);
        layer_matrix = large_depth_matrix(layer_group_indices,:);
        if headers_on == 1
            layer_matrix_wheaders = [discretized_headers;num2cell(layer_matrix)];
            writecell(layer_matrix_wheaders, ...
                [classified_path '\' save_name num2str(layer_depth_km(i)) '_km.csv']);
        else
            writematrix(layer_matrix, ...
                [classified_path '\' save_name num2str(layer_depth_km(i)) '_km.csv']);
        end
        large_depth_matrix(layer_group_indices,:) = 0;
        large_depth_matrix = large_depth_matrix(any(large_depth_matrix,2),:);
        
    end

    % Another step is necessary, and that is now discretizing the data
    % found at each layer in a number of bins decided by the user. After
    % deciding the number of bins, the names of the matrices regarding each
    % layer are read.

    n_bins = 20;
    layer_values = dir([classified_path '\*.csv']);
    layer_values_names = natsortfiles({layer_values.name});

    for i = 1:length(layer_values_names)
        subclassified_path = [classified_path '\subclassification_depth_' ...
            num2str(regexp(layer_values_names{i}, '[\d\.]+', 'match','once')) '_km'];
        if ~exist(subclassified_path, 'dir')
            mkdir (subclassified_path)
        end
        layer_stress_file = readmatrix([classified_path '\' layer_values_names{i}]);
        depth_data = layer_stress_file(:,end-2);
        
        % Generate histograms to show the distribution of centroids or nodes 
        % with depth, based on the number of bins 
        
        figure(cfn+i)
        histogram(depth_data,n_bins);
        xlabel('Depth [m]'); ylabel(['Number of ' histogram_name_part]);
        if i == 1
            title({'Distribution of ' histogram_name_part ' from the ABAQUS model at different depths for',...
                ['layer [0-' num2str(regexp(layer_values_names{i},...
                '[\d\.]+', 'match','once')) '] km']});
            grid on;
            saveas(gcf,[histogram_path '\Layer_[0-' num2str(regexp(layer_values_names{i},...
                '[\d\.]+', 'match','once')) ']_km''_distribution.png'])
        else
            title({'Distribution of ' histogram_name_part ' from the ABAQUS model at different depths for',...
                ['layer [' num2str(regexp(layer_values_names{i-1}, '[\d\.]+', 'match','once'))...
                '-' num2str(regexp(layer_values_names{i}, '[\d\.]+', 'match','once')) '] km']});
            grid on;
            saveas(gcf,[histogram_path '\Layer_[' num2str(regexp(layer_values_names{i-1}, '[\d\.]+', 'match','once'))...
                '-' num2str(regexp(layer_values_names{i}, '[\d\.]+', 'match','once')) ']_km_distribution.png'])
        end

        % Save the matrices based on the binning subdivision 
        [indices,~] = discretize(depth_data,n_bins);
        for j = 1:n_bins
            subdivision_matrix = layer_stress_file(indices == j,:);
            if headers_on == 1
                subdivision_matrix_wheaders = [discretized_headers;num2cell(subdivision_matrix)];
                writecell(subdivision_matrix_wheaders, [subclassified_path '\' matrix_to_save_part ...
                    num2str(regexp(layer_values_names{i}, '[\d\.]+', 'match','once')) ...
                    '_km_bin_' num2str(j) '.csv']);
            else
                writematrix(subdivision_matrix, [subclassified_path '\' matrix_to_save_part ...
                    num2str(regexp(layer_values_names{i}, '[\d\.]+', 'match','once')) ...
                    '_km_bin_' num2str(j) '.csv']);
            end
        end

    end
    
    %----------------------------------------------------------------------

    % We can do the same, which means generating histograms and dividing a
    % larger matrix into binned ones for the entire earth_point part, the
    % one we are interested the most in.
    
    earth_bins = 20;
    [allearth_indices,~] = discretize(earth_data(:,end),earth_bins);
    earth_data_path = [classified_path '\No_layers_subdivision'];
    if ~exist(earth_data_path, 'dir')
        mkdir (earth_data_path)
    end
    if sd_input == 0
        cartesian_coordinates = earth_data(:,9:11);
    else
        cartesian_coordinates = earth_data(:,6:8);
    end
    [lat,lon] = cart2geo(cartesian_coordinates);
    earth_data = [earth_data,lat,lon];
    
    % Save matrices
    
    for j = 1:earth_bins
        subdivision_matrix = earth_data(allearth_indices == j,:);
        if headers_on == 1
            subdivision_matrix_wheaders = [discretized_headers;num2cell(subdivision_matrix)];
            writecell(subdivision_matrix_wheaders, [earth_data_path '\' earth_name_part ...
                'bin_' num2str(j) '.csv']);
        else
            writematrix(subdivision_matrix, [earth_data_path '\' earth_name_part ...
                'bin_' num2str(j) '.csv']);
        end
    end

    cfn = length(findobj('type','figure'));
    % Generating histograms after getting the number of active figures to
    % not overwrite them
    figure(cfn+1)
    histogram(earth_data(:,end-2),earth_bins);
    xlabel('Depth [m]'); ylabel(['Number of ' histogram_name_part]);
    title(['Distribution of ' histogram_name_part ' from the entire ABAQUS model at different depths, n_{bins}='...
        num2str(earth_bins)]);
    grid on;
    saveas(gcf,[histogram_path '\Complete_Earth_distribution.png'])
    
    if headers_on == 1
        earth_matrix_wheaders = [discretized_headers;num2cell(earth_data)];
        writecell(earth_matrix_wheaders, [earth_data_path '\' earth_name_part(1:end-1) '.csv']);
    else
        writematrix(earth_data, [earth_data_path '\' earth_name_part(1:end-1) '.csv']);
    end
end

end

