function [] = B_plots(viscosity_figures_path,alin,a,data_points_indices,part,...
    min_depth,max_depth,min_lat,max_lat,min_lon,max_lon,lat,lon,run)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

variables = [alin,a];
names = {'B_{diff}','B_{disl}'};
fig_counter = findobj('type','figure');
for i=1:size(variables,2)
    figure(i)
    cla
    variable = variables(:,i);
    plot_variable = variable(data_points_indices);
    [Z, refvec] = geoloc2grid(lat(data_points_indices),lon(data_points_indices),...
        plot_variable,0.5);
    load coastlines;
    cmap = colormap('jet');
    alpha 0.7;
    colormap(cmap);
    latlim = [min_lat max_lat];
    lonlim = [min_lon max_lon];
    ax = axesm('stereo','MapLatLimit',latlim,'MapLonLimit',lonlim,'Grid','on','MeridianLabel','on','ParallelLabel','on');
    set(ax,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    geoshow(Z, refvec, 'DisplayType', 'texture');
    plotm(coastlat,coastlon);
    caxis('auto');
    h = colorbar('v');
    if i == 1
        set(get(h,'ylabel'),'string','B_{diff}')
    else
        set(get(h,'ylabel'),'string','B_{disl}')
    end
    title([names{i} ' for part ' part ' with depth range '...
        num2str(min_depth) '-' num2str(max_depth) ' km']);
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    grid on;
    % Save the figures based on whether we are working with
    % stresses or deflections
    saveas(gcf,[viscosity_figures_path '\' names{i} '_' part '_' run...
        '_[' num2str(min_depth) '-' num2str(max_depth) ']_km.png']);
end
end
