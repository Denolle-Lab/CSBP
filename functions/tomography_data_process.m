
clear all;

lat_range = [-90 90];
lon_range = [1 359];
dep_range = [-1500 -100];

load('I:\gloabal tomography by Hilst\Lat.mat');
lat = Lat;
clear Lat;

load('I:\gloabal tomography by Hilst\Long.mat');
lon = Long;
clear Long;

load('I:\gloabal tomography by Hilst\Depth.mat');
dep = Depth;
clear Depth;

load('I:\gloabal tomography by Hilst\dVp.mat');

dep = - dep;

II_lat = lat > lat_range(1) & lat < lat_range(2);
II_lon = lon > lon_range(1) & lon < lon_range(2);
II_dep = dep > dep_range(1) & dep < dep_range(2);
II = II_lat & II_lon & II_dep;

lat = lat(II);
lon = lon(II);
dep = dep(II);
dVp = dVp(II);

lat_axis = unique(lat);
lon_axis = unique(lon);
dep_axis = flipud(unique(dep));

n_lat = numel(lat_axis);
n_lon = numel(lon_axis);
n_dep = numel(dep_axis);

lat_grid = reshape(lat, n_lat, n_lon, n_dep);
lon_grid = reshape(lon, n_lat, n_lon, n_dep);
dep_grid = reshape(dep, n_lat, n_lon, n_dep);
dVp_grid = reshape(dVp, n_lat, n_lon, n_dep);

%% vertical slice

figure;

slice1_x = 151:0.5:156;
slice1_y = 208 - slice1_x;
slice1_z = -300:-20:-900;

slice1_xdata = ones(numel(slice1_z), 1) * slice1_x; 
slice1_ydata = ones(numel(slice1_z), 1) * slice1_y; 
slice1_zdata = slice1_z' * ones(1, numel(slice1_x));

slice(lon_grid, lat_grid, dep_grid, dVp_grid, slice1_xdata, slice1_ydata, slice1_zdata);

cm = colormap_maker(256, 'velocity');
colormap(cm);

%% the justified color range

vi = interp3(lon_grid, lat_grid, dep_grid, dVp_grid, 153.28, 54.874, -600)
caxis([vi - 0.5, vi + 0.5]);

% caxis([-0.5 0.5]);
shading interp;

hold on;
plot3(153.1905, 54.885, -609, '.b', 'MarkerSize', 50);
plot3(153.5032, 55.1008, -609, '.c', 'MarkerSize', 50);
plot3(153.3468, 54.6152, -609, '.g', 'MarkerSize', 50);
plot3(153.4719, 54.1835, -609, '.y', 'MarkerSize', 50);
plot3(153.2530, 54.0396, -609, '.r', 'MarkerSize', 50);

slice2_y = 56:-0.5:53;
slice2_x = 153.28 * ones(size(slice2_y));
slice2_z = -400:-20:-900;

slice2_xdata = ones(numel(slice2_z), 1) * slice2_x; 
slice2_ydata = ones(numel(slice2_z), 1) * slice2_y; 
slice2_zdata = slice2_z' * ones(1, numel(slice2_x));

hold on;
slice(lon_grid, lat_grid, dep_grid, dVp_grid, slice2_xdata, slice2_ydata, slice2_zdata);


cm = colormap_maker(256, 'velocity');
colormap(cm);
caxis([-0.5 0.5]);
shading interp;
caxis([vi - 0.3, vi + 0.3]);


%% horizontal global slice
% slice(lon_grid, lat_grid, dep_grid, dVp_grid, 155, 55, -750);
% cm = colormap_maker(256, 'velocity');
% colormap(cm);
% caxis([-0.8 0.8]);
% shading interp;


%% aftershock.
% epicenter = [153.28, 54.874];
% depth = -609;
% load('I:\BPJ copy\aftershock.mat')
% aftershock_distribution_3D(epicenter, depth, aftershock);

colorbar


