function slab_contour_on
%  slab file needed.
%   Detailed explanation goes here

    load('D:\Program_files\Matlab R2012b\toolbox\Qin\kur_contours.in');
    
    lon = kur_contours(:, 1);
    lat = kur_contours(:, 2);
    dep = kur_contours(:, 3);
      
    xlimit = round(100 * get(gca, 'xlim')) / 100;
    ylimit = round(100 * get(gca, 'ylim')) / 100;
    
    x_length = xlimit(2) - xlimit(1);
    y_length = ylimit(2) - ylimit(1);
    
    xlimit = [xlimit(1) - 2 * x_length  xlimit(2) + 2 * x_length];
    ylimit = [ylimit(1) - 2 * y_length  ylimit(2) + 2 * y_length];
    
    
    x_cor = xlimit(1):0.01:xlimit(2);
    y_cor = ylimit(1):0.01:ylimit(2);
    
    [X Y] = meshgrid(x_cor, y_cor);
    
        
    lon_in_range_index = (lon > xlimit(1) & lon < xlimit(2));
    lat_in_range_index = (lat > ylimit(1) & lat < ylimit(2));
    dep_in_range_index = (dep < 400);
    
    kur_contours_in_range = kur_contours(lon_in_range_index & lat_in_range_index & dep_in_range_index, :);
    
    lon = kur_contours_in_range(:, 1);
    lat = kur_contours_in_range(:, 2);
    dep = kur_contours_in_range(:, 3);
       
    slab = griddata(lon, lat, dep, X, Y, 'linear');
    
    
    
    
    [C, h] = contour(X, Y, slab, -720:20:-500);
    clabel(C, h, 'labelspacing', 600);caxis([-720 -500]);
    
    
    set(gca, 'Xlim', [153.1 153.8]);
    set(gca, 'Ylim', [53.7 55.4]);
    
end

