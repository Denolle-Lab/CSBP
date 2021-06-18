function flag = aftershock_compare
% @function aftershock_compare. Comparing the IBP result and aftershock
%       distribution. Note that the current very is not good enough. Now
%       pre-process is need. I will think about a good version later.
% @param @aftershock. A mtrix containing the aftershock information. 

    load('D:\Program_files\Matlab R2012b\toolbox\Qin\Okhotsk_aftershock.mat');
    epicenter = [153.2530, 54.8850];
    flag = 0;
    
    lon = aftershock(:, 2);
    lat = aftershock(:, 1);
    dep = -aftershock(:, 3);
    mag = aftershock(:, 4);
    
    hold on;
    
    h1 = plot(epicenter(1), epicenter(2), 'pr', 'MarkerFaceColor','r', ...
            'MarkerSize', 20);

    
    mag = mag / 8.3;
    
    cm = colormap(jet(256));
    colormap(jet);
    
    
    for i = 1:size(aftershock, 1)
        color_index = round(256 * (dep(i) - (-720)) / (-(500) - (-720)));
        plot(lon(i), lat(i), 'ok','MarkerSize', mag(i) * 20, 'lineWidth', 2, ...
            'MarkerFaceColor', cm(color_index, :));
            
    end
    colorbar;caxis([-720 -500]);
    

    axis image;
    grid on;


end

