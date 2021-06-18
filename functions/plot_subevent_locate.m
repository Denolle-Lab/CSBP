function plot_subevent_locate(data, n, i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    epicenter = [153.2530, 54.8850];
    strike = 12;
    t_max = data(:, 3);
    t_s = data(:, 4);
    t_e = data(:, 5);
    lon = data(:, 9);
    lat = data(:, 10);
    Mag = data(:, 12);
    Mag = Mag / max(abs(Mag));
    
    % colormap generated 
    cm = my_colormap_maker([1 0 0], [1 1 1], [0 0 1], n);
    face_color = cm(i, :);    
    
    plot(epicenter(1), epicenter(2), 'pr', 'MarkerSize', 20, 'MarkerFaceColor', 'r');
    
    color = cm(i,:);
    
    for i = 1:size(lon,1)
        plot(lon(i), lat(i), 'ok', 'MarkerSize', Mag(i) * 15, 'LineWidth', 2, 'MarkerFaceColor', face_color);
    end
    
    
    

end

