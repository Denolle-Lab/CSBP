function dist_plot(data, n, i)
% @function dist_plot. plot event_dist figure with Yao's subeventreloc.mat 
% @para @data. the field 'newsubeventdata' data in the mat file.
% @param n. the number of result in different band.                   
% @param i. the order of a subevent in one BP result.
%       n and i are used together to identify subevent color, the color
%       from blue to red means time order.

    epicenter = [153.2530, 54.8850];
    strike = 12;
    cm = jet(n);
    event_num = size(data, 1);
    t_max = data(:, 3);
    t_s = data(:, 4);
    t_e = data(:, 5);
    lon = data(:, 9);
    lat = data(:, 10);
    dist = data(:, 8);
    Mag = data(:, 12);
    Mag = Mag / max(Mag);
    Amp = data(:, 11);
    Amp = Amp / max(Amp);
    
     direction = (- strike - 90) / 180 * pi;
     r = [cos(direction) sin(direction)];   
%  the global arc distance   
%     dist = zeros(event_num, 1);
%     
     X = lon - epicenter(1);
     Y = lat - epicenter(2);
     
     pol = sign([X Y] * r');
%     dist = dist * 2 * pi * earthRadius('km') / 360;
     dist = dist .* pol;
        
     hold on;
     
    for i = 1:event_num
        plot(dist(i), t_max(i), '.r','Markersize', 60 * Amp(i), 'MarkerFaceColor', 'r');
    end
    
    line([dist'; dist'], [t_s'; t_e'],'color','r', 'linewidth', 2); 

  
    
    
end

