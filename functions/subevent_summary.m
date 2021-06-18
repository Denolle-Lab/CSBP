function subevent_summary
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    flag = 0;
    a = dir;
    a(1) = [];
    a(1) = []; % the a(1) and a(2) are hidden files.
    hold on;
    
    for i = 1:size(a,1);
        eval(strcat('plot_data = load(''', a(i).name, ...
            ''', ''newsubeventdata'');'));
           % plot_data = load(a(i),name, 'newsubeventdata');
        
        plot_subevent_locate(plot_data.newsubeventdata, size(a,1), i);
            
        clear plot_data;    
            
    end
    
    axis image;
    
    xlim = get(gca, 'xlim');
    ylim = get(gca, 'ylim');
    
%     slab_contour_on;
    
    zoom_out_coef = 0.4;
    
    set(gca, 'xlim', [xlim(1) - zoom_out_coef * (xlim(2) - xlim(1)), xlim(2) + zoom_out_coef * (xlim(2) - xlim(1))]);
    set(gca, 'ylim', [ylim(1) - zoom_out_coef * (ylim(2) - ylim(1)), ylim(2) + zoom_out_coef * (ylim(2) - ylim(1))]);
    
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'box', 'on');

end

