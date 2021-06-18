function flag = dist_summary
% @function dist_summary. Put the subeventSummary_reloc files together.
% @param flag. To show if the function was run successfully. it's 1 when it
%             is successful.
% @param epicenter. the epicenter of the earthquake. a transversal vector
%                   of size of 1 * 2.
% @param strike. the strike of the fault in degree.

    flag = 0;
    a = dir;
    a(1) = [];
    a(1) = []; % the a(1) and a(2) are hidden files.
    hold on;
    
    for i = 1:size(a,1);
        eval(strcat('plot_data = load(''', a(i).name, ...
            ''', ''newsubeventdata'');'));
           % plot_data = load(a(i),name, 'newsubeventdata');
        
        dist_plot(plot_data.newsubeventdata, size(a,1), i);
            
        clear plot_data;    
            
    end
    
    
    Ylim = get(gca, 'Ylim');
    set(gca, 'Ylim', [0 Ylim(2)]);
    Xlim = get(gca, 'Xlim');
    Ylim = get(gca, 'Ylim');
    set(gca, 'XDir', 'reverse');
    xlabel('Along-Rupture Distance (km)', 'FontSize', 12); 
    ylabel('Time (s)', 'FontSize', 12);
    
    for i = 1:6
        plot([0 Xlim(1)], [0 abs(Xlim(1) / i)], 'g--');
        plot([0 Xlim(2)], [0 abs(Xlim(2) / i)], 'g--');
    end
    
%     %figure;
       
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'box', 'on');
    
        
    flag = 1;
    
    
    
    
end

