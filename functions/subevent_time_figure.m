function subevent_time_figure
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    load('subeventreloc.mat');
    coast = load('coast');
    sub_event_data = newsubeventdata;
    epicenter = [153.2530 54.8850];
    
    
    cmap = colormap(jet);
    colormap(cmap);
    ncol = size(cmap,1);
    nsubsrc = size(sub_event_data, 1);
    
    MaxAmpSrc = max(sub_event_data(:, 11));
    MaxTimeSrc = round(max(sub_event_data(:, 5)));
    
    for i = 1:nsubsrc
        colorindex = round((sub_event_data(i,3)/MaxTimeSrc)*(ncol-1)) + 1;
        plot(sub_event_data(i, 9), sub_event_data(i,10), 'bo', 'MarkerSize', 15*sub_event_data(i,11)/MaxAmpSrc, ...
            'MarkerFaceColor', cmap(colorindex,:));
        hold on;
        text(sub_event_data(i, 9), sub_event_data(i,10), num2str(round(sub_event_data(i,3))));
    end
    colorbar; colormap(cmap); caxis([0 MaxTimeSrc]);hold on;
   
    axis image

    slab_contour_on;
    
    xlim = get(gca, 'xlim');
    ylim = get(gca, 'ylim');
    
    
%     MaxDist2hypo = max(sqrt(sub_event_data(:,6).^2 + sub_event_data(:,7).^2)) / 111;
%     ddist = 10 / 111;
%     numdist = ceil(MaxDist2hypo/ddist);
% 
     hold on; plot(epicenter(1), epicenter(2), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
%     for i = 1:numdist
%         if mod(i,5) == 0
%             hold on; plot(epicenter(1) + i*ddist*cosd(0:2:360), ...
%                  epicenter(2) + i*ddist*sind(0:2:360), 'k');
%         else
%             hold on; plot(epicenter(1) + i*ddist*cosd(0:2:360), ...
%                 epicenter(2) + i*ddist*sind(0:2:360), 'k--');
%         end
%     end
    
    plot(coast.long, coast.lat);
    

    
    set(gca, 'xlim', [xlim(1) - 0.1 * (xlim(2) - xlim(1)), xlim(2) + 0.1 * (xlim(2) - xlim(1))]);
    set(gca, 'ylim', [ylim(1) - 0.1 * (ylim(2) - ylim(1)), ylim(2) + 0.1 * (ylim(2) - ylim(1))]);
    
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'box', 'on');
    
    
end

