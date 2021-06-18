function show_event(event,ystart,yending)
%functions to show the event, new abilities to be added, temporary version.
%param @ystart the first sac to be shown
%param @yending the last sac to be shown


if nargin == 1
    ystart = 1;
    yending = event.number_of_sac;
end;

for i = ystart:yending
    plot(event.sac(i).time_coordinate,(normal_lize(event.sac(i).data) + 1+ 2 * (i-1)));
    hold on;
        
    fb = event.sac(i).first_break;
    line([fb , fb],[(0.2 + 2 * (i - 1)) , (1.8 + 2 * (i - 1))]);
    line([fb + event.p_time(2) , fb + event.p_time(2)],...
         [(0.2 + 2 * (i - 1)) , (1.8 + 2 * (i - 1))]);
    hold on;
end;