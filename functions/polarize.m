function event = polarize(event)
%% To make every sac is positive  
% the tech is using the cross corelation coeffient of each wave and the
% stacked wave.

%% get the corelation coeffient.
% note r is a column 
r = corr_coef(event);

%% find the negative wave
nega_wave_num = find(r < 0);

%% change the polarity of target wave.

for i = nega_wave_num' % row vector is correct.
    event.sac(i).data = -1 * event.sac(i).data;
    event.sac(i).pol = 1;
end;

return
    