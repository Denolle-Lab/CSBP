function event = resample(event,freq) % convert event to specific frequency

if nargin == 1
    freq = 10;
end;

for i = 1:event.number_of_sac
    event.sac(i).dt = 1 / freq;
    event.sac(i).length = round(event.sac(i).NPTS * (event.sac(i).DELTA / event.sac(i).dt));
    newt = linspace(event.sac(i).B,event.sac(i).E,event.sac(i).length);   % new time_coordinate
    newa = interp1(event.sac(i).time_coordinate, event.sac(i).data, newt, 'spline');  %new amplitude
    event.sac(i).time_coordinate = newt;
    event.sac(i).data = newa;
end;

return;