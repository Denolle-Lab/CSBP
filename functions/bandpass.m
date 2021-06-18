function event = bandpass(event,freqlow,freqhigh);

if nargin == 1
    freqlow = 0.05;
    freqhigh = 4;
end;

for i = 1:event.number_of_sac
    dt = event.sac(i).DELTA;
    [b,a] = butter(3,[(2 * freqlow * dt),(2 * freqhigh * dt)],'bandpass');
    event.sac(i).data = filtfilt(b,a,event.sac(i).data);

end;

return;
    