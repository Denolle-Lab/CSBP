addpath('./functions');

close all
clear all

output_fold = './Tohoku_HF_AG';
mkdir(output_fold);

A=struct('sac',[]);
A.sac=rd_event_sacs('./data/2011-03-11-mw91-near-east-coast-of-honshu-japanRmResp');
disp('read sac');
event=event_initialize(A);
disp('initialized');



event=bandpass(event,0.05,4);
%event=bandpass(event,0.05,0.2)
disp('bandpass');
event=resample(event,20);
disp('resampled');
event=travel_time_table(event);
disp('travel time table added');
event=field_fill(event,'first_break',60);



figure(22)
show_P_matrix(event);
caxis([-1 1]);

%event=time_cut(event,20,300);
event.p_time=[-2 12];

%%
% % setting the original time and arrival time for each station
% for LL=1:nsta
%     event.sac(LL).O=synsrcTime(JJ);
%     event.sac(LL).A=synsrcTime(JJ)+src2sta_time(LL,JJ);
%     event.sac(LL).STLA=staLat(LL);
%     event.sac(LL).STLO=staLon(LL);   
% end
% 
% 

figure(1);
show_P_matrix(event);
caxis([-1 1]);
title('before polarization');


event=polarize(event);

for KK=1:10
    event=adaptive_stack(event);
end

disp('Polarized and adaptive stacked');

figure(2)
show_P_matrix(event);
caxis([-1 1]);

title('after polarization and stacking');

% 20 and 240s is the time range for the waveforms to be saved. The time
% here is relative to the starting of the waveform (X on the shown Figure 2)
save_aligned(event,20,240,output_fold,'Tohoku_event_data')

