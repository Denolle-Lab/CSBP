% close all
% clear all



A=struct('sac',[]);
%A.sac=rd_event_sacs('/Users/Yin9xun/Work/Earthquake_data/Hinet_data/2016NZ/2016NZ_hinet/Hinet_SACRmResp/Event');
A.sac=rd_event_sacs('.');


disp('read sac');
event=event_initialize(A);
disp('initialized');



%event=bandpass(event,1,2);
event=bandpass(event,0.01,4.5)
disp('bandpass');
event=resample(event,20);
disp('resampled');
event=travel_time_table(event);
disp('travel time table added');
event=field_fill(event,'first_break',30);



figure(22)
show_P_matrix(event);
caxis([-1 1]);

%event=time_cut(event,20,300);
event.p_time=[-2 6];

%%
% setting the original time and arrival time for each station
for LL=1:nsta
    event.sac(LL).O=synsrcTime(JJ);
    event.sac(LL).A=synsrcTime(JJ)+src2sta_time(LL,JJ);
    event.sac(LL).STLA=staLat(LL);
    event.sac(LL).STLO=staLon(LL);   
end



figure(1);
show_P_matrix(event);
caxis([-1 1]);
title('before polarization');


%event=polarize(event);

% for KK=1:10
%     event=adaptive_stack(event);
% end
% 
% disp('Polarized and adaptive stacked');

figure(2)
show_P_matrix(event);
caxis([-1 1]);

title('after polarization and stacking');

%% For practical use, no depth information
%save_aligned(event,10,100,[pwd],['KineSource_L_200_R_0.5_syndata'])
%% Save depth also, For kinematic source only !
save_aligned_depth(event,10,180,[pwd],['KineSource_L_200_R_0.5_syndata'])


figure
Nsac=length(event.sac);
for i=1:Nsac
plot(event.sac(i).time_coordinate,event.sac(i).data/max(abs(event.sac(i).data)),'LineWidth',0.2,'Color',[Nsac+1-i 0 0]/Nsac);
hold on
end