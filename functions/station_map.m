function station_map(event)
%draw station map
% 

nos = event.number_of_sac;

figure;
hold on;

for i = 1:nos
    plot(event.sac(i).STLO,event.sac(i).STLA,'o');
    display([event.sac(i).STLO,event.sac(i).STLA]);
end

plot(event.sac(1).EVLO,event.sac(1).EVLA,'p');

end

