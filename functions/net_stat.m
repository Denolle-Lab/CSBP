function net = net_stat(event)
%the stat of how the data was composed
%i.e. how many sacs from which net

net(1).n = event.sac(1).KNETWK; %name
net(1).q = 1; %quanlity


nos = event.number_of_sac;
m = 0;      % m is the loop control;
k = 0;      % l is the monitor of if the sac is added to current net

for i = 2 : nos
    if k == 0
        m = m + 1;
    else
        k = 0;
    end;
    
    for j = 1 : m;
        if (strcmp(event.sac(i).KNETWK , net(j).n))
            net(j).q = net(j).q + 1;
            k = 1;
        end;
    end;
    
    if k == 0
        net(m + 1).n = event.sac(i).KNETWK;
        net(m + 1).q = 1;
    end;

end

return;
