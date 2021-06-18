%% travel time table calculation using table of tt91.P.
%%

function event = travel_time_table(event)

nos = event.number_of_sac;
TravelTimeFile = 'tt91_o.P';

% read travel time table
ftb = fopen(TravelTimeFile, 'r');
if ftb == -1
    display('no travel time table');
end;
temp = fgetl(ftb); % read first line
nepidist = fscanf(ftb, '%d', 1); % number of epicenter distance in the table, second line
nsrcdep = fscanf(ftb, '%d', 1); % number of depth in the table, second line
temp = fgetl(ftb); 
srcdep = fscanf(ftb, '%f', nsrcdep); % read third line: source depth vector
tmptable = fscanf(ftb, '%f', [nsrcdep+1, nepidist]); % read travel time table
ttable = tmptable(2:(nsrcdep+1),:);
epidist = tmptable(1, :);
fclose(ftb);

for i = 1:nos
    edep = event.sac(i).EVDP; % event depth
    
    if isnan(event.sac(i).STDP)
        event.sac(i).STDP = 0;
    end
    
    sel = (event.sac(i).STEL - event.sac(i).STDP) / 1000; % station depth
    dep = edep + sel;
    deg = event.sac(i).GCARC   ; % the distance between event and station in degrees.

    
    if (deg > min(epidist)) & (deg < max(epidist)) & (dep > min(srcdep)) & (dep < max(srcdep))              
        temp = interp2(epidist, srcdep, ttable,deg,dep,'linear'); % interpolation in the table in minute.
    else
        temp = 1;
    end;
    
    if temp > 0
        event.sac(i).tt_time = temp * 60;%minute to second
        event.sac(i).first_break = event.sac(i).tt_time + event.sac(i).O;  %traveltime table time
        end;
end;

event.range = 40;

fb = first_break(event);

for i = 1:nos
    if fb(i) < 0
        fb(i) = mean(fb);
    end
end

event = set_pick(event,fb);

return;
    
    
    