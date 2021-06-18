%% adaptive stack

function event = adaptive_stack(event)

nos = event.number_of_sac;
range = event.range;
p_length = round((event.p_time(2) - event.p_time(1)) / event.sac(1).dt);  %second to point
sd = zeros(1,p_length);

NOS=nos;
i=1;
while(i<=NOS)
    pill = round((event.sac(i).first_break + event.p_time(1)) / event.sac(i).dt);%the lower limit of the ith SAC
    piul = round((event.sac(i).first_break + event.p_time(2)) / event.sac(i).dt); % the upper limit
    p_i = normal_lize(event.sac(i).data(pill:(piul - 1))); % the stacking
    
    if ~isnan(p_i)
        sd = sd + p_i;
    else
        NOS=NOS-1;
        event.sac(i)=[];
    end
    i=i+1;    
end
nos=NOS;
sd = sd / nos;

for i = 1:nos
    pill = round((event.sac(i).first_break + event.p_time(1)) / event.sac(i).dt); % the lower limit.
    piul = round((event.sac(i).first_break + event.p_time(2)) / event.sac(i).dt);%the upper limit
    p_i = normal_lize(event.sac(i).data(pill:piul));
    [c,lag] = xcorr(p_i,sd);
    offset = lag(maxxc(c));
    event.sac(i).backup.fb = event.sac(i).first_break;
    event.sac(i).first_break = event.sac(i).first_break + offset * event.sac(i).dt;
end;

event.backup.range = event.range;
event.range = 5;

event.number_of_sac=nos;
return;


