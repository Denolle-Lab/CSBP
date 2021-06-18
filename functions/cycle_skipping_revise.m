function event = cycle_skipping_revise(event)
%revise the cycle-skipping

nos = event.number_of_sac;
threshold = 0.5; %the threshold controls the standard of cycle-skipping

r = residual(event);

aver = mean(r);

for i = 1:nos
    r1 = abs((r(i) - aver) / event.p_time(2));
    if r1 > threshold
        if r(i) > aver;
            event.sac(i).first_break = event.sac(i).first_break ...
                                       - event.p_time(2);
        else
            event.sac(i).first_break = event.sac(i).first_break ...
                                       + event.p_time(2);
        end;
    end;
end;

event.range = round(threshold * event.p_time(2));

return;