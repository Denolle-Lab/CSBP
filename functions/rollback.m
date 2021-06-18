function event = rollback(event)
%when one step of the operation is not good, using backup to recover the
%data. only first_break and range can be recovered, can be added.

nos = event.number_of_sac;

event.range = event.backup.range;
for i = 1 : nos
    event.sac(i).first_break = event.sac(i).backup.fb;
end;

end

