function event = net_extract(event,netwk)
%get the event of the same net, others are deleted

nos = event.number_of_sac;

for i = nos:-1:1
    if ~(strcmp(event.sac(i).KNETWK,netwk));
        event = delete_sac(event,i);
    end;    
end;

return;

