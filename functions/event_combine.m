function event_out = event_combine(event1, event2)
%% combine two events to do conventional backprojection
% The combination criteria is to make the onset of different waves aligned.
% This function is not well disigned.
% The user need to ensure the two events are for the same earthquake, and
% the time, depth , p wave length and Mag are same. 
% The speed of this function need to be improved.
 

    nos1 = event1.number_of_sac;
    nos2 = event2.number_of_sac;

    event_out = event1;
    
    data_sac = event1.sac;
    sac2 = event2.sac;
    sac2 = rmfield(sac2, 'EVEl');
    
    for i = (nos1 + 1):(nos1 + nos2)
        data_sac(i) = sac2(i - nos1);
    end
    event_out.sac = data_sac;
    event_out.number_of_sac = nos1 + nos2;
    





end

