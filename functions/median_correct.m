function event = median_correct(event)
%revise the first break time base on 
%the median of the time of stations 500kms or nearer.


    nos = event.number_of_sac;
    D = zeros(nos,nos);
    dl = 100; %distance limit.
    circle = earthRadius('km') * pi / 180;
    
    for i = 1:nos
        for j = (i+1):nos
            s1.x = event.sac(i).STLA;
            s1.y = event.sac(i).STLO;
            s2.x = event.sac(j).STLA;
            s2.y = event.sac(j).STLO;
            d = distance(s1.x,s1.y,s2.x,s2.y);
            d = d * circle;
            D(i,j) = d;
        end
    end
    
    D = D + D';
    
    to = first_break(event); % old first break time
    tn = to; % new time;
    
    for i = 1:nos
        index = find(D(i,:) < dl);
        tn(i) = median(tn(index));
    end
    
    event = set_pick(event,tn);
    
    
    
        