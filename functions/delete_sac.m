function event = delete_sac(event , N)
% delete some traces, two stynax.
% event = delete_sac(event)
% delete the trace whose corrcoef of P_wave with the stacked P_wave is 
% less than the threshold.
% event = delete_sac(event , N)
% delete the sacs specified by elements in N

if nargin == 1
    r = corr_coef(event);
    threshold = 0.5;
    for i = event.number_of_sac:-1:1
        if (abs(r(i) < threshold))
            event.sac(i) = [];
            event.number_of_sac = event.number_of_sac - 1; 
        end;
    end;
else
    
    if N < 1
        r = corr_coef(event);
        threshold = N;
        for i = event.number_of_sac:-1:1
            if (abs(r(i) < threshold))
                event.sac(i) = [];
                event.number_of_sac = event.number_of_sac - 1; 
            end;
        end;
    else    
        N = sort(N,'descend');
        
        for i = N
            event.sac(i) = [];            
        end;
        
        event.number_of_sac = event.number_of_sac - numel(N);
    end
end;
        
return;

