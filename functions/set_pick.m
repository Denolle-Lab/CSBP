function event = set_pick(event , T , N)
% set the time of pick of event directly.there are some syntaxes.
% e = set_pick(e,t)
% set all traces in e to time t
% e = set_pick(e , T)
% set the ith trace to the ith element in vector T.size of T and e must be
% the same
% e = set_pick(e , t , n)
% set the n th trace to time t
% e = set_pick(e , T , n)
% set the traces in specified by vector N to time t
% e = set_pick(e , T , N)
% set the tract N to time T. number of elements must be the same.

nos = event.number_of_sac;

if nargin <= 1
    display('not enough input');
    return;
elseif nargin == 2
    [m1 , n1] = size(T);
    l1 = m1 * n1;
    N_temp = 1:nos;
    if (l1 == 1)
        T = T * ones(1 , nos);
        N_temp = 1:nos;
    elseif (l1 ~= nos)
        display('the element of T and number_of_sac must be the same');
        return;
    end;
elseif nargin > 3
    display('too many inputs');
    return;
else
    [m1 , n1] = size(N);
    l1 = m1 * n1;
    N_temp = N;
    [m2 , n2] = size(T);
    l2 = m2 * n2;
    
    if (l2 == 1)
        T = T * ones(1 , l1);
    elseif (l1 ~= l2)
        display('the number of elements of N and T must be the same');
        return;
    end;
end;

for i = N_temp
    event.sac(i).backup.fb = event.sac(i).first_break;
    event.sac(i).first_break = T(i);
end;
        
return;

end

