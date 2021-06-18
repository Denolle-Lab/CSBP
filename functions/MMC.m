% the algorithm of MMC

function event = MMC(event)

nos = event.number_of_sac;
delta_t = zeros(nos * (nos - 1) / 2 + 1, 1);
k = 1;
range = event.range;

for i = 1:(nos - 1)
     pill = round((event.sac(i).first_break + event.p_time(1))  / event.sac(i).dt); %the lower limit of the ith SAC
     piul = round((event.sac(i).first_break + event.p_time(2)) / event.sac(i).dt); % the upper limit
     p_i = event.sac(i).data(pill:piul); % the data of P wave in sac(i)
     for j = (i+1):nos
        pjll = round((event.sac(j).first_break - range) / event.sac(j).dt) ;%the lower limit of the jth SAC
        pjul = round((event.sac(j).first_break + range + event.p_time(2)) / event.sac(j).dt) ; % the data of P wave in sac(j)
        p_j = event.sac(j).data(pjll:pjul);
        [c,lag] = xcorr(p_i , p_j);
        delta_t(k) = lag(maxxc(c));
        delta_t(k) = delta_t(k) + (range + event.p_time(1)) / event.sac(i).dt;
        k = k + 1;
    end;
end;

delta_t(k) = 0;
A = zeros(k,nos);
k = 1;

for i = 1:(nos - 1)
    for j = (i + 1):nos
        A(k,i) = 1;
        A(k,j) = -1;
        k = k + 1;
    end;
end;
A(k,:) = 1; 

%revise by Qin Weize
temp = mean(delta_t);
delta_t = delta_t - temp; 

t_est = A' * delta_t / nos; %the solution of least square.
%t_est = t_est + temp;


for i = 1:nos
    event.sac(i).first_break = event.sac(i).first_break + t_est(i) * event.sac(i).dt;
end;

return

