function show_P(event,start,ending) %need updating

if nargin == 1  
    start = 1;
    ending = event.number_of_sac;
end;
 
range = 5 * (event.p_time(2) - event.p_time(1)); %display range
nos = event.number_of_sac;
sp = 1;%sac separation distance
Amp = 2;

fb = first_break(event);


aver = mean(fb);
fb = fb - aver;


max0 = min(fb);

% fb = round(fb / event.sac(1).dt);


for i = start:ending

    ll = round((event.sac(i).first_break - range) / event.sac(i).dt);
    ll = round(max(ll , -max0) + 1);
    ul = round((event.sac(i).first_break + range) / event.sac(i).dt);
    plot(event.sac(i).time_coordinate(ll:ul) - fb(i),...
         (Amp * normal_lize(event.sac(i).data(ll:ul)) + 1 + i - start));
    hold on;

end;

axis([aver - range , aver + range , -1 , 2 + sp * (ending - start - 1)] )

line([aver aver],[0.2  (ending - start + 1) * sp - 0.2]);
line([aver aver] + event.p_time(1),[0.2  (ending - start + 1) * sp - 0.2]);
line([aver aver] + event.p_time(2),[0.2  (ending - start + 1) * sp - 0.2]);

return

