function [event,true_pick,tp2] = test2(event)
%generate data for p-search testing
  
nos = event.number_of_sac;
t_range = 5;
b = linspace(0 , 2 * pi , 50);
a = sin(b) + sin(3 * b) + sin(4 * b);  %the P wave need to be pick-up
b = linspace(0 , 2 * pi , 100);
b = cos(2 * b)+ cos(6 * b) + cos(7 * b); 
tp = zeros(nos,1);
tp2 = zeros(nos,1);

for i = 1:nos
    ra = 0.2 * (2 * rand(1 , event.sac(i).length) - 1); %the noise
    rt1 = (t_range * (2 * rand(1) - 1)); %the time shift
    rt2 = (t_range * (2 * rand(1) - 1));
    r1 = round((150 - 2 * t_range + rt1) / event.sac(i).dt) ;%the location of fisrt P
    r2 = round((150 + t_range + rt2) / event.sac(i).dt) ;%the location of main P
    event.sac(i).data = zeros(1 , event.sac(i).length);
    event.sac(i).data(r1 : (r1 + 49)) = a;
    event.sac(i).data(r2 : (r2 + 99)) = 3*b;
    
    event.sac(i).data = event.sac(i).data + ra;
    tp(i) = r1;
    tp2(i) = r2;
end;

true_pick = tp * event.sac(1).dt;
tp2 = tp2 * event.sac(1).dt;

return;

