function event = test1(event)
%generate noised sin wave for test.

nos = event.number_of_sac;
b = linspace(0 , 2 * pi , 100);
a = sin(b) + sin(3 * b);  %the P wave need to be pick-up


for i = 1:nos
    ra = 0.4 * (2 * rand(1 , event.sac(i).length) - 1); %the noise
    rt = (50 * (2 * rand(1) - 1)); %the time shift
  
    r = round(150 / event.sac(i).dt + rt) ;
    
    event.sac(i).data = zeros(1 , event.sac(i).length);
    event.sac(i).data(r : (r + 99)) = a;
    event.sac(i).data = event.sac(i).data + ra;
end

return;