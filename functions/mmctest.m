function event = mmctest(event)
%generate data for mmc test


nos = event.number_of_sac;
b = linspace(0 , 2 * pi , 100);
b = sin(b);

for i = 1:nos
    
    event.sac(i).data = zeros(1 , event.sac(i).length);
    p = round((160 + 5 * (0 * rand(1) - 1))/ event.sac(i).dt) + i;
    event.sac(i).data(p : (p + 99)) = b;
    
end

event = event;

return;

