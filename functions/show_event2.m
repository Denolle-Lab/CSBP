function show_event2(event)
%UNTITLED show the event in the way that firstbreak is aligned
%   Detailed explanation goes here

if nargin == 1
    start = 1;
    ending = event.number_of_sac;
end;

nos = event.number_of_sac;
fb = zeros(1,nos);

for i = 1:nos
    fb(i) = event.sac(i).first_break;
end;

aver = mean(fb);

for i = 1:nos
    fb(i) = event.sac(i).first_break - aver;
end;


for i = start:ending
    plot(event.sac(i).time_coordinate - fb(i),(normalize(event.sac(i).data) + 1+ 2 * (i-1)));
    hold on;
end;

line([aver,aver],[0.2,2 * nos - 0.2]);

end

