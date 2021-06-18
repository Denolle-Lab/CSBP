function P = p_wave(event)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nos = event.number_of_sac;

p_length = round((event.p_time(2) - event.p_time(1)) / event.sac(1).dt);  %second to point
P = zeros(nos,p_length);

for i = 1:nos
    pill = round((event.sac(i).first_break + event.p_time(1)) / event.sac(i).dt);%the lower limit of the ith SAC
    piul = round((event.sac(i).first_break + event.p_time(2)) / event.sac(i).dt); % the upper limit
    P(i,:) = normal_lize(event.sac(i).data(pill:(piul - 1))); % the stacking
end;


return

