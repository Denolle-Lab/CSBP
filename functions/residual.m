function r = residual(event)
%return the risiduals in seconds and in vector
%  

nos = event.number_of_sac;
r = zeros(nos,1);

for i = 1:nos
    r(i) = event.sac(i).first_break - event.sac(i).tt_time;
end;

end

