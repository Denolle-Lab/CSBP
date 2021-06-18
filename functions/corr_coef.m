function r = corr_coef(event)
% calculate the corrcoef of P_wave with the average of correlation
% coeffient of other waves.
%@param @r a vector in column.

nos = event.number_of_sac;

P = p_wave(event);
S = sum(P);
r = zeros(nos , 1);%the output vector

for i = 1:nos
    
    temp = corrcoef(P(i,:), S);
    r(i) = temp(2);
end;

return 

