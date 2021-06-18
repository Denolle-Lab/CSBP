function fb = first_break( event )
%extract the first_break as an vector in column

nos = event.number_of_sac;
fb = zeros(nos,1);

for i = 1:nos
    fb(i) = event.sac(i).first_break;
end;

end

