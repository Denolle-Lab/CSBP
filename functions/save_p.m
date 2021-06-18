function save_p(event,filename) 
%save the data of P wave;
% @filename  @ the name of the data file
%
if nargin == 1
    filename = 'untitled1';
end

filename = strcat(filename,'.mat');

nos = event.number_of_sac;
p_length = (event.p_time(2) - event.p_time(1)) / event.sac(1).dt; % the length of p wave in points

data = struct('nos',nos,...
              'p_time',event.p_time,...
              'P',zeros(nos,p_length),...
              'fb',zeros(nos,1));
for i = 1:nos
    data.fb(i) = event.sac(i).first_break;
    ll = round((event.sac(i).first_break + event.p_time(1)) / event.sac(i).dt + 1);
    ul = round((event.sac(i).first_break + event.p_time(2)) / event.sac(i).dt);
    data.P(i,:) = event.sac(i).data(ll:ul);
    fb(i) = event.sac(i).first_break;
end;

save(filename,'data','-mat')

end

