function event = time_cut(event,start_time, end_time )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   event = time_cut(event,start_time, end_time )

% Revised by Yin

   for i = 1:event.number_of_sac
       

       start_numb = round(start_time / event.sac(i).dt);
       end_numb = round(end_time / event.sac(i).dt);
       
       event.sac(i).data(1:start_numb) = [];
       event.sac(i).data(end_numb+1:end) = [];
       
       
       event.sac(i).time_coordinate(1:start_numb) = [];
       event.sac(i).time_coordinate(end_numb+1:end) = [];
       
       event.sac(i).length = size(event.sac(i).data, 2);
       
   end
       





end

