function event = field_fill( event, field, value )
%%contray to the field_extract function, fill in the specify field with
% given vector or number. 
% @param @event. the event need to be filled
% @param @field. the field you want to change.
% @param @value. the value to be filled in. the parameter can be a 
%                vector or number

    nos = event.number_of_sac;
    number = size(value, 1) * size(value, 2);
    
    
    if number == 1
        
        for i = 1:nos
            eval(strcat('event.sac(i).', field, '= value;'));
        end
        
    elseif number == nos
        
        for i = 1:nos
            eval(strcat('event.sac(i).', field, '= value(i);'));
        end
    
    else
        
        display('wrong input!');
        return;
        
    end
    
    
   

