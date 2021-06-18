function event_out = event_sort( event, field )
% @function. Sorting the event according to the specified field in ascend.
%   syntax
%              event = event(event, field)
%
% @param @event. The event to be process.
% @param @field. Specifying the field for sorting

    % check whether the field is a number.

    eval(['temp = event.sac(1).' field ';']); 
        % temp = event.sac(1).field;
       
    if ~(isnumeric(temp))

        display('The specified field is not numeric.');
        return;
    
    end
    
    if strcmp(field, event.order_field)
        event_out = event;
        return;
    else
        event.order_field = field;
    end
     
    % sort
    
    sort_data = field_extract(event, field);
    [sort_data, index] = sort(sort_data);
    
    
    % generate new event.
    
    event_out = struct(event);
    
    % event_out.sac = [];
    
    for i = 1:event.number_of_sac
        
        event_out.sac(i) = struct(event.sac(index(i)));
        
    end
    
end

