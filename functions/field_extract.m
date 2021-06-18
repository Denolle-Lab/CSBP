function out = field_extract(event, field )
%@function extract the specified field of event.
%if the field is a number, the out is a matrix
%if the field is a string, the out is a cell;
%@param @out. the output, it's zero if failed.
%@param @event. the event data ;
%@param @field. to specify the field;




    eval(strcat('temp = event.sac(1).', field, ';')); 
    nos = event.number_of_sac;
    
    if ischar(temp)
        out = cell(nos, 1);
        
        for i = 1:nos
            eval(strcat('temp = event.sac(i).', field, ';')); 
                % temp = event.sac(i).field;
            out{i} = temp;
        end
        
    elseif isnumeric(temp)
        out = zeros(nos, 1);
        
        for i = 1:nos
            eval(strcat('temp = event.sac(i).', field, ';'));
                % temp = event.sac(i).field;
            out(i) = temp;
        end
    end
    
end

