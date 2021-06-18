function event = duration_check( event)
% check the duration of waveform
%  

    length = field_extract(event, 'length');
    del_index = find(length < 0.5 * mean(length));
    
    event = delete_sac(event, del_index);
    
    
end

