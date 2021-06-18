function [varargout] = matrix_build(event, varargin)
%%  pick the data in sacs to form a matrix
% Syntax
%
%       [M, time_range, P_wave] = matrix_build(event, 'aligned', 'on', 'number', num)
%       [M, time_range, P_wave] = matrix_build(event, 'aligned', 'off', 'number', num, 'time', time_range)
%
% @param @event. the structure need to be pick-up.
% @param @M. the out data matrix. 
% @param @time_range.The time range of data.
% @param @P_wave. The P_wave time range. If aligned is on, the P_wave is a
%                 two-elements vector. If aligned is off, the P_wave is a 
%                 number * 2 vector.
%@param @align. A string, 'on' or 'off'. a parameter decide the data was 
%                in time aligned or first_break aligned.
%@param @num_range. Number of sacs, a column vector.

%% data initialization
    nos = event.number_of_sac;
    dt = event.sac(1).dt;
    maxlength = max(field_extract(event,'length'));
    time_range = [0 maxlength * dt];
    num_range = [1 nos];
    aligned = 'on';
    time_range_flag = 0;% 0 for no time range inputted
  
%% data input    
    if nargin == 1
        aligned = 'on';
    else
        for i = 1:2:size(varargin, 2)
            if strcmp(varargin{i}, 'aligned')
                aligned = varargin{i + 1};
            elseif strcmp(varargin{i}, 'number')
                num_range = varargin{i + 1};
            elseif strcmp(varargin{i}, 'time')
                time_range = varargin{i + 1};
                time_range_flag = 1;
            end
        end
    end
    
    sac_num = num_range(2) - num_range(1) + 1;
    
%% data transform when align is off

    if strcmp(aligned, 'off')
        
        M = nan(sac_num, maxlength);
        
        for i = num_range(1):num_range(2)
            
            M(i, 1:size(event.sac(i).data, 2)) = event.sac(i).data;
            
        end
        fb = round(field_extract(event, 'first_break'));
        P_wave = [fb + event.p_time(1), fb + event.p_time(2)];
       
        M = M(:, round((time_range(1) / dt + 1):(time_range(2) / dt)));
        
    end
    
    
    
    
 
%% data transform when align is on        
    if strcmp(aligned, 'on')        
        
        if time_range_flag  % time range inputted   
            
            fb = field_extract(event, 'first_break'); % first break in seconds.
            
            left_range = time_range(1) - mean(fb);
            
            if (left_range + min(fb)) < 0 % the range is negative in the sac whose fb is minimum.
                left_range = 1 - min(fb);
                display('The left limit of time range is not proper.');
                display(['The left limit is modified to ' num2str(mean(fb) + left_range) ' seconds']);
            end
            
            right_range = time_range(2) - mean(fb);
            max_right_length = min(field_extract(event, 'length') * dt - fb);
                        
            if right_range > max_right_length
                right_range = max_right_length - 1;
                display('The right limit of time range is not proper.');
                display(['The right limit is modified to ' num2str(mean(fb) + right_range) ' seconds']);
            end
            
            time_range = [fb + left_range, fb + right_range]; %data extract time range.
        
        else
            
            P_length = event.p_time(2) - event.p_time(1); % the length of P wave in seconds.
            fb = field_extract(event, 'first_break'); % first break in seconds.
            fb = round(fb);
            left_range = min(2 * P_length, min(fb)); % to guarantee the left_range is above 0;
            max_right_length = min(field_extract(event, 'length') * dt - fb);
            right_range = min(10 * P_length, max_right_length);
            time_range = [fb - left_range, fb + right_range]; %the data extract time range for each sac.
        
        end
        
        data_range = round(time_range / dt) ;% the pick range in points.
        
        M = zeros(sac_num, data_range(1,2) - data_range(1, 1));
        data_range(:, 1) = data_range(:, 1) + 1;
        
        for i = 1:sac_num
            M(i, :) = event.sac(num_range(1) + i - 1).data((data_range(i, 1)) : (data_range(i, 2)));
        end
        
        time_range = mean(time_range);
        P_wave = mean([fb + event.p_time(1), fb + event.p_time(2)]);
        
    end
    
 %% data output
 
    varargout{1} = M;
    varargout{2} = time_range;
    varargout{3} = P_wave;
     
    
end
