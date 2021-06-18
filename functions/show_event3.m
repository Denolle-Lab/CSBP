function show_event3(event , varargin)
%% show the event in RGB color 
% Syntax
%
%      show_event3(event, 'number', num, 'time', time_range)
%
% @param @event. the structure need to be pick-up.
% @param @time_range.The time range of data.
%@param @num. Number of sacs, a column vector.   

%% data initialization

    num_range = [1 event.number_of_sac];
    dt = event.sac(1).dt;
    maxlength = max(field_extract(event,'length'))
    time_range = [0 maxlength * dt];

%% data input

    for i = 1:2:size(varargin, 2)
        if strcmp(varargin{i}, 'number')
            num_range = varargin{i + 1};
        elseif strcmp(varargin{i}, 'time')
            time_range = varargin{i + 1};
        end
    end

%% data process 

    [M, time_range, P_wave] = matrix_build(event, 'aligned', 'off', 'number', num_range, 'time', time_range);
    
    % normalization
    
    normal_coef = max(abs(M'));
    M = M ./ (normal_coef' * ones(1, size(M, 2)));
    
%% data display
    
    % colormap generated 
    rd=[(0:31)/31,ones(1,32)];
    gn=[(0:31)/31,(31:-1:0)/31];
    bl=[ones(1,32),(31:-1:0)/31];
    rwb=[rd',gn',bl'];
    cmaprwb=flipud(rwb);

    imagesc(time_range , num_range , M);
    colormap(rwb);
    colorbar;

end

