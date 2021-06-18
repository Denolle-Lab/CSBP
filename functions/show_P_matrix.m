function show_P_matrix( event, varargin)
%%show the data in RGB color
%  show_P_matrix( event, varargin)
%   varargin:  'time',[t1,t2]
%

%% data initialization
    nos = event.number_of_sac;
    time_range_flag = 0; % no time range inputted as default.
    
%% data input
    for i = 1:2:size(varargin, 2)
        if strcmp(varargin{i}, 'time')
            time_range = varargin{i + 1};
            time_range_flag = 1;
        end
    end
        
%% data process 
    if time_range_flag % time range inputted
        [M time_range p_range] = matrix_build(event, 'time', time_range); % matrix build.
    else        
        [M time_range p_range] = matrix_build(event); % matrix build.
    end
    
    % normalization
    
    normal_coef = max(abs(M'));
    M = M ./ (normal_coef' * ones(1, size(M, 2)));
    stack_wave = mean(M);
    
%% data display

    % colormap generated 
    rd=[(0:31)/31,ones(1,32)];
    gn=[(0:31)/31,(31:-1:0)/31];
    bl=[ones(1,32),(31:-1:0)/31];
    rwb=[rd',gn',bl'];
    cmaprwb=flipud(rwb);

    % display the wave in blue and red color
    clf;
    imagesc(time_range, [1 nos], M);
    colormap(cmaprwb);
    colorbar;
    set(gca,'ydir','normal','FontSize',15);
    ylim = get(gca,'YLim');
    line([p_range(1) p_range(1)], ylim, 'color', [0 0 0], 'Linewidth', 1);%[0 0 0] means black.
    line([p_range(2) p_range(2)], ylim, 'color', [0 0 0], 'Linewidth', 1);
    xlabel('time/s','FontSize',15);
    ylabel('Trace Number','FontSize',15);
    
    
    hold on;
    
    % display the stack wave
    x_coor = linspace(time_range(1), time_range(2), size(M(1, :), 2)); % x coordinate
    y_coor = mean(ylim) + (ylim(2) - mean(ylim)) * stack_wave / 2;
    plot(x_coor, y_coor, 'k-', 'LineWidth', 2);

    
end

    

    




