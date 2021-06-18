function show_P_matrix( event )
%%show the data in RGB color
%

%% data initialization
    nos = event.number_of_sac;
    
    
%% data process 

    [M time_range p_range] = matrix_build(event); % matrix build.
    
    % normalization
    
    normal_coef = max(abs(M'));
    M = M ./ (normal_coef' * ones(1, size(M, 2)));
    
    M = M * 128 + 128;% projected to RGB sapce
    
    
%% data display

    image(time_range, [1 nos], M);
    colormap(jet(256));
    set(gca,'ydir','reverse');
    ylim = get(gca,'YLim');
    line([p_range(1) p_range(1)], ylim, 'marker', '>');
    line([p_range(2) p_range(2)], ylim, 'marker', '<');
    
end

    

    




