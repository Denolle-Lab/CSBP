function aftershock_distribution_3D(epicenter, depth, aftershock )
%UNTITLED2 Summary of this function goes here
%  

    hold on;
    plot3(epicenter(1), epicenter(2), depth, 'pr', 'MarkerSize', 20);
    
    for i = 1:size(aftershock, 1)
        
        plot3(aftershock(i, 2), aftershock(i, 1), - aftershock(i, 3), '.k',...
            'MarkerSize', aftershock(i, 4) / 8.3 * 30);

    end

%     set(gca, 'Zdir', 'reverse');
    %set(gca, 'Ydir', 'reverse');
    
end

