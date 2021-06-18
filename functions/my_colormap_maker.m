function cm = colormap_maker(n, varargin)
%% gernerate color map using linear interpretation.
% sytax
%   
%   cm = m
%
% Now the color must be a RGB color like [0 0 0];
% the N decide how many color make from color1 to color3 passing color2.

    
    
%% data input 

    % the normal colormap used for velocity model.
    if strcmp(varargin{1}, 'velocity')
        varargin = {'r', 'y', 'w', 'c', 'b'};
    end
    
    color_num = size(varargin, 2);
    color_sequence = zeros(color_num, 3);
    
    for i = 1:color_num
        if isnumeric(varargin{i})
            color_sequence(i, :) = varargin{i};
        elseif ischar(varargin{i})
            color_sequence(i, :) = str2rgb(varargin{i});
        end
    end
    
 %% data process
 % the process will be written later.
    
    line = lcm(color_num - 1, n - 1) + 1; % the lines of the large color table
    color_large_table = zeros(line, 3);
    
    color_interval = (line - 1) / (color_num - 1);
    colormap_interval = (line - 1) / (n - 1);
    
    color_index = 1:color_interval:line;
    colormap_index = 1:colormap_interval:line;
    
    color_large_table(color_index, :) = color_sequence;
    
    for i = 1:3
        color_large_table(:, i) = interp1(color_index, color_sequence(:, i), 1:line);
    end
    
    cm = color_large_table(colormap_index, :);
    
    
    

      
end

