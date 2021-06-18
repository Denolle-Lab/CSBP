function cm = colormap_maker(n, varargin)
%% gernerate color map using linear interpolation.
% sytax
%   
%   cm = colormap_maker(n, color1, color2, color3...colorn)
%
% the color is a 1 * 3 rgb vector.
% n is the number of colors in colormap. 

    
    
%% data input 

    % the normal colormap used for velocity model.
    % more colormaps can be added here.
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
%   the process is linear interpolation.
%   the basic process is 
%   1. generate a large color table with lines of least common multiple.
%   2. set value of the lines specified by color serie to the color value
%   3. interpolate the lager table to the right value.
%   4. extract the colormap in the large table.
    
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

