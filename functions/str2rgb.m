function rgb = str2rgb( str )
%% transform a str to a 1*3 vector standing for rgb color  
%   now the color in the color table is 
%
%     color                string
%     red                  'r' or 'red'
%     green                'g' or 'green'
%     blue                 'b' or 'blue'
%     cyan                 'c' or 'cyan'
%     magenta              'm', 'mag' or 'magenta'
%     yellow               'y' or 'yellow'
%     white                'w' or 'white'
%     black                'k' or 'black'
%  
% more colors can be added here later.
    
    
    

    if (strcmp(str, 'r') | strcmp(str, 'red'))
        rgb = [1 0 0];
    elseif (strcmp(str, 'g') | strcmp(str, 'green'))
        rgb = [0 1 0];
    elseif (strcmp(str, 'b') | strcmp(str, 'blue'))
        rgb = [0 0 1];    
    elseif (strcmp(str, 'c') | strcmp(str, 'cyan'))
        rgb = [0 1 1];
    elseif (strcmp(str, 'm') | strcmp(str, 'mag') | strcmp(str, 'magenta'))
        rgb = [1 0 1];
    elseif (strcmp(str, 'y') | strcmp(str, 'yellow'))
        rgb = [1 1 0];
    elseif (strcmp(str, 'w') | strcmp(str, 'white'))
        rgb = [1 1 1];
    elseif (strcmp(str, 'k') | strcmp(str, 'black'))
        rgb = [0 0 0];
    else
        rgb = 0;
        display('Unkown color');
    end
            
        
        
end

