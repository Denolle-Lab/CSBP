function event = sac_selected(event,varargin)
%% Select the sacs with different locations.
% 
%   Syntax 
%
%       event_out = sac_selected(event, property, range, ...)
%
%   Description
%
%       @para @property. The parameters to select the SACs. There are five 
%                        properties now i.e. AZ, GCARC, LAT, LON.
%       @para @range.The range four each parameter, a vector with 2
%                    elements.

%%   data initialization

    if mod(size(varargin, 2), 2) == 1
        display('invalid input');
        return;
    end

    AZ_range = [0 360];
    BAZ_range = [0 360];
    GCARC_range = [0 360];
    LAT_range = [0 360];
    LON_range = [0 360];   
    select = ones(event.number_of_sac, 1);    
    

    
 %% data selection   
    
    for i = 1:2:size(varargin, 2)
        
        if strcmp(varargin{i}, 'AZ')
            
            AZ_range = varargin{i + 1};
            AZ = field_extract(event, 'AZ');
            AZ_select = (AZ >= AZ_range(1)) & (AZ <= AZ_range(2));
            select = select & AZ_select;            
            
        elseif strcmp(varargin{i}, 'GCARC')
        
            GCARC_range = varargin{i + 1};
            GCARC = field_extract(event, 'GCARC');
            GCARC_select = (GCARC >= GCARC_range(1)) & (AZ <= GCARC_range(2));
            select = select & GCARC_select;            
            
        elseif strcmp(varargin{i}, 'BAZ')
        
            BAZ_range = varargin{i + 1};
            BAZ = field_extract(event, 'BAZ');
            BAZ_select = (BAZ >= BAZ_range(1)) & (BAZ <= BAZ_range(2));
            select = select & BAZ_select;       
            
        elseif strcmp(varargin{i}, 'LAN')
        
            LAT_range = varargin{i + 1}; 
            LAT = field_extract(event, 'STLA');
            LAT_select = (LAT >= LAT_range(1)) & (LAT <= LAT_range(2));
            select = select & LAT_select;      
            
        elseif strcmp(varargin{i}, 'LON')
        
            LON_range = varargin{i + 1};
            LON = field_extract(event, 'STLO');
            LON_select = (LON >= LAT_range(1)) & (LON <= LAT_range(2));
            select = select & LON_select;      
        
        end
        
    end
    
    delete = ~select; 
    data_number = (1:event.number_of_sac)';
    delete_number = data_number(delete);
    
    event = delete_sac(event, delete_number);


end

