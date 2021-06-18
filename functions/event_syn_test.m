function event_syn = event_syn_test( event, source, wigle)
%% replae the original SAC with synthetic data
% the funtion is to use the travel time table to get the predicted time.
% the replace the data with wigle to get the synthetic data.
%@ param event. The event structure with travel time calculated.
%@ param source. A matrix containing the source information.
%               The matrix is a chart in the format of:
%                number      time(sec)    x_location(km)      y_location(km)      duration(sec)      
%                  1           10           2                     -2                  3 
%               number is the subevent order.
%               time is base on the predicted time.
%               x_location is east positive and west negative.
%               y_location is north positive and south negative.
%               duration is the duration of each subevent.
%@ param wigle. different wigles can be used to stand for different source.
%               BUT NOW, ONLY triangle wave is used.

    
%% set the sacs to zero

    nos = event.number_of_sac

    for i = 1:nos
        event.sac(i).data = event.sac(i).data * 0;
    end
    
%% the main part.

    source_num = size(source, 2)
        
    for i = 1:source_num
%% subevent information initialization       
        subevent_time_related = source(i, 2);
        subevent_x_location_related = source(i, 3);
        subevent_y_location_related = source(i, 4);
        subevent_duration = source(i, 5);
        
%% calculte the travel time.

        for j = 1:nos
            
            subevent_lon = event.sac(j).EVLO + subevent_x_location_related / 111;
            subevent_lat = event.sac(j).EVLA + subevent_y_location_related;
            
            TravelTimeFile = 'tt91.P';

            % read travel time table
            ftb = fopen(TravelTimeFile, 'r');
            if ftb == -1
                 display('no travel time table');
            end;
            temp = fgetl(ftb); % read first line
            nepidist = fscanf(ftb, '%d', 1); % number of epicenter distance in the table, second line
            nsrcdep = fscanf(ftb, '%d', 1); % number of depth in the table, second line
            temp = fgetl(ftb); 
            srcdep = fscanf(ftb, '%f', nsrcdep); % read third line: source depth vector
            tmptable = fscanf(ftb, '%f', [nsrcdep+1, nepidist]); % read travel time table
            ttable = tmptable(2:(nsrcdep+1),:);
            epidist = tmptable(1, :);
            fclose(ftb);
            %travel time table read end.
            
            % horizontal distance
            ellip = almanac('earth', 'ellipsoid');
            GCARC = dist(subevent_lat, subevent_lon, event.sac(j).STLA, event.sac(j).STLO, ellip, 'degrees');
            
            % vetical distance
            if isnan(event.sac(j).STDP)
                event.sac(j).STDP = 0;
            end
            
            
            
            
            
        
    end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        











end

