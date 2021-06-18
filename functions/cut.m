function event = cut(event, evdata, sac_time)
%@function. Cut the origin sac data base on the specify event data.
%           The cutting range was 300s before the P arrival and 600s
%           after the P arrival.
%           only origin about one-hour-long data can be used. 
%           travel time table is needed.
%@out @event. The event after cutting.
%@out @n. The number of traces that contain the right data of the earthquake.
%@param @event. The event is to be cut.
%@param @evdata. The structure contain the data of the earthquake.
%                The fields of the evdata are: latitude, longtude, elevation,
%                depth, year, julian_day, hour, minute, second.
%                If the evdata is 0, the earthquake data is in the event already. 
%@param @sac_time. The original starting time of the sac, in form of yyyymmddhh.
%                  For example. 2011052413 means 13:00 on May 24th in 2011.

 
        
    if isstruct(evdata)        
        event = field_fill(event, 'EVLA', evdata.latitude);
        event = field_fill(event, 'EVLO', evdata.longtitude);
        event = field_fill(event, 'EVEL', evdata.elevation);
        event = field_fill(event, 'EVDP', evdata.depth);
        event = field_fill(event, 'EVLO', evdata.latitude);
        event = field_fill(event, 'NZYEAR', evdata.year);
        event = field_fill(event, 'NZJDAY', evdata.julian_day);
        event = field_fill(event, 'NZHOUR', evdata.hour);
        event = field_fill(event, 'NZMIN', evdata.minute);
        event = field_fill(event, 'NZSEC', evdata.second);        
    end
    
    sac_year = floor(sac_time / 1000000);
    sac_jday = julian_day(floor(sac_time / 100));
    sac_hour = mod(sac_time, 100);
    
    
    if (event.sac(1).NZYEAR ~= sac_year)        
        display('the earthquake is not in the data');
        return;    
    elseif (event.sac(1).NZJDAY ~= sac_jday)
        display('the earthquake is not in the data');
        return;   
    elseif (event.sac(1).NZHOUR ~= sac_hour)
        display('the earthquake is not in the data');
        return;  
    end
    
    length = field_extract(event, 'E') - field_extract(event, 'B');
    delete_num = find(length < 3600);
    event = delete_sac(event, delete_num);
    
    
    
    event = travel_time_table(event);
    tt_time = field_extract(event, 'tt_time');
    nzmin = field_extract(event, 'NZMIN');
    nzsec = field_extract(event, 'NZSEC');
    dt = field_extract(event, 'dt');
    
    arival_time = tt_time + nzmin * 60 + nzsec;
    aver = mean(arival_time);
    
    r = abs((arival_time - aver) / aver);
    delete_num = find(r > 0.1);
    event = delete_sac(event, delete_num);
    
  
    event = travel_time_table(event);
    tt_time = field_extract(event, 'tt_time');
    nzmin = field_extract(event, 'NZMIN');
    nzsec = field_extract(event, 'NZSEC');
    dt = field_extract(event, 'dt');    
    arival_time = tt_time + nzmin * 60 + nzsec;
    
    time_ll = arival_time - 300; % lower limit of cutting time in seconds. 
    time_ul = arival_time + 600; % upper limit of cutting time in seconds.
    
    ll = round(time_ll ./ dt); % translate seconds to points;
    ul = round(time_ul ./ dt); % translate seconds to points;

    for i = 1:event.number_of_sac
        temp = zeros(size(event.sac(i).DATA1(ll(i):ul(i))));
        temp = event.sac(i).DATA1(ll(i):ul(i));
        event.sac(i).DATA1 = [];
        event.sac(i).DATA1 = temp;
        event.sac(i).B = 0;
        event.sac(i).O = 300 - tt_time(i);
        event.sac(i).E = 900;
        event.sac(i).NPTS = 900 / event.sac(i).dt + 1;
    end
    
    return;
        
        













