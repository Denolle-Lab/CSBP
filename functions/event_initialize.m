function event = event_initialize(event)
% some information useful for next caculation 
if isempty(event)
    display(['***** Null event']);
    return
end;

if ~isstruct(event)
    display(['***** Wrong type']);
    return
end;

% waiting for adding
[m,n] = size(event.sac);
backup.range = 0;
event = struct('number_of_sac',m,...
               'p_time',[0 20],...  % the range of P wave, in seconds
               'range',50,...        % search range in seconds
               'backup',backup,...  %bauk up data 
               'order_field','',...
               'sac',event.sac);

           
           
for i = 1:event.number_of_sac  ;
    
    if isnan(event.sac(i).E)
        event.sac(i).E = event.sac(i).B + (event.sac(i).NPTS - 1) * event.sac(i).DELTA;
    end
    
    if isnan(event.sac(i).AZ)
        
        x = event.sac(i).STLO - event.sac(i).EVLO;
        y = event.sac(i).STLA - event.sac(i).EVLA;
        thita = rad2deg(acos([x y] * [0 ; 1] / distance(0, 0, y, x)));
        
        if x < 0
            thita = 360 - thita;
        end
        
        event.sac(i).AZ = thita;
    else
        thita=event.sac(i).AZ; % Yin revised
        
    end    
    
    if isnan(event.sac(i).BAZ)
        event.sac(i).BAZ = mod(thita + 180, 360);
    end
    
    if isnan(event.sac(i).GCARC)
        event.sac(i).GCARC = distance(event.sac(i).STLA, event.sac(i).STLO, ...
                                    event.sac(i).EVLA, event.sac(i).EVLO);
    end
       
    if isnan(event.sac(i).DIST)
        event.sac(i).DIST = event.sac(i).GCARC / 180 * pi * 6387;
    end
    
    if isnan(event.sac(i).DEPMIN)
        event.sac(i).DEPMIN = min(event.sac(i).DATA1);
    end
    
    if isnan(event.sac(i).DEPMAX)
        event.sac(i).DEPMAX = max(event.sac(i).DATA1);
    end
    
    event.sac(i).data = event.sac(i).DATA1;
    event.sac(i).dt = event.sac(i).DELTA;
    event.sac(i).first_break = 150;
    event.sac(i).length = event.sac(i).NPTS;    % length of each sac
    event.sac(i).data = reshape(event.sac(i).data,1,[]);
    event.sac(i).time_coordinate = linspace(event.sac(i).B,event.sac(i).E,event.sac(i).length);% time_coordinate
    event.sac(i).backup.fb = 0;%backup data fisrt_break;
    event.sac(i).logscale = 0;
    event.sac(i).pol = 1; % polarity;
end;
  
return;



 


