function  map( event, perproty)
% draw a map of the distribution of the stations and events.
% @param @perproty. a perproty string used to modify the perproty of plot. 

    nos = event.number_of_sac;
    
    if nargin == 1
        perproty = '+';
    end
      
    
    figure;
    hold on;
    
    plot(event.sac(1).EVLO, event.sac(1).EVLA,'o');
    
    for i = 1:nos
        plot(event.sac(i).STLO, event.sac(i).STLA, perproty);
    end
    
       


