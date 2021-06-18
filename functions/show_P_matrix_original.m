function show_P_matrix(varargin)
%show the data in RGB color
%

event = varargin{1};

if nargin == 1
    start = 1;
    ending = event.number_of_sac;
    display = '';
elseif nargin == 2
    start = 1;
    ending = event.number_of_sac;
    display = varargin{2};
elseif nargin == 3
    start = varargin{2};
    ending = varargin{3};
    display = varargin{4};
else 
    diplay('too many inputs');
    return;
end;




range = 5 * (event.p_time(2) - event.p_time(1)); %display range in second
N = matrix_build(event);
M = N;
nos = event.number_of_sac;
fb = zeros(1,nos);   

for i = 1:nos
    fb(i) = event.sac(i).first_break;
end;

aver = mean(fb);
lowerline = round(aver / event.sac(1).dt); % the line to mark the lower limit of P
upperline = lowerline + ...
            round((event.p_time(2) - event.p_time(1)) / event.sac(1).dt);

for i = 1:nos
    fb(i) = event.sac(i).first_break - aver; % fb is now the time-shift vector
end;

fb = round(fb / event.sac(1).dt); %transform fb from seconds to points.



max0 = min(fb);

minll = 6000;  
maxul = 0;   %the display range in points in X;
 
for i = start:ending
    ll = round((event.sac(i).first_break - range) / event.sac(i).dt);  % lower limit in point
    ll = round(max(ll , -max0)) + 1; % to guarantee the * will above zero
    ul = round((event.sac(i).first_break + range) / event.sac(i).dt);  % upper limit in point
    N(i , ll : ul) = M(i , ((ll : ul) + fb(i)));   %time shift *
    
    if ll < minll
        minll = ll;
    end;
    
    if ul > maxul
        maxul = ul;
    end; 

end;   

for i = start:ending
    N(i, minll : maxul) = N(i, minll : maxul) / max(abs(N(i, minll : maxul)));
end;

N = N * 128 + 128;

if strcmp(display, 'order_field_on') 
    axes('position',[.9, 0.05, .05, 0.9]);
    set(gca, 'XTick', [], 'YTick', []);
    set(gca,'ydir','reverse');
    color_data = field_extract(event, event.order_field);
    color_data = normal_lize(color_data) * 12;
    image(color_data);
    colormap(gray(12));
    
    
    
    axes('position',[.1, 0.05, 0.75, 0.9]);

    image(([minll maxul] * event.sac(1).dt) , [1 nos], N(: , (minll : maxul)));
    colormap(jet(256));
    set(gca,'ydir','reverse');
    ylim = get(gca,'YLim');
    line([aver,aver] + event.p_time(1),ylim,'marker','>');
    line([aver,aver] + event.p_time(2),ylim,'marker','<');

else 
    
    image(([minll maxul] * event.sac(1).dt) , [1 nos], N(: , (minll : maxul)));
    colormap(jet(256));
    set(gca,'ydir','reverse');
    ylim = get(gca,'YLim');
    line([aver,aver] + event.p_time(1),ylim,'marker','>');
    line([aver,aver] + event.p_time(2),ylim,'marker','<');

end






end




