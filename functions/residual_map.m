function residual_map(event)
%show the residual map of ak135.

nos = event.number_of_sac;
X = zeros(1,nos);
Y = X;
H = Y;

for i = 1:nos
    X(i) = event.sac(i).STLO;%longtitude of the station of the ith sac
    Y(i) = event.sac(i).STLA;%latitude of the station of the ith sac
    H(i) = event.sac(i).first_break - event.sac(i).tt_time;
end;

[x,y]=meshgrid(min(X):0.1:max(X),min(Y):0.1:max(Y));

z=griddata(X,Y,H,x,y,'linear');
% figure(1)
% 
% 
% meshc(x,y,z),rotate3d
% xlabel('X'),ylabel('Y'),zlabel('H'),title('H')


figure(4)

[c,h] = contour(x,y,z,50);xlabel('X'),ylabel('Y'),title('H.linear');
clabel(c,h);

hold on,grid on,grid minor,plot(X,Y,'*'),hold off

% figure(3)
% 
% contour3(x,y,z,50);xlabel('X'),ylabel('Y'),zlabel('H'),title('')


