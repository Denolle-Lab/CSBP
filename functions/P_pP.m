function t = P_pP(event)
%inversion homework 2

event.sac(1).first_break = 300;
dt = event.sac(1).DELTA;
p1 = event.sac(1).data(round((300 / dt) : (303 / dt)));
p2 = event.sac(1).data(round((303 / dt) : (308 / dt)));

[c,lag] = xcorr(p1,p2);
offset = lag(maxxc(c));
t = (offset) * dt + 3;

end

