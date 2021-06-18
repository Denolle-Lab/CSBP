%% auto pick using STA/LTA

function event = auto_pick1(event)

range = event.range ; % searching range in seconds
ns = 40 ; % short-term in points
nl = 400 ;% long-term in points

for i = 1:event.number_of_sac
    bl = round((event.sac(i).first_break - range) / event.sac(i).dt);
    bh = round((event.sac(i).first_break + range) / event.sac(i).dt);
    r = zeros(1,bh);
    d = zeros(1,bh);
    for j = bl:1:bh
        a = event.sac(i).data((j - ns):j);  
        STA = a * a' / ns;
        b = event.sac(i).data((j - nl):j);
        LTA = b * b' / nl;
        if (LTA - 0) < 0.0001
            r(j) = 0;
        else
            r(j) = STA / LTA;
        end
        d(j - 1) = r(j) - r(j - 1);
    end;
    c = d((bl + 2):(bh - 2));
    maxx = maxxc(c) + bl + 1;
    event.sac(i).backup.fb = event.sac(i).first_break;
    event.sac(i).first_break = maxx * event.sac(i).dt;
end;

event.backup.range = event.range;
event.range = 20;

end