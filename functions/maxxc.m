%% calculate the x coordinate of first maxium in data
%% maxxc the output result, double
%% data , vetor in line

function maxxc = maxxc(data)

[m,n] = size(data);
max = -inf;
xc = 0;
for i = 1:n
    if (data(i) - max) > 0.001
        max = data(i);
        xc = i;
    end;
end;

maxxc = xc;

return;