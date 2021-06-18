function N = normal_lize(data)
% normalize data to [0 1];

N = (data - min(data)) / (max(data) - min(data)); 

