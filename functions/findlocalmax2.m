function [idx, locmax, count] = findlocalmax2(data, nmax, neighborsize, irowcirc, icolcirc)

% function to find local maximum values of a input matrix of at least 3x3
%% input:
% data: input data matrix
% nmax: number of local maximums to be found
% neighborsize: [m n] vector: mute the local maximums around a larger local
%               maximum within m points on the row and n points on the
%               column
% irowcirc: index for circular rows (e.g., angles 0 and 359 are close to
%           each other). = 1, for circular rows; = 0, for non-circular rows
% icolcirc: index for circular columns. = 1: circular; = 0: non-circular

%% output:
% idx: 2 x count size array storing the position of local maximums. the
%      first row for row index of loc. max., the second row for column
%      index of loc. max.
% locmax: value of local maximums corresponding to idx position
% count: number of return local maximums. count may be less than nmax, the
%        requested number of local maximums, if the number of loc.max. is
%        less than nmax

%%
if nmax <= 0
    nmax = 1;
end

m = size(data,1);
n = size(data,2);

if m < 3 || n < 3  % matrix size less than 3x3
    idx = [NaN; NaN];
    locmax = NaN;
    count = NaN;
else
    dm0np1 = zeros(m,n);
    dm0np1(1:m,1:n-1) = data(1:m,2:n) - data(1:m,1:n-1);
    dm0ns1 = zeros(m,n);
    dm0ns1(1:m,2:n) = -dm0np1(1:m,1:n-1);
    
    dmp1n0 = zeros(m,n);
    dmp1n0(1:m-1,1:n) = data(2:m,1:n) - data(1:m-1,1:n);
    dms1n0 = zeros(m,n);
    dms1n0(2:m,1:n) = -(dmp1n0(1:m-1,1:n));
    
    dmp1np1 = zeros(m,n);
    dmp1np1(1:m-1,1:n-1) = data(2:m,2:n) - data(1:m-1,1:n-1);
    dms1ns1 = zeros(m,n);
    dms1ns1(2:m,2:n) = -dmp1np1(1:m-1,1:n-1);
    
    dmp1ns1 = zeros(m,n);
    dmp1ns1(1:m-1,2:n) = data(2:m,1:n-1) - data(1:m-1,2:n);
    dms1np1 = zeros(m,n);
    dms1np1(2:m,1:n-1) = -dmp1ns1(1:m-1,2:n);
    
    pattern = sign(dm0np1) + sign(dm0ns1) + sign(dmp1n0) + sign(dms1n0) + ...
        sign(dmp1np1) + sign(dms1ns1) + sign(dmp1ns1) + sign(dms1np1);
    
    IIlocmax = (pattern == -8);
    
    if icolcirc == 1
        II = find(pattern(2:m-1,1) == -5);
        IIlocmax(II+1,1) = 1; clear II;
        II = find(pattern(2:m-1,n) == -5);
        IIlocmax(II+1,n) = 1; clear II;
    end
    
    if irowcirc == 1
        II = find(pattern(1,2:n-1) == -5);
        IIlocmax(1,II+1) = 1; clear II;
        II = find(pattern(m,2:n-1) == -5);
        IIlocmax(m,II+1) = 1; clear II;
    end
    
%     if pattern(1,1) == -3, IIlocmax(1,1) = 1; end
%     if pattern(1,n) == -3, IIlocmax(1,n) = 1; end
%     if pattern(m,1) == -3, IIlocmax(m,1) = 1; end
%     if pattern(m,n) == -3, IIlocmax(m,n) = 1; end
    
    
    % the rows are circularly connected (e.g., angles)
    if irowcirc == 1
        II = find(IIlocmax(1,1:n) == 1); % for first row local max
        for i = 1:length(II)
            if pattern(1,II(i)) == -5 % for first row non-corner local max
                if length(find(data(m,II(i)-1:II(i)+1) > data(1,II(i)))) > 0
                    IIlocmax(1,II(i)) = 0;
                else
                    IIlocmax(m,II(i)-1:II(i)+1) = 0;
                end
            else  % for first row two corner local max points
                if II(i) == 1  % for (1,1) point
                    if data(1,1) < data(m,1) || data(1,1) < data(m,2)
                        IIlocmax(1,1) = 0;
                    end
                elseif II(i) == n % for (1,n) point
                    if data(1,n) < data(m,n) || data(1,n) < data(m,n-1)
                        IIlocmax(1,n) = 0;
                    end
                end
            end
        end
        
        II = find(IIlocmax(m,1:n) == 1); % for last row local max
        for i = 1:length(II)
            if pattern(m,II(i)) == -5 % for last row non-corner local max
                if length(find(data(1,II(i)-1:II(i)+1) > data(m,II(i)))) > 0
                    IIlocmax(m,II(i)) = 0;
                else
                    IIlocmax(1,II(i)-1:II(i)+1) = 0;
                end
            else  % for last row two corner local max points
                if II(i) == 1  % for (m,1) point
                    if data(m,1) < data(1,1) || data(m,1) < data(1,2)
                        IIlocmax(m,1) = 0;
                    end
                elseif II(i) == n % for (m,n) point
                    if data(m,n) < data(1,n) || data(m,n) < data(1,n-1)
                        IIlocmax(m,n) = 0;
                    end
                end
            end
        end
        
    end
    

    % the columns are circularly connected (e.g., angles)
    if icolcirc == 1
        II = find(IIlocmax(1:m,1) == 1); % for first column local max 
        for i = 1:length(II)
            if pattern(II(i),1) == -5 % for first column non-corner local max
                if length(find(data(II(i)-1:II(i)+1,n) > data(II(i),1))) > 0
                    IIlocmax(II(i),1) = 0;
                else
                    IIlocmax(II(i)-1:II(i)+1,n) = 0;
                end
            else  % for first column two corner local max points
                if II(i) == 1  % for (1,1) point
                    if data(1,1) < data(1,n) || data(1,1) < data(2,n)
                        IIlocmax(1,1) = 0;
                    end
                elseif II(i) == m % for (m,1) point
                    if data(m,1) < data(m,n) || data(m,1) < data(m-1,n)
                        IIlocmax(m,1) = 0;
                    end
                end
            end
        end
        
        II = find(IIlocmax(1:m,n) == 1); % for last column local max
        for i = 1:length(II)
            if pattern(II(i),n) == -5 % for last column non-corner local max
                if length(find(data(II(i)-1:II(i)+1,1) > data(II(i),n))) > 0
                    IIlocmax(II(i),n) = 0;
                else
                    IIlocmax(II(i)-1:II(i)+1,1) = 0;
                end
            else  % for last column two corner local max points
                if II(i) == 1  % for (1,n) point
                    if data(1,n) < data(1,1) || data(1,n) < data(2,1)
                        IIlocmax(1,n) = 0;
                    end
                elseif II(i) == m % for (m,n) point
                    if data(m,n) < data(m,1) || data(m,n) < data(m-1,1)
                        IIlocmax(m,n) = 0;
                    end
                end
            end
        end
        
    end
    
    
    locmaxdata = zeros(m,n);
    for i = 1:m
        locmaxdata(i,1:n) = data(i,1:n).*IIlocmax(i,1:n);
    end
%    figure; subplot(1,2,1); imagesc(data); 
%    subplot(1,2,2); imagesc(pattern);  
%    subplot(1,2,2); imagesc(locmaxdata);
    
    midx = [];
    nidx = [];
    alllocmax = [];
    
    % find locations of local maximum
    for i = 1:m
        II = find(IIlocmax(i,:) == 1);
        nII = length(II);
        if nII > 0
            kk = length(midx);
            midx(kk+1:kk+nII) = ones(1,nII)*i;
            nidx(kk+1:kk+nII) = II;
            alllocmax(kk+1:kk+nII) = data(i,II);
        end
        clear II nII
    end
%     check the correctness of the local peaks found    
%     figure; imagesc(data);  hold on; plot(nidx,midx,'*')
    
    nlocalmax = length(midx);
    
    idx = [];
    locmax = [];
    
    if nlocalmax > 0
        minlocmax = min(alllocmax);
        [locmax(1), kk] = max(alllocmax);
        idx(1,1) = midx(kk);
        idx(2,1) = nidx(kk);
        clear kk
        for i = 2:min(nmax, nlocalmax)
            JJ = find(abs(midx - idx(1,i-1)) < neighborsize(1) & abs(nidx - idx(2,i-1)) < neighborsize(2));
            alllocmax(JJ) = minlocmax - 1;
            if irowcirc == 1
                JJ = find(abs(abs(midx - idx(1,i-1)) - m) < neighborsize(1) & abs(nidx - idx(2,i-1)) < neighborsize(2));
                alllocmax(JJ) = minlocmax - 1;
            end
            if icolcirc == 1
                JJ = find(abs(midx - idx(1,i-1)) < neighborsize(1) & abs(abs(nidx - idx(2,i-1)) - n) < neighborsize(2));
                alllocmax(JJ) = minlocmax - 1;
            end
            if irowcirc == 1 && icolcirc == 1
                JJ = find(abs(abs(midx - idx(1,i-1)) - m) < neighborsize(1) & abs(abs(nidx - idx(2,i-1)) - n) < neighborsize(2));
                alllocmax(JJ) = minlocmax - 1;
            end
            clear JJ
            [maxval, kk] = max(alllocmax);
            
            if maxval < locmax(i-1) && maxval >= minlocmax
                locmax(i) = maxval;
                idx(1,i) = midx(kk);
                idx(2,i) = nidx(kk);
            else
                break
            end
            clear kk            
            
        end
        
    end
    
    count = length(locmax);
    
end