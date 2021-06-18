function M = matrix_build(event, start, ending, appear)
%pick the data in all sac to form a matrix
%@M the matrix built
%@param @event. the structure need to be pick-up
%@param @start. the lower limit of the window in seconds
%@param @ending. the upper limit of the window in seconds
%@param @appear. 'normal' or 'aligned' a parameter decide the data was in time aligned or first_break aligned. 
%                



nos = event.number_of_sac;
dt = event.sac(1).dt;
maxlength = max(field_extract(event,'length'));

if nargin == 1
    start = 1;
    ending = maxlength * event.sac(1).dt; 
    appear = 'normal';
end;

if nargin == 3
    appear = 'normal';
end



if strcmp(appear,'aligned')
    
    start = round(start / dt);
    ending = round(ending / dt);
    
    
    fb = first_break(event);
    fb = fb - mean(fb);
    fb = round(fb / dt);
    M = zeros(nos, (ending - start + 1));
    
    for i = 1:nos
        if ((start + fb(i)) < 0) 
             M(i, (- start - fb(i) + 2 ):(ending - start + 1)) = event.sac(i).data(1:(ending + fb(i)));
            
        else
            M(i, :) = (event.sac(i).data((start:ending) + fb(i)))'; % still need proved , not strong enough;
    
        end
    end
else
    
    M = zeros(nos,maxlength);
    
    for i = 1 : nos
        M(i, 1:event.sac(i).length) = event.sac(i).data;
        if event.sac(i).length < maxlength
                 M(i , (event.sac(i).length + 1) : maxlength) = NaN;
        end;
    end;

    start = round(start / event.sac(1).dt);
    ending = round(ending / event.sac(1).dt);

    M = M(: , start : ending);
end

