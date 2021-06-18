function [] = startmatlabpool( size )
%  start matlab pool
%   
    isstart = 0;
    nlabs = matlabpool('size');
    
    if nlabs == 0
        isstart = 1;
    end
    
    if isstart == 1
        if nargin == 0
            matlabpool('open','local');
        else
            try
                matlabpool('open','local',size);
            catch ce
                matlabpool('open','local');
                size = matlabpool('size');
                display(ce.message);
                display(strcat('input size wrong, use default pram , size =',num2str(size)));
            end
        end
    else
        display('matlabpool already lauched');
        if nlabs ~= size
            matlabpool close;
            startmatlabpool(size);
        end
    end
    


