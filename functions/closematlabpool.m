function [] = closematlabpool
% as the title

    nlabs = matlabpool('size');
    
    if nlabs > 0
        matlabpool close;
    end