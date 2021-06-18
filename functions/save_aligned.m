function flag = save_aligned(event, t1, t2, outdir , outname)
%%@function save the aligned data for next program
% the module is designed for different usage in principle.
% this function is designed base on Yao's iterative back projection program
% the using form is 
%
%    flag = save_aligned(event, t1, t2, outdir , outname)
%
% @param @flag. a flage denotes whether the function is run successfully
%               1 for success while 0 for failure.
% @param @event. the event aligned.
% @param @t1. the starting time of the data, the time means the time on
%             aligned time, i.e. time after shifting.
% @param @t2. the end time of the data, the time means the time on
%             aligned time, i.e. time after shifting.
% @param @outdir. the directory that the output file is saved. 
% @param @outname. the name of the output file.

    flag = 0;
    
    
%% produce the data to be saved.
    nos = event.number_of_sac;
    data_matrix = matrix_build(event, 'time', [t1 t2]);
    npara = 14; % I don't know whether it's right or not, but it works well now;
    npts = size(data_matrix,2);
    fb = first_break(event);
    aver = mean(fb);
    e.p_time(1) = t1 - aver;
    e.p_time(2) = t2 - aver;
    r_s = corr_coef(event);    
    
    for i = 1:nos
        event.sac(i).logscale = log(max(abs(data_matrix(i, :))));
        data_matrix(i, :) = normal_lize(data_matrix(i, :));
    end
    
    
%% save the data;
    mkdir(outdir);
    fw = fopen(strcat(outdir, '/', outname), 'w', 'l');
    
    fwrite(fw, 1, 'int'); % header
    fwrite(fw, nos, 'int'); % number of traces
    fwrite(fw, npara, 'int'); % number of parameters for each trace
    fwrite(fw, npts, 'int'); % number of points for each trace
    fwrite(fw, e.p_time(1), 'float32'); % starting time for stacks
    fwrite(fw, e.p_time(2), 'float32'); % ending time for stacks
    fwrite(fw, event.sac(1).dt, 'float32'); % time interval for stacks
    fwrite(fw, 999, 'int'); %trailer
    
    for i = 1:nos
        %% write station parameters
        fwrite(fw, 1, 'int'); %header
        fwrite(fw, event.sac(i).STLA, 'float32'); % station latitude
        fwrite(fw, event.sac(i).STLO, 'float32'); % station longtitude
        fwrite(fw, event.sac(i).GCARC, 'float32'); % Station-to-event great-circle arc length (degree)
        fwrite(fw, event.sac(i).BAZ, 'float32'); % src-station azimuth
        fwrite(fw, event.sac(i).tt_time, 'float32'); % predicted travel time from 1-D reference model
        fwrite(fw, event.sac(i).first_break, 'float32'); % travel time after alignment using cross-correlation
        fwrite(fw, event.sac(i).logscale, 'float32'); % amplitude amplification ratio(log scale) for the trace
        fwrite(fw, event.sac(i).pol, 'float32'); % polarity info(1, 0, -1);
        fwrite(fw, r_s(i), 'float32'); %cross-correlation coefficient between trace and the stacked trace
        fwrite(fw, 999, 'int'); %trailer
        
        %% write data of the corresponding station
        fwrite(fw, 1, 'int'); %header
        fwrite(fw, data_matrix(i, :),'float');
        fwrite(fw, 999, 'int');
    end
    
    fclose(fw);
    
    flag = 1;
        
 end

