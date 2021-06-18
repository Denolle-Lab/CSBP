% remove instrument response (phase shift) for SAC data using
clear all; close all;


sacdatadir = [pwd, '/'];  % SAC data folder
respdatadir = ['./IRISDMC/'];   % RESP or PoleZero response file folder
resptype = 2;  % 1 for RESP response file; 2 for PoleZeros response file (but still convert to velocity)
outsacdir = [pwd , 'RmResp/']; % output sac data folder after response removal
freqrange = [0.006 0.007 19.9 20]; % [f_low_cut f_low_pass f_high_pass f_high_cut]
mintimelength = 300; % minimum length of time duration (s) for each trace kept for further analysis
mkdir(outsacdir);

sacfilelist = dir([sacdatadir '*.SAC']); % obtain all the sac data file names
allsacdata = readsac('*.SAC');


for i = 1:length(sacfilelist)
    
    s = readsac([sacdatadir sacfilelist(i).name]);
    display(['the ' num2str(i)  ' time'])
    indexGood = 1;
    for k = 1:i-1
        if strcmp(s.KSTNM, allsacdata(k).KSTNM)
            display(['repeating station: ' s.KSTNM]);
            indexGood = 0;
        end
    end
    
    if indexGood == 1
        
        if s.DELTA*s.NPTS > mintimelength
            
            s.DATA1 = s.DATA1 - mean(s.DATA1);
            s.DATA1 = detrend(s.DATA1);
            % For different kind of response files
            if resptype == 1
                respfilestr = ['RESP.' s.KNETWK '.' s.KSTNM '.' s.KHOLE '.' s.KCMPNM];
            elseif resptype == 2
                respfilestr = ['*SACPZ*' s.KNETWK '*' s.KSTNM '*' s.KCMPNM];
            end
            respfile = dir([respdatadir respfilestr]);
            respfilenew = [respdatadir respfile.name];
            
            if strcmpi(s.KSTNM, 'TARA')
                a = 0;
            end
            
            
            if length(respfile) == 1
                
                bfwave = rmResp_bpfilter(s.DATA1, 1/s.DELTA, freqrange, respfilenew, resptype);
                
                if isempty(bfwave) % something wrong with the instrumental response file
                    delete(s.FILENAME);
                else
                    display(s.FILENAME);
                    %         figure(11); hold off; plot((1:length(s.DATA1))*s.DELTA, s.DATA1, 'b');
                    %         hold on; plot((1:length(s.DATA1))*s.DELTA, bfwave, 'r');
                    % waitforbuttonpress
                    s.DATA1 = bfwave;
                    s.FILENAME = [outsacdir s.FILENAME];
                    
                    writesac(s);
                    
                    clear s;
                end
            else
                display(['***** No response file for station   ' s.KSTNM]);
            end
        else
            display(['!!!!! waveform duration too short for station   ' s.KSTNM]);
        end
    end
    
end