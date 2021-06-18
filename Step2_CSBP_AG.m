clear all
close all

%% ----- Input Parameters ------ %%

%% get the parameters from the input parameter text
[dfile, dpath] = uigetfile({'*.txt'},'Open Parameter File:');

[ParaName,ParaValue,Comments]=textread([dpath dfile],'%s%f%q',-1);
NP=length(ParaName);

CurrentDataPath = cell2mat(ParaName(1));

copyfile('Step2_CSBP_AG.m', CurrentDataPath);
copyfile([dpath dfile], CurrentDataPath);


FigFolder0 = cd;
tracefile = cell2mat(ParaName(2));
FigFolder1 = cell2mat(ParaName(3));
FigFolder = ['.' cell2mat(ParaName(3))];
TraveTimeFile = cell2mat(ParaName(4));     % travel time table file
copyfile(TraveTimeFile, CurrentDataPath);

eval([ 'estdatamisftRatio=' cell2mat(ParaName(5)) ';']);


for np=6:NP
    eval([cell2mat(ParaName(np)) '=' num2str(ParaValue(np)) ';']);
end

% set the grids in the source region
xlocmin = -Ngrid_X * ddx2 * dxloc/2;
xlocmax = Ngrid_X * ddx2 * dxloc/2;
ylocmin = -Ngrid_Y * ddy2 * dyloc/2;
ylocmax = Ngrid_Y * ddy2 * dyloc/2;

twinhypo = [twinhypo_LB twinhypo_UB]; % start and end time for selecting the first hypocenter subevent 

freqband = [freqband_LB freqband_UB]; % frequency band (Hz) for band-pass filtering the traces
freqbandcvx = [freqbandcvx_LB freqbandcvx_UB]; % frequency band (Hz) for CVX and beamforming

azirange = [azirange_LB azirange_UB];
distrange = [distrange_LB distrange_UB];
lonrange = [lonrange_LB lonrange_UB];
latrange = [latrange_LB latrange_UB];	


faultxpt = faultxpt_LB:faultxpt_dif:faultxpt_UB; % x-coordinates of fault line
faultypt = faultxpt*tand(90-strike); % y-coordinates of fault line

%% change current path to the data file!
CS_computationFile=cd;

cd(CurrentDataPath); 
mkdir(FigFolder);

%% ----- End of Input Parameters ------ %%

% obtain lon and lat grids for 2-D backprojection
xrange = xlocmin:dxloc:xlocmax;
yrange = ylocmin:dyloc:ylocmax;
latloc = qlat + yrange/111.19;
lonloc = qlon + xrange/111.19/cosd(qlat);
nlonpt = length(lonloc);
nlatpt = length(latloc);

% find epicenter location index (iepilon, iepilat)
[A, iepilon] = min(abs(xrange));
[B, iepilat] = min(abs(yrange));

if xrange(iepilon) ~= 0 || yrange(iepilat) ~= 0
    errordlg('Epicenter not at (x=0, y=0)');
end


%% read original realigned trace data from realign_yao_efs codes
fs = fopen(tracefile, 'r', 'l'); % for genereated data
% read ns, npara, nt
header = fread(fs,1,'int');
ns = fread(fs, 1, 'int'); % num of traces
npara = fread(fs, 1, 'int'); % number of parameters for each trace
nt = fread(fs, 1, 'int'); % number of points for each trace
t1 = fread(fs, 1, 'float32') % starting time for stacks
t2 = fread(fs, 1, 'float32') % ending time for stacks
dt = fread(fs, 1, 'float32') % time interval for stacks
trailer = fread(fs,1, 'int');

slat0 = zeros(ns,1); % station latitude
slon0 = zeros(ns,1); % station longitude
del0 = zeros(ns,1);  % epicenter distance
azi0 = zeros(ns,1);  % src-station azimuth
tpred0 = zeros(ns,1); % predicted travel time from 1-D reference model
talign0 = zeros(ns,1); % travel time after alignment using cross-correlation
logscale0 = zeros(ns,1); % amplitude amplification ratio(log scale)for the trace
pol0 = zeros(ns,1); % polarity info (+1, 0, -1 for the trace)
rr0 = zeros(ns,1);  % cross-correlation coefficient between trace and stacked trace
data0 = zeros(ns, nt); % data matrix (ntrace x ntime_point)
ntr0 = ns;
dt = single(dt); % get rid of very small decimal numbers

for i = 1:ns
    % read statin parameters
    header = fread(fs,1,'int');
    slat0(i) = fread(fs, 1, 'float32');
    slon0(i) = fread(fs, 1, 'float32');
    del0(i) = fread(fs, 1, 'float32');
    azi0(i) = fread(fs, 1, 'float32');
    tpred0(i) = fread(fs, 1, 'float32');
    talign0(i) = fread(fs, 1, 'float32');
    logscale0(i) = fread(fs, 1, 'float32');
    pol0(i) = fread(fs, 1, 'float32');
    rr0(i) = fread(fs, 1, 'float32');
    trailer = fread(fs,1, 'int');
    
    % read data of corresponding station
    header = fread(fs,1,'int');
    data0(i, 1:nt) = fread(fs, nt, 'float32');
    trailer = fread(fs,1, 'int');
    
end

fclose(fs);



% II0=find(slon0<0);
% slon0(II0)=slon0(II0)+360;


% specific correction due to the missing the first event at hypocenter
if thopyshift ~= 0	
    display(['Starting and ending time with respect to hypocenter time are changed to:']);													  
    t1 = t1+thopyshift
    t2 = t2+thopyshift
end


%% sort trace according to azimuth
[temp,JJ] = sort(azi0);
slat0 = slat0(JJ); % station latitude
slon0 = slon0(JJ); % station longitude
del0 = del0(JJ);  % epicenter distance
azi0 = azi0(JJ);  % src-station azimuth
tpred0 = tpred0(JJ); % predicted travel time from 1-D reference model
talign0 = talign0(JJ); % travel time after alignment using cross-correlation
logscale0 = logscale0(JJ); % amplitude amplification ratio(log scale)for the trace
pol0 = pol0(JJ); % polarity info (+1, 0, -1 for the trace)
rr0 = rr0(JJ);  % cross-correlation coefficient between trace and stacked trace
data0 = data0(JJ, :); % data matrix (ntrace x ntime_point)


%% read travel time table
ftb = fopen(TraveTimeFile, 'r');
temp = fgetl(ftb); % read first line
nepidist = fscanf(ftb, '%d', 1); % number of epicenter distance in the table, second line
nsrcdep = fscanf(ftb, '%d', 1); % number of depth in the table, second line
temp = fgetl(ftb); 
srcdep = fscanf(ftb, '%f', nsrcdep); % read third line: source depth vector
tmptable = fscanf(ftb, '%f', [nsrcdep+1, nepidist]); % read travel time table
ttable = tmptable(2:(nsrcdep+1),:);
epidist = tmptable(1, :);
fclose(ftb);


%% calculate predicted travel time from epicenter to stations
IIeffsta = ones(1,ntr0);
tepi0 = zeros(ntr0,1);
ellip = almanac('earth','ellipsoid');
epistadist0 = km2deg(distance(qlat*ones(ntr0,1), qlon*ones(ntr0,1), slat0, slon0, ellip));
tepi0 = interp2(epidist, srcdep, ttable, epistadist0, qdep*ones(ntr0,1), 'linear');
tepi0 = tepi0*60;  % minute to second
II = isnan(tepi0);
IIeffsta(II) = 0; % for distance/depth outside the travel time table


%% calculate time shifts of the traces at each grid location with respect 
%% to the epicenter location

% predicted travel time (s) from grid points to stations, -1 value for
% distance/depth outside of travel time table
tlocgrids0 = zeros(nlonpt, nlatpt, ntr0); 
% predicted travel time shifts (in points) between the grid point and the
% epicenter location round((tlocgrids - tepi)/dt)
tpredshiftgrids0 = zeros(nlonpt, nlatpt, ntr0);
ntpredshiftgrids0 = zeros(nlonpt, nlatpt, ntr0); 
                                              

%% CAUTIONS: if the grid information in the input_file has been changed, travel_time.mat should be recalculated!!
ellip = almanac('earth','ellipsoid');
tic

display('calculating travel time shift from travel time table ......');
for ilon = 1:nlonpt
    for jlat = 1:nlatpt
        gridlon = lonloc(ilon);
        gridlat = latloc(jlat);
        locstadist = km2deg(distance(gridlat*ones(ntr0,1), gridlon*ones(ntr0,1), slat0, slon0, ellip));
        tlocgrids0(ilon,jlat,1:ntr0) = interp2(epidist, srcdep, ttable, locstadist, qdep*ones(ntr0,1), 'linear');
        tlocgrids0(ilon,jlat,1:ntr0) = tlocgrids0(ilon,jlat,1:ntr0)*60; % minute to seconds
        tpredshiftgrids0(ilon,jlat,1:ntr0) = tlocgrids0(ilon,jlat,1:ntr0) - reshape(tepi0,1,1,ntr0); % travel time different between grids and epicenter to stations
        ntpredshiftgrids0(ilon,jlat,1:ntr0) = round(tpredshiftgrids0(ilon,jlat,1:ntr0)/dt); % shift points for the grid location
        IIeffsta(isnan(tlocgrids0(ilon,jlat,1:ntr0))) = 0;
        tlocgrids0(ilon,jlat,isnan(tlocgrids0(ilon,jlat,1:ntr0))) = -1; % set travel time to be -1 for distance/depth outside the travel time table
        ntpredshiftgrids0(ilon,jlat,isnan(ntpredshiftgrids0(ilon,jlat,1:ntr0))) = 0; % set shift points be 0 for distance/depth outside the travel time table
    end
end

toc

save('travel_time.mat','gridlat','gridlon','locstadist','tlocgrids0','tpredshiftgrids0','ntpredshiftgrids0','IIeffsta');


%% calculate SNR of each trace and only keep SNR >= minSNR trace with
%% non-zero amplitude and non-zero polarity for next step backprojection
maxamp = max(abs(data0), [], 2);
noiseamp = zeros(1,ntr0);
sigamp = zeros(1,ntr0);
traceSNR = zeros(1,ntr0);
% iptmin = round((iptmin - t1)/dt) + 1;
iptmin = 1;
iptzero = round((0 - t1)/dt) + 1;
iptmax = round((ptmax - t1)/dt) + 1;
SampleF = double(1/dt);
[b,a] = butter(2, [freqband(1)/(SampleF/2) freqband(2)/(SampleF/2)]);
bpdata = zeros(size(data0,1),size(data0,2));
for i = 1:ntr0
    if maxamp(i) == 0 || pol0(i) == 0
        traceSNR(i) = 0;
    else
        bpdata(i,:) = filtfilt(b,a,data0(i,:));
        tracegroup = abs(hilbert(bpdata(i,:)));
        noiseamp(i) = mean(tracegroup(iptmin:iptzero));
        sigamp(i) = max(tracegroup(iptzero:iptmax));
        traceSNR(i) = sigamp(i)/noiseamp(i);
%         figure(1); hold off; plot(data0(i,iptmin:iptmax)); 
%         hold on; plot(bpdata(i,iptmin:iptmax), 'k');
%         hold on; plot(tracegroup(iptmin:iptmax), 'r');
%         title(['SNR = ', num2str(traceSNR(i))]); 
%         waitforbuttonpress
    end
    
end
    
IIgoodSNR = (traceSNR > minSNR); % find trace index with SNR > minSNR

%% find stations within given azimuthal ranges 
nrange = size(azirange,1);
IIazifit = zeros(nrange,ns);
IIazi = zeros(1,ns);
for i = 1:nrange
    IIazifit(i,:) = (azi0 >= azirange(i,1) & azi0 <= azirange(i,2));
    IIazi = IIazi | IIazifit(i,:);
end

%% find stations within given distance ranges
nrange = size(distrange,1);
IIdistfit = zeros(nrange,ns);
IIdist = zeros(1,ns);
for i = 1:nrange
    IIdistfit(i,:) = (del0 >= distrange(i,1) & del0 <= distrange(i,2));
    IIdist = IIdist | IIdistfit(i,:);
end

%% find stations within given longitude range
nrange = size(lonrange,1);
IIlonfit = zeros(nrange,ns);
IIlon = zeros(1,ns);
for i = 1:nrange
    IIlonfit(i,:) = (slon0 >= lonrange(i,1) & slon0 <= lonrange(i,2));
    IIlon = IIlon | IIlonfit(i,:);
end

%% find stations within given latitude range
nrange = size(latrange,1);
IIlatfit = zeros(nrange,ns);
IIlat = zeros(1,ns);
for i = 1:nrange
    IIlatfit(i,:) = (slat0 >= latrange(i,1) & slat0 <= latrange(i,2));
    IIlat = IIlat | IIlatfit(i,:);
end

%% Use the random selected stations if needed
if IndexRandStation == 1
    RandTraceIndexALL=randperm(ntr0);
    NRandTraceUSE=round(ntr0*randPercent);
    RandTraceIndexUSE=RandTraceIndexALL(1:NRandTraceUSE);
    IIrandsta=zeros(1,ntr0);
    IIrandsta(RandTraceIndexUSE)=1;
else
    IIrandsta = ones(1,ntr0);
end
%% obtain new data matrix by removing stations outside of travel time table
%% or station without good SNR data or or without given azimuthal ranges 
%% and distance ranges. 
II = find(IIeffsta == 1 & IIgoodSNR == 1 & IIazi & IIdist & IIlon & IIlat & IIrandsta); 
data = bpdata(II,:); % data selected with SNR > minSNR and bandpass filtered
slat = slat0(II);  % station latitude
slon = slon0(II); % station longitude
del = del0(II);  % epicenter distance
azi = azi0(II);  % src-station azimuth
tpred = tpred0(II); % predicted travel time from 1-D reference model
talign = talign0(II); % travel time after alignment using cross-correlation
logscale = logscale0(II); % amplitude amplification ratio(log scale)for the trace
pol = pol0(II); % polarity info (+1, 0, -1 for the trace)
rr = rr0(II);  % cross-correlation coefficient between trace and stacked trace
ntr = length(II);

data = 10*data/(max(max(abs(data))));


tepi = tepi0(II);
epistadist = epistadist0(II);
tlocgrids = tlocgrids0(:,:,II);
tpredshiftgrids = tpredshiftgrids0(:,:,II);
ntpredshiftgrids = ntpredshiftgrids0(:,:,II);

save([FigFolder 'DiffTravtbl.mat'], 'tpredshiftgrids', 'lonloc', 'latloc', 'qlon', 'qlat');

clear bpdata
clear tepi0 epistadist0 tlocgrids0 tpredshiftgrids0 ntpredshiftgrids0 

%% calculate new amplitude scaling of each trace within [ptmin ptmax]
%% window to check the correctness of logscale from realign_yao_efs

% obtain start and end points of the trace for processing
iptmin = round((ptmin - t1)/dt) + 1;
iptmax = round((ptmax - t1)/dt) + 1;
maxabsamptrace = max(abs(data(:,iptmin:iptmax)), [], 2);
logscale = log10(maxabsamptrace);
% normalize the trace amplitude using the max abs amplitude
for i = 1:ntr
    data(i,:) = data(i,:)/maxabsamptrace(i);
end
% data = 10*data; 
logscale = zeros(1,ntr);


%% calculate ampliude weighting of each station (with correlation coefficient 
%%  larger than minrr) from station density 
% the weighting is the inverse of the number of stations with respect to 
% center station within the statweightdist circle range.

if IndexpospolOnly == 1
    IIstackGood = find(rr > minrr & pol > 0);
else
    IIstackGood = find(rr > minrr);
end
ntrGood = length(IIstackGood);

staAmpWeight = zeros(ntr,1); % station amplitude weight vector
for kk = 1:ntrGood
    i = IIstackGood(kk);   
    dist = deg2km(distance(slat(i)*ones(ntrGood,1), slon(i)*ones(ntrGood,1), slat(IIstackGood), slon(IIstackGood)));
    staAmpWeight(i) = 1/length(find(dist < staweightdist));
end

%% obtain and plot aligned trace with respect to epicentral location and
%% the master (linear) stack trace for cross-correlation and realignment
iptmin = round((ptmin - t1)/dt) + 1;
iptmax = round((ptmax - t1)/dt) + 1;
nptnew = iptmax - iptmin + 1;
trange = twinhypo(1):dt:twinhypo(2);
itrange = round((trange - ptmin)/dt)+1;

newdata = zeros(ntr, nptnew); % shifted traces for epicenter location
epistackdata = zeros(1, nptnew); % stack of the shifted traces for epicenter location
for kk = 1:ntrGood
    i = IIstackGood(kk);
    newdata(i,:) = data(i,iptmin:iptmax)/10^logscale(i);
    epistackdata = epistackdata + staAmpWeight(i)*newdata(i,:)*pol(i);
end
 
stackweight = sum(staAmpWeight(IIstackGood));
epistackdata = epistackdata/stackweight;
epistackdata0 = epistackdata;

%% perform cross-correlation between the stacked trace within the given
%% time window (twinhypo) and given frequency range with all individual 
%% traces to obtain the new additional travel time shifts

tmaster = (twinhypo(1) + twinhypo(2))/2;

ipt1 = round((twinhypo(1) - ptmin)/dt) + 1; % start point (time) of twinhypo
ipt2 = round((twinhypo(2) - ptmin)/dt) + 1; % end point (time) of twinhypo
nptwinhypo = ipt2 - ipt1 + 1; % number of points in twinhypo
if mod(nptwinhypo,2) == 0  % to make sure numwinpt is an odd number
    nptwinhypo = nptwinhypo + 1;
    twinhypo(2) = twinhypo(2) + dt;
    ipt2 = ipt2 + 1;
end
npthalfwinhypo = round((nptwinhypo-1)/2); % number of points for half of the subevent signal window  

maxshiftpt = round(maxshifttime/dt); % max shift points for cross-correlation of subevent signals
totshiftpt = 2*maxshiftpt + 1; % total shiftpoints ( = size of (-maxshiftpt):maxshiftpt) )
epistackdata = epistackdata0;

maxiter = 3;
trshiftpt = zeros(ntr,maxiter); 
% figure;
% color ={'r','g','b','k','m'};
for iter = 1:maxiter 
    mastertrace = epistackdata(ipt1:ipt2);
%     hold on; plot(twinhypo(1):dt:twinhypo(2), mastertrace, color{iter});
    sqmastertr = sum(mastertrace.^2); % sum of amplitude squared of the mastertrace
    statrace = zeros(ntr, nptnew); % traces after alignment with respect to sub-event location
    xcorstack = zeros(ntr, totshiftpt); % cross-correlation function between individual trace 
    maxxcorcoef = zeros(ntr,1); % max abs value of the cross-correlation coef. 
    polxcorcoef = zeros(ntr,1);

    wintrace = zeros(ntr, nptwinhypo);
    iptwinsta = iptmin + (ipt1:ipt2) - 1;
    for k = 1:ntr
        
        wintrace(k,:) = data(k, iptwinsta);
        % cross-correlation
        for ishift = (-maxshiftpt):maxshiftpt
            ic = ishift + maxshiftpt + 1;
            iptxcor = iptwinsta + ishift;
            xcorstack(k,ic) = dot(mastertrace, data(k, iptxcor));
            xcorstack(k,ic) = xcorstack(k,ic)/sqrt(sqmastertr*sum(data(k, iptxcor).^2));
        end
        [maxxcorcoef(k), II] = max(abs(xcorstack(k,:)));
        polxcorcoef(k) = sign(xcorstack(k,II));
        trshiftpt(k,iter) = II - maxshiftpt - 1;
    end
    
    % update stacked trace
    epistackdata = zeros(1,nptnew);
    for i = 1:ntr
        if maxxcorcoef(i) > minrr
            epistackdata = epistackdata + data(i,(iptmin:iptmax)+trshiftpt(i,iter))*polxcorcoef(i);
        end
    end
    stackweight = sum(maxxcorcoef > minrr);
    epistackdata = epistackdata/stackweight;       
    
end

newtrshiftpt = trshiftpt(:,maxiter);

% recalculate the index of good traces for stacking
if IndexpospolOnly == 1
    IIstackGood = find(maxxcorcoef > minrr & polxcorcoef > 0);
    IIGoodIndex = (maxxcorcoef > minrr & polxcorcoef > 0);
else
    IIstackGood = find(maxxcorcoef > minrr);
    IIGoodIndex = (maxxcorcoef > minrr);
end
ntrGood = length(IIstackGood);

% recalculate the station weight for stacking
staAmpWeight = zeros(ntr,1); % station amplitude weight vector
for kk = 1:ntrGood
    i = IIstackGood(kk);   
    dist = deg2km(distance(slat(i)*ones(ntrGood,1), slon(i)*ones(ntrGood,1), slat(IIstackGood), slon(IIstackGood)));
    staAmpWeight(i) = 1/length(find(dist < staweightdist));
end

newdata = zeros(ntr, nptnew);
epistackdata = zeros(1, nptnew); % stack of the shifted traces for epicenter location
for i = 1:ntr 
    newdata(i,:) = data(i,(iptmin:iptmax)+newtrshiftpt(i))*polxcorcoef(i)/(10^logscale(i));
    if IIGoodIndex(i);
        epistackdata = epistackdata + newdata(i,:)*staAmpWeight(i);
    end
end
 
stackweight = sum(staAmpWeight(IIstackGood));
epistackdata = epistackdata/stackweight;
[maxabsamp, id] = max(abs(epistackdata(itrange)));
tmaxampepi = trange(id);  % time corresponding to the abs. max. of the nth root stack

% re-assign c.c. coef and polarity value
rr = maxxcorcoef;
pol = polxcorcoef;

%%
figure(1);

tracemaxamp = max(max(newdata(IIstackGood,:)));
climnewdata = tracemaxamp*[-1 1];
% plot the aligned waveform with Good correlation
subplot(3,6,1:6);
hold on; plot(ptmin:dt:ptmax,ntr0/2*epistackdata/(1.2*max(abs(epistackdata))),'-k','LineWidth',0.8);
%hold on; plot(ptmin:dt:ptmax,ntrGood/2*epistackdata/(1.2*max(abs(epistackdata))),'-b','LineWidth',0.8);
xlim([ptmin,ptmax])
ylim([-1/2*ntr0,0.5*ntr0])
set(gca,'ytick',[]);
%legend('before corr','good corr (>0.6)','Location','SouthEast');



subplot(3,6,7:12);
imagesc(ptmin:dt:ptmax, 1:ntrGood, newdata(IIstackGood,:), climnewdata); colorbar;
hold on; plot(tmaxampepi*[1 1], [1 ntrGood], 'k--', 'LineWidth', 1);
hold on; plot(twinhypo(1)*[1 1], [1 ntrGood], 'k', 'LineWidth', 1);
hold on; plot(twinhypo(2)*[1 1], [1 ntrGood], 'k', 'LineWidth', 1);
newdata2=newdata(IIstackGood,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only Used to Make Movie
newdata2=newdata(IIstackGood,:);
save('waveform.mat','newdata2','ntrGood','climnewdata','dt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Only Used to Make Movie


traceplotratio = 1/maxabsamp;

xlabel('Time (s)', 'FontSize', 14);
ylabel('Trace Number', 'FontSize', 14);
set(gca,'FontSize', 14, 'YDir', 'Normal');

subplot(3,6,13:14);
tempts = -5;
tempte = 15;
ipttempts = round((tempts - ptmin)/dt) + 1;
ipttempte = round((tempte - ptmin)/dt) + 1;
tempmedian = median(max(newdata(IIstackGood,ipttempts:ipttempte),[],2));
imagesc(tempts:dt:tempte, 1:ntrGood, newdata(IIstackGood,ipttempts:ipttempte), tempmedian*[-1 1]); colorbar;
hold on; plot(tmaxampepi*[1 1], [1 ntrGood], 'k--', 'LineWidth', 1);
hold on; plot(twinhypo(1)*[1 1], [1 ntrGood], 'k', 'LineWidth', 1);
hold on; plot(twinhypo(2)*[1 1], [1 ntrGood], 'k', 'LineWidth', 1);

ylabel('Trace Number', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14);
set(gca, 'FontSize', 14, 'YDir', 'Normal');

subplot(3,6,16:17);
plot(maxxcorcoef(IIstackGood),1:ntrGood,'r.');
hold on; plot([mean(maxxcorcoef(IIstackGood)) mean(maxxcorcoef(IIstackGood))],[1 ntrGood],'k--','LineWidth',1); 
set(gca,'FontSize', 14,'Ydir','Normal');
xlim([minrr-0.1 1]);
ylim([1 ntrGood]);
ylabel('Trace Number', 'FontSize', 14);
xlabel('cross-correlation coefficient', 'FontSize', 8);

print('-dpdf', [FigFolder 'hypocenter_align_Goodtrace.pdf']);

clf

%% plot locations of only selected traces for final stack
figure(1); hold off
coast = load('coastlines');
axesm('eqdazim','origin',[qlat qlon])
axis off; framem on;
setm(gca,'MLabelParallel',0);
plotm(coast.coastlat, coast.coastlon, 'Color', [0.6 0.6 0.6]);
nIIGoodPos = 0;
nIIGoodNeg = 0;
for k = 1:ntr
    if (rr(k) > minrr) && (pol(k) > 0)  %positive polarity
        plotm(slat(k),slon(k),'^','MarkerSize', 2);
        hold on
        nIIGoodPos = nIIGoodPos + 1;
    end
end
if IndexpospolOnly ~= 1
    for k = 1:ntr
        if (rr(k) > minrr) && (pol(k) < 0) %negative polarity
            plotm(slat(k),slon(k),'m^','MarkerSize', 2);
            hold on
            nIIGoodNeg = nIIGoodNeg + 1;
        end
    end
end

hold on; plotm(qlat, qlon, 'gp', 'MarkerSize',12, 'MarkerFaceColor','g');
title(['P=', num2str(nIIGoodPos),', N=', num2str(nIIGoodNeg)],'FontSize', 10);

print('-dpdf', [FigFolder 'GoodTraceDistribution.pdf']);

clf

figure(1);
hold off
latlim=[min(slat)-2,max(slat)+2];
lonlim=[min(slon)-2,max(slon)+2];
worldmap(latlim,lonlim);
hold on
plotm(coast.coastlat, coast.coastlon, 'Color', [0 0 0]);

for k = 1:ntr
    if (rr(k) > minrr) && (pol(k) > 0)  %positive polarity
        plotm(slat(k),slon(k),'^','MarkerSize', 3);
        hold on
    end
end

if IndexpospolOnly ~= 1
    for k = 1:ntr
        if (rr(k) > minrr) && (pol(k) < 0) %negative polarity
            plotm(slat(k),slon(k),'m^','MarkerSize', 3);
            hold on
        end
    end
end

title(['P=', num2str(nIIGoodPos),', N=', num2str(nIIGoodNeg)],'FontSize', 10);

print('-dpdf', [FigFolder 'GoodTraceDistribution2.pdf']);
hold off

%% save selected stations 
SLON=slon(IIstackGood);
SLAT=slat(IIstackGood);
STATIONS=[SLON,SLAT];
save station_distribution.dat STATIONS -ascii 

%% source imaging using a convex optimization method (CVX package)

iptminAll = round((ptmin - t1)/dt) + 1; % start point (time) for imaging
iptmaxAll = round((ptmax - t1)/dt) + 1; % end point (time) for imaging
nptnewAll = iptmaxAll - iptminAll; % number of points in the waveform for bpj

iptmin = round((ptmin - t1)/dt) + 1;  % initiation of the first time window start point
iptmax = round((ptmin + segwint - t1)/dt) + 1; % initiation of the first time window end point
nptnew = iptmax - iptmin + 1; % number of points segment window for beamforming

nsegwin = floor((ptmax - ptmin - segwint)/dtsegwin) + 1; % number of window segment for imaging
ntseg = round(segwint/dt) + 1; % number of data points in each time window segment

nfft = 2^(nextpow2(ntseg)); % number of points for trace data for imaging
tukeytaper = tukeywin(ntseg, 0.2)'; % time domain taper window for the data
fs = 1/dt; % data sampling frequency
df = fs/nfft; % frequency interval of the data spectrum
ff = (0:(nfft/2))*df; % frequency vector of the data spectrum
nff = length(ff); % number of f in the spectrum

IIeff = find(ff > freqbandcvx(1) & ff < freqbandcvx(2)); % effective freq index within freqband
neff = length(IIeff); % number of effective freq.
ffeff = ff(IIeff); % freq for inverse imaging



%% Dimensionality reduction matrix (to accelerate the computation)


if IndexDimenRedu ==1    % multiplying a random Dimension-reduction matrix
    ntrDR=floor(ntrGood*DRpercentage);
    DRmatrix=randn(ntrDR,ntrGood)/sqrt(ntrGood);
else
    
    DRmatrix=eye(ntrGood);
    ntrDR=ntrGood;
    
end




Nsrc = nlonpt*nlatpt; % number of grids (src region)
Atime = reshape(tpredshiftgrids(:,:,IIstackGood), Nsrc, ntrGood); % with size Nsrc x ntrGood  
srcGridSpec = zeros(nlonpt, nlatpt, neff, nsegwin); % time-dependent src grids spectrum data (to be inverted for)
srcGridBeam = zeros(nlonpt, nlatpt, neff, nsegwin); % time-dependent beamforming results (forward calculation)
coast = load('coastlines');
cmap = colormap(pink);
cmap = cmap(end:-1:1,:);

mkdir('./align_waveform');
waveform_align_file='./align_waveform';

for iw = 1:nsegwin
    wint1 = ptmin + (iw - 1)*dtsegwin;  % window starting time
    wint2 = ptmin + (iw - 1)*dtsegwin + segwint; % window ending time
    display([num2str(wint1) '--' num2str(wint2) ' s']);
    
    iptwint1 = round((wint1 - t1)/dt) + 1;
    iptwint2 = round((wint2 - t1)/dt) + 1;
    windata = zeros(ntrGood, nfft); % time-domain data
    fftdata = zeros(ntrGood, nfft); % spectrum data
    specdata = zeros(ntrGood, neff); % spectrum data
    
    % window and taper data traces
    for i = 1:ntrGood
        j = IIstackGood(i); % obtain good trace index
        iptwin = (iptwint1:iptwint2) + newtrshiftpt(j); % applying addtional time shifts
        windata(i,1:ntseg) = tukeytaper.*data(j,iptwin)*polxcorcoef(j); % taper the time domain data
    end
    maxdataAmp = max(abs(windata),[],2); % max amplitude of each trace
    maxAmpAll = max(maxdataAmp); % max amplitude of all traces
    
    figure(777)
    clf
    subplot(1,3,1)
    imagesc(wint1:dtsegwin:wint2,1:ntrGood,windata);
    caxis([-1 1])
    colorbar;
    align_file_name=[waveform_align_file '/' num2str(wint1) 's.pdf'];
    print('-dpdf',align_file_name);    
    
    % normalize the trace amplitude by the max amplitude of all traces
    for i = 1:ntrGood
        windata(i,:) = windata(i,:)*maxAmpAll/abs(maxdataAmp(i));
    end
    
    windatanew = windata;
    

    % do FFT for all the trace windows 
    for i = 1:ntrGood
        fftdata(i,1:nfft) = fft(windatanew(i,1:nfft), nfft); % fft --> spectrum data
    end
    specdata = fftdata(:,IIeff); % spectrum data for inverse imaging using CVX   

    dataweight = double(diag(staAmpWeight(IIstackGood))); % station weight matrix
    
    sumweightL1 = sum(staAmpWeight(IIstackGood));
    sumweightL2 = sqrt(sum(staAmpWeight(IIstackGood).^2));
    
    I_row=zeros(1,neff);
    I_column=zeros(1,neff);
    Peak_Amp=zeros(1,neff);
        
    for ifreq = 1:neff % loop over freq
        
        %% sparse source localization inversion using CVX 
        display([' ...... Frequency = ' num2str(ffeff(ifreq)) ' Hz ......']);
        freq = ffeff(ifreq); 
        omega = freq*2*pi;
        Arep = double(exp(-1j*omega*Atime.')); % note: Atime.' means array transpose, but not involve conjugation!
        yobs = double(specdata(:,ifreq)); % data spectrum at certain freq.
        yobsAmp = abs(yobs);
        II = find(yobsAmp > 0);
        yobs(II) = yobs(II)*median(yobsAmp)./yobsAmp(II); % normalize the spectrum applitude to median amplitude

        
        estdatamisftRatio = sort(estdatamisftRatio);        
        estdatamisftRatio = double(estdatamisftRatio);
        
        
        if IndexWavelet == 1  % using wavelet transform to do the CS inversion
            Arep2=Arep;
            for i=1:ntrGood
                Arep_temp=reshape(Arep(i,:),nlonpt,nlatpt);
                Arep_temp_wlt=lwt2(Arep_temp,'db4');
                Arep(i,:)=reshape(Arep_temp_wlt,1,Nsrc);
            end
        end
            
        
        
        %% CVX solution to source grids spectrum
        l1norm=zeros(size(estdatamisftRatio)); % L1norm(x)
        lnorm_dif=zeros(size(estdatamisftRatio)); % L?norm(Ax-b) --difference 
        linearFittingParameter=zeros(size(estdatamisftRatio));%linerfittingParameter=l?norm(Ax-b)/l?norm(b),which depends on the level of random noise. 
        wMatrix=zeros(Nsrc,length(estdatamisftRatio));
        %% $$ to calculate the l1norm and lnorm_dif
        for Id=1:length(estdatamisftRatio)
            
            Dx=ddx2; 
            Dy=ddy2;
            X=struct('X0',[],'index',[],'level',[]); % structure keeping the information of grid: X0--results; index: index(:,2) row and column index
            % in A matrix; level: refining
            % level
            ixy = 1;
            for ix = 1:Ngrid_X
                for iy = 1:Ngrid_Y
                    X.index(ixy,1)=1+(ix-1)*Dx;  % x index
                    X.index(ixy,2)=1+(iy-1)*Dy;  % y index
                    ixy=ixy+1;
                end
            end
            
            X.level=ones(Ngrid_X*Ngrid_Y,1); % refining level (1-no refinement, 2-refined once, 3-refined twice,...)
            X.X0=X.level;
            
            % to get the submatrix of A needed in the calculation
            A_column = (X.index(:,2)-1)*nlonpt+X.index(:,1);  % transform index in the matrix to the column index in A matrix
            A_calculate=Arep(:,A_column); % submatrix of A, used in the computation of
            
            N_cal=length(X.level);  % initial coarse grid number
            
            
            % calculate results using CVX package
            cvx_clear
            cvx_begin
            variable w(N_cal) complex
            minimize(norm(DRmatrix*A_calculate*w-DRmatrix*yobs,Lnorm)+estdatamisftRatio(Id)*ntrDR*norm(w,1));
            cvx_end
            
            X.X0=abs(w); % record the results on the coarse grids
            
         
            % begin to refine grid
            for n_denser=1:log2(ddx2)  %  densify the grid for KK times: Dx=2^(KK)
                disp(['Refining the grid ' num2str(n_denser) ' times']);
                
                % grid densify begin
                [I,J]=find(X.X0>10^(-4)*max(X.X0) & X.X0>10^(-7));
                
                NK=length(I);
                for n_k=1:NK
                    index_temp=X.index(I(n_k),:);
                    
                    D_grid_X=Dx/2^(X.level(I(n_k))); % refined grid (left and right)
                    D_grid_Y=Dy/2^(X.level(I(n_k))); % refined grid (up and down)
                    X.level(I(n_k))=X.level(I(n_k))+1; % add the level of grid points that have been refined
                    
                    index_new = zeros(8,2);
                    level_new = zeros(8,1);
                    
                    ii=1;
                    
                    if index_temp(1)-D_grid_X > 0
                        index_new(ii,:) = index_temp+[-D_grid_X,0]; %left
                        ii=ii+1;
                    end
                    if index_temp(2)-D_grid_Y > 0
                        index_new(ii,:) = index_temp+[0,-D_grid_Y]; %up
                        ii=ii+1;
                    end
                    
                    if index_temp(2)-D_grid_Y > 0 && index_temp(1)-D_grid_X > 0
                        index_new(ii,:) = index_temp+[-D_grid_X,-D_grid_Y]; %up-left
                        ii=ii+1;
                    end
                    
                    if index_temp(2)+D_grid_Y <= nlatpt
                        index_new(ii,:) = index_temp+[0,D_grid_Y]; %down
                        ii=ii+1;
                    end
                    
                    if index_temp(2)+D_grid_Y <= nlatpt && index_temp(1)-D_grid_X > 0
                        index_new(ii,:) = index_temp+[-D_grid_X,D_grid_Y]; %down-left
                        ii=ii+1;
                    end
                    
                    if index_temp(1)+D_grid_X <=nlonpt
                        index_new(ii,:) = index_temp+[D_grid_X,0]; %right
                        ii=ii+1;
                    end
                    
                    if index_temp(2)-D_grid_Y >0 && index_temp(1)+D_grid_X <=nlonpt
                        index_new(ii,:) = index_temp+[D_grid_X,-D_grid_Y]; %up-right
                        ii=ii+1;
                    end
                    if index_temp(2)+D_grid_Y < nlatpt && index_temp(1)+D_grid_X <=nlonpt
                        index_new(ii,:) = index_temp+[D_grid_X,D_grid_Y]; %down-right
                        ii=ii+1;
                    end
                    
                    if ii <= 8
                        index_new((ii):8,:)=[];
                        level_new((ii):8,:)=[];
                    end
                    
                    if ii ~= 1
                        ii=ii-1;
                    end
                    
                    
                    level_new(:) = X.level(I(n_k));
                    
                    
                    X.index=[X.index; index_new];
                    X.X0=[X.X0; zeros(ii,1)];
                    X.level=[X.level; level_new];
                    
                    
                end
                [X.index,IC,IA]=unique(X.index,'rows');
                X.X0=X.X0(IC);
                X.level=X.level(IC);
                
                
                A_column = (X.index(:,2)-1)*nlonpt+X.index(:,1);  % new A matrix corresponding to the grids
                A_calculate = Arep(:,A_column);
                N_cal=length(X.level);
                
                % calculate results using CVX package
                cvx_clear
                cvx_begin
                variable w(N_cal) complex
                minimize(norm(DRmatrix*A_calculate*w-DRmatrix*yobs,Lnorm)+estdatamisftRatio(Id)*ntrDR*norm(w,1));
                cvx_end
                
                X.X0=abs(w);
                
            end
            %
            XXX=zeros(nlatpt,nlonpt);
            for i=1:length(X.X0)
                XXX(X.index(i,2),X.index(i,1))=w(i);
            end
            w=reshape(XXX',Nsrc,1);

%%            
            wMatrix(:,Id)=w;        
            l1norm(Id)=norm(w, 1)*ntrDR;        
            lnorm_dif(Id)=norm(DRmatrix*Arep * w - DRmatrix*yobs, 2);
            linearFittingParameter(Id)=lnorm_dif(Id)/norm(yobs,2);
        end 
        
        %% $$  to remove the unsteady solution points and record the steady solution points
        jj=1;
        kk=length(estdatamisftRatio);
        for ii=1:length(estdatamisftRatio)
            if abs(lnorm_dif(ii)) < 10^(-7)
                jj=jj+1;
            else
                break
            end               
        end
        for ii=length(estdatamisftRatio):(-1):1
            if abs(l1norm(ii)) < 10^(-10)
                kk=kk-1;
            else
                break
            end
        end
        SteadyNum=kk-jj+1;
        l1normSteady=l1norm(jj:kk);
        lnormSteady=lnorm_dif(jj:kk);
        SteadymisftRatio=estdatamisftRatio(jj:kk);
        
        %% scale the steady solution to the (maximum-minimum)
        l1normScaled=zeros(1,SteadyNum);
        lnormScaled=zeros(1,SteadyNum);
        ScaledSum=zeros(1,SteadyNum);
        l1norMax=l1normSteady(1);
        l1norMin=l1normSteady(SteadyNum);
        lnormMax=lnormSteady(SteadyNum);
        lnormMin=lnormSteady(1);
        for ii=1:SteadyNum
            l1normScaled(ii)=(l1normSteady(ii)-l1norMin)/(l1norMax-l1norMin);
            lnormScaled(ii)=(lnormSteady(ii)-lnormMin)/(lnormMax-lnormMin);
            ScaledSum(ii)=sqrt(l1normScaled(ii)^2+lnormScaled(ii)^2);
        end
        
        SteadyData=zeros(6,SteadyNum);
        SteadyData(1,:)=SteadymisftRatio;
        SteadyData(2,:)=l1normSteady;
        SteadyData(3,:)=l1normScaled;
        SteadyData(4,:)=lnormSteady;
        SteadyData(5,:)=lnormScaled;
        SteadyData(6,:)=ScaledSum;
        
        [SUMminimum,II]=min(ScaledSum);
        I=II+jj-1; %to get the real index in the unscaled array               
        w=wMatrix(:,I);
        bestRatio=estdatamisftRatio(I);
        
        
        if IndexWavelet == 1  % do the inverse wavelet for the results 
            w_in=reshape(w,nlonpt,nlatpt);
            w_iwlt=ilwt2(w_in,'db4');
            srcAmpGrids=abs(w_iwlt);
        else
            srcAmpGrids = reshape(abs(w), nlonpt, nlatpt);
        end
        
        [tempMax,I_temp]=max(srcAmpGrids,[],1);
        [Peak_Amp(ifreq),I_column(ifreq)]=max(tempMax);
        I_row(ifreq)=I_temp(I_column(ifreq));
        
        srcGridSpec(:,:,ifreq,iw) = reshape(srcAmpGrids, nlonpt, nlatpt, 1, 1);
        
        
        %% beamforming of the data
        
        if IndexWavelet ==1
            beamfreq=sqrt(abs(Arep2'*(yobs.*staAmpWeight(IIstackGood))).^2)/sum(staAmpWeight(IIstackGood));
        else
            
            beamfreq = sqrt(abs(Arep'*(yobs.*staAmpWeight(IIstackGood))).^2)/sum(staAmpWeight(IIstackGood));
        end
        beamfreq = reshape(beamfreq, nlonpt, nlatpt);
        srcGridBeam(:,:,ifreq,iw) = reshape(beamfreq, nlonpt, nlatpt, 1, 1);
        
               
    end
    
  
    
    matfolder = [FigFolder 'mat/']; 
    mkdir(matfolder);
    matfile = [matfolder 'srcCVXresult_t' num2str((wint1+wint2)/2) 's.mat'];
    srcGridSpecWin = squeeze(srcGridSpec(:,:,:,iw));
    srcGridBeamWin = squeeze(srcGridBeam(:,:,:,iw));
    save(matfile, 'srcGridSpecWin', 'srcGridBeamWin', 'lonloc', 'latloc', 'qlon', 'qlat', 'xrange', 'yrange', 'ffeff', ...
            'ptmin', 'ptmax', 'dtsegwin', 'segwint', 'nsegwin', 'iw');

   
end


%% save output data to mat file

for iw = 1:nsegwin
    wint1 = ptmin + (iw - 1)*dtsegwin;  % window starting time 
    wint2 = ptmin + (iw - 1)*dtsegwin + segwint; % window ending time    
    display([num2str(wint1) '--' num2str(wint2) ' s']);
    matfilewin = [matfolder 'srcCVXresult_t' num2str((wint1+wint2)/2) 's.mat'];
    load(matfilewin);
    srcGridSpec(:,:,:,iw) = reshape(srcGridSpecWin, nlonpt, nlatpt, neff, 1);
    srcGridBeam(:,:,:,iw) = reshape(srcGridBeamWin, nlonpt, nlatpt, neff, 1);
end
    
matfile = [FigFolder 'srcCVXresult.mat'];
save(matfile, 'srcGridSpec', 'srcGridBeam', 'lonloc', 'latloc', 'qlon', 'qlat', 'xrange', 'yrange', 'ffeff', ...
    'ptmin', 'ptmax', 'dtsegwin', 'segwint','newdata2','ntrGood','climnewdata','dt');


%% make a plot input file

fid=fopen('Plotfile_in','w');

fprintf(fid,'%s\n',FigFolder);
fprintf(fid,'%s\n',FigFolder0);
fprintf(fid,'%s\n',[FigFolder0 FigFolder1 'DiffTravtbl.mat']);

fclose(fid);


%% back to the CS computation file
cd(CS_computationFile);
toc
