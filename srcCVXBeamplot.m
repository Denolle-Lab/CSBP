% plot srcimagingCVXBeam.m output results

clear all;
close all;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUST FOR CONVENIENCE!!!
% ftest=0.05;
% while (ftest<=1)
%     close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUST FOR CONVENIENCE!!!

%% some parameters for plotting
DataPath='./Tohoku_HF_AG';

PlotTxt='Plotfile_in';

IndexTopoBackground =0; % =0 not using Topography background; =1 use topography background
% be careful, a topography data file must exist in
% the data result file

IndexSettingEnd =1 ; % =1 set ending time ; =0 not
beginTime = 0;
endingTime = 180;


PlotFile=cd;

IndexSettingRange =1; % =1 set ploting range; =0 not
Xrange=[140 145];
Yrange=[37 40];

% single frequencies to be studies: ftest = [1 0.75 0.5 0.4 0.3 0.2 0.15 0.1 0.075]
ftest = 0.05;  % test freq for plotting, effective only when indexSumAllFreq = 0;
indexSumAllFreq = 1; % = 1, plot result after summing all freqs, =0, only plot test freq results
sumfreqrange = [0.5 1]; % frequency range to be summed (only work if indexSumAllFreq = 1);
indexSpatialSmooth = 0; % apply spatial smoothing for the CVX results
SmoothRadius = 20; % smoothing radius using a Gaussian smoothing function A exp(-d^2/r^2)
ntwinAverage = 1; % number of time windows for averaging the power, must be odd number,
% if ntwinAverage = 1, do not average,
% if ntwinAverage = 3, 5, 7, ..., perform averaging

indexPlotwin = 0; % = 1: plot source imaging results of each time window and make the movie; otherwise: do not plot
Index_Snapshot = 0;%indexPlotwin; % =1 plot snapshot
T_snapshot=10:10:120;



% set strike
strike = 50; % strike of the fault, from USGS WPhase solution
strikevector = [cosd(90-strike);sind(90-strike)];
if strike==0||strike==180
    faultypt = -400:2:400;
    faultxpt = zeros(1,length(faultypt));
else
    faultxpt = -400:2:400; % x-coordinates of fault line
    faultypt = faultxpt*tand(90-strike); % y-coordinates of fault line
end

%% input parameters END

%% move some necessary files to the data file
copyfile('world_borders.shp',DataPath);
copyfile('world_borders.dbf',DataPath);
copyfile('world_borders.shx',DataPath);
copyfile('plate_boundaries',DataPath);

%% change to the data path

cd(DataPath);

if IndexTopoBackground == 1
    load(ETOPOmat);
end


[FilePath]=textread(PlotTxt,'%s',-1);
CSresultdir_temp = cell2mat(FilePath(1));
FigFolder = CSresultdir_temp;
datafile = [CSresultdir_temp 'srcCVXresult.mat'];
travtbl = load([FigFolder 'DiffTravtbl.mat']); % differential travel time table of grids with respect to hypocenter



%%
load(datafile); % load data

NgridsX=length(xrange);
NgridsY=length(yrange);

II = find(lonloc > 180); lonloc(II) = lonloc(II) - 360;
if qlon > 180; qlon = qlon - 360; end

bbox = [lonloc(1) latloc(1);  lonloc(end) latloc(end)];
borders = shaperead('world_borders.shp','BoundingBox',bbox);
% coastline = shaperead('./Coastline/30arcsecond/FullWorldLines_30.shp','BoundingBox',bbox);

nlonpt = length(lonloc);
nlatpt = length(latloc);
Nsrc = nlonpt*nlatpt;
% lon & lat of the strike from monent tensor inversion
faultlat = qlat + faultypt/111.19;
faultlon = qlon + faultxpt/111.19/cosd(qlat);

neff = length(ffeff); % number of freq.

%% setting ending time in case of inaccurate ending time
if IndexSettingEnd ==1
    ptmin=beginTime;
    ptmax=endingTime;
end
%%

nsegwin = floor((ptmax - ptmin - segwint)/dtsegwin) + 1; % number of window segment for imaging
[temp, IIf] = min(abs(ffeff - ftest)); % find index of the closest freq to ftest
ftestnew = ffeff(IIf); % closest freq to ftest

% create folder to save CS and beamforming images of each time step
if indexPlotwin == 1
    load('waveform.mat') ;
    if indexSumAllFreq == 0
        CVXBeamFigfolder = [FigFolder 'CVXBeamFig_' num2str(ftestnew,3) 'Hz/'];
        mkdir(CVXBeamFigfolder);
    else
        CVXBeamFigfolder = [FigFolder 'CVXBeamFig_' num2str(sumfreqrange(1)) '-' num2str(sumfreqrange(2)) 'Hz/'];
        mkdir(CVXBeamFigfolder);
    end
end

%% read plate boundary data
fbdr = fopen('plate_boundaries', 'r');
platebdr = [];
temp = fscanf(fbdr, '%s', 1);
i = 0;
while ~feof(fbdr)
    i = i + 1;
    if strcmp(temp,'>')
        platebdr(i,1:2) = [NaN NaN];
        fgetl(fbdr);
    else
        platebdr(i,1) = str2double(temp);
        platebdr(i,2) = fscanf(fbdr, '%f', 1);
        fgetl(fbdr);
    end
    
    if i > 1 && (abs(platebdr(i,1) - platebdr(i-1,1)) > 180)
        platebdr(i+1,:) = platebdr(i,:);
        platebdr(i,:) = [NaN NaN];
        i = i + 1;
    end
    
    temp = fscanf(fbdr, '%s', 1);
end
fclose(fbdr);


%% plot CVX and beamforming results

% build up (x,y) grids and distance matrix for spatial smoothing
tic
if indexSpatialSmooth == 1
    
    xgrids = zeros(nlonpt, nlatpt);
    ygrids = zeros(nlonpt, nlatpt);
    xlocnew = zeros(Nsrc,1);
    ylocnew = zeros(Nsrc,1);
    distlocal = zeros(Nsrc, Nsrc);
    for i = 1:nlonpt
        xgrids(i,:) = xrange(i);
        ygrids(i,:) = yrange;
    end
    xlocnew = reshape(xgrids, Nsrc, 1);
    ylocnew = reshape(ygrids, Nsrc, 1);
    
    for i = 1:Nsrc
        distlocal(i,:) = sqrt((xlocnew - xlocnew(i)*ones(Nsrc,1)).^2 + (ylocnew - ylocnew(i)*ones(Nsrc,1)).^2);
    end
    
    dxgrid = xrange(2) - xrange(1);
    dygrid = yrange(2) - yrange(1);
    WeightSmooth = (dxgrid*dygrid)/(pi*SmoothRadius)*exp(-distlocal.^2/SmoothRadius^2);
    
end
toc


srcAmpGridsAll = zeros(nlonpt, nlatpt, nsegwin);
beamfreqAll = zeros(nlonpt, nlatpt, nsegwin);

SumPower = zeros(nlonpt, nlatpt);

if indexSumAllFreq == 1
    IIsumfreq = find(ffeff >= sumfreqrange(1) & ffeff <= sumfreqrange(2));
    freqrangestr = [num2str(ffeff(IIsumfreq(1)),3) '-' num2str(ffeff(IIsumfreq(end)),3)];
end


for iw = 1:nsegwin
    wint1 = ptmin + (iw - 1)*dtsegwin;  % window starting time
    wint2 = ptmin + (iw - 1)*dtsegwin + segwint; % window ending time
    display([num2str(wint1) '--' num2str(wint2) ' s']);
    
    if indexSumAllFreq == 1
        if indexSpatialSmooth ~= 1 % without spatial smoothing before summation over freqs.
            % plot src image by stacking results over all frequencies
            srcAmpGrids = squeeze(sum(srcGridSpec(:,:,IIsumfreq,iw).^2, 3));  % for CVX
            srcAmpGrids = sqrt(srcAmpGrids);
            beamfreq = squeeze(sum(srcGridBeam(:,:,IIsumfreq,iw), 3).^2);  % for beamforming
            beamfreq = sqrt(beamfreq);
        else   % with spatial smoothing before summation over freqs.
            srcAmpGridsReshape = zeros(Nsrc,1);
            for ifreq = IIsumfreq(1):IIsumfreq(end)
                srcAmpGridsReshape = srcAmpGridsReshape + WeightSmooth*reshape(srcGridSpec(:,:,ifreq,iw).^2, Nsrc, 1);
            end
            srcAmpGrids = reshape(srcAmpGridsReshape, nlonpt, nlatpt);
            srcAmpGrids = sqrt(srcAmpGrids);
            beamfreq = squeeze(sum(srcGridBeam(:,:,IIsumfreq,iw), 3).^2);  % for beamforming
            beamfreq = sqrt(beamfreq);
            
            SumPower = SumPower + srcAmpGrids.^2;
            
        end
        
    else % only for single freq results at ftestnew
        if indexSpatialSmooth ~= 1
            srcAmpGrids = squeeze(srcGridSpec(:,:,IIf,iw));
            beamfreq = squeeze(srcGridBeam(:,:,IIf,iw));
        else
            srcAmpGridsReshape = WeightSmooth*reshape(srcGridSpec(:,:,IIf,iw).^2, Nsrc, 1);
            srcAmpGrids = reshape(srcAmpGridsReshape, nlonpt, nlatpt);
            srcAmpGrids = sqrt(srcAmpGrids);
            beamfreq = squeeze(srcGridBeam(:,:,IIf,iw));  % for beamforming
            
            SumPower = SumPower + srcAmpGrids.^2;
        end
    end
    
    
    srcAmpGridsAll(:,:,iw) = reshape(srcAmpGrids, nlonpt, nlatpt, 1);
    beamfreqAll(:,:,iw) = reshape(beamfreq, nlonpt, nlatpt, 1);
    
end

%% perform time averaging for the power over nearby time windows and find
%% the local peaks for each time-averaged power

beampeakinfo = zeros(nsegwin, 5);
cvxpeakinfo = zeros(nsegwin, 5); % peak of source amplitude
cvx2ndpeakinfo = zeros(nsegwin, 5); % 2nd peak of source amplitude

cvxpeakpossition = zeros(nsegwin, 2);% possiton vectors of the peak points, used for estimating of the rupture velocity.

% perform averaging
if ntwinAverage >= 3 % average in time
    display('average power in nearby windows ...');
    srcAmpGridsAllNew = zeros(nlonpt, nlatpt, nsegwin);
    beamfreqAllNew = zeros(nlonpt, nlatpt, nsegwin);
    runwinfilter = ones(1,ntwinAverage)/ntwinAverage;
    nhalfwin = (ntwinAverage - 1)/2;
    for ilon = 1:nlonpt
        for jlat = 1:nlatpt
            temp = filter(runwinfilter, 1, [reshape(srcAmpGridsAll(ilon,jlat,:).^2,nsegwin,1); zeros(nhalfwin,1)]);
            srcAmpGridsAllNew(ilon,jlat,:) = sqrt(temp((1+nhalfwin):end)); % power to amplitude
            temp = filter(runwinfilter, 1, [reshape(beamfreqAll(ilon,jlat,:).^2,nsegwin,1); zeros(nhalfwin,1)]);
            beamfreqAllNew(ilon,jlat,:) = sqrt(temp((1+nhalfwin):end)); % power to amplitude
        end
    end
end

% find local peaks for CS and beamforming results
for iw = 1:nsegwin
    wint1 = ptmin + (iw - 1)*dtsegwin;  % window starting time
    wint2 = ptmin + (iw - 1)*dtsegwin + segwint; % window ending time
    display([num2str(wint1) '--' num2str(wint2) ' s']);
    
    if ntwinAverage >= 3
        beamfreq = squeeze(beamfreqAllNew(:,:,iw));
        srcAmpGrids = squeeze(srcAmpGridsAllNew(:,:,iw));
    elseif ntwinAverage == 1
        beamfreq = squeeze(beamfreqAll(:,:,iw));
        srcAmpGrids = squeeze(srcAmpGridsAll(:,:,iw));
    end
    
    
    nmax = 2;
    neighborsize = [3 3];
    irowcirc = 0;
    icolcirc = 0;
    [ind, locmax, count] = findlocalmax2(srcAmpGrids, nmax, neighborsize, irowcirc, icolcirc);
    % note: locmax has been ordered from large to small
    
    if count > 0
        cvxpeakinfo(iw,1) = (wint1 + wint2)/2;  % window center time
        cvxpeakinfo(iw,3) = ind(1,1);  % lon grid point index
        cvxpeakinfo(iw,4) = ind(2,1);  % lat grid point index
        cvxpeakinfo(iw,5) = locmax(1);       % cvx solution max amplitude
    end
    
    if count == 2 && locmax(2)/locmax(1) > 0 % 0.5
        cvx2ndpeakinfo(iw,1) = (wint1 + wint2)/2;  % window center time
        cvx2ndpeakinfo(iw,3) = ind(1,2);  % lon grid point index
        cvx2ndpeakinfo(iw,4) = ind(2,2);  % lat grid point index
        cvx2ndpeakinfo(iw,5) = locmax(2);       % cvx solution max amplitude
    end
    
    [maxBeam, IIlonMaxBeam] = max(max(beamfreq,[],2));
    [maxBeam, IIlatMaxBeam] = max(max(beamfreq,[],1));
    beampeakinfo(iw,1) = (wint1 + wint2)/2;  % window center time
    beampeakinfo(iw,3) = IIlonMaxBeam;  % lon grid point index
    beampeakinfo(iw,4) = IIlatMaxBeam;  % lat grid point index
    beampeakinfo(iw,5) = maxBeam;       % beamforming peak amplitude
    
    cvxpeakpossition(iw,1) = (cvxpeakinfo(iw,3)-37)*(xrange(2)-xrange(1));
    cvxpeakpossition(iw,2) = (cvxpeakinfo(iw,4)-37)*(yrange(2)-yrange(1));
    MaxAmpSrc = max(cvxpeakinfo(:, 5));
end

%% estimating the rupture velocity
projectDistance = cvxpeakpossition*strikevector;
hold off
figure(99)
clf
for iw = 1:nsegwin;
    T = (iw-1) * dtsegwin + (segwint + ptmin)/2;
    figure(99)
    if cvxpeakinfo(iw,5) > 0
        hold on
        plot(T,projectDistance(iw),'o','MarkerSize',25*cvxpeakinfo(iw,5)/MaxAmpSrc,'MarkerEdgeColor','r');
    end
end

clear EpiDistance;
for Velocity=1.6:0.5:1.6
    Time=10:180;
    EpiDistance=(Time-10)*Velocity+25;
    hold on
    plot(Time,EpiDistance,'--b');
end
clear EpiDistance;
clear Time
EpiDistance=zeros(181,1);
Time=0:1:180;
plot(Time,EpiDistance,'-k','LineWidth',2);
axis([0 90 0 150]);
daspect([1,3,1])
ylabel('Distance from epicenter');
xlabel('time(s)');


if indexSumAllFreq ==1
    title(['Distance from epicenter at f=['  num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) 'Hz]'])
    grid on
    hold off
    filename = [FigFolder 'Distance from epicenter at f=['  num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) 'Hz].pdf'];
else
    title(['Distance from epicenter at f=[' num2str(ftestnew,3) 'Hz]']);
    grid on
    hold off
    filename = [FigFolder 'Distance from epicenter at f=[' num2str(ftestnew,3) 'Hz.pdf]'];
end
print('-dpdf', filename);




%% plot the frequency-summed power of the sources
figure;
if indexSpatialSmooth == 1
    %cmap = colormap(pink);
    cmap = colormap(bone);
    cmap = cmap(end:-1:1,:);
    maxpower = max(max(SumPower));
    imagesc(lonloc,latloc,SumPower'/maxpower), colorbar; colormap(cmap);
    hold on;
    v = [0.2 0.4 0.6 0.9];
    [C, h] = contour(lonloc, latloc, SumPower'/maxpower, v,'b-','LineWidth',2); axis equal
    text_handle = clabel(C,h,'FontSize',15,'FontWeight','bold','Color','b');
    
    for i = 1:size(borders)
        hold on; plot(borders(i).X, borders(i).Y, 'k');
    end
    hold on; plot(platebdr(:,1), platebdr(:,2),'b-', 'LineWidth',2); % plot plate boundaries
    hold on; plot(faultlon, faultlat, 'b--', 'LineWidth',1);
    hold on; plot(qlon, qlat, 'm+', 'LineWidth',2,'MarkerSize',18);
    %     for i = 1:length(IIafter)
    %         k = IIafter(i);
    %         hold on; plot(aftershockUSGS(k,8), aftershockUSGS(k,7), 'ko', 'MarkerSize', aftershockUSGS(k,end));
    %     end
    xlim([lonloc(1) lonloc(end)]);
    ylim([latloc(1) latloc(end)]);
    
    if IndexSettingRange ==1
        ylim(Yrange);
        xlim(Xrange);
    end
    
    daspect([1 cosd(qlat) 1]); 
    set(gca,'YDir','normal','FontSize',15); grid on; box on;
    
    if indexSumAllFreq == 0  % single frequency
        title(['CVX: f = ' num2str(ftestnew,3) 'Hz, Max Power = ', num2str(maxpower)]);
        print('-dpdf', [FigFolder 'CVX_SumPower_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_smooth' num2str(SmoothRadius) '.pdf']);
    else % in a frequency range
        title(['CVX: f = [' num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) '] Hz, Max Power = ', num2str(maxpower)]);
        print('-dpdf', [FigFolder 'CVX_SumPower_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win_smooth' num2str(SmoothRadius) '.pdf']);
    end
end

%% calculate the adjusted source time for each time window considering the
%% travel time difference between hypocenter to stations and each grid
%% point to stations (using input differential travel time table ...)

for iw = 1:nsegwin
    if cvxpeakinfo(iw,5) > 0
        cvxmeddifftime = median(travtbl.tpredshiftgrids(cvxpeakinfo(iw,3), cvxpeakinfo(iw,4), :));
        cvxpeakinfo(iw,2) = cvxpeakinfo(iw,1) - cvxmeddifftime;
    end
    
    if beampeakinfo(iw,5) > 0
        beammeddifftime = median(travtbl.tpredshiftgrids(beampeakinfo(iw,3), beampeakinfo(iw,4), :));
        beampeakinfo(iw,2) = beampeakinfo(iw,1) - beammeddifftime;
    end
    
    
    if cvx2ndpeakinfo(iw,5) > 0
        cvx2ndmeddifftime = median(travtbl.tpredshiftgrids(cvx2ndpeakinfo(iw,3), cvx2ndpeakinfo(iw,4), :));
        cvx2ndpeakinfo(iw,2) = cvx2ndpeakinfo(iw,1) - cvx2ndmeddifftime;
    end
end


%% plot time-space distribution of sources
figure;
cmap = colormap(jet);
% cmap = cmap(end:-1:1,:);
ncol = size(cmap,1);

maxmidtime = 600;
IItime = find(cvxpeakinfo(:, 1) <= maxmidtime);

figure(100); hold off;
MaxAmpSrc = max(cvxpeakinfo(:, 5));
MaxTimeSrc = ceil(max(cvxpeakinfo(IItime, 2)));

% subplot(1,2,1);
firstsrcinfo = [];
secondsrcinfo = [];
kk1 = 0;
kk2 = 0;



for i = 1:size(borders)
    hold on; plot(borders(i).X, borders(i).Y, 'k');
end
hold on; plot(platebdr(:,1), platebdr(:,2),'b-', 'LineWidth',2); % plot plate boundaries
hold on; plot(faultlon, faultlat, 'b--', 'LineWidth',1);



for i = 1:nsegwin
    % plot largest peak (not on the boundaries)
    isbdr = (cvxpeakinfo(i,3) == 1 || cvxpeakinfo(i,3) == nlonpt || cvxpeakinfo(i,4) == 1 ||  cvxpeakinfo(i,4) == nlatpt);
    if cvxpeakinfo(i,3) > 0 && cvxpeakinfo(i,1) <= maxmidtime && cvxpeakinfo(i,5)/MaxAmpSrc > 0.0 && (~isbdr)
        colorindex = round((cvxpeakinfo(i,2)/MaxTimeSrc)*(ncol-1)) + 1;
        if colorindex <= 0
            colorindex = 1;
        elseif colorindex > ncol
            colorindex = ncol;
        end
%         plot(lonloc(cvxpeakinfo(i,3)), latloc(cvxpeakinfo(i,4)), 'o', 'MarkerEdgeColor', cmap(colorindex,:), ...
%              'LineWidth',3,'MarkerSize', 30*cvxpeakinfo(i,5)/MaxAmpSrc );
        plot(lonloc(cvxpeakinfo(i,3)), latloc(cvxpeakinfo(i,4)), 'o', 'MarkerEdgeColor', 'k','MarkerFaceColor', ...
            cmap(colorindex,:),'LineWidth', 1, 'MarkerSize', 30*cvxpeakinfo(i,5)/MaxAmpSrc );
        hold on;
        kk1 = kk1 + 1;
        firstsrcinfo(kk1, 1:4) = [lonloc(cvxpeakinfo(i,3)) latloc(cvxpeakinfo(i,4)) cvxpeakinfo(i,2) cvxpeakinfo(i,5)];
    end
    
        % plot second largest peak (not on the boundaries)
        isbdr = (cvx2ndpeakinfo(i,3) == 1 || cvx2ndpeakinfo(i,3) == nlonpt || cvx2ndpeakinfo(i,4) == 1 ||  cvx2ndpeakinfo(i,4) == nlatpt);
        if cvx2ndpeakinfo(i,3) > 0 && cvx2ndpeakinfo(i,1) <= maxmidtime && cvx2ndpeakinfo(i,5)/MaxAmpSrc > 0.0 && (~isbdr)
            colorindex = round((cvx2ndpeakinfo(i,2)/MaxTimeSrc)*(ncol-1)) + 1;
            if colorindex <= 0
                colorindex = 1;
            elseif colorindex > ncol
                colorindex = ncol;
            end
            plot(lonloc(cvx2ndpeakinfo(i,3)), latloc(cvx2ndpeakinfo(i,4)), 'o', 'MarkerEdgeColor', cmap(colorindex,:), ...
                'MarkerFaceColor', cmap(colorindex,:), 'MarkerSize', 16*cvx2ndpeakinfo(i,5)/MaxAmpSrc );
            hold on;
            kk2 = kk2 + 1;
            secondsrcinfo(kk2, 1:4) = [lonloc(cvxpeakinfo(i,3)) latloc(cvxpeakinfo(i,4)) cvxpeakinfo(i,2) cvxpeakinfo(i,5)];
        end
    
end

hold on; plot(qlon, qlat, 'm+', 'LineWidth',2,'MarkerSize',18);


hold on;

set(gca,'FontSize',15);


colorbar; colormap(cmap); caxis([0 MaxTimeSrc]);
xlim([lonloc(1) lonloc(end)]);
ylim([latloc(1) latloc(end)]);

if IndexSettingRange ==1
    ylim(Yrange);
    xlim(Xrange);
end

daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal'); grid on;

if indexSumAllFreq==1
    title(['CVX Peak Time: f = [' num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) '] Hz']);
else
    title(['CVX Peak Time: f = [' num2str(ftestnew,3) '] Hz']);
end
box on

if indexSumAllFreq == 0
    if indexSpatialSmooth ~= 1
        filename = [FigFolder 'CVX_Peaktime_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win'];
        srcmatfile = [FigFolder 'CVX_peakinfo_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win.mat'];
    else
        filename = [FigFolder 'CVX_Peaktime_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius)];
        srcmatfile = [FigFolder 'CVX_peakinfo_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '.mat'];
    end
elseif indexSumAllFreq == 1
    if indexSpatialSmooth ~= 1
        filename = [FigFolder 'CVX_Peaktime_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win'];
        srcmatfile = [FigFolder 'CVX_peakinfo_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win.mat'];
    else
        filename = [FigFolder 'CVX_Peaktime_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '_peaks'];
        srcmatfile = [FigFolder 'CVX_peakinfo_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '_peaks.mat'];
    end
end
save(srcmatfile, 'ftestnew', 'firstsrcinfo', 'secondsrcinfo'); % save peaks information to mat file
print('-dpdf', [filename '.pdf']);
if indexSumAllFreq == 1
    save([FigFolder 'CVX_peakinfo_f' freqrangestr 'Hz.txt'],'-ascii','firstsrcinfo');
end

%% Plot the CVX results with the topography background
%--------------------------------------------------------------------------------------------------------------------------------------
if IndexTopoBackground==1
    figure(101); hold off;
    
    firstsrcinfo = [];
    secondsrcinfo = [];
    kk1 = 0;
    kk2 = 0;
    
    imagesc(maplon,maplat,mapTopo);
    axis xy;
    hold on
    colormap(colormap_etopo);
    %colormap(gray)
    
    for i = 1:size(borders)
        hold on; plot(borders(i).X, borders(i).Y, 'k');
    end
    hold on; plot(platebdr(:,1), platebdr(:,2),'k-', 'LineWidth',2); % plot plate boundaries
    hold on; plot(faultlon, faultlat, 'k--', 'LineWidth',1);
    
    
    for i = 1:nsegwin
        % plot largest peak (not on the boundaries)
        isbdr = (cvxpeakinfo(i,3) == 1 || cvxpeakinfo(i,3) == nlonpt || cvxpeakinfo(i,4) == 1 ||  cvxpeakinfo(i,4) == nlatpt);
        if cvxpeakinfo(i,3) > 0 && cvxpeakinfo(i,1) <= maxmidtime && cvxpeakinfo(i,5)/MaxAmpSrc > 0.0 && (~isbdr)
            colorindex = round((cvxpeakinfo(i,2)/MaxTimeSrc)*(ncol-1)) + 1;
            if colorindex <= 0
                colorindex = 1;
            elseif colorindex > ncol
                colorindex = ncol;
            end
            hold on
            %             plot(lonloc(cvxpeakinfo(i,3)), latloc(cvxpeakinfo(i,4)), 'o', 'MarkerEdgeColor', cmap(colorindex,:), ...
            %                 'MarkerFaceColor', cmap(colorindex,:), 'MarkerSize', 16*cvxpeakinfo(i,5)/MaxAmpSrc );
            plot(lonloc(cvxpeakinfo(i,3)), latloc(cvxpeakinfo(i,4)), 'o', 'MarkerEdgeColor','k','MarkerFaceColor', cmap(colorindex,:), ...
                'LineWidth', 1, 'MarkerSize', 25*cvxpeakinfo(i,5)/MaxAmpSrc );
            hold on;
            kk1 = kk1 + 1;
            firstsrcinfo(kk1, 1:4) = [lonloc(cvxpeakinfo(i,3)) latloc(cvxpeakinfo(i,4)) cvxpeakinfo(i,2) cvxpeakinfo(i,5)];
        end
        
        
    end
    
    hold on; plot(qlon, qlat, 'm+', 'LineWidth',2,'MarkerSize',18);
    
    hold on;
    
    
    
    colorbar;% colormap(cmap); caxis([0 MaxTimeSrc]);
    xlim([lonloc(1) lonloc(end)]);
    ylim([latloc(1) latloc(end)]);
    daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal'); grid on;
    
    if IndexSettingRange ==1
        ylim(Yrange);
        xlim(Xrange);
    end
    
    if indexSumAllFreq==1
        title(['CVX Peak Time: f = [' num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) '] Hz']);
    else
        title(['CVX Peak Time: f = [' num2str(ftestnew,3) '] Hz']);
    end
    
    
    
    %    save(srcmatfile, 'ftestnew', 'firstsrcinfo', 'secondsrcinfo'); % save peaks information to mat file
    print('-dpdf', [filename '_topo.pdf']);
end
%% ----------------------------------------------------------------------------------------------------------------------------------------------


%% Beamforming results ploting
figure(200); hold off;

MaxAmpSrc = max(beampeakinfo(:, 5));
% MaxTimeSrc = ceil(max(beampeakinfo(:, 2)));
MaxTimeSrc = ceil(max(beampeakinfo(IItime, 2)));

for i = 1:size(borders)
    hold on; plot(borders(i).X, borders(i).Y, 'k');
end
hold on; plot(platebdr(:,1), platebdr(:,2),'b-', 'LineWidth',2); % plot plate boundaries
hold on; plot(faultlon, faultlat, 'b--', 'LineWidth',1);


for i = 1:nsegwin
    isbdr = (beampeakinfo(i,3) == 1 || beampeakinfo(i,3) == nlonpt || beampeakinfo(i,4) == 1 ||  beampeakinfo(i,4) == nlatpt);
    if beampeakinfo(i,3) > 0 && beampeakinfo(i,1) <= maxmidtime && (~isbdr)
        colorindex = round((beampeakinfo(i,2)/MaxTimeSrc)*(ncol-1)) + 1;
        if colorindex <= 0
            colorindex = 1;
        elseif colorindex > ncol
            colorindex = ncol;
        end
        %         plot(lonloc(beampeakinfo(i,3)), latloc(beampeakinfo(i,4)), 'o', 'MarkerEdgeColor', cmap(colorindex,:), ...
        %             'MarkerFaceColor', cmap(colorindex,:), 'MarkerSize', 16*beampeakinfo(i,5)/MaxAmpSrc );
        plot(lonloc(beampeakinfo(i,3)), latloc(beampeakinfo(i,4)), 'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',cmap(colorindex,:), ...
            'LineWidth', 1, 'MarkerSize', 25*beampeakinfo(i,5)/MaxAmpSrc );
        hold on;
    end
end


hold on; plot(qlon, qlat, 'm+', 'LineWidth',2,'MarkerSize',18);
hold on;



colorbar; colormap(cmap); caxis([0 MaxTimeSrc]);
xlim([lonloc(1) lonloc(end)]);
ylim([latloc(1) latloc(end)]);
daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal'); grid on; box on;

if IndexSettingRange ==1
    ylim(Yrange);
    xlim(Xrange);
end


if indexSumAllFreq==1
    title(['Beamforming Peak Time: f = [' num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) '] Hz']);
else
    title(['Beamforming Peak Time: f = [' num2str(ftestnew,3) '] Hz']);
end

if indexSumAllFreq == 0
    filename = [FigFolder 'Beam_Peaktime_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius)];
elseif indexSumAllFreq == 1
    filename = [FigFolder 'Beam_Peaktime_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius)];
end
print('-dpdf', [filename '.pdf']);

%% Plot the beamforming results with topograpy background
%--------------------------------------------------------------------------------------------------------------------------------
if IndexTopoBackground == 1
    figure(201); hold off;
    
    imagesc(maplon,maplat,mapTopo);
    axis xy;
    hold on
    colormap(colormap_etopo);
    %colormap(gray);
    
    
    for i = 1:size(borders)
        hold on; plot(borders(i).X, borders(i).Y, 'k');
    end
    hold on; plot(platebdr(:,1), platebdr(:,2),'k-', 'LineWidth',2); % plot plate boundaries
    hold on; plot(faultlon, faultlat, 'k--', 'LineWidth',1);
    
    
    for i = 1:nsegwin
        isbdr = (beampeakinfo(i,3) == 1 || beampeakinfo(i,3) == nlonpt || beampeakinfo(i,4) == 1 ||  beampeakinfo(i,4) == nlatpt);
        if beampeakinfo(i,3) > 0 && beampeakinfo(i,1) <= maxmidtime && (~isbdr)
            colorindex = round((beampeakinfo(i,2)/MaxTimeSrc)*(ncol-1)) + 1;
            if colorindex <= 0
                colorindex = 1;
            elseif colorindex > ncol
                colorindex = ncol;
            end
            hold on;
            %         plot(lonloc(beampeakinfo(i,3)), latloc(beampeakinfo(i,4)), 'o', 'MarkerEdgeColor', cmap(colorindex,:), ...
            %             'MarkerFaceColor', cmap(colorindex,:), 'MarkerSize', 16*beampeakinfo(i,5)/MaxAmpSrc );
            plot(lonloc(beampeakinfo(i,3)), latloc(beampeakinfo(i,4)), 'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',cmap(colorindex,:), ...
                'LineWidth', 1.5, 'MarkerSize', 25*beampeakinfo(i,5)/MaxAmpSrc );
            hold on;
        end
    end
    
    hold on; plot(qlon, qlat, 'm+', 'LineWidth',2,'MarkerSize',18);
    hold on;
    
    
    
    colorbar; %colormap(cmap); caxis([0 MaxTimeSrc]);
    xlim([lonloc(1) lonloc(end)]);
    ylim([latloc(1) latloc(end)]);
    
    if IndexSettingRange ==1
        ylim(Yrange);
        xlim(Xrange);
    end
    
    daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal'); grid on;
    
    if indexSumAllFreq==1
        title(['Beamforming Peak Time: f = [' num2str(ffeff(IIsumfreq(1)),3) '  ' num2str(ffeff(IIsumfreq(end)),3) '] Hz']);
    else
        title(['Beamforming Peak Time: f = [' num2str(ftestnew,3) '] Hz']);
    end
    
    
    print('-dpdf', [filename '_topo.pdf']);
    
end
%% ------------------------------------------------------------------------------------------------------------------------------



%% analyze the source power at different frequencies
figure(300);
srcpowfreqtime = zeros(nsegwin, neff);
Nsrc = nlonpt*nlatpt;
for iw = 1:nsegwin
    for mf = 1:neff
        srcpowfreqtime(iw,mf) = sum(reshape(srcGridSpec(:,:,mf,iw),Nsrc,1).^2);
    end
end
wint = ptmin + ((1:nsegwin) - 0.5)*dtsegwin;
imagesc(ffeff, wint, 10*log10(srcpowfreqtime)); colorbar;
maxdb = 10*log10(max(max(srcpowfreqtime)));
mindb = maxdb - 40;
caxis([mindb maxdb]);
xlabel('f (Hz)');
ylabel('Time (s)');
title('Source power (dB)');
freqrangestr = [num2str(ffeff(1),3) '-' num2str(ffeff(end),3)];
filename = [FigFolder 'SrcPower_dB_f' freqrangestr 'Hz.pdf'];
print('-dpdf', filename);



if indexPlotwin == 1
    h = figure (88);
    cmap = colormap(bone);
    cmap = cmap(end:-1:1,:);
    cmap=colormap(jet);
    PaperWidLength = round([1.0 0.6]*800);
    set(h,'Position',[100 100 PaperWidLength]);
    
    mov(1:nsegwin) = struct('cdata', [],...
        'colormap', []);
    axis tight
    set(gca,'nextplot','replacechildren');
    hold off;
    
    if ntwinAverage >= 3
        beamfreqMax = max(max(max(abs(beamfreqAllNew))));
        srcAmpGridsMax = max(max(max(abs(srcAmpGridsAllNew))));
    elseif ntwinAverage == 1
        beamfreqMax = max(max(max(abs(beamfreqAll))));
        srcAmpGridsMax = max(max(max(abs(srcAmpGridsAll))));
    end
    
    
    for iw = 1:nsegwin
        wint1 = ptmin + (iw - 1)*dtsegwin;  % window starting time
        wint2 = ptmin + (iw - 1)*dtsegwin + segwint; % window ending time
        display([num2str(wint1) '--' num2str(wint2) ' s']);
        
        if ntwinAverage >= 3
            beamfreq = squeeze(beamfreqAllNew(:,:,iw));
            srcAmpGrids = squeeze(srcAmpGridsAllNew(:,:,iw));
            
            beamfreq1=beamfreq;%/beamfreqMax;
            srcAmpGrids1=srcAmpGrids;%/srcAmpGridsMax;
            
            
        elseif ntwinAverage == 1
            beamfreq = squeeze(beamfreqAll(:,:,iw));
            srcAmpGrids = squeeze(srcAmpGridsAll(:,:,iw));
            
            beamfreq1=beamfreq;%/beamfreqMax;
            srcAmpGrids1=srcAmpGrids;%/srcAmpGridsMax;
            
        end
        
        
        
        clf
        subplot(3,2,1:4);
        imagesc(lonloc,latloc,srcAmpGrids1'), colorbar; colormap(cmap);
        %caxis([0 max(max(srcAmpGrids1))]);
        caxis([0 srcAmpGridsMax]);
        %     if indexSumAllFreq == 1
        %         caxis([0 max(cvxpeakinfo(:,5))]);
        %     end
        % hold on; plot(coast.long, coast.lat, 'w', 'LineWidth',1);
        for i = 1:size(borders) % plot border
            hold on; plot(borders(i).X, borders(i).Y, 'w');
        end
        hold on; plot(platebdr(:,1), platebdr(:,2),'w--', 'LineWidth',2); % plot plate boundaries
        %hold on; plot(faultlon, faultlat, 'y--', 'LineWidth',0.5);
        hold on; plot(qlon, qlat, 'yp', 'MarkerSize',12,'MarkerFaceColor','y');
%         if cvxpeakinfo(iw,5) > 0
%             hold on; plot(lonloc(cvxpeakinfo(iw,3)), latloc(cvxpeakinfo(iw,4)), 'm+', 'MarkerSize',9);
%         end
        %         if cvx2ndpeakinfo(iw,5) > 0
        %             hold on; plot(lonloc(cvx2ndpeakinfo(iw,3)), latloc(cvx2ndpeakinfo(iw,4)), 'm+', 'MarkerSize',6);
        %         end
        if IndexSettingRange ==1
            ylim(Yrange);
            xlim(Xrange);
        end
        
        text(-71.3,-30.3,['t=' num2str(cvxpeakinfo(iw,2),3)  ' s'],'Color','w','FontSize',15);
        
        daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal'); grid on;

                
        subplot(3,2,5:6)
        imagesc(ptmin:dt:ptmax, 1:ntrGood, newdata2, climnewdata); colorbar;
        hold on;
        plot([80 80],[0.5 ntrGood-1],'-m','LineWidth',3);
        
        twinBegin=beampeakinfo(iw,1)-segwint/2;
        twinEnd=beampeakinfo(iw,1)+segwint/2;
        
        
        hold on
        plot([twinBegin twinBegin],[1 ntrGood],'-','LineWidth',2);
        hold on
        plot([twinEnd twinEnd],[1 ntrGood],'-','LineWidth',2);
        colormap(jet);
        set(gca,'YDir','normal');
        
        if indexSumAllFreq == 0
            if indexSpatialSmooth ~= 1
                print('-dpdf', [CVXBeamFigfolder 'srcImage_t_' num2str(wint1) '-' num2str(wint2) 's_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win.pdf']);
            else
                print('-dpdf', [CVXBeamFigfolder 'srcImage_t_' num2str(wint1) '-' num2str(wint2) 's_f' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '.pdf']);
            end
        elseif indexSumAllFreq == 1
            if indexSpatialSmooth ~= 1
                print('-dpdf', [CVXBeamFigfolder 'srcImage_t_' num2str(wint1) '-' num2str(wint2) 's_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win.pdf']);
            else
                print('-dpdf', [CVXBeamFigfolder 'srcImage_t_' num2str(wint1) '-' num2str(wint2) 's_f' freqrangestr 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '.pdf']);
            end
        end
        % waitforbuttonpress
        % pause(0.5)
        mov(iw) = getframe(h);
    end
    
    close(88);
    
    numFrames = nsegwin;
    
    animated(1,1,1,numFrames) = 0;
    for k=1:numFrames
        if k == 1
            [animated, cmap] = rgb2ind(mov(k).cdata, 256);
        else
            animated(:,:,1,k) = ...
                rgb2ind(mov(k).cdata, cmap);
        end
    end
    
    if indexSumAllFreq == 0
        if indexSpatialSmooth ~= 1
            filename = [FigFolder 'CVXBeam_Animated_' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win.gif'];
        else
            filename = [FigFolder 'CVXBeam_Animated_' num2str(ftestnew,3) 'Hz_Ave' num2str(ntwinAverage) 'win_Rs' num2str(SmoothRadius) '.gif'];
        end
    elseif indexSumAllFreq == 1
        if indexSpatialSmooth ~= 1
            filename = [FigFolder 'CVXBeam_Animated_' freqrangestr 'Hz.gif'];
        else
            filename = [FigFolder 'CVXBeam_Animated_' freqrangestr 'Hz_Rs' num2str(SmoothRadius) '.gif'];
        end
    end
    
    
    imwrite(animated, cmap, filename, 'DelayTime', 0.5, ...
        'LoopCount', inf);
    movie2avi(mov,'./power_movie.avi','quality',100,'fps',2)
end


%% plot snapshot
if Index_Snapshot == 1
    
    h = figure (1088);
    cmap2=colormap(jet);
    PaperWidLength = round([1.0 0.6]*1500);
    set(h,'Position',[100 100 PaperWidLength]);
    NT=length(T_snapshot);
    II_t=zeros(NT);
    for IT=1:NT
%        [tmp0,II_t(IT)] = min(abs(T_snapshot(IT)-cvxpeakinfo(:,2)));
        [tmp0,II_t(IT)] = min(abs(T_snapshot(IT)-cvxpeakinfo(:,1)));
    end
    
    for IT=1:NT
        if ntwinAverage >= 3
            
            srcAmpGrids = squeeze(srcAmpGridsAllNew(:,:,II_t(IT)));
            srcAmpGrids1=srcAmpGrids;%/srcAmpGridsMax;
            
            
        elseif ntwinAverage == 1
            srcAmpGrids = squeeze(srcAmpGridsAll(:,:,II_t(IT)));
            srcAmpGrids1=srcAmpGrids;%/srcAmpGridsMax;           
        end
        
        subplot(3,4,IT)
        
        imagesc(lonloc,latloc,srcAmpGrids1'); colormap(cmap2);
        caxis([0 srcAmpGridsMax]);

        for i = 1:size(borders) % plot border
            hold on; plot(borders(i).X, borders(i).Y, 'w');
        end
        hold on; plot(platebdr(:,1), platebdr(:,2),'w--', 'LineWidth',2); % plot plate boundaries
        %hold on; plot(faultlon, faultlat, 'y--', 'LineWidth',0.5);
        hold on; plot(qlon, qlat, 'mp', 'MarkerSize',12,'MarkerFaceColor','m');
        
        if IndexSettingRange ==1
            ylim(Yrange);
            xlim(Xrange);
        end
        
%        text(-71.8,-32.2,['t=' num2str(cvxpeakinfo(II_t(IT),2),3)  ' s'],'Color','w','FontSize',10);
        text(-71.8,-32.2,['t=' num2str(cvxpeakinfo(II_t(IT),1))  ' s'],'Color','w','FontSize',10);
        
        daspect([1 cosd(qlat) 1]); set(gca,'YDir','normal','xtick',[],'ytick',[]); grid on;
        
    end
                
        

print('-dpdf','snapshots.pdf')
end


%% get back to the PlotFile
cd(PlotFile);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUST FOR CONVENIENCE!!!
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUST FOR CONVENIENCE!!!