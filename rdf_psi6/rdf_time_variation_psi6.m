%% Compare RDFs over time for a single experiment
% Joe Prisaznuk 12/16/2024

%% Load data to compare early to later video

clc; clear variables; close all

particleFolder = fullfile('..','data','particlePositionRoiColorNN'); % particle data
rdfFolder = fullfile('..','data','rdfRoi'); % pre-saved RDF data

expCode = '241211_P6W05X3'; % experiment name

load(fullfile('..','data','imageMetadata','MDall.mat')); % load all metadata

mdIdx = find(strcmp({MD.expCode}, expCode)); % find MDall index
tElap = MD(mdIdx).T.time - MD(mdIdx).T.time(1);

times = [1 10 16 19]; % which times to plot

frames = sort([1000:1000:18000, 16250, 16500]);
for ii = 1:length(frames)
    if ismember(ii,times)
        fprintf('<strong>%02d: %05d</strong> t = %s hh:mm:ss\n',ii,frames(ii),tElap(frames(ii)));
    else
        fprintf('%02d: %05d t = %s hh:mm:ss\n',ii,frames(ii),tElap(frames(ii)));
    end
end

%% Plot the points as psi_6

figure('Position',[11 61 1480 600])    

T = tiledlayout(1,length(times),'TileSpacing','tight');

cutOff = 80; % cutOff [µm]
numBins = 1000; % how many histogram bins for RDF
sbUM = 50; % scale bar length [µm]
sbHt = 8; % scale bar height [µm]
sbShift = -5; % scale bar offset from corner [µm]

frameTimes = MD(mdIdx).T.time(frames(times)) - MD(mdIdx).T.time(1);
psi6_global = zeros(length(times),1);
psi6_mag_all = cell(1,length(times));

for ii = 1:length(times)
    nexttile; hold on

    ptsFileName = sprintf('%s_%05d_colorPts.mat',expCode,frames(times(ii)));
    load(fullfile(particleFolder,ptsFileName),'pts','DT','d','NN','roiSz','roi','roiPts')

    DTroi = delaunayTriangulation(roiPts(:,1),roiPts(:,2));
    thresholdAngDeg = 120; % acute triangle threshold
    [DTroiNew, theta] = filterAcuteTriangles(DTroi, thresholdAngDeg);

    neighbors = vertexAttachments(DTroiNew); % cell array of triangle indices per point
    psi6_local = zeros(size(roiPts,1),1);
    psi6_mag = zeros(size(roiPts,1),1);
    for jj = 1:size(roiPts,1)
        triIndices = neighbors{jj};
        vertexIDs = unique(DTroiNew.ConnectivityList(triIndices,:));
        vertexIDs(vertexIDs == jj) = []; % exclude self

        angles = atan2(roiPts(vertexIDs,2) - roiPts(jj,2), ...
                       roiPts(vertexIDs,1) - roiPts(jj,1));
        psi6_local(jj) = mean(exp(1i * 6 * angles)); % local order parameter
        psi6_mag(jj) = abs(psi6_local(jj));
    end
    psi6_mag_all{ii} = psi6_mag;
    % psi6_global(ii) = mean(psi6_mag(jj)); % global sixfold order

    % set particle colors
    if ii == length(times) % for legend to be correct
        colormap viridis
        cb = colorbar; cb.Label.String = '\psi_6';
        cb.FontSize = 24;
        cb.Label.FontSize = 30;
        cb.Position = cb.Position + [0.015 0 0.005 0];
        % cb.Label.Position = cb.Label.Position + [0.01 0 0];
        cb.Label.Rotation = 0;
    end
    roiCrop = 250; % how much to shrink ROI [µm]
    R = roi + roiCrop*[1 -1 0.85 -0.85]; % shrink the ROI
    
    % Logical mask for points inside ROI
    inROI = roiPts(:,1) >= R(1) & roiPts(:,1) <= R(2) & ...
            roiPts(:,2) >= R(3) & roiPts(:,2) <= R(4);

    % Compute average ψ₆ for just those points
    psi6_global(ii) = mean(abs(psi6_local(inROI)));

    Rdim = roiSz - roiCrop*[2 2];
    ptsz = 150; % point size in plot
    triplot(DTroiNew,'LineWidth',1,'Color',[0.6 0.6 0.6 0.5])
    
    scatter(DTroiNew.Points(:,1),DTroiNew.Points(:,2), ptsz, psi6_mag','filled')
    set(gca,'YDir','reverse')
    clim([0 1])

    sbPadx = Rdim(2)/10; sbPady = Rdim(2)/4.2; 
    if ii == length(times)
        rectangle('Position',[R(2)-sbUM+sbShift, R(4)-sbHt+sbShift, sbUM, sbHt],...
            'FaceColor','k','EdgeColor','none')
    end
    
    timeText = sprintf('t = %.f min',minutes(frameTimes(ii)));
    title(timeText,'FontSize',28,'FontWeight','normal')
    set(gca,'YDir','reverse','Clipping','on')
    axis equal
    axis(R)
    xticks([]); yticks([])
    box on
    axis off
end

%% Plot psi_6 histograms

figure('Position',[11 61 1200 800])
legString = strings(4,1);
for ii = 1:length(times)
    hold on
    histogram(psi6_mag_all{ii},20,'normalization','probability','BinLimits',[0 1])
    psi6_std = std(psi6_mag_all{ii},0);
    legString(ii) = sprintf('%.f min, \\sigma = %.3f',minutes(frameTimes(ii)),psi6_std);
end
set(gca,'FontSize',24,'LineWidth',1.5)
legend(legString)
xlabel('\psi_6')
ylabel('P(\psi_6)')

%% Plot RDF

figure('Position',[11 61 920 600])

numPks = 4; % number of peaks to find
pksX = zeros(length(times),numPks); % x-value at peak locations
pksY = zeros(size(pksX)); % y-value of peak locations
legendArray = strings(1,length(times));
color = ["#8040E6";"#1AA640";"#E68000";"#A3373C"];

hold on
for ii = 1:length(times)
    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expCode,frames(times(ii)),cutOff,numBins,roiSz);
    load(fullfile(rdfFolder,rdfFile)); % load the RDF roi data
    gR = gRroi;
    % plot(gR.values,gR.histo,'color',color(ii,:),'linewidth',0.5) % raw data
    smoothHist = movmean(gR.histo,30); % apply moving average
    plot(gR.values,smoothHist,'color',color(ii,:),'linewidth',3) % plot smooth curve
    pks = findpeaks(smoothHist,'MinPeakWidth',10,'MinPeakHeight',0.5);
    for jj = 1:numPks
        pksIdx = find(smoothHist == pks(jj)); % find index of peaks
        pksX(ii,jj) = gR.values(pksIdx);  % find x-value of peaks
        pksY(ii,jj) = smoothHist(pksIdx); % find y-value of peaks
    end
    ratios = pksX(ii,:)/pksX(ii,1); % get peak ratios
    legendArray(ii) = sprintf('t = %.f min',minutes(frameTimes(ii)));
end

for ii = 1:length(times)
    for jj = 1 % for jj = 1:numPks
        plot([pksX(ii,jj),pksX(ii,jj)],[0 pksY(ii,jj)], ...
            '--','color',color(ii,:),'linewidth',1) % plot x-value
    end
    legendArray(ii+length(times)) = sprintf('%.1f µm',pksX(ii,jj));
end

axis([0 50 0 2.2])
set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on')
xlh = xlabel('r (μm)');
ylh = ylabel('g(r)');
grid off; box on
legend(legendArray,'FontSize',22,'NumColumns',2,'Box','off')
title('Radial Distribution Function')

%% plot RDFs raw

figure('Position',[11 61 1200 800])

ss = string(NaN(1,numPks+1)); % spacing string
legendArray = strings(1,length(times)*length(ss));

hold on
for ii = 1:length(times)
    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expCode,frames(times(ii)),cutOff,numBins,roiSz);
    load(fullfile(rdfFolder,rdfFile)); % load the RDF roi data
    gR = gRroi;
    plot(gR.values,gR.histo,'color',[hex2rgb(color(ii,:)) 0.5],'linewidth',0.5) % raw data
    smoothHist = movmean(gR.histo,30); % apply moving average
    plot(gR.values,smoothHist,'color',color(ii,:),'LineWidth',2) % plot smooth curve
    pks = findpeaks(smoothHist,'MinPeakWidth',10,'MinPeakHeight',0.5);
    for jj = 1:numPks
        pksIdx = find(smoothHist == pks(jj)); % find index of peaks
        pksX(ii,jj) = gR.values(pksIdx); % find x-value of peaks
        plot([pksX(ii,jj),pksX(ii,jj)],[0 smoothHist(pksIdx)], ...
            '--','color',color(ii,:),'linewidth',1) % plot x-value
    end
    ratios = pksX(ii,:)/pksX(ii,1); % get peak ratios
    label = sprintf('t_{%d} [%d, %.2f, %.2f, %.2f]',times(ii),ratios); % create legend label
    legendArray((numPks+2)*ii-numPks-1:(numPks+2)*ii) = [ss(1), label, ss(2:end)];
end

set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on')
xlim([0 50])
legend(legendArray)
title('Radial Distribution Function')
xlabel('r (µm)');
ylabel('g(r)');
grid on

%% plot more RDFs

figure('Position',[11 61 1200 800])

time = MD(mdIdx).T.time(1:max(frames)) - MD(mdIdx).T.time(1);
adjFrames = cumsum([0 MD(mdIdx).numFrames(1:end-1)'])';
timeIdx = MD(mdIdx).T.frame(frames) - MD(mdIdx).T.frame(1) + adjFrames(MD(mdIdx).T.vidIdx(frames));
colorAll = turbo(length(time));

hold on
for ii = 1:length(frames)
    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expCode,frames(ii),cutOff,numBins,roiSz);
    load(fullfile(rdfFolder,rdfFile)); % load the RDF roi data
    gR = gRroi;
    smoothHist = movmean(gR.histo,25); % apply moving average
    plot(gR.values,smoothHist,'-','color',[colorAll(timeIdx(ii),:) 1],'LineWidth',1.5) % plot smooth curve
    label = sprintf('$t_{%d}$',ii); % create legend label
end
colormap(turbo)
tFinal = hours(max(time(timeIdx)));
cH = colorbar; clim([0 tFinal])
cH.Label.String = 'Time (hours)';

axis([0 40 0 2.2])
xlabel('r (µm)');
ylabel('g(r)');
grid off; box on;
set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on')

%% plot more RDFs (3D)

figure('Position',[11 61 1200 800])

time = MD(mdIdx).T.time(1:max(frames)) - MD(mdIdx).T.time(1);
adjFrames = cumsum([0 MD(mdIdx).numFrames(1:end-1)'])';
timeIdx = MD(mdIdx).T.frame(frames) - MD(mdIdx).T.frame(1) + adjFrames(MD(mdIdx).T.vidIdx(frames));
colorAll = turbo(length(time));

hold on
for ii = 1:length(frames)
    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expCode,frames(ii),cutOff,numBins,roiSz);
    load(fullfile(rdfFolder,rdfFile)); % load the RDF roi data
    gR = gRroi;
    smoothHist = movmean(gR.histo,30); % apply moving average
    currentTimeArray = minutes(time(timeIdx(ii))*ones(size(smoothHist)));
    if ismember(ii,times)
        plot3(gR.values,currentTimeArray,smoothHist,'-',...
            'color',[colorAll(timeIdx(ii),:) 1],'LineWidth',2.5) % plot smooth curve
    else
        plot3(gR.values,currentTimeArray,smoothHist,'-',...
            'color',[colorAll(timeIdx(ii),:) 1],'LineWidth',0.5) % plot smooth curve
    end
    label = sprintf('$t_{%d}$',ii); % create legend label
end
colormap(turbo)
tFinal = hours(max(time(timeIdx)));
cH = colorbar; clim([0 tFinal])
cH.Label.String = 'Time (hours)';

pbaspect([1.5 1.5 1])
xlim([0 40])
view(3)
xlabel('r (µm)');
ylh = ylabel('time (minutes)'); ylh.Rotation = -37.5;
yticks([1 10 100]);
zlabel('g(r)');
grid off; box on;
set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on','ZMinorTick','on')
set(gca,'YScale','log')