%% Compare RDFs vs. salt concentration
% Joe Prisaznuk 4/16/2025

%% Load data

clc; clear variables; close all

particleFolder = fullfile('..','data','particlePositionRoiColorNN'); % particle data
rdfFolder = fullfile('..','data','rdfRoi'); % pre-saved RDF data
load(fullfile('..','data','imageMetadata','MDall.mat')); % load all metadata

expList = ["241205_P6W00X3", "240501_P6W01X1", "241211_P6W05X3"];
saltConc = double(regexp(expList, '(?<=W)\d{2}', 'match', 'once')); % NaCl concentration [mM]
targetFrame = [1300 1300 5000];

imScale = 2; % image resize factor
scale = 1.848; % image scale factor [µm/px]
sbUM = 50; % scale bar length [µm]
sbHt = 10; % scale bar height [µm]
sbShift = -5; % scale bar offset from corner [µm]

%% Plot the points as psi_6

figure('Position',[11 61 1480 470])    

T = tiledlayout(1,length(targetFrame),'TileSpacing','tight');

cutOff = 80; % cutOff [µm]
numBins = 1000; % how many histogram bins for RDF

psi6_global = zeros(length(targetFrame),1);

for ii = 1:length(targetFrame)
    nexttile; hold on
    
    mdIdx = find(strcmp({MD.expCode}, expList(ii))); % find MDall index
    tElap = MD(mdIdx).T.time(targetFrame(ii)) - MD(mdIdx).T.time(1);

    ptsFileName = sprintf('%s_%05d_colorPts.mat',expList(ii),targetFrame(ii));
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
    % psi6_global(ii) = mean(psi6_mag(jj)); % global sixfold order

    % set particle colors
    if ii == length(targetFrame) % for legend to be correct
        colormap viridis
        cb = colorbar; cb.Label.String = '\psi_6';
        cb.FontSize = 24;
        cb.Label.FontSize = 30;
        cb.Position = cb.Position + [0.015 0 0.005 0];
        % cb.Label.Position = cb.Label.Position + [0.01 0 0];
        cb.Label.Rotation = 0;
    end
    newRoiSz = [200 200]; % width, height
    xMid = (roi(2) + roi(1))/2; yMid = (roi(4) + roi(3))/2;
    R = [xMid-newRoiSz(1)/2, xMid+newRoiSz(1)/2,...
                 yMid-newRoiSz(2)/2, yMid+newRoiSz(2)/2];
    
    % Logical mask for points inside ROI
    inROI = roiPts(:,1) >= R(1) & roiPts(:,1) <= R(2) & ...
            roiPts(:,2) >= R(3) & roiPts(:,2) <= R(4);

    % Compute average ψ₆ for just those points
    psi6_global(ii) = mean(abs(psi6_local(inROI)));

    ptsz = 125; % point size in plot
    triplot(DTroiNew,'LineWidth',1,'Color',[0.6 0.6 0.6 0.5])
    
    scatter(DTroiNew.Points(:,1),DTroiNew.Points(:,2), ptsz, psi6_mag','filled')
    set(gca,'YDir','reverse')
    clim([0 1])

    if ii == length(targetFrame)
        rectangle('Position',[R(2)-sbUM+sbShift, R(4)-sbHt+sbShift, sbUM, sbHt],...
            'FaceColor','k','EdgeColor','none')
    end
    
    saltText = sprintf('%d mM',saltConc(ii));
    psiText = sprintf('$\\bar{\\psi_6} = %.3f$',psi6_global(ii));
    text(R(2)-1.5,R(3)-1.5,psiText,'FontSize',30,'Color','k','BackgroundColor','w',...
        'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','latex')

    set(gca,'YDir','reverse','Clipping','on')
    axis equal
    axis(R)
    xticks([]); yticks([])
    box on
    axis off
end

%% Plot RDF

figure('Position',[11 61 920 650])

numPks = 4; % number of peaks to find
pksX = zeros(length(targetFrame),numPks); % x-value at peak locations
pksY = zeros(size(pksX)); % y-value of peak locations
legendArray = strings(1,length(targetFrame));
% color = lines(length(targetFrame)); % colors to plot
color = ["#8040E6";"#1AA640";"#E68000";"#A3373C"];

hold on
for ii = 1:length(targetFrame)
    ptsFileName = sprintf('%s_%05d_colorPts.mat',expList(ii),targetFrame(ii));
    load(fullfile(particleFolder,ptsFileName),'pts','DT','d','NN','roiSz','roi','roiPts')

    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expList(ii),targetFrame(ii),cutOff,numBins,roiSz);
    load(fullfile(rdfFolder,rdfFile)); % load the RDF roi data
    gR = gRroi;
    % plot(gR.values,gR.histo,'color',color(ii,:),'linewidth',0.5) % raw data
    smoothHist = movmean(gR.histo,30); % apply moving average
    plot(gR.values,smoothHist,'color',color(ii,:),'LineWidth',3) % plot smooth curve
    pks = findpeaks(smoothHist,'MinPeakWidth',10,'MinPeakHeight',0.5);
    for jj = 1:numPks
        pksIdx = find(smoothHist == pks(jj)); % find index of peaks
        pksX(ii,jj) = gR.values(pksIdx);  % find x-value of peaks
        pksY(ii,jj) = smoothHist(pksIdx); % find y-value of peaks
    end
    ratios = pksX(ii,:)/pksX(ii,1); % get peak ratios
    legendArray(ii) = sprintf('%d mM',saltConc(ii));
end

for ii = 1:length(targetFrame)
    for jj = 1 % for jj = 1:numPks
        plot([pksX(ii,jj),pksX(ii,jj)],[0 pksY(ii,jj)], ...
            '--','color',color(ii,:),'linewidth',1) % plot x-value
    end
    legendArray(ii+length(targetFrame)) = sprintf('%.1f µm',pksX(ii,jj));
end

axis([0 50 0 2])
% yticks([0 1 2 3])
set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on')
xlh = xlabel('r (μm)');
ylh = ylabel('g(r)');
% xlh.Position = xlh.Position + [0 0.1 0]; % adjust x, y, z
% ylh.Position = ylh.Position + [0.06 0 0]; % adjust x, y, z
grid off; box on
legend(legendArray,'FontSize',22,'NumColumns',2,'Box','off')
title('Radial Distribution Function')


%% plot RDFs raw

figure('Position',[11 61 1200 800])

ss = string(NaN(1,numPks+1)); % spacing string
legendArray = strings(1,length(expList)*length(ss));

hold on
for ii = 1:length(expList)
    ptsFileName = sprintf('%s_%05d_colorPts.mat',expList(ii),targetFrame(ii));
    load(fullfile(particleFolder,ptsFileName),'pts','DT','d','NN','roiSz','roi','roiPts')

    rdfFile = sprintf('%s_%05d_RDF_%d_%d_ROI_%d_%d_median.mat',...
        expList(ii),targetFrame(ii),cutOff,numBins,roiSz);
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
    label = sprintf('t_{%d} [%d, %.2f, %.2f, %.2f]',targetFrame(ii),ratios); % create legend label
    legendArray((numPks+2)*ii-numPks-1:(numPks+2)*ii) = [ss(1), label, ss(2:end)];
end

set(gca,'FontSize',28,'XMinorTick','on','YMinorTick','on','LineWidth',1.5)
xlim([0 50])
legend(legendArray)
title('Radial Distribution Function')
xlabel('r (µm)');
ylabel('g(r)');
grid on