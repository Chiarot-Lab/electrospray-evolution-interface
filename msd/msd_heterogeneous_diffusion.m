%% Tracking results and MSD of particles in PLA v6 DI water run 3
% Joe Prisaznuk 12/16/2024

%% Load images and metadata

clc; clear variables; close all;

expCode = '241205_P6W00X3'; % experiment name

load(fullfile('..','data','imageMetadata','MDall.mat')); % load all metadata

mdIdx = find(strcmp({MD.expCode}, expCode)); % find MDall index

%% Read xml track file

startFrame = 10801; % frame to start at
endFrame = 18000; % frame to end at
calcFrames = startFrame:endFrame; % frames to calculate position at

imScale = 2; % image resize factor
scale = 1.848; % image scale factor [µm/px]
roi = [25 25 125 125]; % image crop dimensions in pixels [x y w h]
roiLim = [roi(1) roi(1) + roi(3) roi(2) roi(2) + roi(4)] + 1;
roiW = roi(3); roiH = roi(4); % ROI size
pad = 5; % show region outside of ROI
W = MD(mdIdx).T.width(startFrame); H = MD(mdIdx).T.height(startFrame);

trackFile = sprintf('%s_%05d_%05d_scale_%dX_crop_%dx_%dy_%dw_%dh_Tracks_%d.mat',...
    expCode,startFrame,endFrame,imScale,roi,length(calcFrames)); % file name

load(trackFile,'X','Y','fullTrackIdx','metadata')

nTracks = metadata.nTracksAttribute; % number of tracks

%% Remove tracks with duplicate X or Y

validIdx = false(size(X,1), 1);

for ii = 1:size(X, 1)
    particle_data = [X(ii,:)' Y(ii,:)']; % Get x, y for a particle across all time
    unique_data = unique(particle_data, 'rows'); % Check if (x, y) pairs are unique
    if size(unique_data, 1) == size(X, 2) % # unique (x, y) pairs = # time points?
        validIdx(ii) = true;
    end
end

X = X(validIdx,:);          Y = Y(validIdx,:); % remove tracks with duplicate values
Xpx = X/scale + roi(1) + 1; Ypx = Y/scale + roi(2) + 1;
nP = size(X, 1); % number of particles after reduction
nT = size(X, 2); % number of time points

%% Calculate MSD and segment particles

[Dt_r, MSD] = fastMSDfft(X, Y);

if strcmp(metadata.timeUnitsAttribute,"ms")
    frameIntSec = metadata.frameIntervalAttribute/1000; % convert to seconds
else
    frameIntSec = metadata.frameIntervalAttribute; % assume seconds by default
end
delay = Dt_r * frameIntSec; 

MSDsegmentIdx = 1:300;  % use early delays for segmentation
MSDcutoff = 0.5;        % choose based on trial and error

avgMSD = mean(MSD(:,MSDsegmentIdx),2);
avgMSD = [avgMSD,(1:length(avgMSD))'];
avgMSD = sortrows(avgMSD,1);
ogSort = avgMSD(:,2); % "original" sorting, to differentiate from below
MSD = MSD(avgMSD(:,2),:); % sort by average MSD
MSDcutIdx = find(avgMSD(:,1) < MSDcutoff,1,'last');

colorMSD = [abyss(MSDcutIdx); spring(nP-MSDcutIdx)];

figure('Position',[11 61 1200 800])

hold on
plot(NaN,NaN,'LineWidth',3,'Color',[colorMSD(MSDcutIdx/2,:),0.5])
plot(NaN,NaN,'LineWidth',3,'Color',[colorMSD(MSDcutIdx+(nP-MSDcutIdx)/2,:),0.5])
for ii = 1:nP
    plot(delay,MSD(ii,:),'color',[colorMSD(ii,:) 0.6],'LineWidth',0.5)
end
plot([delay(2) delay(end)],[MSDcutoff MSDcutoff],'g','LineWidth',3)

axis([delay(1) delay(end) 0.001 250])
set(gca,'XScale','log','YScale','log')
set(gca,'MinorGridLineWidth',0.5,'GridLineWidth',1,'FontSize',24,'LineWidth',1.0)

grid on; box on;
legend(sprintf('Lower group = %d particles',MSDcutIdx),...
    sprintf('Upper group = %d particles',nP-MSDcutIdx),'Location','northwest')
xlabel('Delay Time \tau (s)')
ylabel('MSD ⟨x^2(\tau)⟩ (µm^2)')
title(sprintf('MSD segmentation using range 0 - %.1f s',delay(MSDsegmentIdx(end))))

%% Plot tracks (color by MSD)

figure('Position',[11 61 1000 1000])

DT = delaunayTriangulation(X(:,1),Y(:,1));
thresholdAngDeg = 120; % acute triangle threshold
[DT, ~] = filterAcuteTriangles(DT, thresholdAngDeg);
attachedTriangles = vertexAttachments(DT);

NN = cell(size(X(:,1)));
for ii = 1:nP
    attachedVertices = DT.ConnectivityList(attachedTriangles{ii},:);
    NN{ii} = setdiff(unique(attachedVertices), ii);
end

numNN = cellfun(@length, NN); % array giving # of neighbors for coloring
% set particle colors
color = [0,0,0; lines(7); 1,0,1; 1,0,0; 0,1,0; 0,0,1; 1,1,0]; % particle color list
ptCol = zeros(nP,3); % initialize color array
for k = 1:nP; ptCol(k,:) = color(numNN(k)+1,:); end % assign colors

hold on
for ii = 1:nP
    plot(X(avgMSD(ii,2),:),Y(avgMSD(ii,2),:),'color',[colorMSD(ii,:) 0.6],'LineWidth',0.5)
end
triplot(DT,'k','LineWidth',1,'color',[0 0 0 0.2])
scatter(X(:,1),Y(:,1),100,ptCol)

set(gca,'YDir','reverse','FontSize',24,'LineWidth',1.0)
xlabel('x (µm)'); ylabel('y (µm)')
axis equal; box on
axis([0 roiW*scale 0 roiH*scale] + pad*[-1 1 -1 1])
% axis([-4 66 164 234])

title(sprintf('%s tracks (%d particles)',MD(mdIdx).expName,size(X,1)))
subtitle(sprintf('Frames %d - %d (original)',startFrame,endFrame))

%% Get the mean velocity in each cell
% note: drift correction code was partially generating using ChatGPT o3;
% the output was fully reviewed by the author before implementation.

nCells = 3;  % number of cells along each dimension
ptsIdxDrift = ogSort; % consider all particles

% Compute displacements and midpoints (for each consecutive time step)
dX = diff(X(ptsIdxDrift,:),1,2); dY = diff(Y(ptsIdxDrift,:),1,2);
X_mid = (X(ptsIdxDrift,2:end) + X(ptsIdxDrift,1:end-1)) / 2;
Y_mid = (Y(ptsIdxDrift,2:end) + Y(ptsIdxDrift,1:end-1)) / 2;
% nTime = size(dX, 2) - 100;  % number of displacement time steps
nTime = size(dX, 2);  % number of displacement time steps

x_edges = linspace(min(X(:)), max(X(:)), nCells+1);
y_edges = linspace(min(Y(:)), max(Y(:)), nCells+1);

% Preallocate 3D arrays for drift (cell-averaged velocity) for each time step.
driftX_grid = nan(nCells, nCells, nTime);
driftY_grid = nan(nCells, nCells, nTime);

for t = 1:nTime
    % Get the midpoints and displacements for this time step (as column vectors)
    X_mid_t = X_mid(:, t);
    Y_mid_t = Y_mid(:, t);
    dX_t    = dX(:, t);
    dY_t    = dY(:, t);
    
    % Loop over cells
    c = 1;
    for i = 1:nCells
        for j = 1:nCells
            % Determine the cell boundaries for cell (i,j)
            x_min = x_edges(i); x_max = x_edges(i+1);
            y_min = y_edges(j); y_max = y_edges(j+1);
            
            % Find indices of points in this cell
            cellIdx = X_mid_t >= x_min & X_mid_t < x_max & ...
                      Y_mid_t >= y_min & Y_mid_t < y_max;
            
            if any(cellIdx)
                % Compute the average displacement in the cell
                driftX_grid(i,j,t) = mean(dX_t(cellIdx));
                driftY_grid(i,j,t) = mean(dY_t(cellIdx));
            else
                % If no points fall in the cell, set drift to zero (or NaN if preferred)
                driftX_grid(i,j,t) = 0;
                driftY_grid(i,j,t) = 0;
            end
            nPinCell(c) = sum(cellIdx);
            if t == 1
                fprintf('%d,%d: %d particles\n',i,j,nPinCell(c));
            end
            c = c + 1;
        end
    end
end

fprintf('\nMean # particles per cell: %.2f\n',mean(nPinCell))

%% Integrate the drift velocity to get a drift track per cell

% Use cumtrapz along the 3rd dimension (time) to integrate velocity in each cell.
% The result has the same dimensions as driftX_grid.
drift_trackX = cumtrapz(1, driftX_grid, 3);
drift_trackY = cumtrapz(1, driftY_grid, 3);

% Compute the centers of the grid cells for plotting
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
[X_center, Y_center] = ndgrid(x_centers, y_centers);

figure('Position',[11 61 1000 1000])

vXoa = X(:,end) - X(:,1); vYoa = Y(:,end) - Y(:,1);
quiver(X(:,1),Y(:,1),vXoa,vYoa,'k','AutoScale','off')
hold on;
for k = 1:length(x_edges)
    plot([x_edges(k), x_edges(k)], [min(y_edges), max(y_edges)], 'k--', 'LineWidth', 0.5);
    plot([min(x_edges), max(x_edges)], [y_edges(k), y_edges(k)], 'k--', 'LineWidth', 0.5);
end

for i = 1:nCells
    for j = 1:nCells
        dTx = squeeze(drift_trackX(i,j,:));
        dTy = squeeze(drift_trackY(i,j,:));
        plot(5*dTx+X_center(i,j)+20,5*dTy+Y_center(i,j)+15,'LineWidth',1)
    end
end
title(sprintf('Drift tracks with %d × %d grid',nCells,nCells))
legend('v_{overall}')
xlabel('x (µm)'); ylabel('y (µm)')
axis equal; box on
axis([0 roiW*scale 0 roiH*scale] + pad*[-1 1 -1 1])
set(gca,'YDir','reverse','FontSize',24,'LineWidth',1.0)

%% Assign each particle to a cell based on its initial position

% For each particle, use its initial (time=1) position to determine its cell.
cellIdx_x = discretize(X(:,1), x_edges);  % indices along x (nParticles x 1)
cellIdx_y = discretize(Y(:,1), y_edges);  % indices along y

%% Loop over particles and subtract the drift corresponding to the cell

X_corr = X; Y_corr = Y; % corrected trajectory, this ensures X_corr(:,1) = X(:,1)
for p = 1:nP
    % Determine the cell for particle p (based on its initial position)
    iCell = cellIdx_x(p);
    jCell = cellIdx_y(p);
    
    % If a particle falls outside the grid (NaN), default to no drift correction.
    if isnan(iCell) || isnan(jCell)
        driftX_this = zeros(1, nT);
        driftY_this = zeros(1, nT);
    else
        % Get the drift track for the cell as a row vector.
        % Note: squeeze removes singleton dimensions.
        driftX_this = squeeze(drift_trackX(iCell, jCell, :))';  % 1 x nTime
        driftY_this = squeeze(drift_trackY(iCell, jCell, :))';
    end
    
    % Subtract the cell drift track from the original trajectory.
    % (Assuming X and Y are defined for the same time steps as drift_track arrays.)
    X_corr(p,2:end) = X(p,2:end) - driftX_this;
    Y_corr(p,2:end) = Y(p,2:end) - driftY_this;
end

%% Calculate MSD_DC_velNew - drift correction using velocity cells

[~, MSD_DC_velNew] = fastMSDfft(X_corr, Y_corr);

MSDsegmentIdx = 1:300;  % use early delays for segmentation
MSDcutoff = 0.5;        % choose based on trial and error

avgMSD = mean(MSD_DC_velNew(:,MSDsegmentIdx),2);
avgMSD = [avgMSD,(1:length(avgMSD))'];
avgMSD = sortrows(avgMSD,1);
MSD_DC_velNew = MSD_DC_velNew(avgMSD(:,2),:); % sort by "new" average MSD
% MSD_DC_NN = MSD_DC_NN(ogSort,:); % sort by average MSD
MSDcutIdx = find(avgMSD(:,1) < MSDcutoff,1,'last');

lowMSD = [mean(MSD_DC_velNew(1:MSDcutIdx,:))' std(MSD_DC_velNew(1:MSDcutIdx,:),1,1)'];
highMSD = [mean(MSD_DC_velNew(MSDcutIdx+1:end,:))' std(MSD_DC_velNew(MSDcutIdx+1:end,:),1,1)'];

qLow(1,:) = quantile(MSD_DC_velNew(1:MSDcutIdx,:), 0.1);
qLow(2,:) = quantile(MSD_DC_velNew(1:MSDcutIdx,:), 0.9);
qHigh(1,:) = quantile(MSD_DC_velNew(MSDcutIdx+1:end,:), 0.1);
qHigh(2,:) = quantile(MSD_DC_velNew(MSDcutIdx+1:end,:), 0.9);
lowMSDiqr = [mean(MSD_DC_velNew(1:MSDcutIdx,:))' qLow(1,:)' qLow(2,:)'];
highMSDiqr = [mean(MSD_DC_velNew(MSDcutIdx+1:end,:))' qHigh(1,:)' qHigh(2,:)'];

colorMSD = [abyss(MSDcutIdx); spring(nP-MSDcutIdx)];

figure('Position',[11 61 1200 800])

hold on
plot(NaN,NaN,'LineWidth',3,'Color',[colorMSD(MSDcutIdx+(nP-MSDcutIdx)/2,:),0.5])
plot(NaN,NaN,'LineWidth',3,'Color',[colorMSD(MSDcutIdx/2,:),0.5])
for ii = 1:nP
    plot(delay,MSD_DC_velNew(ii,:),'color',[colorMSD(ii,:) 0.6],'LineWidth',0.5)
end
% plot([delay(2) delay(end)],[MSDcutoff MSDcutoff],'g','LineWidth',3)

axis([delay(1) delay(end) 0.001 250])
set(gca,'XScale','log','YScale','log')
set(gca,'MinorGridLineWidth',0.5,'GridLineWidth',1,'FontSize',24,'LineWidth',1.0)

grid on; box on;
legend(sprintf('Upper group = %d particles',nP-MSDcutIdx),...
    sprintf('Lower group = %d particles',MSDcutIdx),'Location','northwest')
xlabel('Delay Time \tau (s)')
ylabel('MSD ⟨x^2(\tau)⟩ (µm^2)')
title(sprintf('MSD with %d × %d cell drift correction',nCells,nCells))

%% Plot drift-corrected tracks (color by MSD)

figure('Position',[11 61 1000 1000])

DT = delaunayTriangulation(X(:,1),Y(:,1));
thresholdAngDeg = 120; % acute triangle threshold
[DT, ~] = filterAcuteTriangles(DT, thresholdAngDeg);
attachedTriangles = vertexAttachments(DT);

NN = cell(size(X(:,1)));
for ii = 1:nP
    attachedVertices = DT.ConnectivityList(attachedTriangles{ii},:);
    NN{ii} = setdiff(unique(attachedVertices), ii);
end

numNN = cellfun(@length, NN); % array giving # of neighbors for coloring
% set particle colors
color = [0,0,0; lines(7); 1,0,1; 1,0,0; 0,1,0; 0,0,1; 1,1,0]; % particle color list
ptCol = zeros(nP,3); % initialize color array
for k = 1:nP; ptCol(k,:) = color(numNN(k)+1,:); end % assign colors

hold on
for ii = 1:nP
    plot(X_corr(avgMSD(ii,2),:),Y_corr(avgMSD(ii,2),:),'color',[colorMSD(ii,:) 0.6],'LineWidth',0.5)
end
triplot(DT,'k','LineWidth',1,'color',[0 0 0 0.2])
scatter(X_corr(:,1),Y_corr(:,1),100,ptCol)
sc1 = scatter(NaN,NaN,120,'go','MarkerFaceColor','w','Linewidth',4);
sc2 = scatter(NaN,NaN,120,'ro','MarkerFaceColor','w','Linewidth',4);

set(gca,'YDir','reverse','FontSize',24,'LineWidth',1.0)
xlabel('x (µm)'); ylabel('y (µm)')
axis equal; box on
axis([0 roiW*scale 0 roiH*scale] + pad*[-1 1 -1 1])

title(sprintf('%s tracks (%d particles)',MD(mdIdx).expName,size(X,1)))
subtitle(sprintf('Frames %d - %d (drift corrected)',startFrame,endFrame))

%% Calculate theoretical D

kb = 1.380649e-23;    % boltzmann [J/K]
T = 273 + 20;         % temperature [K]
eta = [1e-3 1.81e-5]; % dynamic viscosity water, air [J*s/m3]
weights = [1 10];     % for weighted average
etaAvg = [sum(eta)/2 sum(eta.*weights)/sum(weights)]; % average dynamic viscosity [J*s/m3]
r = 1e-6;             % particle radius [m]

D = kb*T./(6*pi*eta*r); % diffusivity [m^2/s]
D = D*(1e6)^2;          % diffusivity [µm^2/s]

Davg = kb*T./(6*pi*etaAvg*r); % diffusivity [m^2/s]
Davg = Davg*(1e6)^2;         % diffusivity [µm^2/s]

fprintf(['\nEinstein diffusivity D = k*T/(6*π*η*r)' ...
    '\nD_w = %.3g µm²/s,   D_a = %.3g µm²/s\n' ...
    'D_avg = %.3g µm²/s, D_avgWeighted = %.3g µm²/s\n'],D(1),D(2),Davg(1),Davg(2))

%% Plot mean MSD and fit different parts of the curves

figure('Position',[11 61 1200 800])

% lineColor = lines(2);
% lineColor = [mean(colorMSD(1:MSDcutIdx,:)); mean(colorMSD(MSDcutIdx:end,:))];
lineColor = [0.1 0.25 0.8; 0.8 0.25 0.1];

patchIdx = 7199;
lowMSDband = [lowMSDiqr(2:patchIdx,2)', fliplr(lowMSDiqr(2:patchIdx,3)')];
% lowMSDband(lowMSDband <= 0) = 1e-6;
highMSDband = [highMSDiqr(2:patchIdx,2)', fliplr(highMSDiqr(2:patchIdx,3)')];
% highMSDband(highMSDband <= 0) = 1e-6;

hold on
patch([delay(2:patchIdx) fliplr(delay(2:patchIdx))],lowMSDband,...
    lineColor(1,:),'FaceAlpha',0.25,'EdgeColor',lineColor(1,:),...
    'LineWidth',0.5,'EdgeAlpha',0.5)
patch([delay(2:patchIdx) fliplr(delay(2:patchIdx))],highMSDband,...
    lineColor(2,:),'FaceAlpha',0.25,'EdgeColor',lineColor(2,:),...
    'LineWidth',0.5,'EdgeAlpha',0.5)

grid off; box on;
% axis([delay(1) delay(100) 0 3])
axis([delay(1) delay(end) 1e-3 250])
set(gca,'XScale','log','YScale','log')
set(gca,'MinorGridLineWidth',0.5,'GridLineWidth',1)
xlabel('\tau (s)')
ylabel('MSD (µm^2)')
% title('Mixed Mobility on DI Water/Air Interface'); % no title for paper

% fitting manually
lowFit1Idx = 1:20;
lowFit1 = fittype('a*x^b'); % linear 'fittype' variable of the form f(a,x) = a*x
low1 = fit(delay(lowFit1Idx)',lowMSDiqr(lowFit1Idx)',lowFit1,'StartPoint',[0 0]);
low1Curve = low1.a * delay.^low1.b;
low1D = low1.a / 4; low1C = confint(low1); low1DC = low1C(:,1) / 4;

highFit1Idx = 1:5;
highFit1 = fittype('a*x^b'); % linear 'fittype' variable of the form f(a,x) = a*x
high1 = fit(delay(highFit1Idx)',highMSDiqr(highFit1Idx)', highFit1,'StartPoint',[0 0]);
high1Curve = high1.a * delay.^high1.b;
high1D = high1.a / 4; high1C = confint(high1); high1DC = high1C(:,1) / 4;

plot(delay,lowMSDiqr(:,1),'o','Color',lineColor(1,:),'MarkerSize',7,'LineWidth',1)
plot(delay,highMSDiqr(:,1),'o','Color',lineColor(2,:),'MarkerSize',7,'LineWidth',1)

plot(delay,low1Curve, 'b-','LineWidth',1)
plot(delay,high1Curve,'r-','LineWidth',1)

plot(delay(8:40),D(1)*4*delay(8:40),'k-.','LineWidth',2)
% plot(delay(10:30),Davg(1)*4*delay(10:30),'c-','LineWidth',2)

fprintf('\nLinear fit of the average MSD_low1 curve with 95%% confidence interval:\n')
fprintf('D = %.3g [ %.3g - %.3g ] %s \n', low1D,  low1DC(1),  low1DC(2),  'µm²/s');
fprintf('α = %.3g [ %.3g - %.3g ]\n', low1.b,  low1C(1,2),  low1C(2,2));

fprintf('\nLinear fit of the average MSD_high1 curve with 95%% confidence interval:\n')
fprintf('D = %.3g [ %.3g - %.3g ] %s \n', high1D, high1DC(1), high1DC(2), 'µm²/s');
fprintf('α = %.3g [ %.3g - %.3g ]\n', high1.b,  high1C(1,2),  high1C(2,2));

plot(delay(8:40),0.1*delay(8:40),'k-','LineWidth',2)
text(delay(18),0.1*delay(18),'1',...
    'HorizontalAlignment','left','VerticalAlignment','top')
set(gca,'FontSize',28,'LineWidth',1.0)
xticks([0.1 1 10 100])
yticks([0.01 0.1 1 10 100])
xlim([0 430])

lh = legend('','','⟨x^2(\tau)_{low}⟩','⟨x^2(\tau)_{high}⟩', ...
    sprintf('%.3f\\tau^{%.3f}',low1.a,low1.b),...
    sprintf('%.3f\\tau^{%.3f}',high1.a,high1.b),...
    'D_w','Location','northwest','FontSize',24,'NumColumns',1,'box','off');
lh.Position = lh.Position + [0.01 0 0 0];
