%% Methanol Supernatants UV-vis Spectra Data Processing and Grouped Plotting
% Written by Joe Prisaznuk 4/23/2025

%% Load consolidated .mat file

clc; clear variables; close all;

dataFile = 'spectra_data.mat';
load(dataFile)

SD = squeeze(struct2cell(spectraData)); % for group plot configuration
for ii = 1:size(SD,2); fprintf('%02d: %s\n',ii,SD{1,ii}); end

%% Compare absorbance of methanol supernatant to NO CUVETTE case - zoom axes

figure('Position',[11 61 1200 800]); hold on

dataIdx = [24 8 12 16];
legendEntries = ["New (HPLC)","Old","PS1300"];
plotTitle = "Methanol Supernatants Measured 5/15"; % define simplified title

W = spectraData(dataIdx(1)).wavelength;  % same wavelength for all data
I0 = spectraData(dataIdx(1)).intensity;  % no cuvette intensity
win = 5; % smoothing window
[~,xIdx1] = min(abs(W - 390)); % define search range for adjustment factor
[~,xIdx2] = min(abs(W - 410)); % define search range for adjustment factor

for ii = 2:length(dataIdx)
    I = spectraData(dataIdx(ii)).intensity; % raw data
    T = I ./ I0;                      % transmission
    T(T < 0) = NaN;                   % negative is meaningless
    T(T > 1) = NaN;                   % big is meaningless
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    Asmooth = movmean(A, win);    % smooth absorbance
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(W, AsmoothNorm, 'LineWidth', 1.5);
end

xlim([200 600]); ylim([-0.1 1.1])
xlabel('Wavelength (nm)'); % Axes labels
ylabel('Normalized Absorbance');

set(gca,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',24,'LineWidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','GridLineWidth',1)
title(sprintf('%s',plotTitle));
legend(legendEntries);
subtitle('Absorbance relative to no cuvette signal')

insetPos = [0.45 0.3 0.42 0.35]; % [x y width height] in normalized units
insetAx = axes('Position', insetPos);

hold(insetAx, 'on')
for ii = 2:length(dataIdx)
    I = spectraData(dataIdx(ii)).intensity; % raw data
    T = I ./ I0;                      % transmission
    T(T < 0) = NaN;                   % negative is meaningless
    T(T > 1) = NaN;                   % big is meaningless
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    Asmooth = movmean(A, win);    % smooth absorbance
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(insetAx, W, AsmoothNorm, 'LineWidth', 1.5);
end
plot(insetAx, [W(1) W(end)],[0 0],'k--','LineWidth',0.5)
xlim(insetAx, [380 520]); ylim(insetAx, [-0.006 0.02]);
set(insetAx,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',20)
yticks(insetAx, [0 0.01 0.02])
yticklabels(insetAx, ["0","0.01","0.02"])

%% Compare absorbance of methanol supernatant to RESPECTIVE SOLVENT - zoom axes

figure('Position',[11 61 1200 800]); hold on

dataIdx = [8 12 16; 2 5 2];
legendEntries = ["HPLC","Old","PS1300"];
plotTitle = "Methanol Supernatants Measured 5/15"; % define simplified title

I0idx = 24;
W = spectraData(I0idx).wavelength;  % same wavelength for all data
win = 5; % smoothing window
[~,xIdx1] = min(abs(W - 390)); % define search range for adjustment factor
[~,xIdx2] = min(abs(W - 410)); % define search range for adjustment factor

for ii = 1:size(dataIdx,2)
    I = spectraData(dataIdx(1,ii)).intensity; % raw data
    I0 = spectraData(dataIdx(2,ii)).intensity; % respective solvent intensity
    T = I ./ I0;                      % transmission
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    Asmooth = movmean(A, win);    % smooth absorbance
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(W, AsmoothNorm, 'LineWidth', 1.5);
end

xlim([200 600]); ylim([-0.2 1.1])
xlabel('Wavelength (nm)'); % Axes labels
ylabel('Normalized Absorbance');

set(gca,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',24,'LineWidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','GridLineWidth',1)
title(sprintf('%s',plotTitle));
legend(legendEntries);
subtitle('Absorbance relative to respective solvent')

insetPos = [0.45 0.32 0.42 0.35]; % [x y width height] in normalized units
insetAx = axes('Position', insetPos);

hold(insetAx, 'on')
for ii = 1:size(dataIdx,2)
    I = spectraData(dataIdx(1,ii)).intensity; % raw data
    I0 = spectraData(dataIdx(2,ii)).intensity; % respective solvent intensity
    T = I ./ I0;                      % transmission
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    Asmooth = movmean(A, win);    % smooth absorbance
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(insetAx, W, AsmoothNorm, 'LineWidth', 1.5);
end
plot(insetAx, [W(1) W(end)],[0 0],'k--','LineWidth',0.5)
xlim(insetAx, [380 520]);
ylim(insetAx, [-0.006 0.02]);
set(insetAx,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',20)
yticks(insetAx, [0 0.01 0.02])
yticklabels(insetAx, ["0","0.01","0.02"])

%% Compare absorbance of methanol supernatant to HPLC SOLVENT - zoom axes

figure('Position',[11 61 1200 800]); hold on

dataIdx = [16 8; 2 2];
legendEntries = ["PS1300","PS 2 μm"];
plotTitle = "Methanol Supernatants Measured 5/15"; % define simplified title

I0idx = 24;
W = spectraData(I0idx).wavelength;  % same wavelength for all data
win = 5; % smoothing window
[~,xIdx1] = min(abs(W - 600)); % define search range for adjustment factor
[~,xIdx2] = min(abs(W - 650)); % define search range for adjustment factor
colors = [0.3 0.3 0.3; 0.15 0.85 0.15]; % for overlaid lines

for ii = 1:size(dataIdx,2)
    I = spectraData(dataIdx(1,ii)).intensity; % raw data
    I0 = spectraData(dataIdx(2,ii)).intensity; % respective solvent intensity
    T = I ./ I0;                      % transmission
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    Asmooth = movmean(A, win);    % smooth absorbance
    max(abs(Asmooth))
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(W, AsmoothNorm, 'Color', [colors(ii,:)], 'LineWidth', 2);
end

xlim([205 550]); ylim([-0 1.05])
xlabel('Wavelength (nm)'); % Axes labels
ylabel('Normalized Absorbance');

set(gca,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',24,'LineWidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','ClippingStyle','rectangle')
title(sprintf('%s',plotTitle));
legend(legendEntries);
subtitle('Absorbance relative to HPLC solvent')

insetPos = [0.46 0.35 0.4 0.38]; % [x y width height] in normalized units
insetAx = axes('Position', insetPos);

hold(insetAx, 'on')
for ii = 1:size(dataIdx,2)
    I = spectraData(dataIdx(1,ii)).intensity; % raw data
    I0 = spectraData(dataIdx(2,ii)).intensity; % respective solvent intensity
    T = I ./ I0;                      % transmission
    A = log10(1 ./ T);                % absorbance
    aZeroAvg = mean(A(xIdx1:xIdx2));  % mean absorbance in a range
    A = A - aZeroAvg;                 % shift down
    ANorm = A / max(abs(A)); % normalized absorbance
    Asmooth = movmean(A, 50);    % smooth absorbance
    AsmoothNorm = Asmooth / max(abs(Asmooth)); % normalized absorbance
    plot(insetAx, W, ANorm, 'Color', [colors(ii,:) 0.25], 'LineWidth', 0.5);
    plot(insetAx, W, AsmoothNorm, 'Color', [colors(ii,:)], 'LineWidth', 2);
end
xlim(insetAx, [380 520]); ylim(insetAx, [-0.002 0.024]);
set(insetAx,'Box','on','XMinorTick','on','YMinorTick','on',...
    'FontSize',20,'ClippingStyle','rectangle')
yticks(insetAx, [0 0.01 0.02])
yticklabels(insetAx, ["0","0.01","0.02"])