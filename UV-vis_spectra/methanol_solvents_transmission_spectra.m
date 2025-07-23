%% Methanol Solvents UV-vis Spectra Data Processing and Grouped Plotting
% Written by Joe Prisaznuk 5/15/2025

%% Load consolidated .mat file

clc; clear variables; close all;

dataFile = 'spectra_data.mat';
load(dataFile)

SD = squeeze(struct2cell(spectraData)); % for group plot configuration
for ii = 1:size(SD,2); fprintf('%02d: %s\n',ii,SD{1,ii}); end

%% Compare methanol raw signals including 4/23/25 AND 5/7/25

figure('Position',[11 61 1200 650]); hold on

dataIdx = [24 2 5 27 26 25];
legendEntries = ["no cuvette","HPLC 5/15","Old 5/15",...
    "HPLC 5/7","Old 5/7","Old 4/23"];
plotTitle = "Methanol Solvents over Time"; % define simplified title

W = spectraData(dataIdx(1)).wavelength; % same wavelength for all data

for ii = 1:length(dataIdx)
    if dataIdx(ii) == 25
        plot(W, 0.8*spectraData(dataIdx(ii)).intensity, 'LineWidth', 1.5);
    else
        plot(W, spectraData(dataIdx(ii)).intensity, 'LineWidth', 1.5);
    end
end

xlim([200 W(end)]); ylim([0 2^16])
xlabel('Wavelength (nm)'); % Axes labels
ylabel('Intensity (a.u.)');

set(gca,'Box','on','XMinorTick','on','YMinorTick','on')
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',24,'LineWidth',1.5)
grid on; set(gca,'GridLineWidth',1)
title(sprintf('%s',plotTitle));
legend(legendEntries);
grid on

%% Compare transmission of methanol grades to no cuvette case including 4/23/25 AND 5/7/25

figure('Position',[11 61 1200 650]); hold on

dataIdx = [24 2 5 27 26 25];
legendEntries = ["HPLC spec","HPLC 5/15 (30 min spin)",...
    "Old 5/15 (30 min spin)","HPLC 5/7","Old 5/7","Old 4/23"];
ls = ["-","-","--","-","--","--"]; % line styles
C = lines(7);
colArray = [0 0 0; C(1,:); C(1,:); C(2,:); C(2,:); C(3,:)];
plotTitle = "Methanol Solvents over Time"; % define simplified title

W = spectraData(dataIdx(1)).wavelength;  % same wavelength for all data
I0 = spectraData(dataIdx(1)).intensity;  % no cuvette intensity
win = 1; % smoothing window
[~,xIdx1] = min(abs(W - 380)); % define search range for adjustment factor
[~,xIdx2] = min(abs(W - 420)); % define search range for adjustment factor

ref_W = [210, 220, 225, 230, 235, 240, 250, 260, 280];
ref_I = [45, 65, 70, 85, 90, 95, 95, 98, 98];
plot(ref_W,ref_I,'k-o','MarkerSize',3,'LineWidth',1.5);

for ii = 2:length(dataIdx)
    I = spectraData(dataIdx(ii)).intensity; % raw data
    T = 100 * I ./ I0;                      % transmission
    tMaxAvg = mean(T(xIdx1:xIdx2));         % mean transmission 400 - 800 nm
    Tadj = T * 100 / tMaxAvg;               % adjusted transmission
    Tsmooth = movmean(Tadj, win);           % smooth transmission
    plot(W, Tsmooth, ls(ii), 'Color', colArray(ii,:), 'LineWidth', 1.5);
end

xlim([200 420]); ylim([0 120])
xlabel('Wavelength (nm)'); % Axes labels
ylabel('Transmission (%)');

set(gca,'Box','on','XMinorTick','on','YMinorTick','on')
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',24,'LineWidth',1.5)
grid on; set(gca,'GridLineWidth',1)
title(sprintf('%s',plotTitle));
legend(legendEntries,'Location','southeast');
grid on