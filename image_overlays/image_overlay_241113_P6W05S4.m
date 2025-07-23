%% Export image with time and scale bar overlaid
% Joe Prisaznuk 4/15/2025

%% Load images and metadata

clc; clear variables; close all;
    
particleFolder = fullfile('..','data','particlePositionRoiColorNN'); % particle data

expCode = '241113_P6W05S4'; % experiment name

I = imageDatastore(fullfile('..','data','imageData8bit',expCode),...
    'IncludeSubfolders',true); % Create imageDatastore

load(fullfile('..','data','imageMetadata','MDall.mat')); % load all metadata

mdIdx = find(strcmp({MD.expCode}, expCode)); % find MDall index
tElap = MD(mdIdx).T.time - MD(mdIdx).T.time(1);

times = 1:2; % which times to plot

ionVid = 3; % ionizer source video
frames = [800 1150] + sum(MD(mdIdx).numFrames(1:ionVid-1));

for ii = 1:length(frames)
    if ismember(ii,times)
        fprintf('<strong>%02d: %05d t = %s hh:mm:ss (vid %d frame %d)</strong>\n',...
            ii,frames(ii),tElap(frames(ii)),MD(mdIdx).T.vidIdx(frames(ii)),...
            MD(mdIdx).T.frame(frames(ii)));
    else
        fprintf('%02d: %05d t = %s hh:mm:ss (vid %d frame %d)\n',ii,frames(ii),...
            tElap(frames(ii)),MD(mdIdx).T.vidIdx(frames(ii)),MD(mdIdx).T.frame(frames(ii)));
    end
end

%% Show the images with text overlaid

imScale = 2; % image resize factor
scale = 1.848; % image scale factor [µm/px]
sbUM = 50; % scale bar length [µm]
sbHt = 10; % scale bar height [µm]
sbShift = 10; % scale bar offset from corner [µm]
frameTimes = MD(mdIdx).T.time(frames(times)) - MD(mdIdx).T.time(1);
frameTimes.Format = 'hh:mm';
numPts = zeros(length(times),1);

for ii = 1:length(times)

    ptsFileName = sprintf('%s_%05d_colorPts.mat',expCode,frames(times(ii)));
    load(fullfile(particleFolder,ptsFileName),'pts','DT','d','NN','roiSz','roi','roiPts')

    numPts(ii) = length(roiPts);
    W = MD(mdIdx).T.width(frames(times(ii))); H = MD(mdIdx).T.height(frames(times(ii)));
    figure('Position',[11 61 imScale*roiSz(1)/scale imScale*roiSz(2)/scale],'MenuBar','none')
    axes('Position',[0 0 1 1])

    roiCrop = 20; % how much to shrink ROI [µm]
    R = imScale*roi/scale + roiCrop*[1 -1 1 -1]; % shrink the ROI
    cropIdx = round(R);
    im = readimage(I,ii); % read image
    imUp = imresize(im,imScale);         % upscale image
    imCrop = imUp(cropIdx(3):cropIdx(4),cropIdx(1):cropIdx(2));
    imCropAdj = imadjust(imCrop);
    imshow(imCropAdj);
    
    if ii == 1
        timeText = sprintf('t_1 = %.1f min',minutes(frameTimes(ii)));
    else
        timeText = sprintf('t_2 - t_1 = %.1f sec',seconds(frameTimes(2)-frameTimes(1)));
    end
    text(cropIdx(2)-cropIdx(1)-1,1,timeText,'FontSize',32,'Color','w','BackgroundColor','k',...
        'HorizontalAlignment','right','VerticalAlignment','top')

    rectangle('Position',[cropIdx(2)-cropIdx(1)-(sbUM+sbShift)*imScale/scale,...
        cropIdx(4)-cropIdx(3)-(sbHt+sbShift)*imScale/scale, sbUM*imScale/scale,...
        sbHt*imScale/scale], 'FaceColor','w','EdgeColor','none')
end
    
aggregatateSize = numPts(1)/numPts(2);
noAgg = 600; % assume 500 particles didn't aggregate
noAgg = 2*numPts(2) - numPts(1);
aggregatateSizeEstimate = (numPts(1)-noAgg)/(numPts(2)-noAgg);
fprintf('\nRatio of np_i/np_f = %d/%d = %.2f\n',...
    numPts(1),numPts(2),aggregatateSize)
fprintf('Ratio of (np_i-nA)/(np_f-nA) = %d/%d = %.2f\n',...
    numPts(1)-noAgg,numPts(2)-noAgg,aggregatateSizeEstimate)
