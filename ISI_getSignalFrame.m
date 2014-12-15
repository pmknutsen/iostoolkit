function signalFrame = ISI_getSignalFrame(ISIdata, prmts, figureflag)
%
%
% Usage:
%   signalFrame=ISI_getSignalFrame(ISIdata,stiminterval,figureflag)
%
%  given ISIdata, isolate signal frame within given stiminterval.
%  plot the stimulus frame and the signal frame if flag is true
%

stimi = prmts.stimInterval;

defaultstimi=[0.5 1.5];
if nargin==1
    figureflag=1; %default to plot figure
    stimi=defaultstimi;
elseif nargin==2
    stimi=defaultstimi;
elseif nargin==3
    %do nothing
else
    error('incorrect number of input arguments');
end

% make sure stimi is rounded down to valid frame timepoints
% e.g.  if frame bins are 0.2 sec, then [0.5 1.5] should be
%       rounded to [0.4 1.4]
binpersec = ISIdata.frame_rate / ISIdata.bin_duration;
stimi = floor(stimi ./ (1/binpersec)) .* (1/binpersec); % added, Per 012312

% Average deltaR across baseline, stim, and "poststim",
stimFrameix = (ISIdata.nPreStimFrames+1):(ISIdata.nPreStimFrames+ISIdata.nStimFrames);
stimFrame = mean(ISIdata.deltaSignal(:,:,stimFrameix), 3);
maxval1 = max(max(stimFrame));
minval1 = min(min(stimFrame));

% Calculate the frames for x1 - x2 sec after stim
signalFrameix = (ISIdata.nPreStimFrames+binpersec*stimi(1)):(ISIdata.nPreStimFrames+binpersec*stimi(2));
signalFrame = mean(ISIdata.deltaSignal(:,:,signalFrameix), 3);
maxval2=max(max(signalFrame)); minval2=min(min(signalFrame));

if figureflag
    %%
    clim = prmts.climAll ./ 1000;

    % Smooth evoked maps
    if ~isnan(prmts.smoothSigma) && prmts.smoothSigma > 0
        F = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
        stimFrameSmooth = filter2(F, stimFrame, 'valid');
        signalFrameSmooth = filter2(F, signalFrame, 'valid');
    else
        stimFrameSmooth = stimFrame;
        signalFrameSmooth = signalFrame;
    end
    
    hFig = figure('name',['Evoked map - ' prmts.name], 'visible', 'off');
    
    set(hFig, 'position', [0 0 800 350], 'numbertitle', 'off')
    centerfig(hFig)
    set(hFig, 'visible', 'on')
    colormap(gray(256))

    hAx = subplot(1, 2, 1);
    hIm = imagesc(stimFrameSmooth, 'parent', hAx);
    title('Stimulation frames');

    hAxR = subplot(1, 2, 2);
    hImR = imagesc(signalFrameSmooth, 'parent', hAxR);

    axis([hAx hAxR], 'image', 'off');

    set([hAx hAxR], 'clim', clim);
    

    title(sprintf('Post stim window (%.1f - %0.1f s)', stimi(1), stimi(2) ));
    %%
end

return