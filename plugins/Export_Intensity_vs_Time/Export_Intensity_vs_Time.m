function Export_Intensity_vs_Time(hObject, eventdata, handles)
% Export_Intensity_vs_Time
% Export intensity within ROI for individual frames in every frame and
% export to .mat file.
%
% The exported .mat file is suitable for import into Spiky.
%

% Get data structure
hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
tUserData = get(hGUI, 'UserData');
ISIdata = [];
if isempty(tUserData)
    warndlg('No data found. You must load data before running this function.', mfilename)
else
    if isfield(tUserData, 'frameStack')
        ISIdata = tUserData;
    end
end

% Get parameters
global ISI;
tParams = ISI.GetParams(hObject, eventdata, handles);

% Get ROIs
ISIdata = ISI_selectSignalROI(ISIdata, tParams.filesQueue(1), tParams.sColMap);
tROI = ISIdata.analysisSignalROI;
mROI = double(tROI(1).mROI);

% Iterate over trials
mSignal = nan(size(ISIdata.frameStack));
for ti = 1:size(ISIdata.frameStack, 1)
    
    % Average of all frames in this trial
    baseframeix = 1:size(ISIdata.frameStack, 2);
    baselineFrames = single(cell2mat(ISIdata.frameStack(ti, baseframeix)));
    baselineFrames = reshape(baselineFrames(:), [ISIdata.frameSizeYX length(baseframeix)]);
    baselineTrial = mean(baselineFrames, 3)'; 

    % Median intensity correction
    baselineTrialMedianCorr = baselineTrial - median(baselineTrial(:));

    % Iterate over all frames in trial
    for fi = 1:size(ISIdata.frameStack, 2)
        if isempty(ISIdata.frameStack{ti, fi})
            warning('Missing frames in ISIdata.frameStack')
        end
        mFrameSignal = single((ISIdata.frameStack{ti, fi}' - baselineTrial) ./ baselineTrial);

        % Mask signal by ROI
        mFrameSignal = mFrameSignal .* mROI;
        
        % Get average signal intensity within ROI
        mSignal(ti, fi) = nanmean(mFrameSignal(:)); % average intensity of all ROI pixels
    end
end

% Plot
nT = (1 / (ISIdata.frame_rate / ISIdata.bin_duration));
hFig = figure('NumberTitle', 'off', 'Name', 'IOS Intensity vs Time');
hAx = subplot(2,2,1);
plot(hAx, mSignal', '-')
set(hAx(1), 'ylim', [min(mSignal(:)) max(mSignal(:))].*1.2, 'xlim', [0 size(mSignal, 2)*nT])
xlabel(hAx(1), 'Trial duration (s)')
ylabel(hAx(1), 'dF/F')
title('All trials')

hAx(2) = subplot(2,2,2);
vMean = mean(mSignal);
vErr = std(mSignal) ./ sqrt(size(mSignal, 2));
hold on
plot(hAx(2), (1:size(mSignal, 2))*nT, vMean, 'k-', 'LineWidth', 2)
plot(hAx(2), (1:size(mSignal, 2))*nT, vMean-vErr, 'k--', 'LineWidth', 1)
plot(hAx(2), (1:size(mSignal, 2))*nT, vMean+vErr, 'k--', 'LineWidth', 1)
axis tight
set(hAx(2), 'xlim', [0 size(mSignal, 2)*nT])
title('Average dF/F (mean +/i std err)')
xlabel(hAx(2), 'Trial duration (s)')
ylabel(hAx(2), 'dF/F')

% Plot all lines concatenated
hAx(3) = subplot(2,1,2);
cla
hold on
for ti = 1:size(mSignal, 1)
    vTime = (1:size(mSignal, 2)) + (size(mSignal, 2) * (ti - 1)); % frames
    vTime = vTime .* nT;
    % Plot a vertical line to denote new trial
    plot(hAx(3), [vTime(1) vTime(1)], [-1 1], 'g-')
    % Plot trial intensity
    plot(hAx(3), vTime, mSignal(ti, :), 'k')
end
set(hAx(3), 'ylim', [min(mSignal(:)) max(mSignal(:))].*1.2, 'xlim', [0 length(mSignal(:))*nT])
xlabel(hAx(3), 'Cumulative trial duration (s)')
ylabel(hAx(3), 'dF/F')

% Export traces to .mat file
% Each trial is saved to a separate variable

% Root path
sPath = tParams.filesQueue(1).path2dir;

% Create a directory to store exported data
sNewDir = [tParams.sAppName ' Export'];
mkdir(sPath, sNewDir);
sPath = [sPath sNewDir];
sFile = 'IOS_dF_All_Trials.mat';
sPath = [sPath filesep sFile];

% Create variables in local workspace
for ti = 1:size(ISIdata.frameStack, 1)
    eval(sprintf('dF_Trial_%d = mSignal(ti, :);', ti));
end
dF_Fs = ISIdata.frame_rate / ISIdata.bin_duration;

% Save trial variables
save(sPath, '-regexp', 'dF_*');

return
