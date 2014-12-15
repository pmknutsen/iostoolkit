function Streamer_QuickLook(hObject, eventdata, handles)
% Streamer_QuickLook
%
%
% Per M Knutsen <pmknutsen@gmail.com>
% June 2014
%
%   TODO
%       Add switch to re-use data in MAT file
%           OR
%       Option to continue analysing file after loading incomplete MAT
%           Add indicator in T saying which frame was reached?
%       Recycle baseline frame from loaded MAT file
%       Add white area around axes (give some space)
%
%

T(1).Date = date();

% This is the default _conf file that will be loaded in the absence of
% a video specific *.mat or *_conf.m
sDefaultConfFile = 'DefaultStreamerConfFile';

% Get iostoolkit GUI handles
handles = guidata(hObject);

% Get filename and path from GUI
T.sPath = get(handles.path, 'string');
T.sFileName = get(handles.data_filename, 'string');
T.sFileName = T.sFileName(1:end-4);

% Load settings from disk
sPwd = pwd;
sConfFile = [T.sPath T.sFileName '_conf.m'];
if exist(sConfFile, 'file')
    cd(T.sPath)
    sConfFile = [T.sFileName '_conf'];
else
    sConfFile = sDefaultConfFile;
end
if exist(sConfFile, 'file')
    eval(sConfFile);
else
    error('Could not locate any configuration file for:\n%s', [T.sPath T.sFileName '.bin'])
end
T.sConfFile = sConfFile;
cd(sPwd)

% Load existing T from disk
if T.bLoad
    if exist([T.sPath T.sFileName '.mat'], 'file')
        R = load([T.sPath T.sFileName '.mat'], 'T');
        disp(sprintf('Results loaded from %s', [T.sPath T.sFileName '.mat']))
        % Copy R.T into T without replacement
        csFields = fieldnames(R.T);
        
        for f = 1:length(csFields)
            if ~isfield(T, csFields{f})
                T.(csFields{f}) = R.T.(csFields{f});
            end
            
        end
        if isfield(T, 'AA_proc')    % CM 20140723 allows to not load AA_proc...
            ...                     (done to allow proper initialisation of the graph)
                T= rmfield(T,'AA_proc');
        end
    end
end

% Set first value in colormap to black
T.mColMap = padarray(T.mColMap, 1, 'pre');

% Get file info (resolution etc) from text file associated with .bin
% Settings file string

% Open handle to associated .txt file
hFID = fopen([T.sPath T.sFileName '.txt']);
if (hFID == -1), error('Failed opening settings file.'); end
[cSettings, ~] = textscan(hFID, '%s');
fclose(hFID)

% Get acquisition settings
mSettings = cSettings{1};
T.vVideoRes = str2num(mSettings{2, :}); % Video resolution (height * width * bitdepth)
vFrames = str2num(mSettings{3, :});
T.nVideoFPS = vFrames(1); % frame/sT
T.nVideoNumFrames = vFrames(2); % total# of frames in file

% Create smoothing filter
T.nSigma = str2num(get(handles.smooth_sigma, 'string'));
if ~isempty(T.nSigma)
    T.mWin = fspecial('gaussian', T.nSigma*3, T.nSigma);
end

% Open file handle
FID = fopen([T.sPath T.sFileName '.bin']);
if (hFID == -1), error('Failed opening video file.'); end

% Check that not too many baseline frames were set
if (T.nBaselineFrames * T.nNumColChans) >= T.nVideoNumFrames
    error('Set a lower number of baseline frames');
end

% Calculate baselines
hWait = waitbar(0, 'Calculating baselines...');
mBaselineSums = zeros([T.vVideoRes([2 1]) T.nNumColChans]);

% Get background frames of 1st channel
vBgFrames = (T.nSkipInitialFrames+1):T.nNumColChans:(T.nVideoNumFrames-T.nNumColChans);
vBgFrames = vBgFrames(unique(round(linspace(1, length(vBgFrames), T.nBaselineFrames))));

for fi = 1:length(vBgFrames)
    f = vBgFrames(fi); % abs frame#
    for c = 0:(T.nNumColChans-1)
        fc = f + c; % abs frame# of c'th channel
        if ishandle(hWait)
            waitbar(fc/T.nVideoNumFrames, hWait)
        else
            break
        end
        if isfield(T, 'vROI') % read ROI of frame
            mFrame = double(ISI_readStreamer(FID, fc, T.vVideoRes, T.vROI));
        else
            mFrame = double(ISI_readStreamer(FID, fc, T.vVideoRes));
        end
        mBaselineSums(:, :, c+1) = mBaselineSums(:, :, c+1) + mFrame;
    end
    if ~ishandle(hWait), break; end
end
if ishandle(hWait), close(hWait); end
T.mBaseline = mBaselineSums ./ fi;

% Calculate radial profiles of channel baselines
% Note: It is assumed frame has equal width and height
hWait = waitbar(0, 'Calculating radial intensity profiles...');

T.nRadIntensityPxBin = 8;
T.vRadIntensity = [];

mD = zeros(T.vVideoRes(1:2) / T.nRadIntensityPxBin);
nO = T.vVideoRes(1) / 2; % origin

for c = 1:T.nNumColChans
    % Create distance matrix
    vI = 1:T.nRadIntensityPxBin:T.vVideoRes(1);
    vJ = 1:T.nRadIntensityPxBin:T.vVideoRes(2);
    for i = 1:length(vI)
        for j = 1:length(vJ)
            mD(i, j) = sqrt((nO-vI(i))^2 + (nO-vJ(j))^2);
        end
    end
    mD = mD - min(mD(:));
    vD = unique(mD(:))';
    
    % Reduce resolution further by limiting radial steps
    nSteps = 100;
    vD = vD(unique(round(linspace(1, length(vD), nSteps))));
    T.nRadIntensityDist = vD;
    
    % Resize baseline image to match size of distance matrix
    mBaseline = imresize(T.mBaseline(:, :, c), size(mD));
    
    % Compute average intensity value at each radial distance
    for i = 1:length(vD)
        if ~ishandle(hWait), break; end
        waitbar((i * c)/(length(vD) * T.nNumColChans), hWait)
        T.vRadIntensity(c, i) = mean(mBaseline(mD == vD(i)));
    end
    
    % Fit gaussians to radial intensity profiles
    tOptions = statset('Display', 'off', 'MaxIter', 10, 'TolFun', 1e-2, 'TolX', 1e-4);
    fGaussFun = inline('b(1) .* (1/sqrt(2*pi)/b(3)*exp(-(x-b(2)).^2/2/b(3)/b(3)))', 'b', 'x');
    
    vX = [-vD vD];
    vY = [T.vRadIntensity(c, :) T.vRadIntensity(c, :)];
    vY = vY - min(vY);
    vY = vY ./ max(vY);
    
    vB = nlinfit(vX, vY, fGaussFun, [500 0 200], tOptions);
    %vYi = fGaussFun(vB, vX);
    
    % Create radial mask
    mG = fspecial('gaussian', T.vVideoRes(1:2), vB(3));
    mG = mG ./ max(mG(:));
    vY = T.vRadIntensity(c, :);
    mG = mG .* (max(vY) - min(vY));
    mG = mG + min(vY);
    T.mRadIntensityMask(:, :, c) = mG;
    T.vRadIntensityMaskSigma(c) = vB(3);
end
if ishandle(hWait), close(hWait); end

% Initialize circular frame buffer
T.mCircBuffer = nan(T.nNumColChans, T.vVideoRes(2), T.vVideoRes(1), T.nCircBufferSize);
T.vCircBufferIndx = nan(T.nNumColChans, T.nCircBufferSize);

% Get order of frames based on standard deviation
T.vChOrder = 1:T.nNumColChans;
if length(T.nNumColChans) == 2
    mA = double(ISI_readStreamer(FID, T.nSkipInitialFrames+1, T.vVideoRes));
    mA = mA - T.mRadIntensityMask(:, :, 1);
    mB = double(ISI_readStreamer(FID, T.nSkipInitialFrames+2, T.vVideoRes));
    mB = mB - T.mRadIntensityMask(:, :, 2);
    if std(mA(:)) > std(mB(:)) % NORMAL :TO USE WHEN THE BLUE IMAGE IS MORE VARIABLE CELINE COMMENTED 20140802 TO CORRECT FOR cm_142 THAT HAD A VERY VARIANT RED IMAGE
        %if std(mA(:)) < std(mB(:)) % EXOTIC CELINE TO USE WHEN THE BLUE IMAGE IS LESS VARIABLE CELINE COMMENTED 20140802 TO CORRECT FOR cm_142 THAT HAD A VERY VARIANT RED IMAGE
        T.vChOrder = fliplr(T.vChOrder); % order is [red blue ...]
    end
end
T.mBaseline = T.mBaseline(:, :, T.vChOrder);

% Initialize figure
[T.hFig, T.hAx, T.hImg, T.hGraphs, T.hFrameSlider, T.hStopButton] = InitFig(T);

for i = 1:length(T.hGraphs)
    set(T.hGraphs(i), 'xdata', 1:T.nVideoNumFrames, 'ydata', nan(T.nVideoNumFrames, 1))
end
T.mGraphs = nan(T.nVideoNumFrames, T.nNumColChans);

set(T.hFrameSlider, 'value', (T.nSkipInitialFrames+1), ...
    'min', (T.nSkipInitialFrames+1), 'max', T.nVideoNumFrames)

set(T.hAx(end), 'xlim', [T.nSkipInitialFrames+1 T.nVideoNumFrames] ./ T.nVideoFPS);

T.hFrameLegend = text(1, 1, '000 s', 'parent', T.hAx(1));
set(T.hFrameLegend, 'units', 'normalized', 'position', [.02 .96 0],...
    'Color', [1 1 1]);

% Iterate over video frames
vFrames = (T.nSkipInitialFrames+1):T.nVideoNumFrames;
T.f = get(T.hFrameSlider, 'value');
f = 0;
T.vCLim = [0 1];

%% Initialize AVI container
if exist('bSaveAVI', 'var')
    if bSaveAVI
        % Use VideoWriter if available, else avifile (deprecated)
        sVideoFile = fullfile(T.sPath, 'QuickLookMovie.avi');
        if exist('VideoWriter', 'file')
            oVideo = VideoWriter(sVideoFile);
            oVideo.FrameRate = nAVIFPS;
            open(oVideo);
        else
            oVideo = avifile(sVideoFile);
            oVideo.Compression = 'none';
            oVideo.FPS = nAVIFPS;
        end
    end
end

%%
nCounter = 0;
while f <= round( (length(vFrames) - T.nFramestep) / 2) * 2
    % Quit loop if window was closed
    if ~ishandle(T.hFig), break; end
    if get(T.hStopButton, 'value'), break; end
    
    % Check if slider was moved
    if get(T.hFrameSlider, 'value') ~= T.f && (f ~= 0)
        f = find(vFrames == round(get(T.hFrameSlider, 'value')));
    else
        f = f + T.nFramestep;
    end
    if f > length(vFrames), break; end
    T.f = vFrames(f);
    
    % Check if frame is in list of ignored frames
    if isfield(T, 'vIgnoreFrames')
        if any(T.vIgnoreFrames == T.f)
            set(T.hFrameSlider, 'value', T.f)
            continue
        end
    end
    
    T.mFrame = double(ISI_readStreamer(FID, T.f, T.vVideoRes));
    
    % Determine channel
    T.nCh = mod(T.f - 1, T.nNumColChans) + 1;
    
    % Update circular frame buffer
    nCircIndx = mod(f-1, T.nCircBufferSize*T.nNumColChans) ./ T.nNumColChans + 1 ...
        - mod(f-1, T.nNumColChans) ./ T.nNumColChans;
    T.vCircBufferIndx(T.nCh, nCircIndx) = T.f;
    T.mCircBuffer(T.nCh, :, :, nCircIndx) = T.mFrame;
    
    % Run sandbox function
    T = T.SandboxCallback(T);
    
    % Update frame display
    if T.bView && ishandle(T.hFig)
        % Update slider position
        set(T.hFrameSlider, 'value', T.f) % <- important to not skip this in the loop!
        
        if (mod(T.f, T.nDisplaystep) == 0)
            set(T.hImg(T.vChOrder(T.nCh)), 'cdata', T.mFrame)

            % Set axis limits
            vXlim = get(T.hAx(T.vChOrder(T.nCh)), 'xlim');
            if any(vXlim > size(T.mFrame))
                set(T.hAx(T.vChOrder(T.nCh)), 'clim', T.vCLim, ...
                    'xlim', [0 size(T.mFrame, 1)], ...
                    'ylim', [0 size(T.mFrame, 2)] );
            else
                % Update only CLim
                set(T.hAx(T.vChOrder(T.nCh)), 'clim', T.vCLim)
            end
            
            % Update time stamp
            set(T.hFrameLegend, 'string', sprintf('%4.1f s', (T.f / T.nVideoFPS)))
            
            % Update graphs
            if isfield(T, 'mGraphs')
                nStart = (T.nSkipInitialFrames) + T.nCh;
                mYY = T.mGraphs(nStart:end, T.nCh);
                mXX = (nStart:size(T.mGraphs, 1)) ./ T.nVideoFPS; % s
                iRem = isnan(mYY);
                mYY(iRem) = [];
                mXX(iRem) = [];
                set(T.hGraphs(T.nCh), 'xdata', mXX, 'ydata', mYY);
                set(T.hAx(T.nNumColChans + T.nCh), 'ylim', T.vCLim);
            end
            
            % Grab frame for AVI file
            if exist('bSaveAVI', 'var')
                if bSaveAVI
                    set([T.hFrameSlider T.hStopButton], 'visible', 'off')
                    try mFrame = getframe(T.hFig);
                    catch, break; end
                    try     writeVideo(oVideo, mFrame);
                    catch,  oVideo = addframe(oVideo, mFrame); end
                end
            end
            
        end
    else
        % Display progress in command window
        nCounter = nCounter + 1;
        if nCounter > 1000
            disp(sprintf('Frame# %d', T.f))
            nCounter = 0;
        end
    end
    drawnow
end
%%
fclose(FID);

% Write frames to disk and close AVI file
if exist('bSaveAVI', 'var')
    if bSaveAVI
        try     close(oVideo);
        catch,  oVideo = close(oVideo); end
    end
end

% Save T to disk
if T.bSave
    sOutputFile = fullfile(T.sPath, [T.sFileName '.mat']);
    % TODO Always ask before overwriting an existing MAT file
    if exist(sOutputFile, 'file')
        sAns = questdlg('Do you want to overwrite the existing data file?', ...
            'Streamer Quick Look', ...
            'Do not overwrite', 'Overwrite', 'Backup old, and save', 'Do not overwrite');
        switch sAns
            case 'Do not overwrite'
                sOutputFile = [];
            case 'Backup old, and save'
                % Backup old file
                sBackupFile = fullfile(T.sPath, [T.sFileName datestr(now, '_HHMM_mmddyy') '.mat']);
                copyfile(sOutputFile, sBackupFile);
        end
    end
    if ~isempty(sOutputFile)
        save(sOutputFile, 'T');
        fprintf('Results saved to %s', [T.sPath T.sFileName '.mat'])
    end
end

return

% Initialize figure with one axis for each color channel
function [hFig, hAx, hImg, hGraphs, hSlider, hToggle] = InitFig(T)
hFig = figure;
set(hFig, 'toolbar', 'figure', 'renderer', 'painters', 'color', 'w')

for i = 1:T.nNumColChans
    hAx(i) = axes('position', [(i-1)*(1/T.nNumColChans) 0.3 (1/T.nNumColChans) .7]);
    hImg(i) = imagesc(zeros(T.vVideoRes(1:2)), 'parent', hAx(i));

    if size(T.mColMap, 3) > 1
        colormap(hAx(i), T.mColMap(:,:,i))
    else
        colormap(hAx(i), T.mColMap)
    end
end
axis(hAx, 'off')
axis(hAx, 'image')

% Plotting axes
mCols = [1 0 0; 0 0 1; 0 1 0; 1 0 1];
for i = 1:T.nNumColChans
    hAx(length(hAx)+1) = axes('position', [0 .05 1 .25]);
    hGraphs(i) = plot(hAx(end), nan, '-', 'color', mCols(i, :));
    %set(hAx(end), 'fontsize', 8, 'xcolor', mCols(i, :), 'ycolor', mCols(i, :), 'color', 'none');
    set(hAx(end), 'fontsize', 8, 'ycolor', mCols(i, :), 'color', 'none');
end
set(hAx(end-i+1), 'color', [1 1 1]);
linkaxes(hAx((end-i+1):end), 'x')
h = zoom(hFig);
setAxesZoomMotion(h, hAx((end-i+1):end), 'horizontal')

% Frame selector
hSlider = uicontrol(hFig, 'style', 'slider');
set(hSlider, 'units', 'normalized' ...
    , 'Position', [0 0 .95 .05] ...
    , 'Style', 'slider');

hToggle = uicontrol(hFig, 'style', 'togglebutton');
set(hToggle, 'units', 'normalized' ...
    , 'Position', [0.95 0 .05 .05] ...
    , 'Style', 'togglebutton', 'background', 'r');

return


