function Export_Streamer_file_to_AVI(hObject, eventdata, handles)
% Export selected data file to AVI file
%
% This plug-in is different from the built-in method 'Save AVI movie' as it
% does not require the datafile to be pre-loaded into memory. As such, it
% can export any arbitrary size datafile.
% 
% Plugin only supports .bin data files from the Streamer applications.
%
% Available customization options:
% 
%   Start frame
%   End frame
%   Skip every N'th frame
%   Adjust output framerate
%   Text overlays (frame# and absolute time)
%

% Get filename and path
T.sPath = get(handles.path, 'string');
T.sFileName = get(handles.data_filename, 'string');

% Check that file exists
if ~exist(fullfile(T.sPath, T.sFileName), 'file')
    warndlg('The selected data file could not be found on the filesystem.')
    return
end

% Get data type (.dat or .bin)
T.sType = T.sFileName(end-2:end);
T.sBasefile = T.sFileName(1:end-4);

% Get data info (framerate, resolution etc)
switch T.sType
    case 'bin'

    % Read settings file
    T.sSettingsFile = [T.sPath T.sBasefile '.txt'];
    hFID = fopen(T.sSettingsFile);
    if (hFID == -1)
        warndlg('Failed opening video settings file.');
    end
    [cSettings, ~] = textscan(hFID, '%s');
    mSettings = cSettings{1};
    fclose(hFID);

    % Video resolution (px * px * bitdepth)
    T.vResolution = str2num(mSettings{2, :});
    
    % Framerate and number of frames
    vFrames = str2num(mSettings{3, :});
    T.nFPS = vFrames(1);
    T.nNumFrames = vFrames(2);
    
    case 'dat'
        warndlg('Sorry, this plugin only support .bin files (Streamer).');
        return
end

% TODO Get options from user.
persistent p_nStartFrame p_nEndFrame p_nSkipFrame p_nFPSOut p_bOverlays p_nRescale
if isempty(p_nStartFrame), p_nStartFrame = 1; end
if isempty(p_nEndFrame), p_nEndFrame = T.nNumFrames; end
if isempty(p_nSkipFrame), p_nSkipFrame = 1; end
if isempty(p_nFPSOut), p_nFPSOut = T.nFPS; end
if isempty(p_bOverlays), p_bOverlays = 1; end
if isempty(p_nRescale), p_nRescale = 2; end

cPrompt = {'Start (frame)', ...
    'End (frame)', ...
    'Interval (frames)', ...
    'Framerate (frames/s)', ...
    'Show overlays (1 is yes)', ...
    'Image rescale factor' };

cAns = inputdlg(cPrompt, 'Export Streamer to AVI', 1, ...
    {num2str(p_nStartFrame), num2str(p_nEndFrame), num2str(p_nSkipFrame), ...
    num2str(p_nFPSOut), num2str(p_bOverlays), num2str(p_nRescale)});
if isempty(cAns), return, end

p_nStartFrame = str2num(cAns{1});
p_nEndFrame = str2num(cAns{2});
p_nSkipFrame = str2num(cAns{3});
p_nFPSOut = str2num(cAns{4});
p_bOverlays = str2num(cAns{5});
p_nRescale = str2num(cAns{6});

T.nFPS = p_nFPSOut;
if p_nStartFrame < 1, p_nStartFrame = 1; end
if p_nStartFrame > T.nNumFrames, p_nStartFrame = T.nNumFrames; end
if p_nEndFrame < 1, p_nEndFrame = 1; end
if p_nEndFrame > T.nNumFrames, p_nEndFrame = T.nNumFrames; end

% Initialize figure for displaying video frames
hFig = figure;
hAx = axes();
hImg = imagesc(ones(T.vResolution(1:2)), 'parent', hAx);
truesize(hFig, T.vResolution(1:2) ./ p_nRescale)

colormap(hAx, gray(256))
axis(hAx, 'image', 'off')
if p_bOverlays
    hTxt = text(10, 30, 'Frame# 0', 'parent', hAx, 'color', 'r', 'fontsize', 12);
end

% Initiate object for storing AVI
% Note: Use VideoWrite if it is available, otherwise use the older
% (deprecated) AVI write function avifile
T.sVideoFile = fullfile(T.sPath, [T.sBasefile '.avi']);
if exist('VideoWriter', 'file')
    oVideo = VideoWriter(T.sVideoFile);
    oVideo.FrameRate = T.nNumFrames;
    open(oVideo);
else
    oVideo = avifile(T.sVideoFile);
    oVideo.Compression = 'none';
    oVideo.FPS = T.nFPS;
end

% Open handle to data file
FID = fopen(fullfile(T.sPath, T.sFileName));

% Iterate over frames, display, grab & save
for f = p_nStartFrame:p_nSkipFrame:p_nEndFrame
    % Read frame
    mImg = ISI_readStreamer(FID, f, T.vResolution);
    
    % Check if figure still exists
    if ishandle(hFig)
        % Display frame
        set(hImg, 'cdata', mImg)

        % Text overlays
        if p_bOverlays
            nSec = f / T.nFPS; % absolute time in seconds
            set(hTxt, 'string', sprintf('Frame# %d, %.2f s', f, nSec));
        end

        drawnow
    end
    
    % Get video frame
    try
        mFrame = getframe(hAx);
    catch MExcept
        break
    end
    
    if exist('VideoWriter', 'file')
        writeVideo(oVideo, mFrame);
    else
        try
            oVideo = addframe(oVideo, mFrame);
        catch MExcept
            warndlg('An error occurred when adding a frame.')
            break
        end
    end
    
end

% Close open file handles
oVideo = close(oVideo);
fclose(FID);

% Close figure
if ishandle(hFig)
    close(hFig)
end

return



%% OLD CODE BELOW

if ~isfield(ISIdata,'deltaSignal')
    error('deltaSignal does not exist');
end

hmovie = figure('color',[.7 .7 .7],'name',['Trial averaged moving average ' prmts.name]);
set(hmovie,'DoubleBuffer','on');
set(gca,'xlim',[1 ISIdata.frameSizeYX(2)],'ylim',[1 ISIdata.frameSizeYX(1)],...
    'NextPlot','replace','Visible','off');
colormap(gray);

mvi = 0; % frame counter for movie
clear MOVIE;
secperbin=ISIdata.bin_duration/ISIdata.frame_rate;

sFilename = fullfile(prmts.path2dir,prmts.name);
[path, name, ext] = fileparts(sFilename);
if ~strcmpi(ext,'.avi')
    k=strfind(sFilename,ext);
    if ~isempty(k)
        sFilename=strcat(sFilename(1:k(end)-1),'.avi');
    end
end
if isempty(ext)
    sFilename = strcat(sFilename,'.avi');
end

mSize = [];
for fi = 1:ISIdata.nFramesPerTrial
    img = ISIdata.deltaSignal(:,:,fi);
    
    % Smooth image
    if ~isnan(prmts.smoothSigma)
        mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
        img = single(filter2(mWin, img, 'same'));
    end
    
    % Mask frame
    if prmts.useManualMask
        tMask = prmts.manualMask;
        if ~isempty(tMask)
            img(~tMask.mROI) = NaN;
            % Crop image
            [vI, vJ] = find(~isnan(img));
            vI = [min(vI) max(vI)];
            vJ = [min(vJ) max(vJ)];
            img = imcrop(img, [vJ(1) vI(1) diff(vJ) diff(vI)]);
        end
    end
    
    % Replace NaN's with zeros
    img(isnan(img)) = 0;
    
    if ~ishandle(hmovie), return, end
    figure(hmovie)
    imagesc(img);
    
    title(sprintf('Frame # %d', fi));
    text(4,ISIdata.frameSizeYX(1)-10, sprintf('%5.2f s', fi*secperbin),'color','k','backgroundcolor','w','fontweight','bold')
    set(gca,'clim',ISIdata.climAll); colorbar ('horizontal');
    if fi > ISIdata.nPreStimFrames  && fi <= (ISIdata.nPreStimFrames+ISIdata.nStimFrames)
        w = floor(ISIdata.frameSizeYX(1)*0.02); % make patch width 2% of frame
        patch([10 10 10+w 10+w],[10 10+w 10+w 10], 'red');
    end
    axis image; axis off
    truesize(gcf, size(img))
    
    mvi = mvi + 1;
    MOVIE(mvi) = getframe(gca);
    mSize(end+1, :) = size(MOVIE(mvi).cdata);
    drawnow
end

% Crop all frames in MOVIE to smallest frame (may vary due to buggy matlab plotting)
vMinSize = min(mSize, [], 1);
for f  = 1:length(MOVIE)
    MOVIE(f).cdata = MOVIE(f).cdata(1:vMinSize(1), 1:vMinSize(2), 1:1:vMinSize(3));
end

try
    movie2avi(MOVIE, sFilename)
catch ME
    avimov = close(avimov);
    error(ME.message);
end

return