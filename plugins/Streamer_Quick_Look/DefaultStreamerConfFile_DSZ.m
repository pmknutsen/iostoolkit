%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User variables
% Note: Don't modify these inside your Sandbox function.
% Save these variables as a Matlab (.m) script

% Pixel resolution (microns)
T.nPixRes = 4;

% Circular buffer size (N last frame to keep in memory)
T.nCircBufferSize = 5;

% Number of baseline frames
T.nBaselineFrames = 500;

% Region of interest to read from frames [l t w h]
T.vROI = [200 300 100 150];

% Number of initial frames to skip
T.nSkipInitialFrames = 200;

% Step between analysed frames
T.nFramestep = 1;

% Step between displayed frames
% If nNumColChans is an even number this number should be t should be an
% odd number
T.nDisplaystep = 1;

% Sandbox function - @YourFunctionName
% Any Matlab function that you want to execute for each frame.
% Your function should accept T as its sole input parameter.
T.SandboxCallback = @IOSAnalysisDSZ;

% Number of color channels (typically 2)
T.nNumColChans = 1;

% Colormap
T.mColMap = gray(1024);

% Pixel binning during estimation of radial intensity
% Note: Higher number reduces execution time.
T.nRadIntensityPxBin = 10;

% Boolean parameters
%  0 = no      1 = yes
% Enable intensity mask normalization
T.bEnableNorm = 1;
% View frames
T.bView = 1;
% Load T if .mat file with same name exists on disk
T.bLoad = 1;
% Save T when last frame is processed
T.bSave = 1;

% Frames to ignore (are neither read nor processed)
T.vIgnoreFrames = [];

% Spawn Spiky and open associated .daq file
sFile = fullfile(T.sPath, [T.sFileName '.spb']);
if exist(sFile, 'file')
    % Launch spiky, load file and get data structure
    spiky
    global Spiky
    Spiky.main.OpenFile([], sFile)
    [FV, ~] = Spiky.main.GetStruct();
    
    % Get camera framerate
    iCamShut = strcmp('CameraShutter', {FV.tChannelDescriptions.sDescription});
    sCamShutCh = FV.tChannelDescriptions(iCamShut).sChannel;
    vCamShutUp = FV.tData.([sCamShutCh '_Up']);
    nCamFrameRate = 1 / median(diff(vCamShutUp)); % frames/s

    % Get camera synchronization time
    iExSyncShut = strcmp('ExSyncIn', {FV.tChannelDescriptions.sDescription});
    sExSyncCh = FV.tChannelDescriptions(iExSyncShut).sChannel;
    vExSyncUp = FV.tData.([sExSyncCh '_Up']); % s
    
    % Get whisking frames relative to synchronization signal
    vWhiskingUp = FV.tData.WhiskingFrames_Up - vExSyncUp(1); % s

    % Compute frame numbers
    vWhiskingFrames = round(vWhiskingUp * nCamFrameRate);
    
    % Get frames *without* whisking (to be ignored)
    T.vIgnoreFrames = setdiff(1:length(vCamShutUp), vWhiskingFrames);    
end

% End of user variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
