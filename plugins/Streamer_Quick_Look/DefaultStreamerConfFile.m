%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User variables
% Note: Don't modify these inside your Sandbox function.
% Save these variables as a Matlab (.m) script

% Pixel resolution (microns)
T.nPixRes = 0.6;

% Circular buffer size (N last frame to keep in memory)
T.nCircBufferSize = 5;

% Number of baseline frames
T.nBaselineFrames = 500;

% Number of initial frames to skip
T.nSkipInitialFrames = 13000;

% Step between analysed frames
T.nFramestep = 1;

% Step between displayed frames
% If nNumColChans is an even number this number should be t should be an
% odd number
T.nDisplaystep = 7;

% Sandbox function - @YourFunctionName
% Any Matlab function that you want to execute for each frame.
% Your function should accept T as its sole input parameter.
T.SandboxCallback = @IOSAnalysis;

% Number of color channels (typically 2)
T.nNumColChans = 2;

% Colormap
% Different colormap can be used for each color channel (requires R2014B)
T.mColMap(:,:,1) = fireice(2^8);
T.mColMap(:,:,2) = gray(2^8);

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
% Save frames to AVI file
bSaveAVI = 1;
nAVIFPS = 25;

% End of user variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
