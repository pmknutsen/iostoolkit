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
T.nSkipInitialFrames = 20;

% Number of frames to step between
T.nFramestep = 1;

% Sandbox function - @YourFunctionName
% Any Matlab function that you want to execute for each frame.
% Your function should accept T as its sole input parameter.
T.SandboxCallback = @IOSAnalysis;

% Number of color channels (typically 2)
T.nNumColChans = 2;

% Colormap
T.mColMap = gray(2^6);

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

% End of user variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
