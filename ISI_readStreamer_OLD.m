function [mFrames, varargout] = readstreamer(sFile, varargin)
% READSTREAMER Read Streamer binary files into Matlab
%
% M = readstreamer(F)
%   open file in string F and read all frames into matrix M
%
% M = readstreamer(F, S, E)
%   open file in string F and read frames from S to E
%
% M = readstreamer(FID, F, [I J B])
%   read frame number F from an existing file handle FID, where the frame
%   has resolution [I J B], where I and J is the size and B is the
%   bit-depth. This syntax does not support reading a range of frames.
%
% M = readstreamer(..., ROI)
%   read only the region of interest specified in ROI, where ROI is a
%   vector [L T W H] denoting the Left/Top corner and Width/Height of the
%   ROI, respectively.
%
% [M T] = readstreamer(...)
%   Returns the structure T with basic information on the data file.
%
% Note:
%   When readstreamer() is used with a file identifier (FID), the file must
%   be opened and closed externally.
%
% Example:
%   FID = fopen('filename.bin');
%   mFrame = ISI_readStreamer(FID, [1024 1024 12]);
%   fclose(FID);
%   imagesc(mFrame)
%
% TODO:
%   Add optional argument to read an ROI only (will speed up reading)
%
%
% Per M Knutsen <pmknutsen@gmail.com>, 2014
%

% Check inputs
if nargin < 1
    error('Too few input parameters.');
end

% Check if first input parameter is a file identifier
% If not, then read resolution etc from associated text file
if isnumeric(sFile)
    if length(varargin) < 2
        error('Too few input parameters.');
    end
    nFrom = varargin{1}; % frame to read
    nTo = varargin{1};
    hFID = sFile;
    vResolution = varargin{2};
else
    % Open settings file
    [sPath, sFileName, ~] = fileparts(sFile);
    
    % Settings file string
    sSettingsFile = [fullfile(sPath, sFileName) '.txt'];
    
    % Read acquisition settings
    if ~exist(sSettingsFile, 'file')
        error('The video settings file cannot be found.')
        return;
    end
    hFID = fopen(sSettingsFile);
    if (hFID == -1), error('An error occurred when opening the settings file.'); end
    
    [cSettings, ~] = textscan(hFID, '%s');
    mSettings = cSettings{1};
    fclose(hFID);

    % Video resolution (height * width * bitdepth)
    vResolution = str2num(mSettings{2, :}); %#ok<*ST2NM>
    
    % Framerate and number of frames
    vFrames = str2num(mSettings{3, :});
    nFPS = vFrames(1);
    nNumFrames = vFrames(2);

    % Missed frames
    nNumMissedFrames = str2num(mSettings{4, :});
    nNumMissedFrames = nNumMissedFrames(1);
    if (nNumMissedFrames > 0)
        error('There are missed frames. Not programmed to deal with these yet.');
    end

    % Get range of frames to read
    if nargin > 1, nFrom = varargin{1}; else nFrom = 1; end
    if nargin > 2, nTo = varargin{2}; else nTo = nNumFrames; end

    % Open data file
    hFID = fopen(sFile);
    if (hFID == -1), error('Failed opening data file.'); end

end

% Get bytesize of full frame
nFullFrameBytes = prod(vResolution(1:2));
if vResolution(3) > 8
    nFullFrameBytes = nFullFrameBytes * 2;
end

% Check if last input is an ROI vector with a length of 4
bOffset = 0;
vROI = [];
if length(varargin{end}) == 4
    vROI = varargin{end};
    vROI(4) = min([vResolution(2) vROI(2)+vROI(4)]); % height
    vResolution(1) = min([vResolution(1) vROI(1)+vROI(3)]); % width
    
    % Compute the frame byte offset (from top only)
    bOffset = vROI(2) - 1;
    if vResolution(3) > 8
        bOffset = bOffset * 2;
    end
end

% Read frames in range
vFrames = nFrom:nTo;
if vResolution(3) > 8
    mFrames = uint16(zeros(vResolution(2), vResolution(1), length(vFrames)));
else
    mFrames = uint8(zeros(vResolution(2), vResolution(1), length(vFrames)));
end

% Get bytesize of frame part to read
nFrameBytes = prod(vResolution(1:2));
if vResolution(3) > 8
    nFrameBytes = nFrameBytes * 2;
end

% Type conversion
if vResolution(3) > 8
    sType = 'uint16=>uint16';
else
    sType = 'uint8=>uint8';
end

for f = vFrames
    fseek(hFID, ((f-1) * nFullFrameBytes) + bOffset, -1); % does resetting origin to start of file slow things down?
    try
        mFrames(:, :, vFrames == f) = fread(hFID, vResolution([2 1]), sType);
    catch mExcep
        error('Failed assigning frame. May have reached end of file...')
    end
end

% Crop mFrames if vROI is set
if ~isempty(vROI)
    mFrames = mFrames(1:(vROI(4)-vROI(2)), vROI(1):end, :);
end

if ~isnumeric(sFile)
    fclose(hFID);
end

% Populate varargout with a structure containing video info (resolution etc)
tInfo(1).vResolution = vResolution(1:2);
tInfo(1).sFile = sFile;
tInfo(1).nBitDepth = vResolution(3);
tInfo(1).nBytesPerFrame = nFrameBytes;
if exist('nNumMissedFrames', 'var')
    tInfo(1).nNumMissedFrames = nNumMissedFrames;
end
tInfo(1).sType = sType;
if exist('nNumFrames', 'var')
    tInfo(1).nNumFrames = nNumFrames;
end
if exist('nFPS', 'var')
    tInfo(1).nFPS = nFPS;
end
tInfo(1).nStartFrame = nFrom;
tInfo(1).nEndFrame = nTo;
varargout{1} = tInfo;

return