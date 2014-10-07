function [mFrames, varargout] = ISI_readStreamer(sFile, varargin)
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
% Written by Per M Knutsen <pmknutsen@gmail.com>
% June 2014
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
    
    % Data file string
    sDataFile = [sPath sFileName '.bin'];
    
    % Settings file string
    sSettingsFile = [sPath sFileName '.txt'];
    
    % Read acquisition settings
    hFID = fopen(sSettingsFile);
    if (hFID == -1), error('Failed opening settings file.'); end
    
    [cSettings, ~] = textscan(hFID, '%s');
    mSettings = cSettings{1};
    fclose(hFID)

    % Video resolution (px * px * bitdepth)
    vResolution = str2num(mSettings{2, :});
    
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

% Read frames in range
vFrames = nFrom:nTo;
if vResolution(3) > 8
    mFrames = uint16(zeros(vResolution(1), vResolution(2), length(vFrames)));
else
    mFrames = uint8(zeros(vResolution(1), vResolution(2), length(vFrames)));
end

% Get size of frame in bytes
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
    fseek(hFID, (f-1) * nFrameBytes, -1); % does this slow things down?
    try
        mFrames(:, :, vFrames == f) = fread(hFID, vResolution(1:2), sType);
    catch mExcep
        error('Failed assigning frame. May have reached end of file...')
    end
end

if ~isnumeric(sFile)
    fclose(hFID);
end

% Populate varargout with a structure containing video info (resolution etc)
tInfo(1).vResolution = vResolution(1:2);
tInfo(1).sFile = sFile;
tInfo(1).nBitDepth = vResolution(3);
tInfo(1).nBytesPerFrame = nFrameBytes;
tInfo(1).nNumMissedFrames = nNumMissedFrames;
tInfo(1).sType = sType;
tInfo(1).nNumFrames = nNumFrames;
tInfo(1).nFPS = nFPS;
tInfo(1).nStartFrame = nFrom;
tInfo(1).nEndFrame = nTo;
varargout{1} = tInfo;

return
