function mMask = ISI_createVesselMask(handles)
% Create and return a vessel mask
% 
% Usage:
%   M = ISI_createVesselMask(P)
% 
%   where M is the vessel mask and P is the parameters structure (tParams)
%
% TODO
%   Improve vessel mask by recycling code from the Streamer_Quick_Look
%   plugin
%   Smooth before thresholding
%


global ISI

% Get current parameters
tParams = ISI.GetParams(handles);

% Load vessel image
mImg = fullfile(tParams.filesQueue.path2dir, tParams.filesQueue.refImage);
vThresh = tParams.filesQueue.maskthresh;

ISIdata = ISI.GetISIData(handles);

% Create mask based on threshold only
if ~tParams.useVesselMask || ~exist(mImg, 'file')
    mMask = ones(ISIdata.frameSizeYX);
else
    mMask = imread(mImg);
    mMask = mat2gray(mMask); % we dont need to normalize, but do it for clarity
    mMask = im2bw(mMask, vThresh);
end

return