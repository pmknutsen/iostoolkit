function ISI_viewMask(handles)
% View the masks over the vessel image
%
% Only the selected masks are displayed, ie. vessel or manual masks or the
% combination of both
%
% Ssed to compare actual coverage, see what areas are being masked out.
%

global ISI

% Get current parameters
tParams = ISI.GetParams(handles);
tParams = tParams.filesQueue(1);

% Load vessel image
sFilename = fullfile(handles.pathstr, get(handles.vessel_filename,'string'));

if ~exist(sFilename, 'file')
    warndlg('Vessel image not found.', 'IOS Toolkit');
    return
end

mVesselImgOrig = imread(sFilename);
mVesselImg = mat2gray(mVesselImgOrig);

% Get vessel mask
mVesselMask = ISI_createVesselMask(handles);
mVesselImg(~mVesselMask) = NaN;
mMask = mVesselMask;

% Get manual mask
if tParams.useManualMask
    tMask = tParams.manualMask;
    if ~isempty(tMask)
        mVesselImg(~tMask.mROI) = NaN;
        mMask = mMask | tMask.mROI;
    end
end

% Remove zeros and ones
mVesselImg(mVesselImg == 0) = NaN;
mVesselImg(mVesselImg == 1) = NaN;

% Initialize figure
hFig = findobj('tag', 'ISI_maskPreview');
if isempty(hFig)
    hFig = figure('visible', 'off');
    set(hFig, 'tag', 'ISI_maskPreview')
else
    figure(hFig)
end
set(hFig, 'position', [1 1 750 400])
centerfig(hFig)
set(hFig, 'visible', 'on')

% Original image
hAx = subplot(1,3,1);
imagesc(mVesselImgOrig);
axis(hAx, 'image', 'off');
title('Original')
colormap(hAx, gray(256))

% Mask
hAx = subplot(1,3,2);
imagesc(mMask);
axis(hAx, 'image', 'off');
title('Mask');

% Masked image
hAx = subplot(1,3,3);
hv = imshow(mVesselImg);  
axis(hAx, 'image', 'off');
hold(hAx, 'on')
hmask = imshow(mMask);
set(hmask,'alphadata', 0.3, 'CData', mMask);
title('Masked image')

return