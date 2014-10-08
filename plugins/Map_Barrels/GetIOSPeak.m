function [mMaskMap, mMap] = GetIOSPeak(sFile, nThresh)
% Display average IOS map and select ROI interactively
%
% Usage:
%   After map loads in figure, click once on the center of the IOS peak and
%   then select the radius of a circular ROI.
%
%

% Load .dat file
load(sFile)

% Check that ISIdata structure was loaded. If not, then exit.
if ~exist('ISIdata', 'var')
    disp(sprintf('GetIOSPeak: %s is not an IOS file.', sFile))
    mMaskMap = [];
    mMap = [];
    return
end

% Invert map (signals are positive)
mMap = double(ISIdata.signalFrame .* -1);

% Get sigma factor for smoothing of maps from GUI
nSigma = str2double(get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'smooth_sigma'), 'string'));

% Smooth map (sigma set in GIU)
if isnumeric(nSigma) && ~isempty(nSigma) && ~isnan(nSigma)
    mWin = fspecial('gaussian', nSigma*3, nSigma);
    mMap = single(filter2(mWin, mMap, 'same'));
end

% Apply manual mask (if set in GUI)
if get(findobj('tag', 'chk_manualmask'), 'value') == 1
    mManualMask = get(findobj('tag', 'btn_setmanualmask'), 'userdata');
    if ~isempty(mManualMask) && all(size(mManualMask.mROI) == size(mMap))
        mMap(~mManualMask.mROI) = NaN;
    end
end

% Determine clim from 1st and 99th percentile of intensity values
% TODO Get CLIM from GUI
vCLim = prctile(mMap(:), [.5 99.5]);

% Initialize figure
sTag = 'map_barrels_mark_ios_region';
hFig = findobj('tag', sTag);
if isempty(hFig)
    hFig = figure;
    set(hFig, 'tag', sTag);
end
clf(hFig)

% Display average IOS map
hAx = axes();
colormap(gray(256))
imagesc(mMap, 'parent', hAx)
hold(hAx, 'on')
set(hAx, 'clim', vCLim)
hTit = title(sFile);
set(hTit, 'interpreter', 'none');
axis(hAx, 'equal', 'off')
set(hAx, 'ydir', 'reverse')

% Get ROI interactively. The ROI should be drawn conservatively around the
% activated region.

% Get circular ROI
title(hAx, sprintf('%s\nSelect center of activated region\nPress Enter to repeat the previous map', sFile))
[X1, Y1] = ginput(1);
if isempty(X1)
    mMaskMap = 1;
    mMap = [];
    return
end

xp = [X1 Y1];
hCirc = plot(X1, Y1, 'r.');
title(hAx, sprintf('%s\nResize circle to include the entire activated region\n', sFile))
set(hFig, 'WindowButtonMotionFcn', {@mousemove, hCirc, xp, hAx});
k = waitforbuttonpress;
set(hFig, 'WindowButtonMotionFcn', '');

% Convert ROI to a mask
vX = get(hCirc, 'xdata');
vY = get(hCirc, 'ydata');
mBW = poly2mask(vX, vY, size(mMap, 1), size(mMap, 2));

%mBW = roipoly; % get polygon ROI

if isempty(mBW)
    mMaskMap = [];
    mMap = [];
    return
end

mMaskMap = mMap .* mBW;
mMaskMap(mMaskMap == 0) = NaN;

return


function mousemove(object, eventdata, hCirc, bp, hAx)
cp = get(hAx, 'CurrentPoint');
r = norm([cp(1,1) - bp(1) cp(1,2) - bp(2)]);
theta = 0:.1:2*pi;
xc = r * cos(theta)+bp(1);
yc = r * sin(theta)+bp(2);
set(hCirc, 'XData', xc);
set(hCirc, 'YData', yc);
return