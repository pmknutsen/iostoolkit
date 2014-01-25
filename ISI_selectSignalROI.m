function ISIdata = ISI_selectSignalROI(ISIdata, prmts, cmap)
% ISI_selectSignalROI
% Manually select one or more regions of interest (ROI) on vessel or signal
% images.
%
%
%

persistent p_tROI

% Display reference image and let user select ROI for analysis.
if ~isfield(ISIdata, 'signalFrame')
    % Get the average frame over stimulus interval, suppressing the comparison figure
    ISIdata.signalFrame = getISIsignalframe(ISIdata, prmts.stimInterval, 0);
end
mFrame = ISIdata.signalFrame;
mWin = ones(prmts.imgBin);

% Smooth image
if ~isnan(prmts.smoothSigma)
    mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
    mFrame = single(filter2(mWin, mFrame, 'same'));
end

% Read vessel image
sFile = fullfile(prmts.path2dir, prmts.refImage);
mVesselImg = imread(sFile);
mVesselImg = mat2gray(mVesselImg);

% Display vessel image
hFig = figure('Name', 'Select Region of Interest', 'NumberTitle', 'off');

vPos = get(hFig, 'position');
set(hFig, 'position', [vPos(1:2) 700 350]);

mCols = get(hFig, 'DefaultAxesColorOrder');
mCols = [mCols; mCols ./ 2];
hAx(1) = subplot(1,2,1);
imshow(mVesselImg)
axis image off; hold on
hTit = title('Vessel Image');
PlotROIs(p_tROI, mCols);

% Display signal image
hAx(2) = subplot(1,2,2);
imagesc(mFrame)
axis image off; hold on
if isfield(ISIdata, 'climAll')
    set(hAx(2), 'Clim', prmts.climAll./1000)
elseif isfield(ISIdata, 'climAll')
    set(hAx(2), 'Clim', ISIdata.climAll)
end
hTit(2) = title('Average Signal Image');
PlotROIs(p_tROI, mCols);


colormap(cmap)

sAns = questdlg('Choose image to draw ROI on', 'ROI Image', ...
    'Vessel Image', 'Average Signal Image', 'Average Signal Image');
if isempty(sAns), return; end

if strcmp(sAns, 'Vessel Image')
    axes(hAx(1))
    set(hAx(1), 'xcolor', 'r', 'ycolor', 'r', 'xtick', [], 'ytick', [], 'visible','on')
else
    axes(hAx(2))
    set(hAx(2), 'xcolor', 'r', 'ycolor', 'r', 'xtick', [], 'ytick', [], 'visible','on')
end

hRemObj = findobj('Tag', 'REMOVE');
if ~isempty(hRemObj)
    xlabel('Press ESC to use last-set ROI(s)', 'color', 'r')
    delete(hRemObj)
end

cROImask = {};
cXi = {};
cYi = {};
while 1
    [cROImask{end+1}, cXi{end+1}, cYi{end+1}] = roipoly;
    axes(hAx(1))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cROImask), :))
    axes(hAx(2))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cROImask), :))

    % Check if ESC was pressed
    if isempty(cROImask{end})
        break
    else
        sAns = questdlg('Do want to draw one more ROI?', 'ROI', ...
            'Yes', 'No', 'No');
        if isempty(sAns) || strcmp(sAns, 'No')
            break
        end
    end
end

% Plot ROIs
axes(hAx(1));
PlotROIs(p_tROI, mCols);

axes(hAx(2));
PlotROIs(p_tROI, mCols);
drawnow

% Store ROI in persistent variable
if ~isempty([cROImask{:}])
    p_tROI = struct('mROI', cROImask, 'vXi', cXi, 'vYi', cYi);
end
ISIdata.analysisSignalROI = p_tROI;

return


function PlotROIs(p_tROI, mCols)
% Plot last-set ROIs
if ~isempty(p_tROI)
    for i = 1:length(p_tROI)
        plot(p_tROI(i).vXi, p_tROI(i).vYi, 'o-', ...
            'color', mCols(i,:), 'Tag', 'REMOVE')
    end
end
return
