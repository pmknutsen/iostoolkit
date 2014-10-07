function ISIdata = ISI_isolateBarrel(ISIdata, mMapRaw, prmts, saveFig)
% ISIdata=isolatebarrel(ISIdata,signalFrame,contourquantile,selectROI)
%
% given signalFrame and desired quantile for contours, return the contours
% in the ISIdata struct.
%
% If the saveFig flag is set, a .fig file is written out with just the
% smoothed signal image
%
%

if ~isfield(ISIdata,'vesselmask')
    warning('isolatebarrel:nomask','No vessel mask in dataset.');
end
if nargin < 4 || isempty(saveFig) % selectROI flag doesn't exist
    saveFig = 0;
end

% Crop evoked map by an ROI
if isfield(ISIdata,'analysisROI')
    mMap = mMapRaw .* (ISIdata.analysisROI.ROI);
else
    mMap = mMapRaw;
end

% Mask evoked map using the vesselmask (if used)
mMap = mMap .* ISIdata.vesselmask;

% Smooth evoked map using a sigma of at least 2
if ~isnan(prmts.smoothSigma) && prmts.smoothSigma > 0
    F = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
else
    F = fspecial('gaussian', 2*3, 2);
end
mMap = filter2(F, mMap, 'same');

% Mask image with a user-define ROI
if prmts.useManualMask
    tMask = prmts.manualMask;
    if ~isempty(tMask)
        mMap(~tMask.mROI) = NaN;
    end
end

% Initialize figure
hFig = figure('name',['Evoked map - ' prmts.name], 'visible', 'off');
set(hFig, 'position', [0 0 800 350], 'numbertitle', 'off')
centerfig(hFig)
set(hFig, 'visible', 'on')
colormap(gray(256))

% Make copy of the cropped, masked and smoothed evoked map
dipIm = mMap;
ISIdata.signalFrameFilt = dipIm;

hAx = subplot(1, 2, 1);
hImg = imagesc(dipIm, 'parent', hAx);
axis(hAx, 'image', 'off')
hold(hAx, 'on');
title(hAx, sprintf('Evoked map - %s', prmts.Whisker{1}), 'interpreter', 'none');
set(hAx, 'clim', prmts.climAll ./ 1000);

% Create second axis for vessel image
hAxRight = subplot(1, 2, 2);

if saveFig
    savefilename = fullfile(prmts.path2dir, prmts.name);
    [p n e] = fileparts(savefilename);
    saveas(hFig, fullfile(p, [n '.png']), 'png');
end

% Get contours
contour_level = prmts.ContourLineVals;
contour_level = sort(contour_level,'ascend');
nContourLevels = numel(contour_level);

if nContourLevels>1
    caxis(contour_level([1 end])-1)
else
    caxis([-4e-4 4e-4]);
end

% Colorrange=[.00010 .00015 .0002 .00025 .0003 .00035];
contour_color={'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r+-' 'g+-' 'b+-'};
nContColor = length(contour_color);

% Ensure contour levels are in ascending order (lower value -> deoxigenated area)
contourLinesLeftPerLevel = zeros(nContourLevels, 1);
gotLowest = 0;
for iLevel = 1:nContourLevels
    [xcontour{1,iLevel}, h] = contour(hAx, dipIm, [-1 -1]+contour_level(iLevel),contour_color{mod(iLevel,nContColor)+1}, ...
        'LineWidth', 2);

    % Delete small objects
    hChildren = get(h,'children');
    for iChild = 1:numel(hChildren)
        np = numel(get(hChildren(iChild),'Faces'));
        if np < prmts.minFaces || np > prmts.maxFaces
            delete(hChildren(iChild));
        end
    end

    % Check how many contour points left for current level
    hChildren = get(h,'children');
    contourLinesLeftPerLevel(iLevel) = numel(hChildren);

    % Keep lowest for later use
    if contourLinesLeftPerLevel(iLevel) > 0
        if ~gotLowest 
            XdataLowest = get(hChildren, 'Xdata');
            YdataLowest = get(hChildren, 'Ydata');
            if isnumeric(XdataLowest) %if only one contour, still display to give chance to skip
                XdataLowest={XdataLowest};
                YdataLowest={YdataLowest};
            end
            if iscell(XdataLowest) %only one contour
                %                 need to ask user to choose
                %                 keyboard
                %mjp 2011.07.29 added ability to select more than one contour
                idx = ISI_selectContour(XdataLowest, YdataLowest, prmts, dipIm, contour_level(iLevel)-1, hAxRight);
                drawnow;
                if idx > 0
                    Xcontour = XdataLowest{idx(1)};
                    Ycontour =YdataLowest{idx(1)};
                    if numel(idx)>1
                        for ix=2:length(idx)
                            Xcontour=[Xcontour; NaN; XdataLowest{idx(ix)};];
                            Ycontour=[Ycontour; NaN; YdataLowest{idx(ix)};];
                        end
                        Xcontour = [Xcontour; NaN];
                        Ycontour = [Ycontour; NaN];
                    end
                    gotLowest = 1;
                    lowestLevel = iLevel;
                end
            else
                gotLowest = 1;
                lowestLevel = iLevel;
                Xcontour = XdataLowest;
                Ycontour = YdataLowest;
            end

        end
    end
end

% Update figure with selected contour
cla(hAx)
cla(hAxRight)

hImg = imagesc(dipIm, 'parent', hAx);

%if nContourLevels > 1
%    caxis(hAx, contour_level([1 end])-1)
%else
%    caxis(hAx, prmts.climAll ./ 1000);
%end

% Display lowest contour
if isfield(prmts, 'refImage')
    if ~isempty(prmts.refImage)
        imRef = imread(fullfile(prmts.path2dir,prmts.refImage));

        % Compute size difference between filtered "dip" image and reference
        [deltaRC] = size(imRef) - size(dipIm);
        deltaR = fix(deltaRC(1)/2);
        deltaC = fix(deltaRC(2)/2);
        
        imshow(imRef, 'parent', hAxRight)
        x = Xcontour + deltaR;
        y = Ycontour + deltaC;
        ISIdata.contourM = [x y];
        
        title([prmts.name '   Contour = ' num2str(contour_level(lowestLevel)-1)],'interpreter','none');
        
        % Plot selected contour into both axis
        plot(hAx, x, y, 'r-', 'LineWidth', 2)
        plot(hAxRight, x, y, 'r-', 'LineWidth', 2)
        
        mx = mean(x(~isnan(x)));
        my = mean(y(~isnan(y)));
        
        %mjp 2011.09.22 should be centering txt from mean; for now it's offset
        text(mx, my, prmts.Whisker, 'color', 'r', 'fontSize', 10, 'parent', hAx);
        text(mx, my, prmts.Whisker, 'color', 'r', 'fontSize', 10, 'parent', hAxRight);
        
        figname = get(gcf,'name');
        saveas(hFig, fullfile(prmts.path2dir, [figname(1:end-4) '.pdf']), 'pdf');
    end
else
    fprintf('\nNo reference image for %s', prmts.name);
end

return