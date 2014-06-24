function [mMask] = GetRegionMask(mFrame, nPixRes, nBinThresh, nMedFilt)
% Find vessels in an image
% Automated steps:
%   1   Median filter at a resolution of 5-by-5 microns
%   2   Binarize image
%   3   Grow mask by local percentile decrease in threshold in 5 loops
%
% Usage:
%   GetRegionMask(mFrame, nPixRes, nBinThresh, nMedFilt)
%       where
%       mFrame      Image to estimate mask form
%       nPixRes     Pixel resolution of the image in microns/pix
%       nBinThresh  Grayscale value to threshold the median filtered image.
%                   If nBinThresh is NAN or empty ([]), a GUI appears so that
%                   the value can be selected interactively with a slider.
%       nMedFilt    Size of the median filter in microns
%
%
% Per M Knutsen <pmknutsen@gmail.com>
% June 2014
%

vVideoRes = size(mFrame);
nMedFilt = round(nMedFilt / nPixRes);

% Median filter and normalize image (0 -> 1)
mImg = mFrame;
mImg = medfilt2(mImg, [nMedFilt nMedFilt]);
mImg = mImg - min(mImg(:));
mImg = mImg ./ max(mImg(:));

% Binarize image (select threshold interactive if missing or NAN)
if isempty(nBinThresh) || isnan(nBinThresh)
    hFig = figure;
    hSlider = uicontrol(hFig, 'style', 'slider');
    set(hSlider, 'units', 'normalized' ...
        , 'Position', [0 0.02 1 .05] ...
        , 'Style', 'slider', 'value', .5 ...
        , 'SliderStep', [0.01 .1]);
    set(hFig, 'toolbar', 'figure')
    
    while ishandle(hFig)
        figure(hFig)
        colormap(gray)
        hAx(1) = subplot(1,2,1);
        imagesc(mImg)
        axis image off
        
        hAx(2) = subplot(1,2,2);
        mMask = im2bw(mImg, get(hSlider, 'value'));
        mMask = FineTuneMask(mMask, vVideoRes, nMedFilt);
        
        imagesc(mMask)
        axis image off

        nBinThresh = get(hSlider, 'value');
        title(hAx(2), sprintf('Threshold = %.2f  Black is ON', nBinThresh));
        title(hAx(1), 'Close window when done');
        
        waitfor(hSlider, 'value')
    end
    %%
end

mMask = im2bw(mImg, nBinThresh);
mMask = FineTuneMask(mMask, vVideoRes, nMedFilt);

return


function mMask = FineTuneMask(mMask, vVideoRes, nMedFilt)


% Move across image in 2^median_filter sized regions
% In each region,   1 - find mean intensity value of vessel
%                   2 - recruit pixels with near/darker intensity values
nS = 2^nMedFilt; % pixel size
nP = 0.01;
for a = 1:5
    for i = 0:(vVideoRes(1)/nS)-1 % one loop, assume image is square
        for j = 0:(vVideoRes(2)/nS)-1 % one loop, assume image is square
            % cutout
            vI = [i*nS+1 (i+1)*nS];
            vJ = [j*nS+1 (j+1)*nS];
            
            mCutout = mImg(vI(1):vI(2), vJ(1):vJ(2));
            mCutoutMask = mMask(vI(1):vI(2), vJ(1):vJ(2));

            % Get mean vessel intensity
            nAvgInt = min(mCutout(mCutoutMask));

            nThresh = nAvgInt + (nAvgInt*nP);
            mCutoutMask(mCutout < nThresh) = 0;
            mMask(vI(1):vI(2), vJ(1):vJ(2)) = mCutoutMask;
        end
    end
end

% Bridge close, non-connected regions
mMask = ~bwmorph(~mMask, 'bridge', 2);

% Fill holes
mMask = ~imfill(~mMask, 'holes');

% Remove connection regions with an aread smaller than X
IL = bwlabel(~mMask);
R = regionprops(~mMask, 'Area');
ind = find([R.Area] > (nMedFilt^2));
mMask = ~ismember(IL, ind);


return

