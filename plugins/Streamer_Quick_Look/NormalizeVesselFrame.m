function T = NormalizeVesselFrame(T)
% Normalize vessel image
%
%
% Per M Knutsen <pmknutsen@gmail.com>
% June 2014
%

% Normalize baseline
% Baseline is only used for estimating image intensity limits below
persistent p_mBaselineNorm p_vCLim;
if isempty(p_mBaselineNorm)
    p_mBaselineNorm = T.mBaseline(:, :, T.nCh) ./ T.mRadIntensityMask(:, :, T.nCh);
    p_vCLim = [min(p_mBaselineNorm(:))*1.2 max(p_mBaselineNorm(:))*1];
    T.mBaselineNorm = p_mBaselineNorm;
end

% Normalize image by the average radial intensity mask
T.mFrame = T.mFrame ./ T.mRadIntensityMask(:, :, T.nCh);

% Create an unsharp mask
persistent p_mUnsharpMask
if isempty(p_mUnsharpMask)
    % Mirror pad frame
    mF = T.mFrame;
    mF = [  flipud(fliplr(mF))  flipud(mF)  flipud(fliplr(mF)); ...
        fliplr(mF)              mF          fliplr(mF); ...
        flipud(fliplr(mF))  flipud(mF)  flipud(fliplr(mF)) ];
    nS = 100;
    mWin = fspecial('gaussian', nS*3, nS);
    p_mUnsharpMask = filter2(mWin, mF, 'same');
    p_mUnsharpMask = p_mUnsharpMask((T.vVideoRes(1)+1):(T.vVideoRes(1)*2), ...
        (T.vVideoRes(2)+1):(T.vVideoRes(2)*2));
end

% Subtract unsharp mask to sharpen image
T.mFrame = T.mFrame - p_mUnsharpMask;
T.vCLim = p_vCLim;

return
