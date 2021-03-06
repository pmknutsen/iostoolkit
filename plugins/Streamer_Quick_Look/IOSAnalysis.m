function T = IOSAnalysis(T)
% IOSAnalysis
%
% Analyse intrinsic signals (IOS) and blood vessel diameter
%
% Run this function via IOSToolkit plug-in Streamer_Quick_Look()
%

% Process channels
if T.nCh == T.vChOrder(1)   % Red channel
    if ~isfield(T, 'cMasks'), return; end

    % Compute intrinsic signal (dR/R)
    T.mFrame = (T.mFrame - T.mBaseline(:, :, T.vChOrder(T.nCh))) ./ T.mBaseline(:, :, T.vChOrder(T.nCh));
    mFrame_orig = T.mFrame;
    
    % Bin image pixels - OPTIONAL
    T.mFrame = BinImagePixels(T.mFrame, 16);
    if ~all(size(T.mFrame) == size(T.cMasks{T.nCh})) % do once
        T.cMasks{T.nCh} = BinImagePixels(T.cMasks{T.nCh}, 16);
        T.cMasks{T.nCh} = round(T.cMasks{T.nCh});
    end
    
    % Set c-limits on displayed images
    T.vCLim = [-.01 .01]; % signal ratio

    % Mask image (show only non-vessel regions)
    if isfield(T, 'cMasks')
        mFrame = T.mFrame .* ~T.cMasks{T.nCh};
        mFrame(mFrame == 0) = nan;
    end
    
    % Get average signal amplitude across all pixels
    T.mGraphs(T.f, 1) = nanmean(mFrame(:)); %CM_20140610 to grow to have the diameter channel
    
    % Smooth image by gaussian (sigma is set in GUI)
    if isfield(T, 'mWin')
        T.mFrame = filter2(T.mWin, mFrame_orig, 'same');
    else
        T.mFrame = mFrame_orig;
    end
    
elseif T.nCh == T.vChOrder(2)   % Blue channel (vessels)

    % Normalize vessel image
    T = NormalizeVesselFrame(T);
    T = CM_IOS_SANDBOX_CELINE(T);

    % Mask image (show non-vessels) - OPTIONAL
    %if isfield(T, 'cMasks')
    %    T.mFrame = T.mFrame .* T.cMasks{T.nCh};
    %    T.mFrame(T.mFrame == 0) = nan;
    %end
    
    % Get average blue signal amplitude
    %T.mGraphs(T.f, 2) = nansum(T.mFrame(:));
    T.mGraphs(T.f, 2) = NaN; % don't plot
    
    % Get intensity range
    T.vCLim = [-.28 .45]; % [min(T.mFrame(:))*1.2 max(T.mFrame(:))*1];

    % Compute a vessel and non-vessel masks from the 'blue' image
    if ~isfield(T, 'cMasks')
        hWait = waitbar(1/3, 'Calculating masks...');

        % Mask of vessels
        T.cMasks{T.vChOrder(2)} = GetRegionMask(T.mFrame, T.nPixRes, [], 10); % vessel
        waitbar(2/3, hWait)

        % Mask of non-vessel regions
        T.cMasks{T.vChOrder(1)} = GetRegionMask(T.mFrame.*-1, T.nPixRes, [], 20);

        close(hWait)
    end
    
else
    error('Unrecognized video channel.');
end

return
