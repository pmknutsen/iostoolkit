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
    
    % Bin image pixels
    T.mFrame = BinImagePixels(T.mFrame, 16);
    if ~all(size(T.mFrame) == size(T.cMasks{T.nCh})) % do once
        T.cMasks{T.nCh} = BinImagePixels(T.cMasks{T.nCh}, 16);
        T.cMasks{T.nCh} = round(T.cMasks{T.nCh});
    end
    
    % Set c-limits on displayed images
    T.vCLim = [-.04 .04]; % signal ratio

    % Mask image (show only non-vessel regions)
    if isfield(T, 'cMasks')
        T.mFrame = T.mFrame .* ~T.cMasks{T.nCh};
        T.mFrame(T.mFrame == 0) = nan;
    end
    
    % Get average signal amplitude across all pixels
    T.mGraphs(T.f, T.nCh) = nanmean(T.mFrame(:));
    
elseif T.nCh == T.vChOrder(2)   % Blue channel (vessels)

    % Normalize vessel image
    T = NormalizeVesselFrame(T);

    % Mask image (show non-vessels)
    if isfield(T, 'cMasks')
        T.mFrame = T.mFrame .* T.cMasks{T.nCh};
        T.mFrame(T.mFrame == 0) = nan;
    end
    
    % Get signal amplitude
    % CELINE: You can calculate vessel diameters here and insert into
    % mGraphs for plotting in the window, like below.
    T.mGraphs(T.f, T.nCh) = nansum(T.mFrame(:));

    % Get intensity range
    T.vCLim = [-.07 .11]; % [min(T.mFrame(:))*1.2 max(T.mFrame(:))*1];

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

% Smooth image by gaussian (sigma is set in GUI)
if isfield(T, 'mWin')
    T.mFrame = filter2(T.mWin, T.mFrame, 'same');
end

return
