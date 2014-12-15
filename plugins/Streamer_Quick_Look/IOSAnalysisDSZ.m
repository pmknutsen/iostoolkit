function T = IOSAnalysisDSZ(T)
% IOSAnalysis
%
% Analyse intrinsic signals (IOS) and blood vessel diameter
%
% Run this function via IOSToolkit plug-in Streamer_Quick_Look()
%

% Process channels
if T.nCh == T.vChOrder(1)   % Red channel
    %if ~isfield(T, 'cMasks'), return; end

    % Compute intrinsic signal (dR/R)
    mFrame = (T.mFrame - T.mBaseline(:, :, T.vChOrder(T.nCh))) ./ T.mBaseline(:, :, T.vChOrder(T.nCh));
    mFrame = mFrame - median(mFrame(:));

    persistent mSumFrame nCount
    if isempty(mSumFrame)
        mSumFrame = mFrame;
        nCount = 1;
    else
        mSumFrame = mSumFrame + mFrame;
        nCount = nCount + 1;
    end
    T.mFrame = mSumFrame ./ nCount;
    
    % Update time marker in Spiky (if open)
    global Spiky
    if ~isempty(Spiky)
        persistent nFrameRate nTimeBegin
        if isempty(nFrameRate)
            [FV, ~] = Spiky.main.GetStruct();
            iCamShut = strcmp('CameraShutter', {FV.tChannelDescriptions.sDescription});
            sCamShutCh = FV.tChannelDescriptions(iCamShut).sChannel;
            nFrameRate = 1 / median(diff(FV.tData.([sCamShutCh '_Up']))); % frames/s
            nTimeBegin = FV.tData.([sCamShutCh '_TimeBegin']);
        end
        nNow = (T.f / nFrameRate) + nTimeBegin;
        Spiky.main.SetTimeMarker(nNow)
    end
    
    % Bin image pixels
    %T.mFrame = BinImagePixels(T.mFrame, 8);
    
    % Average all frames
    
    if isfield(T, 'cMasks')
        if ~all(size(T.mFrame) == size(T.cMasks{T.nCh})) % do once
            T.cMasks{T.nCh} = BinImagePixels(T.cMasks{T.nCh}, 16);
            T.cMasks{T.nCh} = round(T.cMasks{T.nCh});
        end
    else
    end
    
    % Set c-limits on displayed images
    T.vCLim = [-.005 .005]; % signal ratio

    % Mask image (show only non-vessel regions)
    if isfield(T, 'cMasks')
        T.mFrame = T.mFrame .* ~T.cMasks{T.nCh};
        T.mFrame(T.mFrame == 0) = nan;
    end
    
    % Get average signal amplitude across all pixels
    T.mGraphs(T.f, 1) = nanmean(T.mFrame(:));

    
elseif T.nCh == T.vChOrder(2)   % Blue channel (vessels)
    % There is no blue channel in the DSZ dataset so we should not be here
    keyboard
else
    error('Unrecognized video channel.');
end

% Smooth image by gaussian (sigma is set in GUI)
if isfield(T, 'mWin')
    T.mFrame = filter2(T.mWin, T.mFrame, 'same');
end

return
