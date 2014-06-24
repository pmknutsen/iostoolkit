function T = SandboxExample(T)
% Sandbox function example for Streamer_Quick_Look()
%
%

% Analysis mode
%	-1	IOS-type derivatives
%	0	Raw
%	1	Derivatives of consequtive frames
bMakeDiff = 0;

% Get intensity range
mBaseline = T.mBaseline(:, :, T.nCh);
T.vCLim = [min(mBaseline(:))*1.2 max(mBaseline(:))*1];

if bMakeDiff == 1
    T.vCLim = [-.025 .025];
elseif bMakeDiff < 0
    T.vCLim = [-.025 .025];
end

% Compute signal
if bMakeDiff == 1
    % Temporal derivative
    
    % Subtract previous frame
    nPrevIndx = T.vCircBufferIndx(T.nCh, :) == T.f - T.nNumColChans;
    if ~any(nPrevIndx), return; end % we're at 1st frame
    
    T.mFrame = T.mFrame - squeeze(T.mCircBuffer(T.nCh, :, :, nPrevIndx));
    
    % Divide by baseline
    T.mFrame = T.mFrame ./ T.mBaseline(:, :, T.nCh);
    
elseif bMakeDiff == -1
    % IOS-type analysis
    
    % Compute dR/R
    T.mFrame = (T.mFrame - T.mBaseline(:, :, T.nCh)) ./ T.mBaseline(:, :, T.nCh);
end

% TODO Filter
%mFrame = filter2(mWin, mFrame, 'same');


return