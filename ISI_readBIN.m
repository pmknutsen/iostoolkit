function [ISIdata] = ISI_readBIN(prmts)
% ISI_READBIN  Load ISI data from .BIN files
% Load trials from binary file generated by Streamer application.
% 
% Notes on use:
%  If function reports an Out of Memory error, try enabling pixel binning
%  to reduce memory usage when the function runs.
% 
%

if ~prmts.DoLoad, return; end

% Open file
ISIdata = [];
path2file = fullfile(prmts.path2dir, prmts.name);

% If filenname contains a wildcard, then concatenate all files in the
% current directory assuming each file is a single trial
tFiles = dir(path2file);
if isempty(tFiles)
    warning('Cannot locate .bin files');
end

% TODO
%   Dont load excluded trials (if set)
%   Load only trials in the included range (if set)

%if ~isempty(prmts.Trials2Use)
%    trial_range = sprintf('[%d %d]', prmts.Trials2Use(1), prmts.Trials2Use(end));
%else
%    trial_range = sprintf('[0 %d]', ntrials);
%end

% Initialize
ISIdata.frameStack = {};

% Load trials
hWait = waitbar(0,'Loading trials...');
for t = 1:length(tFiles)
    % Update waitbar
    if ishandle(hWait)
        waitbar(t/length(tFiles), hWait, sprintf('Loading trial %d of %d', t, length(tFiles)))
    else
        error('File read aborted.');
    end
    
    % Load all frames of current file
    sPath = fullfile(fileparts(path2file), tFiles(t).name);
    [mFrames, tInfo] = ISI_readStreamer(sPath);
    
    % Cast frames to single and concatenate trials
    for f = 1:size(mFrames, 3)
        if prod(prmts.imgBin) > 1
            ISIdata.frameStack{t, f} = single(downsamp2d(mFrames(:,:, f), prmts.imgBin)); % bin
        else
            ISIdata.frameStack{t, f} = single(mFrames(:,:, f));
        end
    end
end
close(hWait)

% FIX BELOW

% Store parameters relevant for further analysis
ISIdata.ntrials = t;
ISIdata.frame_rate = tInfo.nFPS;
ISIdata.nFramesPerTrial = tInfo.nNumFrames;
ISIdata.bin_duration = 1; % frames that were binned
ISIdata.nsec = tInfo.nNumFrames * ISIdata.bin_duration; % trial duration, s
ISIdata.nPreStimFrames = (ISIdata.frame_rate ./ ISIdata.bin_duration) * prmts.preStimDurSec;
ISIdata.nStimFrames = (ISIdata.frame_rate ./ ISIdata.bin_duration) * prmts.stimDurSec;
ISIdata.nPostStimFrames = ISIdata.nFramesPerTrial - (ISIdata.nPreStimFrames + ISIdata.nStimFrames);
ISIdata.frameSizeYX = size(ISIdata.frameStack{1, 1});

return


% OLD .DAT CODE BELOW

% Print file into to prompt
fprintf('\nFilename:\t\t%s\nTrials:\t\t\t%d\nFrame Rate:\t\t%d frames/s\nBin Duration:\t%d frames / %.2f s\nFrame Size:\t\t%dx%d px\nFrames/trial:\t%d\nBit Depth:\t\t%d\nTrial Duration:\t%d s\nTrials Used:\t%s\n',...
    prmts.name, ntrials, frame_rate, bin_duration, bin_duration_sec, size_x,size_y,nFramesPerTrial,bit_depth1,nsec,trial_range);

% Check that frames/trial estimate is consistent with header information
% Note: These numbers WILL deviate if the trial duration is non-integer.
%       In these cases, we will issue a warning to prompt.
if nFramesPerTrial ~= [nsec*(frame_rate/bin_duration)]
    disp(sprintf('Warning: Frames/trial estimated from filesize (%.0f) is different from that retrieved from file headers (%d)', nFramesPerTrial, nsec*(frame_rate/bin_duration)))
end

h2fig = findobj('Tag', 'ISIanalysisGUI_fig');%ensure function can still be run w/o GUI
if ~isempty(h2fig);centerfig(hWait, h2fig);end
drawnow;

