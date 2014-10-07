function [prmts, ISIdata] = ISI_analysis(prmts)
% [prmts,ISIdata] = ISI_analysis(prmts)   
% Main function for IOS analysis.
% 
%

sColMap = 'gray';

% Set the no light file data, subtracting the DC offset
if prmts.useNoLight
    % Get the nolight data
    if ~exist(fullfile(prmts.Rnolight.path2dir, prmts.Rnolight.name),'file')
        choice=questdlg('Continue without NoLight file?','Continue?','Yes','No','Yes');
        if strcmpi(choice,'No')
            return;
        end
        drawnow
        Rnolight=[];
    else
        darkdata = ISI_read(prmts.Rnolight);
        % make darkfield average
        Rnolightframes=nan([darkdata.frameSizeYX darkdata.nFramesPerTrial]);
        Rnolighttrial=nan([darkdata.frameSizeYX darkdata.ntrials]);
        for i=1:darkdata.ntrials
            for j=1:darkdata.nFramesPerTrial
                Rnolightframes(:,:,j)=cell2mat(darkdata.frameStack(i,j))'; %transpose for correct orientation
            end
            Rnolighttrial(:,:,i) = mean(Rnolightframes,3);
        end
        Rnolight = mean(Rnolighttrial,3);
    end
else
    drawnow
    Rnolight=[];
end

% Cycle through files
for iFile = 1 : numel(prmts.filesQueue)

    % Get data from GUI, if already set
    hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
    tUserData = get(hGUI, 'UserData');
    ISIdata = [];
    if ~isempty(tUserData)
        if isfield(tUserData, 'frameStack')
            ISIdata = tUserData;
        end
    end
    
    % Load data if 1) this is the first time loading this trial, or 2) the
    % Load Data button was pressed (but not the Average Trials or Analyze
    % buttons)
    if isempty(ISIdata) || ...
            (prmts.filesQueue.DoLoad && ~prmts.filesQueue.DoTrialAverage && ~prmts.filesQueue.DoAnalyze)
        % read data
        ISIdata = ISI_read(prmts.filesQueue(iFile));
        
        % store data in GUI
        set(hGUI, 'UserData', ISIdata)
    end

    % Exit here if we should not average trials
    if ~prmts.filesQueue.DoTrialAverage
        return
    end
    
    % Compute trial-averaged frames, if:
    % 1) averages have not already been computed, or
    % 2) the Average Trials button was pressed
    if ~isfield(ISIdata, 'deltaSignal') || ...
            (prmts.filesQueue.DoTrialAverage && ~prmts.filesQueue.DoAnalyze)
        ISIdata = ISI_calc_dRR(ISIdata, prmts.filesQueue(iFile),Rnolight);
        set(hGUI, 'UserData', ISIdata) % store data in GUI
    end

    % Exit now if we should not run analysis
    if ~prmts.filesQueue.DoAnalyze
        return
    end
    
    ISIdata.climAll = prmts.filesQueue(iFile).climAll ./ 1000;
    
    % Get the average frame over stimulus interval
    if get(findobj('tag', 'chk_writesigframe'), 'value')
        % Plot the average signal frame
        signalFrame = ISI_getSignalFrame(ISIdata, prmts.filesQueue(iFile), 1);
    else   
        % Suppress plot
        signalFrame = ISI_getSignalFrame(ISIdata, prmts.filesQueue(iFile), 0);
    end
    
    % Plot and save all frames
    if (prmts.saveFigAll)
        ISI_plotAllFrames(ISIdata, ISIdata.climAll, prmts.filesQueue, sColMap);
        drawnow;
        savefilename = fullfile(prmts.filesQueue.path2dir, prmts.filesQueue.name);
        [p n e] = fileparts(savefilename);
        set(gcf, 'Name', savefilename, 'numberTitle', 'off')
        saveas(gcf, fullfile(p,[n '_AllFrames.png']), 'png');
    end
    
    ISIdata.stimInterval = prmts.filesQueue(iFile).stimInterval; %interval we averaged over
    ISIdata.signalFrame = signalFrame;
    
    % Plot trial-average signal of ROI across time
    if prmts.filesQueue(iFile).selectSignalROI
        ISIdata = ISI_selectSignalROI(ISIdata, prmts.filesQueue(iFile), sColMap);
        if isfield(ISIdata, 'analysisSignalROI')
            ISI_plotSignalByTime(ISIdata, prmts.filesQueue(iFile), sColMap);
        end
    end

    % Plot spatial intensity profile of user-selected line across time
    if prmts.filesQueue(iFile).selectSignalProfile
        ISIdata = ISI_selectSignalProfile(ISIdata, prmts.filesQueue(iFile), sColMap);
        if isfield(ISIdata, 'analysisSignalProfile')
            ISI_plotProfileByTime(ISIdata, prmts.filesQueue(iFile), sColMap);
        end
    end
    
    % Isolate barrel
    if prmts.filesQueue(iFile).runBarrelFinder
        ISIdata.vesselmask = ISI_createVesselMask(prmts);
        
        % Select contour for barrel
        ISIdata = ISI_isolateBarrel(ISIdata, signalFrame, prmts.filesQueue(iFile),prmts.saveFig);
        drawnow;
    end

    % Generate AVI movie from loaded data file
    if prmts.saveMovie
        ISIdata = ISI_writeMovie(ISIdata, prmts.filesQueue(iFile));
    end
    
    % Save mat file  
    drawnow; % update the figures before saving the mat file
    pause(0.01);
    savefilename = fullfile(prmts.filesQueue(iFile).path2dir,prmts.filesQueue(iFile).name);
    [p n e]=fileparts(savefilename);
    
    if prmts.saveToMat
        if ~strcmp(e,'.mat')
            matfilename = fullfile(p, [n '.mat']);
        end
        
        %cut out the frameStack field, since it's huge and unneccessary
        %mjp 2011.10.13
        ISIdata = rmfield(ISIdata, 'frameStack');
        
        save(matfilename, 'ISIdata');
    end

    % store data and analyzed parameters in GUI
    set(hGUI, 'UserData', ISIdata)

    
end %cycling files