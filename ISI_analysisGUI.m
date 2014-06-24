function varargout = ISI_analysisGUI(varargin)
% ISI_analysisGUI Run the IOSTOOLKIT GUI
%
% Usage:
%   < IOS_analysisGUI >
%   Runs the IOS Toolkit GUI
%
%   Alternative syntax for running subroutines:
%     global IOS
%     IOS.SUB(arg1, ..., argn)
%   where SUB is the subroutine to run.
%

% IOS Toolkit - Matlab toolkit for analysis of IOS data
% Copyright (C) 2014 Application Authors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 

% Last Modified by GUIDE v2.5 22-Jan-2014 11:12:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ISI_analysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @ISI_analysisGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to execute just before ISI_analysisGUI is made visible
function ISI_analysisGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok

handles.pathstr = pwd;     % initial file directory
set(handles.path,'string',handles.pathstr);

set(hObject,'Tag', 'ISIanalysisGUI_fig')

% Choose default command line output for ISI_analysisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add plugin folder to path
addpath([fileparts(mfilename('fullpath')) filesep 'plugins']);

% Print copyright message
clc
disp('IOS Toolkit - Matlab toolkit for analysis of IOS data')
disp('Copyright (C) 2014 by Application Authors')
disp('This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are')
disp(sprintf('are welcome to redistribute it under certain conditions; see LICENSE for details.\n'))

% Get list of all internal sub-routines
disp('Initializing function handles...')
csStr = mlintmex('-calls', which(mfilename));
[~,~,~,~,subs] = regexp(csStr, '[S]\d* \d+ \d+ (\w+)\n');
cSubs = [subs{:}]';

% Generate function handles for all sub-routines
global ISI;
IOS = struct([]);
for i = 1:length(cSubs)
    ISI(1).(cSubs{i}) = eval(['@' cSubs{i}]);
end

disp(sprintf('\nNote on use:\nYou can run IOS Toolkit routines directly with syntax ISI.SUB() (where ISI is global).\n'))

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs from this function are returned to the command line.
function varargout = ISI_analysisGUI_OutputFcn(hObject, eventdata, handles) %#ok
varargout{1} = handles.output;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file/path info
function path_Callback(hObject, eventdata, handles) %#ok
curpath=get(handles.path,'string');
if ~isempty(curpath) && ~exist(curpath,'dir') % directory doesnt exist
    warndlg(sprintf('The folder: \n%s\ndoes not exist.',curpath));
    set(handles.path,'string',pwd);
    handles.pathstr=pwd;
else
    handles.pathstr=curpath;
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the root path of this application
function sPath = GetRootPath()
sPath = which(mfilename);
nIndx = findstr([mfilename '.m'], sPath);
sPath = sPath(1:nIndx-1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in getPath.
function getPath_Callback(hObject, eventdata, handles)  %#ok
oldpath=handles.pathstr;
if isempty(oldpath) || ~isstr(oldpath), oldpath=pwd;  end
curpath=uigetdir(oldpath);
if curpath==0, curpath=pwd; end %cancelled out
set(handles.path,'string',curpath);
handles.pathstr=curpath;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_filename_Callback(hObject, eventdata, handles)  %#ok
datafile=get(handles.data_filename,'string');
if ~exist(fullfile(handles.pathstr,datafile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',datafile));
    set(handles.data_filename,'string','');
else
    name=regexp(handles.params.data_filename,'(.+)_\w{6}.dat','tokens');
    set(handles.whiskername,'string',name{1});
end
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function data_filename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in getdatafile.
function getdatafile_Callback(hObject, eventdata, handles)

if isempty(eventdata)
    curdir = pwd;
    cd(handles.pathstr); % move to path directory
    [datafile, path2file] = uigetfile( ...
        {'*.dat','LabView files (*.dat)'; ...
        '*.bin','Streamer files (*.bin)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file');
    
    if datafile == 0 % user pressed cancel
        return;
    end
    
    cd(curdir); % return to original directory
    set(handles.data_filename, 'string', datafile);
    
    % Update path to selected file
    set(handles.path,'string', path2file);
    handles.pathstr =path2file;
else
    % Use passed file name passed as function argument
    datafile = eventdata;
    set(handles.data_filename, 'string', datafile);
end

% Try to get whisker name
sWhiskerID = regexp(datafile, '(.+)_\w{6}.dat', 'tokens');

if ~isempty(sWhiskerID)
    set(handles.whiskername,'string', sWhiskerID{1});
else
    warndlg('Could not find whisker name. Please enter it manually.', 'IOSToolkit');
    set(handles.whiskername, 'string', 'Unknown');
end

guidata(hObject, handles);

hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
set(hGUI, 'UserData', [])

% Since a file was selected manually, uncheck the 'Process all files' checkbox
set(handles.process_all_files, 'value', 0)
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vessel_filename_Callback(hObject, eventdata, handles)  %#ok
vesselfile=get(handles.vessel_filename,'string');
if ~exist(fullfile(handles.pathstr,vesselfile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',vesselfile));
    set(handles.vessel_filename,'string','');
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function vessel_filename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in getvessel.
function getvessel_Callback(hObject, eventdata, handles) %#ok
curdir=pwd;
cd(handles.pathstr); %load to path directory
oldfile=get(handles.vessel_filename,'string');
vesselfile=uigetfile({'*.png'});
if vesselfile==0, vesselfile=oldfile; end %cancelled out
cd(curdir);
set(handles.vessel_filename,'string',vesselfile);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function whiskername_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nolight_filename_Callback(hObject, eventdata, handles) %#ok
nolightfile=get(handles.nolight_filename,'string');
if ~exist(fullfile(handles.pathstr,nolightfile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',nolightfile));
    set(handles.nolight_filename,'string','');
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function nolight_filename_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in getnolight.
function getnolight_Callback(hObject, eventdata, handles) %#ok
nolightfile=uigetfile('*.dat');
set(handles.nolight_filename,'string',nolightfile);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function preStimDurSec_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function stimDurSec_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function contour_min_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function contour_max_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function contour_step_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function stiminterval_lo_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskthresh_Callback(hObject, eventdata, handles) %#ok
value=str2double(get(hObject,'string'));
if isnan(value) || ( value<0 || value>1 )
    errordlg('Mask threshold value must be between 0 and 1','Incorrect mask threshold');
    set(hObject,'string','0');
end
set(handles.slider_maskThresh, 'value', value)
% mask preview window is open, then update it
hFig = findobj('tag', 'ISI_maskPreview');
if ~isempty(hFig)
    btn_viewmask_Callback(hObject, eventdata, handles)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function maskthresh_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Executes during object creation, after setting all properties.
function stiminterval_hi_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a structure of parameters used for analysis
% Most parameters are read from the GUI handles/objects
function tParams = GetParams(hObject, eventdata, handles)

% General parameters
tParams.useVesselMask   = get(handles.chk_vesselmask, 'value');
tParams.sColMap         = 'gray'; % default colormap
tParams.sAppName        = 'IOS Toolkit';

% Manual mask
tParams.filesQueue.useManualMask    = get(handles.chk_manualmask, 'value');
tParams.filesQueue.manualMask       = get(handles.btn_setmanualmask, 'userdata');

% File parameters
tParams.filesQueue.path2dir   = handles.pathstr;
tParams.filesQueue.name       = get(handles.data_filename, 'string');
tParams.filesQueue.refImage   = get(handles.vessel_filename, 'string');
tParams.filesQueue.precision  = 'int16';

% Range of trials to analyse
if ~isempty( [get(handles.trials2use_Start, 'string')] ) && ~isempty( [get(handles.trials2use_End, 'string')] )
    tParams.filesQueue.Trials2Use = str2double(get(handles.trials2use_Start, 'string')):str2double(get(handles.trials2use_End, 'string'));
else
    tParams.filesQueue.Trials2Use = [];
end

% Trials to exclude from analysis
tParams.filesQueue.Trials2Exclude = [];
if ~isempty(get(handles.trials2exclude, 'string'))
    tParams.filesQueue.Trials2Exclude = str2num(get(handles.trials2exclude, 'string'));
end

% Filename
tParams.filesQueue.Whisker    = get(handles.whiskername, 'string');

% Some defaults - Not sure what these do.... TODO: Document!
tParams.filesQueue.minFaces   = 50; % ignores contours that would be too small
tParams.filesQueue.maxFaces   = 1000; % ignores contours that are too big

% Contour lines
contour_min     = str2double(get(handles.contour_min, 'string'));
contour_step    = str2double(get(handles.contour_step, 'string'));
contour_max     = str2double(get(handles.contour_max, 'string'));
tParams.filesQueue.ContourLineVals    = contour_min:contour_step:contour_max;

% Stimulus interval
nInterval_Lo = str2double(get(handles.stiminterval_lo, 'string'));
nInterval_Hi = str2double(get(handles.stiminterval_hi, 'string'));
tParams.filesQueue.stimInterval = [nInterval_Lo nInterval_Hi];

% Smoothing (sigma)
tParams.filesQueue.smoothSigma  = str2double(get(handles.smooth_sigma, 'string'));

% Pre- and post- stimulus durations
tParams.filesQueue.preStimDurSec  = str2double(get(handles.preStimDurSec,'string'));
tParams.filesQueue.stimDurSec     = str2double(get(handles.stimDurSec,'string'));

% Clim on all displayed images
climinterval_lo     = str2double(get(handles.climinterval_lo, 'string'));
climinterval_hi     = str2double(get(handles.climinterval_hi, 'string'));
tParams.filesQueue.climAll = [climinterval_lo climinterval_hi];

% Image binning
imgbin_x = str2double(get(handles.imgbin_x, 'string'));
imgbin_y = str2double(get(handles.imgbin_y, 'string'));
tParams.filesQueue.imgBin = [imgbin_x imgbin_y];

tParams.filesQueue.maskthresh     = str2double(get(handles.maskthresh, 'string'));

% If we are using the nolight normalization, set those parameters
tParams.useNoLight       = get(handles.chk_nolight,'value');
tParams.Rnolight         = tParams.filesQueue;
tParams.Rnolight.name    = get(handles.nolight_filename,'string');
tParams.Rnolight.Whisker ='none';

% Trial averaging options
tParams.filesQueue.useMedianCorrection = get(handles.chk_mediancorrection, 'value');
tParams.filesQueue.useMotionCorrection = get(handles.chk_motioncorrection, 'value');

% Remove axes in the main GUI
delete(findobj(handles.ISIanalysisGUI_fig, 'type', 'axes'))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the actions that are executed when the Load, Average or Analyze
% buttons are pressed, and the processing routines that will be executed.
function tParams = GetActions(hObject, eventdata, handles, tParams)

% By default, Load, Average and Analyse
tParams.filesQueue.DoLoad = true;
tParams.filesQueue.DoTrialAverage = true;
tParams.filesQueue.DoAnalyze = true;

% If this function was called via the Info button, instruct next
% scripts to only display file info, but not load or analyze data
if hObject == handles.btn_explore
    tParams.filesQueue.DoLoad = false;
    tParams.filesQueue.DoTrialAverage = false;
    tParams.filesQueue.DoAnalyze = false;
end

% If this function was called via the Load Data button, instruct next
% scripts to only load data
if hObject == handles.btn_loaddata
    tParams.filesQueue.DoLoad = true;
    tParams.filesQueue.DoTrialAverage = false;
    tParams.filesQueue.DoAnalyze = false;
end

% If this function was called via the Average Trial button, instruct next
% scripts to load data (if not already done so) and then average trials.
if hObject == handles.btn_averagetrials
    tParams.filesQueue.DoLoad = true;
    tParams.filesQueue.DoTrialAverage = true;
    tParams.filesQueue.DoAnalyze = false;
end

% Processing routines to run

% 'Save results to MAT file'
tParams.saveToMat       = get(handles.chk_savemat, 'value');

% 'View average frame' TODO: What does it do??
tParams.saveFig         = get(handles.chk_writesigframe, 'value');

% 'View all frames'
tParams.saveFigAll      = get(handles.chk_writesigframeall, 'value');

% 'Locate evoked region'
tParams.filesQueue.runBarrelFinder     = get(handles.chk_runBarrelFinder, 'value');

% 'Save AVI movie'
tParams.saveMovie       = get(handles.chk_savemovie, 'value');

% 'ROI vs time'
tParams.filesQueue.selectSignalROI     = get(handles.chk_selectSignalROI, 'value');

% ' Profile vs time'
tParams.filesQueue.selectSignalProfile = get(handles.chk_selectSignalProfile, 'value');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Analysis.
% Scrape the form, create the param struct and pass it to the ISI_analysis
% function.
%
% We need to stack the structs, because ISI_analysis expects a struct
% containing an array of file parameter structs. We are only running one
% file at a time, so we just have the filesQueue struct within setParams
function btn_run_Callback(hObject, eventdata, handles)

% Iterate over files
if get(handles.process_all_files, 'value')
    % Get list of all .dat files in directory
    sFileList = dir(fullfile(handles.pathstr,'*.dat'));
    bForceAllOption = 1;
else
    sFileList = get(handles.data_filename, 'string');
    bForceAllOption = 0;
end

if isempty(sFileList)
    warndlg('No file has been selected.', 'ISI'); return;
end

if isstruct(sFileList)
    nLoop = length(sFileList);
else
    nLoop = 1;
end

% Iterate over files
% Usually, only the selected file is processed unless the 'Process all'
% option is checked in the GUI
for nFi = 1:nLoop
    set(handles.process_all_files, 'value', bForceAllOption);
    if bForceAllOption
        % Update file name in GUI
        getdatafile_Callback(hObject, sFileList(nFi).name, handles)
        % Create waitbar
        if ~exist('hWait', 'var')
            hWait = waitbar(nFi/length(sFileList), 'Processing all .dat files. Please wait...');
        end
        if ishandle(hWait)
            waitbar(nFi/length(sFileList), hWait);
        else break; end
    end
    
    % Save all GUI settings and associate _settings.mat file current .dat file
    SaveGUIState(handles)
    
    % Seed the parameters structure with default and GUI values
    setParams = GetParams(hObject, eventdata, handles);

    % Decide actions (Load, Average and/or Analyse)
    setParams = GetActions(hObject, eventdata, handles, setParams);
    
    % Run analysis
    try
        ISI_analysis(setParams); %ISI_analysis returns setParams, but we aren't using it
    catch e
        errordlg({['Error using ==> ', e.stack(1).name, ' at ' num2str(e.stack(1).line)],...
            '',e.message,''},...
            'Oops.');
        break;
    end
    
end

% Close waitbar
if exist('hWait', 'var')
    if ishandle(hWait)
        close(hWait)
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_explore.
function btn_explore_Callback(hObject, eventdata, handles) %#ok
btn_run_Callback(hObject, eventdata, handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_viewmask.
function btn_viewmask_Callback(hObject, eventdata, handles) %#ok
if isempty(get(handles.vessel_filename,'string'))
    warndlg('No vessel filename given');
else
    ISI_viewMask(fullfile(handles.pathstr, get(handles.vessel_filename,'string')),...
        str2double(get(handles.maskthresh,'string')));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function trials2use_Start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function trials2use_End_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function framebinsize_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function framebinsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_loaddata.
function btn_loaddata_Callback(hObject, eventdata, handles)
btn_run_Callback(hObject, eventdata, handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function climinterval_lo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function climinterval_hi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function imgbin_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function imgbin_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function smooth_sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_apply_clim.
function btn_apply_clim_Callback(hObject, eventdata, handles)
% Apply new colormap limits (clim) to all open figures
hFig = findobj(0, 'type', 'figure');
nCLimLo = str2double(get(handles.climinterval_lo, 'string')); % permil
nCLimHi = str2double(get(handles.climinterval_hi, 'string')); % permil
nCLimLo = nCLimLo / 1000; % percent
nCLimHi = nCLimHi / 1000; % percent
% Iterate over figures
for hf = hFig(:)
    set(findobj(hf, 'type', 'axes'), 'clim', [nCLimLo nCLimHi])
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_averagetrials.
% Load data, if not already done, and compute trial averages
function btn_averagetrials_Callback(hObject, eventdata, handles)
btn_run_Callback(hObject, eventdata, handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in process_all_files.
function process_all_files_Callback(hObject, eventdata, handles)
if get(hObject, 'value')
    set(handles.data_filename, 'enable', 'off')
else
    set(handles.data_filename, 'enable', 'on')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run plugin from pop-up menu
function popup_plugins_Callback(hObject, eventdata, handles)
% Run selected plugin
sPlugin = get(hObject, 'string');
sPluginId = get(hObject,'value');
% Get plugin path
sPwd = pwd;
sPath = mfilename('fullpath');
vIndx = findstr(filesep, sPath);
sPath = [sPath(1:vIndx(end)) 'plugins' filesep sPlugin{sPluginId} filesep];

% Add path
addpath(sPath)

% Remove axes that were drawn into the GUI window
delete(findobj(gcf, 'type', 'axes'))

% Run plugin
eval(sprintf('%s(hObject, eventdata, handles);',sPlugin{sPluginId}))

% Remove axes that were drawn into the GUI window
delete(findobj(gcf, 'type', 'axes'))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate pop-up menu with list of plug-ins in /plugins/ folder
function popup_plugins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Get list of plugins
sPath = mfilename('fullpath');
vIndx = findstr(filesep, sPath);
sPath = [sPath(1:vIndx(end)) 'plugins' filesep];
sPwd = pwd;
cd(sPath)
tDir = dir;
cPlugins = {};
entryIsDir = [tDir.isdir]; % keep only directories
tDir = tDir(entryIsDir);
for i = 3:length(tDir)
    if tDir(i).isdir
        cPlugins{i-2} = tDir(i).name;
    end
end
cd(sPwd)

% Update popup menu in GUI
set(hObject, 'string', cPlugins)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_setmanualmask.
function btn_setmanualmask_Callback(hObject, eventdata, handles) %#ok
tMask = ISI_setManualMask(fullfile(handles.pathstr, get(handles.vessel_filename,'string')), hObject);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on slider movement.
function slider_maskThresh_Callback(hObject, eventdata, handles)
set(handles.maskthresh, 'string', num2str(get(hObject, 'value')))
% mask preview window is open, then update it
hFig = findobj('tag', 'ISI_maskPreview');
if ~isempty(hFig)
    btn_viewmask_Callback(hObject, eventdata, handles)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function slider_maskThresh_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
function trials2exclude_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in btn_activitymonitor.
function btn_activitymonitor_Callback(hObject, eventdata, handles)

% Get data from GUI, if already set
hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
tUserData = get(hGUI, 'UserData');
ISIdata = [];
if ~isempty(tUserData)
    if isfield(tUserData, 'frameStack')
        ISIdata = tUserData;
    end
end
if isempty(ISIdata), warndlg('You must first load a .dat file', 'ISI Analysis'); return, end

% Produce an average of all frames in a single movie
% use this to look for motion artefacts within single trials
nTrials = size(ISIdata.frameStack, 1);

hFig = figure;
colormap gray

nSigma = str2double(get(handles.smooth_sigma, 'string'));
if ~isnan(nSigma)
    mWin = fspecial('gaussian', nSigma*3, nSigma);
else
    mWin = NaN;
end

nrow = floor(sqrt(nTrials));
ncol = ceil(nTrials / nrow);

t = 1;
for j = 1:nrow
    for i = 1:ncol
        if t > nTrials, break, end
        cFrames = ISIdata.frameStack(t, :);
        mFrames = reshape(cell2mat(cFrames), [512 512 size(cFrames, 2)]);
        mDiffFrames = diff(mFrames, 1, 3);
                
        mImg = mean(mDiffFrames, 3)';
        
        mCtrlFrame1 = mean(mFrames(:,:,1:(size(mFrames,3)/2)), 3);
        mCtrlFrame2 = mean(mFrames(:,:,(size(mFrames,3)/2):end), 3);
        mImg = (mCtrlFrame1 - mCtrlFrame2)';
        
        axes('position', [(i*(1/ncol))-(1/ncol) 1-(j*(1/nrow))  1/ncol 1/nrow]) % modified Per Jan 10th 2012
        
        % Smooth frame
        if ~isnan(mWin)
            mImg = single(filter2(mWin, mImg, 'same'));
        end
        
        if ~ishandle(hFig), return, end
        figure(hFig)
        imagesc(mImg)
        
        axis image off
        hTxt = text(15, 30, num2str(t), ...
            'color', 'w', 'backgroundcolor', 'k'); % modified Per Jan 10th 2012
        t = t + 1;
        drawnow
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save state of GUI
function SaveGUIState(handles)
vChild = findobj(handles.ISIanalysisGUI_fig);
cExceptions = {'ISIanalysisGUI_fig'}; % tags of handles that should not be saved
tState = struct([]);
for c = 1:length(vChild)
    sTag = get(vChild(c), 'tag');
    if any(strcmp(sTag, cExceptions)), continue, end
    if isempty(sTag), continue, end
    if isprop(vChild(c), 'value')
        tState(1).(sTag).value = get(vChild(c), 'value');
    end
    if isprop(vChild(c), 'string')
        tState(1).(sTag).string = get(vChild(c), 'string');
    end
    if isprop(vChild(c), 'userdata')
        tState(1).(sTag).userdata = get(vChild(c), 'userdata');
    end
end
sSaveAs = strrep(fullfile(handles.pathstr, get(handles.data_filename,'string')), '.dat', '_UIValues.mat');
try
    save(sSaveAs, 'tState')
catch
    disp([mfilename ': Failed to save GUI state.'])
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restore state of GUI form last saved settings
function LoadGUIState(handles)
vChild = findobj(handles.ISIanalysisGUI_fig);
sLoadFile = strrep(fullfile(handles.pathstr, get(handles.data_filename,'string')), '.dat', '_UIValues.mat');
if ~exist(sLoadFile, 'file')
    return
end
try % may fail is file is corrupt, as has happened...
    load(sLoadFile)
catch
    return
end
if ~exist('tState', 'var'), return, end
csFieldnames = fieldnames(tState);
for t = 1:length(csFieldnames)
    if isfield(handles, csFieldnames{t})
        if ishandle(handles.(csFieldnames{t}))
            if isfield(tState.(csFieldnames{t}), 'value')
                if ~isempty(tState.(csFieldnames{t}).value)
                    set(handles.(csFieldnames{t}), 'value', tState.(csFieldnames{t}).value)
                end
            end
            if isfield(tState.(csFieldnames{t}), 'string')
                if ~isempty(tState.(csFieldnames{t}).string)
                    set(handles.(csFieldnames{t}), 'string', tState.(csFieldnames{t}).string)
                end
            end
            if isfield(tState.(csFieldnames{t}), 'userdata')
                if ~isempty(tState.(csFieldnames{t}).userdata)
                    set(handles.(csFieldnames{t}), 'userdata', tState.(csFieldnames{t}).userdata)
                end
            end
        end
    end
end
guidata(handles.ISIanalysisGUI_fig, handles);
return
