function iostoolkit
% Alias function for ISI_analysisGUI()
%
% Contributors:
%
%   Per M Knutsen, UCSD
%   Michael Pesavento, UCSD
%   Pablo Blinder, UCSD
%   Patrick Drew, UCSD
%

ISI_analysisGUI();

% Add all sub-folders recursively to path\
% This will include the /plugin/ folder and all its sub-folders to the
% path.
disp('Initializing paths...')
sPath = which('iostoolkit');
sPath = sPath(1:end-12);
addpath(genpath(sPath), '-end');

return