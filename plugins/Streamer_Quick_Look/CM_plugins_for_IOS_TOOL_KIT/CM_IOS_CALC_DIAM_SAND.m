
function [T]=CM_IOS_CALC_DIAM_SAND (T,thickness)

To_init=~isfield(T,'AA_proc'); % if it doesn't exist you have to create it and start calculating the diameter
CST_MAX=length(T.AA);% AA contains the info necessary to measure the vessel diameter
%%=T.AA(1).thickness;
T.Dia_index=((((T.f-T.nSkipInitialFrames)+T.vChOrder(2))/2)-(T.vChOrder(2)-1));% CM_to fix with logic done


%for CST=1:1:CST_MAX % GOES THROUGH the different line profiles
for CST= 3 % PK Temp fix to plot just one vessel
    box_lim=T.AA(CST).box_lim;
    Xi_S=T.AA(CST).Xi_S;
    Yi_S=T.AA(CST).Yi_S;
    box_im=T.mFrame(box_lim(1):box_lim(2),box_lim(3):box_lim(4));%EXTRACTED THE BOX FROM FRAME
    %[FAT_LINE_AVG,~,~,~]=CM_ROTATE_LINE_ORIENTATION_SAND(Xi_S,Yi_S,box_im,thickness,0);
    FAT_LINE_AVG = improfilefast(box_im,Xi_S,Yi_S,thickness);
    
    %FAT_DET=-detrend(FAT_LINE_AVG);
    FAT_DET=-FAT_LINE_AVG;
    
    [dia_temp p1_vec_temp p2_vec_temp an_profile]= CM_calcFWHM(detrend(FAT_DET),1); % CELINE ADDED THE DETREND 20140803 TO COMPENSATE FOR UNSUFFICIENT IR IREGULAR RADIAL CORRECTION
    
    if (To_init) % preallocate the data and set the data
        T.nb_of_blue_frames=floor((T.nVideoNumFrames-T.nSkipInitialFrames)/2);
        T.first_blue_frame=T.Dia_index;
        T.AA_proc(CST).Vess_Cross_im=nan(size(FAT_DET,2),T.nb_of_blue_frames);
        T.AA_proc(CST).Vess_Cross_im_proc=nan(size(an_profile,2),T.nb_of_blue_frames);
        T.AA_proc(CST).dia=nan(T.nb_of_blue_frames,1);
        T.AA_proc(CST).p1_vec=nan(T.nb_of_blue_frames,1);
        T.AA_proc(CST).p2_vec=nan(T.nb_of_blue_frames,1);
        T.AA_proc(CST).thickness=thickness;
        T.AA_proc(CST).Dia_index=1;
        T.mGraphs(1:T.nVideoNumFrames,T.nNumColChans+CST) =nan;
        T=CM_MODIFY_GRAPH(T);% adds the display of the diameter
    end
    T.AA_proc(CST).Dia_index=T.Dia_index;
    T.AA_proc(CST).Vess_Cross_im(:,T.Dia_index)=FAT_DET; % to keep as a persistent variable when running frames
    T.AA_proc(CST).Vess_Cross_im_proc(:,T.Dia_index)=an_profile; % to keep as a persistent variable when running frames
    T.AA_proc(CST).dia(T.Dia_index)=dia_temp;
    T.AA_proc(CST).p1_vec(T.Dia_index)=p1_vec_temp;
    T.AA_proc(CST).p2_vec(T.Dia_index)=p2_vec_temp;
    T.mGraphs(T.f, T.nNumColChans+CST) =dia_temp;
    
    if T.bView && ishandle(T.hFig)
        nStart = (T.nSkipInitialFrames) + T.nCh; %% NOT CLEAR
        mYY = T.mGraphs(nStart:T.nNumColChans:end,T.nNumColChans+CST);
        mXX = (nStart:T.nNumColChans:size(T.mGraphs, 1)) ./ T.nVideoFPS;
        set(T.hGraphs(T.nNumColChans+CST), 'xdata', mXX, 'ydata', mYY);
        drawnow
    end
    
end

% CM_IOS_PLOT_LINES(AA_proc)
% CM_IOS_PLOT_ALL_PROFILES(AA_proc)
end
%
% function [T]=CM_IOS_PLOT_DIA(T)
%  CST_MAX=length(T.AA_proc);
% for CST=1:1:CST_MAX % GOES THROUGH the different line profiles
%     T.mGraphs(T.f, T.nNumColChans+CST) = nansum(T.mFrame(:));
%     T.AA_proc(CST).dia(T.Dia_index)=dia_temp;
%
% end
% end

function [width point1 point2 data] = CM_calcFWHM(data,smoothing,threshold)
data = double(data);
% smooth data, if appropriate
if nargin < 2 % smoothing not passed in, set to default (none)
    smoothing = 1;
end
if smoothing > 1
    data = conv(data,rectwin(smoothing) ./ smoothing);
end
baseline_to_sub=min(data(smoothing:(length(data)-smoothing)));
data =data-baseline_to_sub; % The minimum is calculated on the second to the avant dernier point
%% GIVE THE PROPER THRESHOLD RATIO
if nargin < 3
    threshold_ratio=2;% changed from 2 to 4 CM
    threshold = max(data)/threshold_ratio;
end
%%
aboveI = find(data > threshold);    % all the indices where the data is above half max
if isempty(aboveI)     % nothing was above threshold!
    width = 0;point1=0; point2=0;
    return
end

firstI = aboveI(1);                 % index of the first point above threshold
lastI = aboveI(end);                % index of the last point above threshold

if (firstI-1 < 1) || (lastI+1) > length(data) % interpolation would result in error, set width to zero and just return ...
    
    width = 0;point1=0;point2=0;
    return
end
% use linear interpolation to get a more accurate picture of where the max was
% find value difference between the point and the threshold value,
% and scale this by the difference between integer points ...
point1offset = (threshold-data(firstI-1)) / (data(firstI)-data(firstI-1));
point2offset = (threshold-data(lastI)) / (data(lastI+1)-data(lastI));
point1 = firstI-1 + point1offset;
point2 = lastI + point2offset;
width = point2-point1;
end

function [T]=CM_MODIFY_GRAPH(T)

%hAx=T.hAx;
% Plotting axes
figure (T.hFig);

for i = 1:T.nNumColChans
    handle_num=T.hAx(T.nNumColChans+i);
    set(handle_num,'position', [0 .2  1 .2]);% sets the placement of the axis for IOS
end
mCols=varycolor(length(T.AA));
%     mCols = [1 0 0; 0 0 1; 0 1 0; 1 0 1];
T.hAx(end+1) = axes('position', [0 0.05 1 0.14]);% adds one axis to put all the diameters
hold(T.hAx(end), 'on')
set(T.hAx(end), 'color','none')

for i = 1:length(T.AA)
    T.hGraphs(end+1) = plot(T.hAx(end), nan, '-', 'color', mCols(i, :));
    set(T.hAx(end), 'fontsize', 8, 'ycolor', mCols(i, :));
    set(T.hGraphs(end), 'xdata', 1:T.nVideoNumFrames, 'ydata', nan(T.nVideoNumFrames, 1));
    set(T.hAx(end), 'xlim', [(T.nSkipInitialFrames+1) T.nVideoNumFrames]);
end
%;

linkaxes(T.hAx(T.nNumColChans+1:end), 'x')
h = zoom(T.hFig);
setAxesZoomMotion(h, T.hAx(T.nNumColChans+1:end), 'horizontal')

return
end

function ColorSet=varycolor(NumberOfPlots)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.
%
%     Yellow and White colors were not used because of their poor
%     translation to presentations.
%
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
%
%         figure
%         hold on;
%
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

%Created by Daniel Helmick 8/12/2008

error(nargchk(1,1,nargin))%correct number of input arguements??
error(nargoutchk(0, 1, nargout))%correct number of output arguements??

%Take care of the anomolies
if NumberOfPlots<1
    ColorSet=[];
elseif NumberOfPlots==1
    ColorSet=[0 1 0];
elseif NumberOfPlots==2
    ColorSet=[0 1 0; 0 1 1];
elseif NumberOfPlots==3
    ColorSet=[0 1 0; 0 1 1; 0 0 1];
elseif NumberOfPlots==4
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1];
elseif NumberOfPlots==5
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0];
elseif NumberOfPlots==6
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0; 0 0 0];
    
else %default and where this function has an actual advantage
    
    %we have 5 segments to distribute the plots
    EachSec=floor(NumberOfPlots/5);
    
    %how many extra lines are there?
    ExtraPlots=mod(NumberOfPlots,5);
    
    %initialize our vector
    ColorSet=zeros(NumberOfPlots,3);
    
    %This is to deal with the extra plots that don't fit nicely into the
    %segments
    Adjust=zeros(1,5);
    for m=1:ExtraPlots
        Adjust(m)=1;
    end
    
    SecOne   =EachSec+Adjust(1);
    SecTwo   =EachSec+Adjust(2);
    SecThree =EachSec+Adjust(3);
    SecFour  =EachSec+Adjust(4);
    SecFive  =EachSec;
    
    for m=1:SecOne
        ColorSet(m,:)=[0 1 (m-1)/(SecOne-1)];
    end
    
    for m=1:SecTwo
        ColorSet(m+SecOne,:)=[0 (SecTwo-m)/(SecTwo) 1];
    end
    
    for m=1:SecThree
        ColorSet(m+SecOne+SecTwo,:)=[(m)/(SecThree) 0 1];
    end
    
    for m=1:SecFour
        ColorSet(m+SecOne+SecTwo+SecThree,:)=[1 0 (SecFour-m)/(SecFour)];
    end
    
    for m=1:SecFive
        ColorSet(m+SecOne+SecTwo+SecThree+SecFour,:)=[(SecFive-m)/(SecFive) 0 0];
    end
    
end
end

