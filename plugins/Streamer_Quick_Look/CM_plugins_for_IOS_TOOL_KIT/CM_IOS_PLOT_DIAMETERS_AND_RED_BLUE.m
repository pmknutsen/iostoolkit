% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
% Xi_N=Xi_S+box_lim(3); Yi_N=Yi_N+box_lim(1);


function CM_IOS_PLOT_DIAMETERS_AND_RED_BLUE (T)

counter_max=size(T.mGraphs,2);
ColorSet = varycolor(counter_max);
if (counter_max<1)
    return
else
    figure
    hAx(1)= subplot (3,1,1);
    red_CH=T.mGraphs(23:2:end,1);
    blue_CH=T.mGraphs(22:2:end,2);
    
    length_min = min(length(blue_CH),length(red_CH));
    
    time_axis=(1:length_min)/(T.nVideoFPS/2)-1/(T.nVideoFPS/2);
    
    red_CH = red_CH(1:length_min);
    blue_CH = blue_CH(1:length_min);
    
    plot ( hAx(1),time_axis,red_CH);
    hAx(2)=     subplot (3,1,2);
    
    plot ( hAx(2),time_axis,blue_CH);
    hAx(3)= subplot (3,1,3);
    
    for counter=1:1:counter_max-2
        Y_DIA= T.mGraphs(22:2:end,counter+2);% +2 because the first 2 are taken by red and blue CH
        Y_DIA = Y_DIA(1:length_min);
        X_DIA=time_axis;
        
        plot ( hAx(3),X_DIA,Y_DIA,'Color',ColorSet(counter,:));
        hold on
    end
    
    
end

linkaxes(hAx(1:end),'x')

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







