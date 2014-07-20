function CM_IOS_PLOT_ALL_PROFILES_SAND(AA_proc)
CST_MAX=length(AA_proc);
%CST_MAX=80;

ColorSet = varycolor(CST_MAX);
figure (250)

clf
subplot (3,1,1)
for CST=1:1:CST_MAX
    if (CST==1)
        set(gca, 'ColorOrder', ColorSet);
        hold all;
    end
    plot  ((AA_proc(CST).dia))
    legend250(CST)={num2str(CST)};
  
end
hold off
hleg250 = legend(legend250,'location','NorthEastOutside');
  title(['DIAMETER FROM ' num2str(1) ' TO LINE ' num2str(CST_MAX) ])


for CST=1:1:CST_MAX
    figure (112+CST)
    clf
    

    subplot(3,1,1)
  
    imagesc (AA_proc(CST).Vess_Cross_im)
    imagesc (AA_proc(CST).Vess_Cross_im_proc)
      title(['LINE NUMBER ' num2str(CST)])
    
    colormap (jet)
    hold on
    plot(AA_proc(CST).p1_vec)
    plot(AA_proc(CST).p2_vec)
    plot ((AA_proc(CST).p2_vec+AA_proc(CST).p1_vec)./2)
    hold off
    
    subplot(3,1,2)
    plot(AA_proc(CST).dia)
    subplot(3,1,3)
    plot(AA_proc(CST).p1_vec)
    hold on
    
    plot(AA_proc(CST).p2_vec)
    hold on
    
    plot ((AA_proc(CST).p2_vec+AA_proc(CST).p1_vec)./2)
    hold off
end
figure (250)
% subplot (3,1,2)
% title(['IOS RED AND BLUE' ])
% plot (AA_proc(1).IOS_RED,'r')
% subplot (3,1,3)
% plot (AA_proc(1).IOS_BLUE,'b')
%%
% count=1
% AA_proc=AA_proc_B;
%  IOS_RED=AA_proc(count,1).IOS_RED ;
%  IOS_BLUE=AA_proc(count,1).IOS_BLUE ;
%  DIA_ART=AA_proc(1,3).dia; 
% 
% IOS_RED_AMP=max(IOS_RED(:))-min(IOS_RED(:));
%   IOS_RED_NORM=IOS_RED-min(IOS_RED(:));
%   IOS_RED_NORM=(IOS_RED_NORM/IOS_RED_AMP)*100;
%  
%   
%   DIA=DIA_ART;
%   DIA_AMP=max(DIA_ART(:))-min(DIA_ART(:));
%   DIA=DIA-min(DIA(:));
%   DIA=(DIA/DIA_AMP)*100;
%   
%    IOS_BLUE_NORM=IOS_BLUE;
%   IOS_BLUE_AMP=max(IOS_BLUE(:))-min(IOS_BLUE(:));
%   IOS_BLUE_NORM=IOS_BLUE_NORM-min(IOS_BLUE_NORM(:));
%   IOS_BLUE_NORM=(IOS_BLUE_NORM/IOS_BLUE_AMP)*100;
%   
%   figure
%   
%    plot(DIA,'g')
%    hold on
%   plot(IOS_RED_NORM,'r')
%   hold on
%     plot(IOS_BLUE_NORM,'b')
% hold off




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

