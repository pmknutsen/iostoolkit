% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
% Xi_N=Xi_S+box_lim(3); Yi_N=Yi_N+box_lim(1);


function CM_IOS_PLOT_LINES_SAND(STR_LINE)
%% PLOT FIGURE AND CALL IMPROFILE TO DEFINE THE BOX, TO DEFINE THE

counter_max=length(STR_LINE);
ColorSet = varycolor(counter_max);
if (counter_max<1)
    return
else
    figure (134) % ALL LINES
    clf    
    [AVG_image_to_display]=CM_SUBTRACT_BACKGROUND(STR_LINE(1).AVG_im,200,400,0);
    imagesc(AVG_image_to_display)
    
    hold on
    axis ('image');colormap (gray)
    
    for counter=1:1:counter_max
        
        CX_T=STR_LINE(counter).Xi_S+STR_LINE(counter).box_lim(3);
        CY_T=STR_LINE(counter).Yi_S+STR_LINE(counter).box_lim(1);
        box_lim=STR_LINE(counter).box_lim;
        [AVG_image_to_display]=CM_SUBTRACT_BACKGROUND(STR_LINE(counter).AVG_im,200,400,0);
       box_im=AVG_image_to_display(box_lim(1):box_lim(2),box_lim(3):box_lim(4));
       figure (134) %IMAGE WITH ALL COLORED LINE
        plot (CX_T,CY_T,'Color',ColorSet(counter,:))
        text(CX_T(1),CY_T(1),['' num2str(counter)],...
            'HorizontalAlignment','right',...
            'FontSize',12,'Color',ColorSet(counter,:))
       
        figure (135+counter)%INDIVIDUAL IMAGE FOR EACH LINE 
        subplot(1,2,1)
        imagesc(box_im)
        axis ('image');colormap (gray)
        hold on
        plot (STR_LINE(counter).Xi_S,STR_LINE(counter).Yi_S,'Color',ColorSet(counter,:))
        hold off
        
        subplot(1,2,2)
        imagesc(AVG_image_to_display)
        axis ('image');colormap (gray)
        hold on
        plot (CX_T,CY_T,'Color',ColorSet(counter,:))
        title(['LINE NUMBER ' num2str(counter)])
        
        hold off
        
    end
    figure(134)
    hold off
    
end
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







