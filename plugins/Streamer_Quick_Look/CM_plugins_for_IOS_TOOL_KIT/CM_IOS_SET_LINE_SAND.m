% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
% Xi_N=Xi_S+box_lim(3); Yi_N=Yi_N+box_lim(1);


function [Xi_S,Yi_S,box_lim,AVG_image,to_delete]=CM_IOS_SET_LINE_SAND(AVG_image,thickness,figure_number,line_number)

%% PLOT FIGURE AND CALL impixel TO DEFINE THE CENTER OF THE BOX
[AVG_image_to_display]=AVG_image;
figure (figure_number);
if (line_number==1)
    clf
    imagesc(AVG_image_to_display)
    axis ('image');colormap (gray)
end
[CX,CY,~] = impixel();
SP=100; % SP is the number of pixel spreading from the center of the box
center_of_the_box_Y=mean(CY); % added CM 2140722 to prevent the box from too cropping rotations by making a square box
center_of_the_box_X=mean(CX); % added CM 2140722 to prevent the box from too cropping rotations

box_lim(1)=ceil(center_of_the_box_Y-SP);
box_lim(2)=floor(center_of_the_box_Y+SP);
box_lim(3)=ceil(center_of_the_box_X-SP);
box_lim(4)=floor(center_of_the_box_X+SP);

if box_lim(1)<=1; box_lim(1)=1;end
if box_lim(2)>=size(AVG_image,1); box_lim(2)=size(AVG_image,1);end
if box_lim(3)<=1; box_lim(3)=1;end
if box_lim(4)>=size(AVG_image,2); box_lim(4)=size(AVG_image,2);end

box_im=AVG_image_to_display(box_lim(1):box_lim(2),box_lim(3):box_lim(4));
close (figure_number)

%% CALL IMPROFILES ON THE CROPPED IMAGE TO DEFINE THE VESSEL PERPENDICULAR TRACE
hfig2=figure (26+figure_number+1);
clf
set (hfig2,'position',  [800   100   500   500]);gcf %[left bottom width height]
imagesc(box_im)
axis ('image');colormap (gray)
[~,~,~,Xi_S,Yi_S]=improfile;
CX_T=Xi_S+box_lim(3);
CY_T=Yi_S+box_lim(1);

hfig2=figure (26+figure_number+1);
set (hfig2,'position',  [100   100   500   500]);gcf; %[left bottom width height]
clf
subplot(2,2,1)
imagesc(box_im)
axis ('image');colormap (gray)
hold on
plot (Xi_S,Yi_S,'r')
hold off
subplot(2,2,2)
imagesc(AVG_image_to_display)
axis ('image');colormap (gray)
hold on
plot (CX_T,CY_T,'r')
hold off

box_im=AVG_image_to_display(box_lim(1):box_lim(2),box_lim(3):box_lim(4));%EXTRACTED THE BOX FROM FRAME
to_delete=CM_QUICK_ROT(Xi_S,Yi_S,box_im,thickness,hfig2);
end

% refer to CM_ROTATE_LINE_ORIENTATION_SAND for explanations


function [to_delete]=CM_QUICK_ROT(CX,CY,Im,thickness,graph_handle)
CX_or=CX;
CY_or=CY;

CX=[CX(1) CX(end)]; % Extremities X of the line in original image
CY=[CY(1) CY(end)]; % Extremities Y of the line in original image

CX=CX-size(Im,2)/2; % SUBTRACTION OF THE ORIGIN
CY=CY-size(Im,1)/2;

thetaR=atan2(CY(2)-CY(1),CX(2)-CX(1)); %  CALCULATION OF THE ROTATION ANGLE OF THE LINE  THE IMAGE
[thetaP,rhoP]=cart2pol(CX,CY);  %  CALCULATION OF THE POLAR COORDINATES before ROTATION AND ORIGIN CENTERING
thetaN=thetaP - thetaR;

[Xnew,Ynew]=pol2cart(thetaN,rhoP); %CALCULATION OF THE LINE COORDINATES IN THE ROTATED IMAGE CENTERED
rotated_Im=imrotate(Im,rad2deg(thetaR),'bilinear','crop'); % CROPPED ROTATED IMAGE

Xnew=Xnew+(size(Im,2)/2);
Ynew=Ynew+(size(Im,1)/2);
Xnew=round(Xnew);
Ynew_ofs(1)=round (Ynew(1)-thickness/2);
Ynew_ofs(2)=round (Ynew(2)+thickness/2);
Ynew=round(Ynew);
%FAT_LINE_AVG=rotated_Im(round(Ynew(1)-thickness/2):round(Ynew(1)+thickness/2),round(Xnew(1)):round(Xnew(2)));
% in case the indices are out of bounds
if (min(Xnew)<1 || max(Xnew)>(size(Im,2)) || (Ynew_ofs(1)<1) || (Ynew_ofs(2)>(size(Im,1))))
    to_delete=1;
else 
    to_delete=0;
end

figure (graph_handle);
subplot(2,2,3);
imagesc(Im);colormap (gray);axis 'equal';
hold on;
plot (CX_or,CY_or);
hold off;
subplot(2,2,4);
imagesc(rotated_Im);colormap (gray);axis 'equal';
hold on;
%plot (Xnew,Ynew)
plot ([Xnew(1) Xnew(1) Xnew(2) Xnew(2) Xnew(1)],[Ynew_ofs(1) Ynew_ofs(2) Ynew_ofs(2) Ynew_ofs(1) Ynew_ofs(1)]);

if (to_delete)
   title('DELETED','FontWeight','bold','FontSize',20,'Color',[1 1 0],...
    'BackgroundColor',[1 0 0]);
end
if (to_delete==0)
  title('RECORDED','FontWeight','bold','FontSize',20,'Color',[0 0 1],...
      'BackgroundColor',[0 1 0]);
end

hold off

end






