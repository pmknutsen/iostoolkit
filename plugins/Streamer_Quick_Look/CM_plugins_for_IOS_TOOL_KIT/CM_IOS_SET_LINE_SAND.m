% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
% Xi_N=Xi_S+box_lim(3); Yi_N=Yi_N+box_lim(1);


function [Xi_S,Yi_S,box_lim,AVG_image]=CM_IOS_SET_LINE_SAND(AVG_image,figure_number)

%% PLOT FIGURE AND CALL IMPROFILE TO DEFINE THE BOX, TO DEFINE THE
if (isempty (figure_number))
    figure_number=-1;
end
[AVG_image_to_display]=AVG_image;
hfig1=figure (25+figure_number+1);
clf
imagesc(AVG_image_to_display)
axis ('image');colormap (gray)

[CX,CY,~]=improfile;
SP=50; % SP is the number of pixel used to extend the box on the sides
box_lim(1)=floor(min(CY))-SP;
box_lim(2)=ceil(max(CY))+SP;
box_lim(3)=floor(min(CX))-SP;
box_lim(4)=ceil(max(CX))+SP;
box_im=AVG_image_to_display(box_lim(1):box_lim(2),box_lim(3):box_lim(4));
close (25+figure_number+1)
%% CALL IMPROFILES ON THE CROPPED IMAGE TO DEFINE THE VESSEL PERPENDICULAR TRACE
figure (26+figure_number+1)
clf
imagesc(box_im)
axis ('image');colormap (gray)
[~,~,~,Xi_S,Yi_S]=improfile;
CX_T=Xi_S+box_lim(3);
CY_T=Yi_S+box_lim(1);

hfig2=figure (26+figure_number+1);
clf
subplot(1,2,1)
imagesc(box_im)
axis ('image');colormap (gray)
hold on
plot (Xi_S,Yi_S,'r')
hold off
subplot(1,2,2)
imagesc(AVG_image_to_display)
axis ('image');colormap (gray)
hold on
plot (CX_T,CY_T,'r')
hold off


end








