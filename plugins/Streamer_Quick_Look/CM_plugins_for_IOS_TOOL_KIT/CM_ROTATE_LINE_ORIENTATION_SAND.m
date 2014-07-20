%CM 20140605
% This code defines the angle of a line traced on an image and rotates the
% image of this angle and the line
% the line is then horizontal and the data can be extracted and averaged in
% a rectangular straight box

% INPUTS
% CX : X Coordinates of the line to profile in Im
% CY: Y Coordinates of the line to profile in Im
% NOTE: Only the extremities of the line will be taken into account
% Im : Image either an average or a raw frame
% thickness : Defines how thick the extracted line profile should be
% NOTE: the extracted line is averaged in the Y dimension to give a cleaner
% profile than with a simple improfile

% Xnew= X Coordinates of the line to profile in rotated_Im
% Ynew = Y Coordinates of the line to profile in rotated_Im
% NOTE: mean(round(Xnew))is the pixel size of the profile
% NOTE: Ynew(1) and Ynew(2) should be equal because the returned line is
% horizontal

% rotated_Im = cropped image after rotation to have the line be horizontal
% FAT_LINE_AVG= wanted input


function [FAT_LINE_AVG,Xnew,Ynew,rotated_Im]=CM_ROTATE_LINE_ORIENTATION_SAND(CX,CY,Im,thickness,display)

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
FAT_LINE_AVG=rotated_Im(Ynew_ofs(1):Ynew_ofs(2),Xnew(1):Xnew(2));
FAT_LINE_AVG=mean(FAT_LINE_AVG,1);

if (display==1)    
    figure (90)
    imagesc(Im);colormap (gray);axis 'equal'
    hold on
    plot (CX_Before_ROT,CY_Before_ROT)
    hold off
    
    figure (91)
    imagesc(rotated_Im);colormap (gray);axis 'equal'
    hold on
    plot (Xnew,Ynew)
    hold off
    
    figure (92)
    plot (FAT_LINE_AVG)
end

end