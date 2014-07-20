%%
% HSIZE=200; works well  20140608
% SIGMA=400;
function [Sub_image]=CM_SUBTRACT_BACKGROUND_SAND(Image,HSIZE,SIGMA,display)
Image=double(Image);
mWin = fspecial('gaussian', HSIZE, SIGMA);
Image_to_sub = imfilter(Image,mWin,'symmetric','conv');
Sub_image=Image-Image_to_sub;

if (display)
    figure(450)
    
    set(gcf, 'Position', [100 100 800 600])
    subplot(2,3,1)
    imagesc(Image)
    title ('\fontsize{10} BEFORE BACK SUBTRACT')
    colormap(gray);axis 'image'
    
    subplot(2,3,2)
    imagesc(Sub_image)
    colormap(gray);axis 'image'
    title ('\fontsize{10} AFTER BACK SUBTRACT')
    
    subplot(2,3,3)
    imagesc(Image_to_sub)
    colormap(gray);axis 'image'
    title ('\fontsize{10} FILTER')
    
    subplot(2,3,4)
    plot (mean(Image,1),'r')
    hold on
    plot (mean(Image,2),'b')
    title(['\fontsize{10} BEFORE {\color{blue}DIM 1 ''\color{red}DIM 2}'])
    hold off
    
    
    subplot(2,3,5)
    plot (mean(Sub_image,1),'r')
    hold on
    plot (mean(Sub_image,2),'b')
    title(['\fontsize{10} AFTER {\color{blue}DIM 1 ''\color{red}DIM 2}'])
    hold off
end





