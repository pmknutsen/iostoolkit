function mImg = BinImagePixels(mImg, nB)
% Bin nB-by-nB pixels in an image (sum of neighbouring pixels)
% Note that the *sum* and not the average of pixels is returned.
%
% Per M Knutsen <pmknutsen@gmail.com>
% June 2014
%

vSize = size(mImg);

% Pad matrix to ensure that both dimensions are divisable by nB
mImg = padarray(mImg, mod(nB - mod(vSize, nB), nB), 'symmetric', 'post');
vSize = size(mImg);

mImg = sum(reshape(mImg, nB, []), 1);
mImg = reshape(mImg, vSize(1)/nB, []).'; % note transpose

mImg = sum(reshape(mImg, nB, []), 1);
mImg = reshape(mImg, vSize(2)/nB, []).'; % note transpose

mImg = mImg ./ (nB^2);

return