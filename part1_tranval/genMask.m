function imMask = genMask(im)

imSize = size(im);

thres = 0.2;
hsize = 7;
sigma = 2;
GauF = fspecial('gaussian', hsize, sigma);

im = double(im)/255;
imMask = imfilter(im, GauF, 'replicate', 'same', 'conv');
imMask(imMask > thres) = 1;
imMask(imMask <= thres) = 0;

SE = strel('rectangle', [3 3]);
imMask = imdilate(imMask, SE);