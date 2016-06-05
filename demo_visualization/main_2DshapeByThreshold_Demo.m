clear
close all
clc

%%
dirDataset = './database';
categoryList = dir(dirDataset);
categoryList = categoryList(3:end);


thres = 0.3;
hsize = 7;
sigma = 2;
GauF = fspecial('gaussian', hsize, sigma);

for categID = 1:1 %numel(categoryList)
    imList = dir( fullfile(dirDataset, categoryList(categID).name , '*.jpg') );
    for imID = 1:1%numel(imList)
        im = imread( fullfile(dirDataset, categoryList(categID).name, imList(imID).name) );
        
        im = double(im)/255;
        subplot(3,3,1); imshow(im); title('original image');
        
        imBinary = im;
        imBinary(imBinary > thres) = 1;
        imBinary(imBinary <= thres) = 0;
        subplot(3,3,2); imshow(imBinary); title('image-binary');
        
        
        imBinaryBlur = imfilter(imBinary, GauF, 'replicate', 'same', 'conv');
        subplot(3,3,3); imshow(imBinaryBlur); title('image-binary-blur');
        
        imBlur = imfilter(im, GauF, 'replicate', 'same', 'conv');
        subplot(3,3,4); imshow(imBlur); title('image-blur');
        
        imBlurBinary = imBlur;
        imBlurBinary(imBlurBinary > thres) = 1;
        imBlurBinary(imBlurBinary <= thres) = 0;
        subplot(3,3,5); imshow(imBlurBinary); title('image-blur-binary');
        
        
        imBlurBinarySmall = imresize(imBlurBinary, [50, 50], 'nearest');
        subplot(3,3,6); imshow(imBlurBinarySmall);
        title('image-blur-binary-small');
        
        [x, y] = find(imBlurBinarySmall==1);
        K = convhull(x, y);
        imMask = 0*imBlurBinarySmall;
        [maskContour, maskArea] = convexHullDraw(imMask, x(K), y(K));
        subplot(3,3,7);
        imshow(maskContour);
        title('convex hull contour');
        subplot(3,3,8);
        imshow(maskArea);
        title('convex hull area');
        
        
        imBlurBinarySmallMask = imBlurBinarySmall;
        imBlurBinarySmallMaskcomplement = 1-imBlurBinarySmallMask;
        
        cntCMPO = bwconncomp(imBlurBinarySmallMaskcomplement, 4);
        CMPOholeList = zeros(numel(cntCMPO.PixelIdxList),1);
        for i = 1:length(cntCMPO.PixelIdxList)
            if isempty(find(cntCMPO.PixelIdxList{i}==1, 1)) &&...
                    isempty(find(cntCMPO.PixelIdxList{i}==size(imBlurBinarySmallMask,1), 1)) && ...
                    isempty(find(cntCMPO.PixelIdxList{i}==numel(imBlurBinarySmallMask(:,1:end-1))+1, 1)) && ...
                    isempty(find(cntCMPO.PixelIdxList{i}==numel(imBlurBinarySmallMask), 1))
                
                CMPOholeList(i) = 1;
                imBlurBinarySmallMask(cntCMPO.PixelIdxList{i}) = 1;
            end
        end
        subplot(3,3,9);
        imshow(imBlurBinarySmallMask);
        title('pollen shape mask');
        
        
        fprintf('%s\n', imList(imID).name);        
    end
end



%%
%{
imBlurBinarySmallMask = imBlurBinarySmall;
imBlurBinarySmallMaskcomplement = 1-imBlurBinarySmallMask;

cntCMPO = bwconncomp(imBlurBinarySmallMaskcomplement, 4);
CMPOholeList = zeros(numel(cntCMPO.PixelIdxList),1);
for i = 1:length(cntCMPO.PixelIdxList)
    if isempty(find(cntCMPO.PixelIdxList{i}==1, 1)) &&...
            isempty(find(cntCMPO.PixelIdxList{i}==size(imBlurBinarySmallMask,1), 1)) && ...
            isempty(find(cntCMPO.PixelIdxList{i}==numel(imBlurBinarySmallMask(:,1:end-1))+1, 1)) && ...
            isempty(find(cntCMPO.PixelIdxList{i}==numel(imBlurBinarySmallMask), 1))
        
        CMPOholeList(i) = 1;
        imBlurBinarySmallMask(cntCMPO.PixelIdxList{i}) = 1;
    end
end
%}



















