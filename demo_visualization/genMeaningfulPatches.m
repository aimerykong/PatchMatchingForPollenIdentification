function [xList,yList] = genMeaningfulPatches(im, patchNum, patchSize, threshVal)
% generating meaning patches and returning their center positions in the
% image.
% 
% see also patchSIFT, main_patchExtraction_CNNfeature.m
%
% Shu Kong
% 04/13/2015

%%
imSize = size(im);
y = randperm(imSize(1)-patchSize);
x = randperm(imSize(2)-patchSize);

candiNum = min([200, length(y), length(x)]);
y = y(1:candiNum);
x = x(1:candiNum);

BW1 = edge(im,'Canny');
%figure(1); subplot(1,2,2); imshow(BW1, 'border', 'tight');
%subplot(1,2,1); imshow(im); hold on;
%patchSet = zeros(patchSize, patchSize, patchNum);

count = 0;
xList = zeros(patchNum,1);
yList = zeros(patchNum,1);
for i = 1:candiNum
    patchTMPbinary = BW1( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    %patchTMPorginal = im( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    
    if sum(patchTMPbinary(:)) > patchSize^2*threshVal
        count = count + 1;        
        %xListUL(count) = x(i);
        %yListUL(count) = y(i);
        
        xList(count) = x(i)+patchSize/2;
        yList(count) = y(i)+patchSize/2;
        %patchSet(:,:,count) = patchTMPorginal;
        %line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
    end
    if count >= patchNum
        break;
    end
end

if count < patchNum
    fprintf('Error occurs!\n');
end


