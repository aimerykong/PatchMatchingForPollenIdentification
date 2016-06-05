% test on how to extract convolutional neural network (CNN) feature for a
% single image.
%
%
%
% Shu Kong
% 04/13/2015

%clear
close all
clc;

% setup MtConvNet in MATLAB
% run ./matconvnet/matlab/vl_setupnn
% if ~exist('net', 'var')
%     net = load('./matconvnet/imagenet-vgg-verydeep-19.mat') ;
% end

%% obtain and preprocess an image
maxLayerFeaMap = 27;

categoryName = 'glauca - modern';
imFileName = 'P glauca B1500 Pos 3.tif_Files.jpg';


im = imread( fullfile( './database', categoryName,imFileName) );
im3Mode = repmat(im, [1, 1, 3]);

% im_ = single(im3Mode) ; % note: 255 range
% im_ = im_ - imresize(net.normalization.averageImage, [size(im_,1), size(im_,2)]) ;
% 
% vl_simplenn_display(net)

%% run the CNN
%res = pollenFeatureByCNN(net, im_, maxLayerFeaMap) ;

%% randomly select a set of meaningful patches

patchSize = 52;
imSize = size(im);

patchNum = 30;

y = randperm(imSize(1)-patchSize);
x = randperm(imSize(2)-patchSize);

candiNum = 200;
y = y(1:candiNum);
x = x(1:candiNum);

figure(1); subplot(2,3,1); imshow(im, 'border', 'tight'); title('random patch')
for i = 1:patchNum
    line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
end


BW1 = edge(im,'Canny');
figure(1); subplot(2,3,2); imshow(BW1, 'border', 'tight'); title('Canny Edges')


threVal = 0.0001;%0.15;
xList = zeros(patchNum,1);
yList = zeros(patchNum,1);
figure(1); subplot(2,3,3); imshow(im); hold on; title(strcat('selective patch edgeTreshold=', num2str(threVal)));
%patchSet = zeros(patchSize, patchSize, patchNum);
count = 0;
for i = 1:candiNum
    patchTMPbinary = BW1( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    patchTMPorginal = im( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    
    if sum(patchTMPbinary(:)) > patchSize^2*threVal
        count = count + 1;
        
        %xListUL(count) = x(i);
        %yListUL(count) = y(i);
        
        xList(count) = x(i)+patchSize/2;
        yList(count) = y(i)+patchSize/2;
        %patchSet(:,:,count) = patchTMPorginal;
        line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
    end
    if count >= patchNum
        break;
    end
end




threVal = 0.1;
count = 0;
xList = zeros(patchNum,1);
yList = zeros(patchNum,1);
figure(1); subplot(2,3,4); imshow(im); hold on; title(strcat('selective patch edgeTreshold=', num2str(threVal)));
for i = 1:candiNum
    patchTMPbinary = BW1( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    patchTMPorginal = im( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    
    if sum(patchTMPbinary(:)) > patchSize^2*threVal
        count = count + 1;
        
        xList(count) = x(i)+patchSize/2;
        yList(count) = y(i)+patchSize/2;

        line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
    end
    if count >= patchNum
        break;
    end
end





threVal = 0.15;
count = 0;
xList = zeros(patchNum,1);
yList = zeros(patchNum,1);
figure(1); subplot(2,3,5); imshow(im); hold on; title(strcat('selective patch edgeTreshold=', num2str(threVal)));
for i = 1:candiNum
    patchTMPbinary = BW1( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    patchTMPorginal = im( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    
    if sum(patchTMPbinary(:)) > patchSize^2*threVal
        count = count + 1;
        
        xList(count) = x(i)+patchSize/2;
        yList(count) = y(i)+patchSize/2;

        line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
    end
    if count >= patchNum
        break;
    end
end




threVal = 0.2;
count = 0;
xList = zeros(patchNum,1);
yList = zeros(patchNum,1);
figure(1); subplot(2,3,6); imshow(im); hold on; title(strcat('selective patch edgeTreshold=', num2str(threVal)));
for i = 1:candiNum
    patchTMPbinary = BW1( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    patchTMPorginal = im( y(i):y(i)+patchSize-1, x(i):x(i)+patchSize-1 );
    
    if sum(patchTMPbinary(:)) > patchSize^2*threVal
        count = count + 1;
        
        xList(count) = x(i)+patchSize/2;
        yList(count) = y(i)+patchSize/2;

        line([x(i), x(i)+patchSize-1, x(i)+patchSize-1, x(i), x(i)],  [y(i), y(i), y(i)+patchSize-1, y(i)+patchSize-1, y(i)], 'linewidth', 1, 'color', 'm' )
    end
    if count >= patchNum
        break;
    end
end
%A = reshape(patchSet, [patchSize, patchSize*size(patchSet,3)]);

%% test iteratively saving features at layers

%{
layerID = 21;
imFeaMap = res(layerID).x;
A = max(imFeaMap, [], 3);

xSpan = linspace(1, size(A,2), size(im,2));
ySpan = linspace(1, size(A,1), size(im,1));

xListNew = round(xSpan(xList));
yListNew = round(ySpan(yList));

a = sub2ind(size(A), yListNew, xListNew);
A(a) = Inf;
figure; imagesc(A); axis image;

destFolderName = '.\database_CNNfeature';


for i = 14:maxLayerFeaMap    
    if strcmp( net.layers{i}.type, 'conv') || strcmp( net.layers{i}.type, 'relu')        
        imFeaMap = res(i).x;
        A = max(imFeaMap, [], 3);
        xSpan = linspace(1, size(A,2), size(im,2));
        ySpan = linspace(1, size(A,1), size(im,1));
        
        xListNew = round(xSpan(xList));
        yListNew = round(ySpan(yList));
        
        patchFeat = zeros( size(res(i).x,3), numel(xListNew) );
        for pp = 1:length(xList)
            patchFeat(:, pp) = squeeze( res(i).x(yListNew(pp), xListNew(pp),:) );
        end
        [junk, imFileNameTMP, imFileNameExt] = fileparts(imFileName);
        
        if ~isdir( fullfile(destFolderName, strcat('layer_', num2str(i), '_', categoryName)) )
            mkdir(fullfile(destFolderName, strcat('layer_', num2str(i), '_', categoryName)));
        end
        
        filename = fullfile(destFolderName, strcat('layer_', num2str(i), '_', categoryName), strcat(imFileNameTMP,'.mat'));
%        save(filename, 'patchFeat');
        fprintf('net layer-%d conv --- res layer-%d -- %s\n', i, i+1, filename);
    end
end
%}




