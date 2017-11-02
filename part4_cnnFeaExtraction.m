% extract the CNN features at different layers for all images in the
% dataset.
%
%
% Shu Kong
% 04/13/2015

%clear
close all
clc;

% setup MtConvNet in MATLAB
run ../fossilPollen/code/matconvnet/matlab/vl_setupnn
if ~exist('net', 'var')
    net = load('../fossilPollen/code/matconvnet/imagenet-vgg-verydeep-19.mat') ;
end
vl_simplenn_display(net);

%% obtain and preprocess an image
patchSize = 52;
maxLayerFeaMap = 24;
patchNum = 50;
overlapBorder = 20;

databaseSourceDir = 'database_binaryMask_canonicalShape_globalContrastNorm';
databaseDestinationDir = './database_binaryMask_CNNfeature_globalContrastNorm';

% databaseSourceDir = './database_binaryMask_canonicalShape';
% databaseDestinationDir = './database_binaryMask_CNNfeature';

% databaseSourceDir = './database_thumbnail_canonicalShape';
% databaseDestinationDir = './database_thumbnail_CNNfeature';

%% extract CNN features
if ~isdir( databaseDestinationDir )
    mkdir(databaseDestinationDir);
end

categoryList = dir( strcat(databaseSourceDir, '/* fossil'));

fprintf('fetch data...\n');
for categoryID = 1:numel(categoryList)
    categoryName = categoryList(categoryID).name;
    fprintf('\ncategory-%s...\n', categoryName);
    imList = dir( fullfile(databaseSourceDir, categoryName,'*.jpg') );
    
    for imIDX = 1:numel(imList)
        fprintf('.');
        imFileName = imList(imIDX).name;
        im = imread( fullfile(databaseSourceDir, categoryList(categoryID).name, imFileName) );
        imBKUP = im;
        imSize = size(im);
        
        % select local patches
        patchLoc = genCanonicalPatches(imBKUP, patchSize, patchNum, overlapBorder);
        
        im3Mode = repmat(im, [1, 1, 3]);
        im_ = single(im3Mode) ; % note: 255 range
        im_ = im_ - imresize(net.normalization.averageImage, [size(im_,1), size(im_,2)]) ;
        res = pollenFeatureByCNN(net, im_, maxLayerFeaMap) ;

                
        for layerID = 24:maxLayerFeaMap
            if strcmp( net.layers{layerID}.type, 'conv') || strcmp( net.layers{layerID}.type, 'relu')
                imFeaMap = res(layerID+1).x;
                A = max(imFeaMap, [], 3);
                
                patchLocLayer = patchLoc-1;
                patchLocLayer = bsxfun(@rdivide, patchLocLayer, reshape(size(imBKUP),[],1));
                patchLocLayer = bsxfun(@times, patchLocLayer, reshape(size(A),[],1));
                patchLocLayer = round(patchLocLayer+1);
                
                feaSize = size(res(layerID+1).x);
                patchFeat = zeros( size(res(layerID+1).x,3)+2, size(patchLocLayer,2) );
                for pp = 1:size(patchLocLayer,2)
                    patchFeat(:, pp) = [squeeze( res(layerID+1).x(patchLocLayer(1,pp), patchLocLayer(2,pp),:)); patchLocLayer(:,pp)];
                end
                [junk, imFileNameTMP, imFileNameExt] = fileparts(imFileName);
                
                if ~isdir( fullfile(databaseDestinationDir, strcat('layer_', num2str(layerID), '_', categoryName)) )
                    mkdir(fullfile(databaseDestinationDir, strcat('layer_', num2str(layerID), '_', categoryName)));
                end
                
                filename = fullfile(databaseDestinationDir, strcat('layer_', num2str(layerID), '_', categoryName), strcat(imFileNameTMP,'.mat'));
                save(filename, 'patchFeat', 'feaSize','imSize');
                %fprintf('net layer-%d conv --- res layer-%d -- %s\n', layerID, layerID+1, filename);
            end
        end
    end
end
fprintf('done!\n');





