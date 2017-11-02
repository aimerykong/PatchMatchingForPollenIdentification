clear
close all
clc

addpath(genpath('../fossilPollen/code/ompbox10'));
addpath(genpath('../fossilPollen/code/spams-matlab'));

%% load dataset
numTrain = 65; % the number of training data in each class

% databaseDIR = './mixedSamples_thumbnail_CNNfeature';
databaseDIR = './mixedSamples_binaryMask_CNNfeature_globalContrastNorm';

setbyLayer = dir(databaseDIR);
a = zeros(numel(setbyLayer)-2,1);
for i = 3:numel(setbyLayer)
    aa = strfind(setbyLayer(i).name, '_');
    a(i-2) = str2double(  setbyLayer(i).name(aa(1)+1:aa(2)-1)  );
    setbyLayer(i).layerID = a(i-2);
end
layerIDall = unique(a);
% accList = zeros(2, length(layerIDall));

shapeType = strfind(databaseDIR, '_');
shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);

%% fetch the training and validation set
%{
for LAYERID = layerIDall(:)' % features extract at a specific layer
        
    fprintf('layer-%d\n', LAYERID);
    categoryList = [];
    for i = 1:numel(setbyLayer)
        if setbyLayer(i).layerID == LAYERID
            categoryList = [categoryList, setbyLayer(i)];
        end
    end
    
    pollenName = cell(1,numel(categoryList));
    for i = 1:numel(categoryList)
        dataSet = {};
        imList = dir( fullfile(databaseDIR, categoryList(i).name, '*.mat') );
        
        a = strfind(categoryList(i).name, '_');
        categoryName = categoryList(i).name(a(end)+1:end);
        
        pollenName{i} = categoryName;
        clusterLabel = zeros(1, numel(imList));
        imID = zeros(1,numel(imList));
        fprintf('\tpollen--%s\n', categoryName);
        for imgID = 1:numel(imList)
            imgName = fullfile(databaseDIR, categoryList(i).name, imList(imgID).name );
            tmpMat = load(imgName);
            dataSet{end+1} = tmpMat;
            imID(imgID) = imgID;
            
            a = strfind(imgName,'_K');
            b = strfind(imgName,'.mat');
            clusterLabel(imgID) = str2num(imgName(a+2:b-1));
        end
%         save( [categoryName, '_layer_', num2str(LAYERID), '_shapeType_', shapeType, '.mat'] );
        save( [pollenName{i}, '_layer_' num2str(LAYERID) '_shapeType_' shapeType '_globalContrastNorm.mat'] );
    end
end
%}

%% split into training and validation set
curLayer = 24; %[22, 23, 24, 25]
pollenName = {'mixedTestSet'};
trainDict = [];
trainDictLoc = [];
DictClassLabel = [];
DictClusterLabel = [];

testData = {};
testDataClassLabel = [];
testDataClusterLabel = [];
testDataLoc = {};

fprintf('fetch test set, features at layer-%d\n', curLayer);
categID = 1;
fprintf('\t%s  ', pollenName{categID});
tmpMat = load( [pollenName{categID}, '_layer_' num2str(curLayer) '_shapeType_' shapeType '_globalContrastNorm.mat'] );
imList = tmpMat.imList;
%% testing set
fprintf('testing set  ');
for teID = 1:length(tmpMat.dataSet)
    imFeat = tmpMat.dataSet{teID}.patchFeat(1:end-2,:);
    tmp = sqrt(sum(imFeat.^2,1));
    tmpIdx = find(tmp>1);
    imFeat(:,tmpIdx) = imFeat(:,tmpIdx) ./ repmat( tmp(tmpIdx), size(imFeat,1), 1 );
    testData{end+1} = imFeat;
    
    patchLoc = tmpMat.dataSet{teID}.patchFeat(end-1:end,:);
    patchLoc = patchLoc-1;
    imSize = tmpMat.dataSet{teID}.imSize;
    feaSize = tmpMat.dataSet{teID}.feaSize(1:2);
    patchLoc = bsxfun(@rdivide, patchLoc, feaSize(:));
    patchLoc = bsxfun(@times, patchLoc, imSize(:));
    patchLoc = patchLoc + 1;
    patchLoc = bsxfun(@minus, patchLoc, imSize(:)/2);
    testDataLoc{end+1} = patchLoc;
    
    testDataClassLabel(end+1) = categID;
    testDataClusterLabel(end+1) = tmpMat.clusterLabel(teID) ;
end
fprintf('done\n');

save(['test_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], ...
    'testData', 'testDataLoc', 'testDataClassLabel', 'testDataClusterLabel', 'imList');


