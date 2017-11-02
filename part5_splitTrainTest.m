clear
close all
clc

addpath(genpath('../fossilPollen/code/ompbox10'));
addpath(genpath('../fossilPollen/code/spams-matlab'));

%% load dataset
numTrain = 65; % the number of training data in each class

% databaseDIR = './database_thumbnail_CNNfeature';
% databaseDIR = './database_binaryMask_CNNfeature';
databaseDIR = './database_binaryMask_CNNfeature_globalContrastNorm';

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
%%{
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
        b = strfind(categoryList(i).name, ' ');
        categoryName = categoryList(i).name(a(end)+1:b-1);
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
        save( [categoryName, '_layer_', num2str(LAYERID), '_shapeType_', shapeType, '_globalContrastNorm.mat'] );
%         save( [categoryName, '_layer_', num2str(LAYERID), '_shapeType_', shapeType, '.mat'] );
    end
end
%}

%% split into training and validation set
curLayer = 24; %[22, 23, 24, 25]
pollenName = {'critchfieldii', 'glauca', 'mariana'};
trainDict = [];
trainDictLoc = [];
DictClassLabel = [];
DictClusterLabel = [];
DictImageLabel = [];

testData = {};
testDataClassLabel = [];
testDataClusterLabel = [];
testDataLoc = {};

fprintf('split into training/validation sets, features at layer-%d\n', curLayer);
imgCount = 1;
for categID = 1:length(pollenName)
    fprintf('\t%s  ', pollenName{categID});
    tmpMat = load( [pollenName{categID}, '_layer_' num2str(curLayer) '_shapeType_' shapeType '_globalContrastNorm.mat'] );
    
    %% training set
    fprintf('training set  ');    
    for trID = 1:numTrain
        imFeat = tmpMat.dataSet{trID}.patchFeat(1:end-2,:);
        tmp = sqrt(sum(imFeat.^2,1));
        tmpIdx = find(tmp>1);
        imFeat(:,tmpIdx) = imFeat(:,tmpIdx) ./ repmat( tmp(tmpIdx), size(imFeat,1), 1 );        
        trainDict = [trainDict imFeat];
        
        patchLoc = tmpMat.dataSet{trID}.patchFeat(end-1:end,:);
        patchLoc = patchLoc-1;
        imSize = tmpMat.dataSet{trID}.imSize;
        feaSize = tmpMat.dataSet{trID}.feaSize(1:2);
        patchLoc = bsxfun(@rdivide, patchLoc, feaSize(:));
        patchLoc = bsxfun(@times, patchLoc, imSize(:));
        patchLoc = patchLoc + 1;
        patchLoc = bsxfun(@minus, patchLoc, imSize(:)/2);
        trainDictLoc = [trainDictLoc, patchLoc]; 
        
        DictClassLabel = [DictClassLabel ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*categID ];
        DictClusterLabel = [DictClusterLabel ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*tmpMat.clusterLabel(trID) ];
        DictImageLabel = [DictImageLabel, ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*imgCount];
        imgCount = imgCount + 1;
    end
    
    %% testing set
    fprintf('testing set  ');
    for teID = numTrain+1:length(tmpMat.dataSet)
        imFeat = tmpMat.dataSet{teID}.patchFeat(1:end-2,:);
        tmp = sqrt(sum(imFeat.^2,1));
        tmpIdx = find(tmp>1);
        imFeat(:,tmpIdx) = imFeat(:,tmpIdx) ./ repmat( tmp(tmpIdx), size(imFeat,1), 1 );          
        testData{end+1} = imFeat;
        
        patchLoc = tmpMat.dataSet{teID}.patchFeat(end-1:end,:);
        patchLoc = patchLoc-1;
        imSize = tmpMat.dataSet{trID}.imSize;
        feaSize = tmpMat.dataSet{trID}.feaSize(1:2);
        patchLoc = bsxfun(@rdivide, patchLoc, feaSize(:));
        patchLoc = bsxfun(@times, patchLoc, imSize(:));
        patchLoc = patchLoc + 1;
        patchLoc = bsxfun(@minus, patchLoc, imSize(:)/2);
        testDataLoc{end+1} = patchLoc;
        
        testDataClassLabel(end+1) = categID;
        testDataClusterLabel(end+1) = tmpMat.clusterLabel(teID) ;
    end
    fprintf('done\n');
end
save(['trval_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], ...
    'trainDict', 'trainDictLoc', 'DictClassLabel', 'DictClusterLabel', 'DictImageLabel', ...
    'testData', 'testDataLoc', 'testDataClassLabel', 'testDataClusterLabel');


