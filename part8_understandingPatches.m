clear 
close all
clc
addpath('../toolbox/exportFig');

%%{
%% load dictionary

K = 300;  % select K desired receptive fields
lambdaL = 0.1;
lambdaD = 0.1; % discriminative term (2) -- make the selected points more "class-pure"
lambdaE = 200; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
lambdaB = 2;

% K = 600;  % select K desired receptive fields
% lambdaL = 0.01;
% lambdaD = 0.01; % discriminative term (2) -- make the selected points more "class-pure"
% lambdaE = 20; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
% lambdaB = 2;

curLayer = 24; %[22, 23, 24, 25]

databaseDIR = './database_binaryMask_CNNfeature_globalContrastNorm';
shapeType = strfind(databaseDIR, '_');
shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);

fprintf('load dataset...\n');
load(['trval_layer' num2str(curLayer) '_' shapeType '.mat']);
pollenName = {'critchfieldii', 'glauca', 'mariana'};

dictName = ['exemplarDict_K' num2str(K) 'L' num2str(lambdaL) '_D' num2str(lambdaD) '_E' num2str(lambdaE) '_B' num2str(lambdaB) '_globalContrastNorm.mat' ];
load(dictName);

% smaller dictionary size?
% exemplarIndex = exemplarIndex(1:150);

dict = trainDict(:,exemplarIndex);
dictLoc = trainDictLoc(:,exemplarIndex);
% dictLoc = trainDictLoc(:,exemplarIndex)*locationScaleFactor;
dictLabel = DictClassLabel(exemplarIndex);


%% load image information of training set
numTrain = 65; % the number of training data in each class

atomLoc_featSpace = [];
atomLoc_orgImgSpace = [];
atomImgLabel = [];
atomClassLabel = [];
atomClusterLabel = [];
imgList = {};
imgOrgSizeList = [];
imgFeaSizeList = [];
imgCount = 1;
for categID = 1:length(pollenName)
    fprintf('\t%s  ', pollenName{categID});
    tmpMat = load( [pollenName{categID}, '_layer_' num2str(curLayer) '_shapeType_' shapeType '_globalContrastNorm.mat'] );
    
    fprintf('training set  ');    
    for trID = 1:numTrain
        curImgName = tmpMat.imList(trID).name;
        [~,curImgName,~] = fileparts(curImgName);
        imgList{end+1} = [curImgName, '.jpg']; % image name
        
        patchLoc = tmpMat.dataSet{trID}.patchFeat(end-1:end,:);
        atomLoc_featSpace = [atomLoc_featSpace, patchLoc]; % location in feature space
        
        patchLoc = patchLoc-1;
        imSize = tmpMat.dataSet{trID}.imSize;
        imgOrgSizeList = [imgOrgSizeList, imSize(:)]; % size of original image
        
        feaSize = tmpMat.dataSet{trID}.feaSize(1:2);
        imgFeaSizeList = [imgFeaSizeList, feaSize(:)]; % size of feature maps
        
        patchLoc = bsxfun(@rdivide, patchLoc, feaSize(:));
        patchLoc = bsxfun(@times, patchLoc, imSize(:));
        patchLoc = patchLoc + 1;
        patchLoc = bsxfun(@minus, patchLoc, imSize(:)/2);        
        atomLoc_orgImgSpace = [atomLoc_orgImgSpace, patchLoc]; % location in image space
        
        atomClassLabel = [atomClassLabel ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*categID ];
        atomClusterLabel = [atomClusterLabel ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*tmpMat.clusterLabel(trID) ];
        atomImgLabel = [atomImgLabel, ones(1,size(tmpMat.dataSet{trID}.patchFeat,2))*imgCount];
        imgCount = imgCount + 1;
    end
    
    fprintf('done\n');
end
clear tmpMat
%}

%%
patchMosaic = zeros(1000,1000,3);
countMosaic = zeros(1000,1000,3);
MosaicCenter = size(countMosaic)/2 ;
MosaicCenter = floor(MosaicCenter(1:2));
patchCountAtClassLevel = zeros(1, 3);


for i = 1:length(exemplarIndex)
    exemplarID = exemplarIndex(i);
    curImgNameID = atomImgLabel(exemplarID);
    curPollenClass = atomClassLabel(exemplarID);
    curImgName = imgList{curImgNameID};
    curImgFeaSize = imgFeaSizeList(:, curImgNameID);
    curImgOrgSize = imgOrgSizeList(:, curImgNameID);
    curPatchLoc_featSpace = atomLoc_featSpace(:, exemplarID);
    curPatchLoc_orgImgSpace = atomLoc_orgImgSpace(:, exemplarID);
    
    patchCountAtClassLevel(curPollenClass) = patchCountAtClassLevel(curPollenClass) + 1;
    im = imread( fullfile(['database_' shapeType '_canonicalShape'], [pollenName{curPollenClass} ' fossil'], curImgName) );
    im = [zeros(size(im,1),2), im];
    im = [im, zeros(size(im,1),2)];
    im = [im; zeros(2, size(im,2))];
    im = [zeros(2, size(im,2)); im];
    
    patchCenterH = floor(size(im,1)/2 + curPatchLoc_orgImgSpace(1));
    patchCenterW = floor(size(im,2)/2 + curPatchLoc_orgImgSpace(2));
    
    curPatch = im( patchCenterH-26:patchCenterH+25, patchCenterW-26:patchCenterW+25);
%     curPatch = im( patchCenterH-50:patchCenterH+49, patchCenterW-50:patchCenterW+49);
    
    locTOBE = floor( MosaicCenter(:) + curPatchLoc_orgImgSpace(:) );
    
    patchMosaic( locTOBE(1)-26:locTOBE(1)+25, locTOBE(2)-26:locTOBE(2)+25, curPollenClass) = ...
        patchMosaic( locTOBE(1)-26:locTOBE(1)+25, locTOBE(2)-26:locTOBE(2)+25, curPollenClass) + double(curPatch);    
    countMosaic( locTOBE(1)-26:locTOBE(1)+25, locTOBE(2)-26:locTOBE(2)+25, curPollenClass) = ...
        countMosaic( locTOBE(1)-26:locTOBE(1)+25, locTOBE(2)-26:locTOBE(2)+25, curPollenClass) + 1;
    
    
%     patchMosaic( locTOBE(1)-50:locTOBE(1)+49, locTOBE(2)-50:locTOBE(2)+49, curPollenClass) = ...
%         patchMosaic( locTOBE(1)-50:locTOBE(1)+49, locTOBE(2)-50:locTOBE(2)+49, curPollenClass) + double(curPatch);    
%     countMosaic( locTOBE(1)-50:locTOBE(1)+49, locTOBE(2)-50:locTOBE(2)+49, curPollenClass) = ...
%         countMosaic( locTOBE(1)-50:locTOBE(1)+49, locTOBE(2)-50:locTOBE(2)+49, curPollenClass) + 1;
end

%% demonstrate
countMosaic(countMosaic==0) = 1;
figure;
for i = 1:length(pollenName)
    
    %patchMosaic(:,:,i) = max(patchMosaic(:,:,i),[],3);
    A = patchMosaic(:,:,i) ./ countMosaic(:,:,i);
    imshow( uint8( A ) );
    title(['patches belonging to ' pollenName{i}, ' (#' num2str(patchCountAtClassLevel(i)) ')' ] );
    
    figureName = ['./patches_' pollenName{i} '_K' num2str(K) 'L' num2str(lambdaL) '_D' num2str(lambdaD) '_E' num2str(lambdaE) '_B' num2str(lambdaB) '_globalContrastNorm'];
    export_fig(figureName); % , '-transparent', '-png'
end










