clear;
close all;
clc;

% K = 300;  % select K desired receptive fields
% lambdaL = 0.1;
% lambdaD = 0.1; % discriminative term (2) -- make the selected points more "class-pure"
% lambdaE = 200; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
% lambdaB = 2;


K = 600;  % select K desired receptive fields
smallK = 512; % extract first smallK atoms as the dictionary rather than the whole
lambdaL = 0.01;
lambdaD = 0.01; % discriminative term (2) -- make the selected points more "class-pure"
lambdaE = 20; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
lambdaB = 2;


numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

locationScaleFactor = 0.001;

T = 35; % nonzero elements
lambda = 0.2; % weight of location-aware penalty

%% load dataset
% databaseDIR = './database_thumbnail_CNNfeature';
% databaseDIR = './database_binaryMask_CNNfeature';
databaseDIR = './database_binaryMask_CNNfeature_globalContrastNorm';

shapeType = strfind(databaseDIR, '_');
shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);

fprintf('load dataset...\n');
load(['trval_layer' num2str(curLayer) '_' shapeType '.mat']);
pollenName = {'critchfieldii', 'glauca', 'mariana'};

dictName = ['exemplarDict_K' num2str(K) 'L' num2str(lambdaL) '_D' num2str(lambdaD) '_E' num2str(lambdaE) '_B' num2str(lambdaB) '_globalContrastNorm.mat' ];
load(dictName);

% smaller dictionary size?
exemplarIndex = exemplarIndex(1:smallK);

dict = trainDict(:,exemplarIndex);
dictLoc = trainDictLoc(:,exemplarIndex)*locationScaleFactor;
dictLabel = DictClassLabel(exemplarIndex);
        
load( ['test_layer' num2str(curLayer) '_' shapeType  '_globalContrastNorm.mat'] );
testDataClassLabel = testDataClassLabel*0-1;

%% location-aware sparse coding on test set
labelpredTest = zeros(1,numel(testData));
errorListTest = zeros(numel(pollenName),numel(testData));
scFeatTest = cell(1,numel(testData));
locFeatTest = cell(1,numel(testData));
for i_img = 1:numel(testData) % 1:numel(testData)
    
    fprintf('\t image-%d/%d ---', i_img, numel(testData));
    imFeat = testData{i_img};
    patchLoc = testDataLoc{i_img} * locationScaleFactor;
    
    % overall area of the patches cover
    boundaryPoints = boundary(patchLoc(1,:)', patchLoc(2,:)');
    xv = [patchLoc(1,boundaryPoints), patchLoc(1,boundaryPoints(1))];
    yv = [patchLoc(2,boundaryPoints), patchLoc(2,boundaryPoints(2))];
    polyArea = polyarea(xv,yv);
    
    % variance-related
    [varFea, score]= pca(patchLoc');
    
    locFeatTest{i_img} = [min(patchLoc,[],2); max(patchLoc,[],2); max(sqrt(sum(patchLoc.^2))); polyArea; varFea(:);]; % max(score',[],2); min(score',[],2) ];
    
    W = zeros(size(dict,2), size(patchLoc,2));
    for i = 1 : size(patchLoc,2)
        a = bsxfun(@minus, dictLoc, patchLoc(:,i) );
        a = sum(a.^2);
        a = sqrt(a(:));
        W(:, i) = a;
    end
    
    param.L = T; % not more than param.L non-zeros coefficients (default: min(size(D,1),size(D,2)))
    param.lambda = lambda;
    param.numThreads = -1; % number of processors/cores to use; the default choice is -1 and uses all the cores of the machine
    param.mode = 2; % penalized formulation
    
    A = mexLassoWeighted(imFeat, dict, W, param);
    
    tmpErrorList = ones(numel(pollenName),1)*Inf;
    for categID = 1:numel(pollenName)
        a = find( dictLabel == categID );
        err = imFeat-dict(:,a)*A(a,:);
        err = sum(err(:).^2);
        tmpErrorList(categID) = err;
    end
    errorListTest(:,i_img) = tmpErrorList(:);
    [valMIN, idxMIN] = min(tmpErrorList);
    labelpredTest(i_img) = idxMIN;
    if idxMIN == testDataClassLabel(i_img)
        flag = 'correct';
    else
        flag = 'wrong';
    end
    
    scFeatTest{i_img} = A;
    
    fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));
end

%% get the average-pooling sparse codes for the final representation
scXTest = zeros(size(scFeatTest{1},1), length(scFeatTest));
lcXtest = zeros(size(locFeatTest{1},1), length(scFeatTest));
for i = 1:length(scFeatTest)
    scXTest(:,i) = mean(scFeatTest{i},2);
    lcXtest(:,i) = locFeatTest{i};
%     scXTest(:,i) = max(scFeatTest{i},[],2);
end
% scXTest = scXTest ./ repmat( sqrt(sum(scXTest.^2,1)), size(scXTest,1), 1 );
scXTest = [scXTest;lcXtest]; %

%% load classifier
classifier = load( 'valscFeatOnExemplar_layer24_binaryMask_T35_lambda0.2_ResultACC0.79708_globalContrastNorm.mat', 'U', 'alphaIdx', 'alphaList' );
alpha = classifier.alphaList(classifier.alphaIdx);
U = classifier.U;
clear classifier

%% learning classifier from training set and test on validation set
% accuracy on validation set
YhatTest = U'*scXTest;
[~,predLabelTest] = max(YhatTest,[],1);

save('testOnMixedSet_globalContrastNorm.mat', 'predLabelTest', 'imList');
fn = fopen('testOnMixedSet_globalContrastNorm.csv', 'w');
fprintf(fn, 'fileName,predLabelIdx,predLabel\n');

for i = 1:length(predLabelTest)
    fileName = imList(i).name;
    a = strfind(fileName, '_K');
    fileName = fileName(1:a-1);
    fprintf(fn, '%s,%d,%s\n', fileName, predLabelTest(i), pollenName{predLabelTest(i)});
end
fclose(fn);


