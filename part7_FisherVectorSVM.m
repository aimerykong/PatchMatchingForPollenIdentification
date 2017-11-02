clear
close all
clc;

run('vlfeat/toolbox/vl_setup');

%% test example
%{
numFeatures = 5000;
dimension = 2;
data = rand(dimension, numFeatures);

numClusters = 30;

tic
[means, covariances, priors] = vl_gmm(data, numClusters);
toc

numDataToBeEncoded = 1000;
dataToBeEncoded = rand(dimension, numDataToBeEncoded);

tic
encoding = vl_fisher(dataToBeEncoded, means, covariances, priors);
toc
%}
% VLAD encoding
%{
tic
centers = vl_kmeans(data, numClusters);
toc

tic
kdtree = vl_kdtreebuild(centers);
nn = vl_kdtreequery(kdtree, centers, dataToBeEncoded);
toc
%}
%%
%{
numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

databaseDIR = './database_binaryMask_CNNfeature_globalContrastNorm';
shapeType = strfind(databaseDIR, '_');
shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);

fprintf('load dataset...\n');
load(['trval_layer' num2str(curLayer) '_' shapeType '.mat']);
pollenName = {'critchfieldii', 'glauca', 'mariana'};

%% GMM, and encoding training set
numClusters = 64;
fprintf('GMM...');
tic
[means, covariances, priors] = vl_gmm(trainDict, numClusters);
toc

fprintf('encoding training images...\n');
numTrainData = length(unique(DictImageLabel));
trainSet = [];
trainDataClassLabel = [];
for i_img = 1:numTrainData
    a = find(DictImageLabel==i_img);
    imFeat = trainDict(:,a);
    encoding = vl_fisher(imFeat, means, covariances, priors);
    trainSet = [trainSet, encoding];
    trainDataClassLabel = [trainDataClassLabel, DictClassLabel(a(1)) ];
end

%% encoding test set
fprintf('encoding testing images...\n');
testSet = [];
for i_img = 1:length(testData)
    encoding = vl_fisher(testData{i_img}, means, covariances, priors);
    testSet = [testSet, encoding];
end
save(sprintf('VGG_FV_SVM_numCluster%d.mat', numClusters), 'trainSet', 'testSet', 'trainDataClassLabel', 'testDataClassLabel' );
%}

%% SVM
fprintf('learning linear SVM...');
addpath(genpath('../toolbox/libsvm-3.20/matlab'));

numClusters = 128;
load(sprintf('VGG_FV_SVM_numCluster%d.mat', numClusters) );


% A = trainSet;
% A(A>=0) = 1;
% A(A<0) = -1;
% trainSet = A.*sqrt(abs(trainSet));
% a = sqrt(sum(trainSet.^2, 1));
% trainSet = trainSet ./ repmat(a, size(trainSet,1), 1);
% 
% A = testSet;
% A(A>=0) = 1;
% A(A<0) = -1;
% testSet = A.*sqrt(abs(testSet));
% a = sqrt(sum(testSet.^2, 1));
% testSet = testSet ./ repmat(a, size(testSet,1), 1);

accList = [];
ccList = [];
for cc = 20:15:200
    model = libsvm_svmtrain(trainDataClassLabel(:), trainSet', sprintf('-c %d -t 0',cc) ); % linear
    [predLabelVal, accuracy, dec_values] = libsvm_svmpredict(testDataClassLabel(:), testSet', model); % test the training data
    acc = mean(testDataClassLabel(:) == predLabelVal(:));
    ccList(end+1) = cc;
    accList(end+1) = acc;
    fprintf('\nSVM on validation set\n\taccuracy=%.4f cc=%d\n', acc, cc );
end

