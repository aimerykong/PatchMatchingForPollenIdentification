clear;
close all;
clc;

addpath(genpath('../fossilPollen/code/ompbox10'));
addpath(genpath('../fossilPollen/code/spams-matlab'));

% K = 300;  % select K desired receptive fields
% lambdaL = 0.1;
% lambdaD = 0.1; % discriminative term (2) -- make the selected points more "class-pure"
% lambdaE = 200; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
% lambdaB = 2;


K = 600;  % select K desired receptive fields
smallK = 300; % extract first smallK atoms as the dictionary rather than the whole
lambdaL = 0.01;
lambdaD = 0.01; % discriminative term (2) -- make the selected points more "class-pure"
lambdaE = 20; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
lambdaB = 2;

numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

locationScaleFactor = 0.0014;

T = 45; % nonzero elements default 40
lambdaLoc = 0; % weight of location-aware penalty default 0.09
lambda1 = 0.05; % weight for sparsity

%{
----> lambda1
80.146% (0.095) 80.5839(0.09) 81.3139%(0.085) 81.6058%(0.08)
81.3139%(0.07) 82.0438%(0.06) 82.9197%(0.055) 83.0657%(0.051)
83.2117%(0.05) 83.0657%(0.049) 82.4818%(0.045) 81.8978%(0.04) 80.146%(0.03) 

---- lambda1=0.05 ----> locationScaleFactor
78.2482%(0.0005) 82.1898%(0.0008) 83.2117%(0.001) 82.4818%(0.0011)
82.9197%(0.0012) 79.1241%(0.002) 76.9343%(0.003)  
%}

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
%% visualize exemplars' locations        
%{
dictLocTMP = trainDictLoc(:,exemplarIndex);
img = imread('./database_binaryMask_canonicalShape/glauca fossil/296_Nelson Lake NE906 Pos 5_K1.jpg');
sz = size(img);

dictLocTMP = bsxfun(@plus, dictLocTMP, sz(:)/2);

imshow(img);
hold on;
for i = 1:size(dictLocTMP,2)
    rectangle('position',[dictLocTMP(2,i)-26 dictLocTMP(1,i)-26 52 52], 'EdgeColor','r')
end
hold off;
%}
%% location-aware sparse coding on training set
numTrainData = length(unique(DictImageLabel));

trainClassLabel = zeros(1, numTrainData);
labelpredTrain = zeros(1, numTrainData);
errorListTrain = zeros(numel(pollenName), numTrainData);
scFeatTrain = cell(1, numTrainData);
locFeatTrain = cell(1, numTrainData);

dictGram = (dict'*dict);
dictMultiplier = (dict'*dict)\dict';
dictInverse = pinv(dict);

for i_img = 1:numTrainData
    
    fprintf('\t image-%d/%d ---', i_img, numTrainData );
    a = find(DictImageLabel==i_img);
    imFeat = trainDict(:,a);
    patchLoc = trainDictLoc(:, a) * locationScaleFactor;
    trainClassLabel(i_img) = DictClassLabel(a(1));
    
    % overall area of the patches cover
    boundaryPoints = boundary(patchLoc(1,:)', patchLoc(2,:)');
    xv = [patchLoc(1,boundaryPoints), patchLoc(1,boundaryPoints(1))];
    yv = [patchLoc(2,boundaryPoints), patchLoc(2,boundaryPoints(2))];
    polyArea = polyarea(xv,yv);
    
    % variance-related
    [varFea, score]= pca(patchLoc');
    
    % leftmost location, right most location, maximum radius, area covered
    locFeatTrain{i_img} = [min(patchLoc,[],2); max(patchLoc,[],2); max(sqrt(sum(patchLoc.^2))); polyArea; varFea(:); ...
         ]; % max(score',[],2); min(score',[],2)   % max(score(:));min(score(:));];
    
    A = zeros(size(dict,2), size(patchLoc,2));
    W = zeros(size(dict,2), size(patchLoc,2));
    for i = 1 : size(patchLoc,2)
        a = bsxfun(@minus, dictLoc, patchLoc(:,i) );
%         a = sum(a.^2); a = sqrt(a(:)); a = a./max(a(:));
        a = sum(abs(a)); a = a.^0.5 ./max(a(:));
        W(:, i) = a(:);        
        %{
        a = a./max(a);
        b = bsxfun(@minus, dict, imFeat(:,i));
%         b = sum(abs(b),2);
        b = sqrt( sum( b.^2 ) );
        b = b(:) ./max(b(:));
        a = b + lambda*a(:);
%         a = b .* (lambda*a(:));
        [~, IDX] = sort(a, 'ascend');
        idx = IDX(1:T);
        subDict = dict(:, idx);
        a = (subDict'*subDict)\subDict'*imFeat(:,i);        
%         a = dictMultiplier*imFeat(:,i);
%         a = a(idx);
        %}        
        %{
        sca = (dictGram+lambda*diag(W(:, i)).^2) \ (dict'*imFeat(:,i));
        a = abs(sca);
        [~, IDX] = sort(a, 'descend');
        idx = IDX(1:T);
        A(idx, i) = sca(idx);
        %}
        %{
        v = (dictGram+lambdaLoc*diag(W(:, i)).^2) \ (dict'*imFeat(:,i));
        a = sign(v) .* max( [abs(v)-2*lambda1, zeros(length(v),1)], [], 2) ;
        A(:, i) = a;
        %}
        
%         v = dictMultiplier*imFeat(:,i);
%         a = sign(v) .* max( [abs(v)-2*lambda1*W(:, i), zeros(length(v),1)], [], 2) ;
%         A(:, i) = a;
    end
    V = dictMultiplier*imFeat;
    tmpV = abs(V) - 2*lambda1*W;
%     tmpV = repmat(tmpV,[1,1,2]);
%     tmpV(:,:,2) = 0;
%     A = sign(V).* max(tmpV,[],3);
    tmpV(tmpV<0)=0;
    A = sign(V).* tmpV;
    
    
    tmpErrorList = ones(numel(pollenName),1)*Inf;
    for categID = 1:numel(pollenName)
        a = find( dictLabel == categID );
        err = imFeat-dict(:,a)*A(a,:);
        err = sum(err(:).^2);
        tmpErrorList(categID) = err;
    end
    errorListTrain(:,i_img) = tmpErrorList(:);
    [valMIN, idxMIN] = min(tmpErrorList);
    labelpredTrain(i_img) = idxMIN;
    if idxMIN == testDataClassLabel(i_img)
        flag = 'correct';
    else
        flag = 'wrong';
    end
    
    scFeatTrain{i_img} = A;
    
    fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));
end

% accuracy
acc = mean(trainClassLabel == labelpredTrain);
fprintf('\naccuracy=%.4f (lambda1=%.4f, lambdaLoc=%.4f, T=%d)\n\n', acc, lambda1, lambdaLoc, T);

for categID = 1:numel(pollenName)
    a = find(trainClassLabel==categID);
    accTMP = mean(trainClassLabel(a) == labelpredTrain(a));
    fprintf('%s acc:%.4f\n', pollenName{categID}, accTMP);
end

%% location-aware sparse coding on validation set
labelpredVal = zeros(1,numel(testData));
errorListVal = zeros(numel(pollenName),numel(testData));
scFeatVal = cell(1,numel(testData));
locFeatVal = cell(1,numel(testData));
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
    
    locFeatVal{i_img} = [min(patchLoc,[],2); max(patchLoc,[],2); max(sqrt(sum(patchLoc.^2))); polyArea; varFea(:);...
         ]; % max(score',[],2); min(score',[],2)   % max(score(:));min(score(:));];
    
    A = zeros(size(dict,2), size(patchLoc,2));
    W = zeros(size(dict,2), size(patchLoc,2));
    for i = 1 : size(patchLoc,2)
        a = bsxfun(@minus, dictLoc, patchLoc(:,i) );
%         a = sum(a.^2); a = sqrt(a(:)); a = a./max(a(:));
        a = sum(abs(a)); a = a.^0.5 ./max(a(:));
        W(:, i) = a(:);   
        %{
        a = a./max(a);
        b = bsxfun(@minus, dict, imFeat(:,i));
%         b = sum(abs(b),2);
        b = sqrt( sum( b.^2 ) );
        b = b(:) ./max(b(:));
        a = b + lambda*a(:);
%         a = b .* (lambda*a(:));
        [~, IDX] = sort(a, 'ascend');
        idx = IDX(1:T);
        subDict = dict(:, idx);
        a = (subDict'*subDict)\subDict'*imFeat(:,i);        
%         a = dictMultiplier*imFeat(:,i);
%         a = a(idx);
        %}
        %{
        sca = (dictGram+lambda*diag(W(:, i)).^2) \ (dict'*imFeat(:,i));
        a = abs(sca);
        [~, IDX] = sort(a, 'descend');
        idx = IDX(1:T);
        A(idx, i) = sca(idx);
        %}
        %{
        v = (dictGram+lambdaLoc*diag(W(:, i)).^2) \ (dict'*imFeat(:,i));
        a = sign(v) .* max( [abs(v)-2*lambda1, zeros(length(v),0)], [], 2) ;
        A(:, i) = a;
        %}
        
%         v = dictMultiplier*imFeat(:,i);
%         a = sign(v) .* max( [abs(v)-2*lambda1*W(:, i), zeros(length(v),1)], [], 2) ;
%         A(:, i) = a;
    end
    V = dictMultiplier*imFeat;
    tmpV = abs(V) - 2*lambda1*W;
%     tmpV = repmat(tmpV,[1,1,2]);
%     tmpV(:,:,2) = 0;
%     A = sign(V).* max(tmpV,[],3);
    tmpV(tmpV<0)=0;
    A = sign(V).* tmpV;
    
    tmpErrorList = ones(numel(pollenName),1)*Inf;
    for categID = 1:numel(pollenName)
        a = find( dictLabel == categID );
        err = imFeat-dict(:,a)*A(a,:);
        err = sum(err(:).^2);
        tmpErrorList(categID) = err;
    end
    errorListVal(:,i_img) = tmpErrorList(:);
    [valMIN, idxMIN] = min(tmpErrorList);
    labelpredVal(i_img) = idxMIN;
    if idxMIN == testDataClassLabel(i_img)
        flag = 'correct';
    else
        flag = 'wrong';
    end
    
    scFeatVal{i_img} = A;
    
    fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));
end

%% sparse-coding based accuracy
acc = mean(testDataClassLabel == labelpredVal);
fprintf('\naccuracy=%.4f (lambda1=%.4f, lambdaLoc=%.4f, T=%d)\n\n', acc, lambda1, lambdaLoc, T);
for categID = 1:numel(pollenName)
    a = find(testDataClassLabel==categID);
    accTMP = mean(testDataClassLabel(a) == labelpredVal(a));
    fprintf('\t%s acc:%.4f\n', pollenName{categID}, accTMP);
end

%% get the average-pooling sparse codes for the final representation
scXTrain = zeros(size(scFeatTrain{1},1), length(scFeatTrain));
YTrain = zeros(numel(pollenName), length(scFeatTrain));
locTrain = zeros(numel(pollenName), length(scFeatTrain));
lcXTrain = zeros(size(locFeatTrain{1},1), length(scFeatTrain));

for i = 1:length(scFeatTrain)
%     TMP = scFeatTrain{i};
%     TMP(TMP<0) = 0;
%     scXTrain(:,i) = mean(TMP,2);
    scXTrain(:,i) = mean(scFeatTrain{i},2);
%     scXTrain(:,i) = mean( abs(scFeatTrain{i}),2);
%     scXTrain(:,i) = max(scFeatTrain{i},[],2);

    lcXTrain(:,i) = locFeatTrain{i};
    YTrain(trainClassLabel(i),i) = 1;
end
% scXTrain = scXTrain ./ repmat( sqrt(sum(scXTrain.^2,1)), size(scXTrain,1), 1 );

scXVal = zeros(size(scFeatVal{1},1), length(scFeatVal));
YVal = zeros(numel(pollenName), length(scFeatVal));
lcXval = zeros(size(locFeatVal{1},1), length(scFeatVal));
for i = 1:length(scFeatVal)
%     TMP = scFeatVal{i};
%     TMP(TMP<0) = 0;
%     scXVal(:,i) = mean(TMP,2);
    scXVal(:,i) = mean(scFeatVal{i},2);
%     scXVal(:,i) = mean(abs(scFeatVal{i}),2);
%     scXVal(:,i) = max(scFeatVal{i},[],2);

    lcXval(:,i) = locFeatVal{i};
    YVal(testDataClassLabel(i),i) = 1;
end
% scXVal = scXVal ./ repmat( sqrt(sum(scXVal.^2,1)), size(scXVal,1), 1 );

% A = scXTrain; B = scXTrain;
% A(A<0) = 0; B(B>0) = 0;
% scXTrain = [A; B;lcXTrain];
% 
% A = scXVal; B = scXVal;
% A(A<0) = 0; B(B>0) = 0;
% scXVal = [A; B;lcXval];

scXTrain = [scXTrain;lcXTrain]; % location statistics for shape
scXVal = [scXVal;lcXval]; %

% scXTrain = scXTrain ./ repmat( sqrt(sum(scXTrain.^2,1)), size(scXTrain,1), 1 );
% scXVal = scXVal ./ repmat( sqrt(sum(scXVal.^2,1)), size(scXVal,1), 1 );

%% svm on the sparse codes
addpath(genpath('../toolbox/libsvm-3.20/matlab'));

model = libsvm_svmtrain(trainClassLabel', scXTrain', '-s 0 -c 25 -t 0'); % linear 
[predLabelVal, accuracy, dec_values] = libsvm_svmpredict(testDataClassLabel(:), scXVal', model); % test the training data
predLabelVal = predLabelVal';
acc = mean(testDataClassLabel(:) == predLabelVal(:));
fprintf('\non validation set\n\taccuracy=%.4f (lambda1=%.4f, lambdaLoc=%.4f, T=%d)\n\n', acc, lambda1, lambdaLoc, T);

for categID = 1:numel(pollenName)
    a = find(testDataClassLabel==categID);
    accTMP = mean(testDataClassLabel(a) == predLabelVal(a));
    fprintf('\t\t%s acc:%.4f (#:%d)\n', pollenName{categID}, accTMP, length(a));
    
end

for categID = 1:numel(pollenName)
    a = find(testDataClassLabel==categID);
    b = find(predLabelVal==categID);
    fprintf('\ton validation Class%d #GT:%d #pred:%d\n', categID, length(a), length(b));
end

%% learning classifier from training set and test on validation set
% alphaList = 0.5:0.01:1;
alphaList = [0.01:0.0001:0.015];

accMatTrain = zeros(numel(pollenName)+1, length(alphaList));
accMatVal = zeros(numel(pollenName)+1, length(alphaList));

for alpha_i = 1:length(alphaList)
    alpha = alphaList(alpha_i);
    
    U = (scXTrain*scXTrain' + alpha*eye(size(scXTrain,1)))\scXTrain*YTrain';
    
    % accuracy on training set
    YhatTrain = U'*scXTrain;
    [~,predLabelTrain] = max(YhatTrain,[],1);
    
    acc = mean(trainClassLabel == predLabelTrain);
%     fprintf('on training set\n\taccuracy=%.4f (lambda=%.4f, T=%d)\n', acc, lambda, T);
    accMatTrain(end,alpha_i) = acc;
    for categID = 1:numel(pollenName)
        a = find(trainClassLabel==categID);
        accTMP = mean(trainClassLabel(a) == predLabelTrain(a));
%         fprintf('\t\t%s acc:%.4f (#:%d)\n', pollenName{categID}, accTMP, length(a));
        accMatTrain(categID,alpha_i) = accTMP;
    end
    
    % accuracy on validation set
    YhatVal = U'*scXVal;
    [~,predLabelVal] = max(YhatVal,[],1);
    
    % accuracy on validation set
    acc = mean(testDataClassLabel == predLabelVal);
%     fprintf('\non validation set\n\taccuracy=%.4f (lambda=%.4f, T=%d)\n', acc, lambda, T);
    accMatVal(end,alpha_i) = acc;
    
    for categID = 1:numel(pollenName)
        a = find(testDataClassLabel==categID);
        accTMP = mean(testDataClassLabel(a) == predLabelVal(a));
%         fprintf('\t\t%s acc:%.4f (#:%d)\n', pollenName{categID}, accTMP, length(a));
        accMatVal(categID,alpha_i) = accTMP;
    end
end

%disp(accMatVal)
[bestAcc, alphaIdx] = max(accMatVal(end,:));
fprintf('\nbest overall acc: %.4f, alpha=%.4f\n', bestAcc, alphaList(alphaIdx));
disp(accMatVal(:,alphaIdx))

for categID = 1:numel(pollenName)
    a = find(trainClassLabel==categID);
    b = find(predLabelTrain==categID);
    fprintf('on train Class%d #GT:%d #pred:%d\n', categID, length(a), length(b));
    
    a = find(testDataClassLabel==categID);
    b = find(predLabelVal==categID);
    fprintf('on validation Class%d #GT:%d #pred:%d\n', categID, length(a), length(b));
end

%% save results
save( ['./result4_nonSparse/valscFeatOnExemplar_layer' num2str(curLayer) '_' shapeType '_T' num2str(T) '_lambda1' num2str(lambda1) '_ResultACC' num2str(bestAcc) '_globalContrastNorm.mat'], ...
    'bestAcc', 'alphaIdx',...       
    'alphaList', 'accMatTrain', 'accMatVal', 'pollenName', 'scFeatTrain', 'trainClassLabel', 'scXTrain', 'scFeatVal', 'testDataClassLabel', 'scXVal', 'U');
