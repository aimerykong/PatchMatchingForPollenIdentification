clear
close all
clc

addpath(genpath('../fossilPollen/code/ompbox10'));
addpath(genpath('../fossilPollen/code/spams-matlab'));

%% setup
numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

locationScaleFactor = 0.001;

Tlist = 40:20:40;
lambdaList = 0.0005; % [0.00005, 0.0001, 0.0005, 0.001];

for Tid = 1:length(Tlist)
    for LambdaID = 1:length(lambdaList)
        T = Tlist(Tid); % nonzero elements
        lambda = lambdaList(LambdaID); % weight of location-aware penalty
        
        % databaseDIR = './database_thumbnail_CNNfeature';
        databaseDIR = './database_binaryMask_CNNfeature';
        shapeType = strfind(databaseDIR, '_');
        shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);
        
        load(['trval_layer' num2str(curLayer) '_' shapeType '.mat'])
        pollenName = {'critchfieldii', 'glauca', 'mariana'};
        
        trainDictLoc = trainDictLoc * locationScaleFactor;
        
        %% location-aware sparse coding for classification
        labelpred = zeros(1,numel(testData));
        errorList = zeros(numel(pollenName),numel(testData));
        count = 1;
        scFeat = cell(1,numel(testData));
%         WList = cell(1,numel(testData));
        for i_img = numel(testData):-1:1 % 1:numel(testData)
            
            fprintf('\t image-%d/%d ---', i_img, numel(testData));
            imFeat = testData{i_img};
            patchLoc = testDataLoc{i_img} * locationScaleFactor;
            
            W = zeros(size(trainDict,2), size(patchLoc,2));
            for i = 1 : size(patchLoc,2)
                a = bsxfun(@minus, trainDictLoc, patchLoc(:,i) );
                a = sum(a.^2);
                a = sqrt(a(:));
                W(:, i) = a;
            end
            
            param.L = T; % not more than param.L non-zeros coefficients (default: min(size(D,1),size(D,2)))
            param.lambda = lambda;
            param.numThreads = -1; % number of processors/cores to use; the default choice is -1 and uses all the cores of the machine
            param.mode = 2; % penalized formulation
            
            A = mexLassoWeighted(imFeat, trainDict, W, param);
            
            tmpErrorList = ones(numel(pollenName),1)*Inf;
            for categID = 1:numel(pollenName)
                a = find( DictClassLabel == categID );
                err = imFeat-trainDict(:,a)*A(a,:);
                err = sum(err(:).^2);
                tmpErrorList(categID) = err;
            end
            errorList(:,i_img) = tmpErrorList(:);
            [valMIN, idxMIN] = min(tmpErrorList);
            labelpred(i_img) = idxMIN;
            if idxMIN == testDataClassLabel(i_img)
                flag = 'correct';
            else
                flag = 'wrong';
            end
            
%             WList{i_img} = W;
            scFeat{i_img} = A;
            
            fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));           
        end
        
        %% accuracy
        acc = mean(testDataClassLabel == labelpred);
        fprintf('\naccuracy=%.4f (lambda=%.4f, T=%d)\n\n', acc, lambda, T);
        
        for categID = 1:numel(pollenName)
            a = find(testDataClassLabel==categID);
            accTMP = mean(testDataClassLabel(a) == labelpred(a));
            fprintf('%s acc:%.4f\n', pollenName{categID}, accTMP);
        end
        
        save( ['trval_layer' num2str(curLayer) '_' shapeType '_T' num2str(T) '_lambda' num2str(lambda) '_ResultACC' num2str(acc) '.mat'], ...
            'errorList', 'labelpred', 'scFeat', 'testDataClassLabel', 'pollenName');
        
    end
end

















