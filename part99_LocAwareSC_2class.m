% Only focus on two pollen - glauca and mariana for easy analysis
% 
%
%
%
%
%

clear
close all
clc

addpath(genpath('../fossilPollen/code/ompbox10'));
addpath(genpath('../fossilPollen/code/spams-matlab'));

%% setup
numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

Tlist = 10:20:10;
lambdaList = 0.00005; % [0.00005, 0.0001, 0.0005, 0.001];

locationScaleFactor = 0.001;

for Tid = 1:length(Tlist)
    for LambdaID = 1:length(lambdaList)
        T = Tlist(Tid); % nonzero elements
        lambda = lambdaList(LambdaID); % weight of location-aware penalty
        
        % databaseDIR = './database_thumbnail_CNNfeature';
        databaseDIR = './database_binaryMask_CNNfeature';
        shapeType = strfind(databaseDIR, '_');
        shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);
        
        load(['trval_layer' num2str(curLayer) '_' shapeType '.mat'])
        pollenName = { 'glauca', 'mariana'};        
%         pollenName = {'critchfieldii', 'glauca', 'mariana'};
        tmp = find(DictClassLabel~=1);
        DictClassLabel = DictClassLabel(tmp)-1;
        trainDict = trainDict(:,tmp);
        trainDictLoc = trainDictLoc(:,tmp);

        %% location-aware sparse coding for classification
        labelpred = [];
        errorList = [];
        count = 1;
        
        a = find(testDataClassLabel~=1);
        testDataClassLabel = testDataClassLabel(a)-1;
        testData = testData(a);
        for i_img = 1:numel(testData)
            fprintf('\t image-%d/%d ---', i_img, numel(testData));
            imFeat = testData{i_img};
            patchLoc = testDataLoc{i_img};
            
            %% cross-class spare coding
            %%{
            SCtype = 'cross';
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
            %}
            %% within-class sparse coding
            %{
            SCtype = 'within';
            
            param.L = T; % not more than param.L non-zeros coefficients (default: min(size(D,1),size(D,2)))
            param.lambda = lambda;
            param.numThreads = -1; % number of processors/cores to use; the default choice is -1 and uses all the cores of the machine
            param.mode = 2; % penalized formulation
            
            tmpErrorList = ones(numel(pollenName),1)*Inf;
            for categID = 1:numel(pollenName)
                a = find( DictClassLabel == categID );
                
                subtrainDict = trainDict(:,a);
                subtrainDictLoc = trainDictLoc(a);
                W = zeros(length(a), size(patchLoc,2));
                for i = 1 : size(patchLoc,2)
                    a = bsxfun(@minus, subtrainDictLoc, patchLoc(:,i) );
                    a = sum(a.^2);
                    a = sqrt(a(:));
                    W(:, i) = a;
                end
                A = mexLassoWeighted(imFeat, subtrainDict, W, param);
                err = imFeat-subtrainDict*A;
                err = sum(err(:).^2);
                tmpErrorList(categID) = err;
            end                   
            %}
            %% recording results
            errorList = [errorList, tmpErrorList(:)];
            [valMIN, idxMIN] = min(tmpErrorList);
            labelpred = [labelpred, idxMIN];
            if idxMIN == testDataClassLabel(i_img)
                flag = 'correct';
            else
                flag = 'wrong';
            end
            fprintf('GT:%d pred:%d %s (%.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2));
%             fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));
            
        end
        
        %% accuracy
        acc = mean(testDataClassLabel == labelpred);
        fprintf('\naccuracy=%.4f (lambda=%.4f, T=%d)\n\n', acc, lambda, T);
        
        accTensor = [];
        numTensor = [];
        for categID = 1:numel(pollenName)
            a = find(testDataClassLabel==categID);
            accTMP = mean(testDataClassLabel(a) == labelpred(a));
            fprintf('%s acc:%.4f (#image:%d)\n', pollenName{categID}, accTMP, length(a));
            accTensor(end+1) = accTMP;
            numTensor(end+1) = length(a);
        end
        
        save( ['junk_trval_layer' num2str(curLayer) '_' shapeType '_T' num2str(T) '_lambda' num2str(lambda) '_SCtype_' SCtype '_ResultACC' num2str(acc) '.mat'], ...
            'errorList', 'labelpred', 'accTensor', 'numTensor');        
    end
end















