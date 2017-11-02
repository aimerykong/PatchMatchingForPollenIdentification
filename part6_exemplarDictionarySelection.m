clear;
close all;
clc;

%% prelude for exemplar selection
%sigma = .2; % sigma controls the Gaussian kernel that convert pair-wise distance to similarity

K = 600;  % select K desired receptive fields

lambdaL = 0.01; % term (1) location coverage power
lambdaD = 0.01; % discriminative term (2) -- make the selected points more "class-pure"
lambdaE = 20; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
lambdaB = 2; % balance term (4) -- make the size of clusters represented by selected point similar

ChoiceGraphBuilding = 'overall_knn'; % which way to build the graph, overall_knn or individual_knn or nothing
knn = 11; % knn-graph to measure RF similarity
eta = .4; % soft-threshold the similarity graph to remove some "dissimilar" pairs

distBound = 0.30;

pointNum = 140;

%% dataset setup
numTrain = 65; % the number of training data in each class
curLayer = 24; %[22, 23, 24, 25]

locationScaleFactor = 0.001;

Tlist = 20:20:20;
lambdaList = 0.0005; % [0.00005, 0.0001, 0.0005, 0.001];

for Tid = 1:length(Tlist)
    for LambdaID = 1:length(lambdaList)
        T = Tlist(Tid); % nonzero elements
        lambda = lambdaList(LambdaID); % weight of location-aware penalty
        
        %% load dataset
        % databaseDIR = './database_thumbnail_CNNfeature';
        databaseDIR = './database_binaryMask_CNNfeature_globalContrastNorm';
        shapeType = strfind(databaseDIR, '_');
        shapeType = databaseDIR(shapeType(1)+1:shapeType(2)-1);
        
        fprintf('load dataset...\n');
        load(['trval_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat']);
        pollenName = {'critchfieldii', 'glauca', 'mariana'};
        
        %% graph construction
        fprintf('graph construction/loading...\n');
        trainDictLoc = trainDictLoc * locationScaleFactor; % rescale locations for later comparable use
        
        if exist(['SimilarityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'file') ...
                && exist(['SimilarityGraphFull_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'file')
            load(['SimilarityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat']);
            load(['SimilarityGraphFull_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'])
        else
            [SimilarityGraph, FullConnectGraph] = graphConstruction(trainDict, DictClassLabel, 'overall_knn', knn);
            save(['SimilarityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'SimilarityGraph');
            save(['SimilarityGraphFull_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'FullConnectGraph');
        end
        if exist(['LocationAffinityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'file')
            load(['LocationAffinityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat']);
        else
            [~,LocationAffinityGraph] = graphConstruction(trainDictLoc, DictClassLabel);
            save(['LocationAffinityGraph_layer' num2str(curLayer) '_' shapeType '_globalContrastNorm.mat'], 'LocationAffinityGraph');
        end
        
        %% exemplar selection
        dictName = ['exemplarDict_K' num2str(K) 'L' num2str(lambdaL) '_D' num2str(lambdaD) '_E' num2str(lambdaE) '_B' num2str(lambdaB) '_globalContrastNorm.mat' ];
        fprintf('exemplar selection...\n');
        if exist(dictName, 'file')
            load(dictName);
        else
            tic
            exemplarIndex = DiscriminativeExemplarSelection_LazyGreedy(SimilarityGraph, LocationAffinityGraph, FullConnectGraph, DictClassLabel, K, lambdaL, lambdaD, lambdaE, lambdaB);
            toc
            save(dictName, 'exemplarIndex');
        end
        
        %% visualize exemplars' locations        
        dict = trainDict(:,exemplarIndex);
        dictLoc = trainDictLoc(:,exemplarIndex);
        dictLabel = DictClassLabel(exemplarIndex);
        
        img = imread('./database_binaryMask_canonicalShape/glauca fossil/296_Nelson Lake NE906 Pos 5_K1.jpg');
        sz = size(img);
        
        dictLoc = bsxfun(@plus, dictLoc/locationScaleFactor, sz(:)/2);
        
        imshow(img);
        hold on;
        for i = 1:size(dictLoc,2)
            rectangle('position',[dictLoc(2,i)-26 dictLoc(1,i)-26 52 52], 'EdgeColor','r')
        end
        hold off;
        
        %{
        %%  quatitatively test exemplars for prediction on the validation set by location-aware sparse coding
        dict = trainDict(:,exemplarIndex);
        dictLoc = trainDictLoc(:,exemplarIndex);
        dictLabel = DictClassLabel(exemplarIndex);
        
        labelpred = zeros(1,numel(testData));
        errorList = zeros(numel(pollenName),numel(testData));
        count = 1;
        scFeat = cell(1,numel(testData));
        for i_img = numel(testData):-1:1 % 1:numel(testData)
            
            fprintf('\t image-%d/%d ---', i_img, numel(testData));
            imFeat = testData{i_img};
            patchLoc = testDataLoc{i_img} * locationScaleFactor;
            
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
            errorList(:,i_img) = tmpErrorList(:);
            [valMIN, idxMIN] = min(tmpErrorList);
            labelpred(i_img) = idxMIN;
            if idxMIN == testDataClassLabel(i_img)
                flag = 'correct';
            else
                flag = 'wrong';
            end
            
            scFeat{i_img} = A;
            
            fprintf('GT:%d pred:%d %s (%.4f, %.4f, %.4f)\n', testDataClassLabel(i_img), idxMIN, flag, tmpErrorList(1), tmpErrorList(2), tmpErrorList(3));
        end
        
        % accuracy
        acc = mean(testDataClassLabel == labelpred);
        fprintf('\naccuracy=%.4f (lambda=%.4f, T=%d)\n\n', acc, lambda, T);
        
        for categID = 1:numel(pollenName)
            a = find(testDataClassLabel==categID);
            accTMP = mean(testDataClassLabel(a) == labelpred(a));
            fprintf('%s acc:%.4f\n', pollenName{categID}, accTMP);
        end
        save( ['valscFeatOnExemplar_layer' num2str(curLayer) '_' shapeType '_T' num2str(T) '_lambda' num2str(lambda) '_ResultACC' num2str(acc) '.mat'], ...
            'errorList', 'labelpred', 'scFeat', 'testDataClassLabel', 'pollenName');
        %}
    end
end





