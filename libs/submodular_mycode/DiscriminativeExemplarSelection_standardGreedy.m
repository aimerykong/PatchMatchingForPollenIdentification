function [exemplarIndex] = DiscriminativeExemplarSelection_standardGreedy(SimilarityGraph, LocationAffinityGraph, FullConnectGraph, labelList, K, lambdaL, lambdaD, lambdaE, lambdaB)
% This code is to select a set of exemplars discriminatively from a
% dataset (training set) by solving a well-defined submodular function.
%
% The submodular function consists of several terms
%   term1 [representational power] -- selected exemplar should be
%       representative;
%   term2 [location Coverage] -- selected exemplar should cover all
%       possible locations in the image;
%   term3 [discriminative power] -- selected exemplar should be
%       discriminative in terms of classification;
%   term4 [class-specific dictionary size equality ] -- selected exemplars
%       that mainly represent data from a specific class should be of
%       purity;
%   term5 [cluster size equality] -- exemplars which data from a specific
%       class are belonging to should cluster together, and the cluster
%       size should be similar for all clusters.
%
% In our pollen project, there are two graphs that record feature
% similarity and physical locations, respectively, such that the selected
% patches as exemplars should not only represent training data, but also
% spread over all possible locations in the image.
%
% In addition, the third graph is also needed to calculate discirmination
% power. This graph is not sparse, so that all training data can find its
% own closest exemplar.
%
%
% input
%   similarityGraph         -- a graph of NxN for term1
%   locationAffinityGraph   -- a graph of NxN for term2
%   fullConnectGraph        -- a graph of NxN for term3
%   labelList               -- record the class label of each data point
%   K                       -- the number of exemplars to select
%   lambdaL                 -- control importance of term 2 (location coverage)
%   lambdaD                 -- control importance of term 3 (discriminative power) 
%   lambdaE                 -- control importance of term 4 (dictionary size) 
%   lambdaB                 -- control importance of term 5 (cluster size) 
%
%
% output
%   exemplarIndex              -- the index of selected exemplars
%
%
% reference:
%       Shu Kong et al., "Collaborative Receptive Field Learning", arXiv, 2014
%       Zhuolin Jinag et al., "Submodular Dictionary Learning for Sparse Coding", CVPR, 2012
%       Ming-Yu Liu et al., "Entropy rate superpixel segmentation", CVPR, 2011
%
% 
% Copyright (c) 2015 Shu Kong <skong2@uci.edu>
%           http://www.ics.uci.edu/~skong2/
%              updated: Nov. 17, 2015

%% greedily optimize the submodular function and select the most desired RF's
fprintf('seek for %d most desirable exemplars...\n', K);

nClass = length(unique(labelList)); % number of classes

SelectedFlag = zeros( 1, size(SimilarityGraph, 2) );
CountNumRFperImg = zeros(1, nClass);
exemplarIndex = zeros(1, K);

nClass = length(unique(labelList));

storeIteration = zeros(size(SimilarityGraph,1), K);
flag = true;
k = 1;
while flag  
    %% current function return    
    %% first term -- representation power
    SelectAlready = find(SelectedFlag == 1);    
    if isempty(SelectAlready)
        term1 = 0;
    else
        tmpRows = SimilarityGraph(:, SelectAlready);
        tmpSimilarity = max(tmpRows,[],2);
        term1 = sum(tmpSimilarity);
    end
    Hcur = term1; % first term
    
    %% second term -- location coverage power
    SelectAlready = find(SelectedFlag == 1);    
    if isempty(SelectAlready)
        term2 = 0;
    else
        tmpRows = LocationAffinityGraph(:, SelectAlready);
        tmpSimilarity = max(tmpRows,[],2);
        term2 = sum(tmpSimilarity);
    end
    Hcur = Hcur + lambdaL*term2; % first term
    
    %% third term -- discriminative term
    tmpColumns = SimilarityGraph(:, SelectAlready);
    classLabels = repmat(labelList(:), 1, size(tmpColumns,2));
    classLabels(classLabels==0) = 0;
    count = zeros(nClass, size(tmpColumns,2) );
    for c = 1:nClass
        count(c,:) = sum(classLabels==c,1);
    end
    term3 = max(count,[],1);
    term3 = sum(term3);
    term3 = term3/nClass - length(SelectAlready) ;
    Hcur = Hcur + lambdaD*term3;
    
    %% plus equal-size term (fourth one)
    term4 = sum(log(1+CountNumRFperImg));
    Hcur = Hcur + lambdaE*term4;        
    
    %% plus balance cluster size term (fifth term)
    term5 = 0;
    tmpColumns = FullConnectGraph(:, SelectAlready);  
    [~, simCluster] = max(tmpColumns,[],2);
    for i = 1:length(SelectAlready)
        tmp = length(find(simCluster==i));
        if tmp == 0
            tmp = 0;
        else
            tmp = tmp/size(simCluster,1);
            tmp = -tmp*log(tmp);
        end
        term5 = term5 + tmp;
    end
    term5 = term5  - length(SelectAlready);
    Hcur = Hcur + lambdaB*term5; 
  
    %% calculate the benefit gain
    SelectNotYet = find(SelectedFlag ~= 1);
    predictH = zeros(1, length(SelectNotYet));
    PredictSelectSet = [repmat( SelectAlready(:), 1, length(SelectNotYet) ); SelectNotYet];
    
    %% first term
    term1List = zeros(1, length(SelectNotYet));
    for i = 1:length(SelectNotYet)
        tmpRows = SimilarityGraph(:, PredictSelectSet(:,i));
        tmpSimilarity = max(tmpRows,[],2);
        term1 = sum(tmpSimilarity);
        term1List(i) = term1;
    end
    predictH = term1List; % first term    
    
    %% second term
    term2List = zeros(1, length(SelectNotYet));
    for i = 1:length(SelectNotYet)
        tmpRows = LocationAffinityGraph(:, PredictSelectSet(:,i));
        tmpSimilarity = max(tmpRows,[],2);
        term2 = sum(tmpSimilarity);
        term2List(i) = term2;
    end
    predictH = predictH + lambdaL*term2List; % second term
    
    %% third term -- discriminative term   
    term3ListTMP = zeros(1, length(SelectNotYet));
    for j = 1:length(SelectNotYet)
        tmpColumns = SimilarityGraph(:, PredictSelectSet(:,j));
        classLabels = repmat(labelList(:), 1, size(tmpColumns,2));
        classLabels(classLabels==0) = 0;
        count = zeros(nClass, size(tmpColumns,2) );
        for c = 1:nClass
            count(c,:) = sum(classLabels==c);
        end
        term3 = max(count,[],1);
        term3 = sum(term3);        
        term3 = term3/nClass - length(PredictSelectSet(:,j)) ;
        term3ListTMP(j) = term3;
    end
    predictH = predictH + lambdaD*term3ListTMP;
    term3List = term3ListTMP;    
        
    %% plus balance term (fourth one)
    g_balance = repmat( CountNumRFperImg(:), 1, size(predictH,2) );
    a = ( 0:length(SelectNotYet)-1 ) * size( g_balance, 1 );
    g_balance(a+labelList(SelectNotYet)) = g_balance(a+labelList(SelectNotYet)) + 1;
    term4List = sum( log(1+g_balance), 1 );
    predictH = predictH + lambdaE*term4List;

    %% plus balance cluster size term (fifth term)
    term5List = zeros(1, length(SelectNotYet));
    for j = 1:length(SelectNotYet)
        tmpColumns = FullConnectGraph(:, PredictSelectSet(:,j));
        [~, simCluster] = max(tmpColumns,[],2);
        term5 = 0;
        for i = 1:length(PredictSelectSet(:,j))
            tmp = length(find(simCluster==i));
            if tmp == 0
%                 fprintf('!\n');
                tmp = 0;
            else
                tmp = tmp/size(simCluster,1);
                tmp = -tmp*log(tmp);
            end
            term5 = term5 + tmp;
        end
        term5 = term5  - length(PredictSelectSet(:,j));
        term5List(j) = term5;
    end
    predictH = predictH + lambdaB*term5List;
        
%     fprintf('\titer-%d, [term1=%.4f(rep) term2=%.4f(loc) term3=%.4f(disc) \n\t\tterm4=%.4f(dicSize) term5=%.4f(clusterSize)]\n', k, max(term1List), lambdaL*max(term2List), lambdaD*max(term3List), max(lambdaE*term4List), max(lambdaB*term5List) );
    
    %% calculate the marginal gain and find the most desirable exemplar
    predictGain = predictH(:) - Hcur; % calculate info gain by openning a specific facility
    [y, idx] = max(predictGain); % find the most informative facility which brings the most gain
    storeIteration(SelectNotYet, k) = predictGain;
%     fprintf('\t\tmax gain=%.4f\n\n', y );
    
    if k >= K % jump out loop when no gains or fixed iterations exceeded
%     if y <= 0 || k >= K % jump out loop when no gains or fixed iterations exceeded
        SelectedFlag(SelectNotYet(idx)) = 1; % select the most informative one
        CountNumRFperImg(labelList(SelectNotYet(idx))) = CountNumRFperImg(labelList(SelectNotYet(idx))) + 1;
        exemplarIndex(k) = SelectNotYet(idx);
        
        flag = false;
        break;
    end

    SelectedFlag(SelectNotYet(idx)) = 1; % select the most informative one
    CountNumRFperImg(labelList(SelectNotYet(idx))) = CountNumRFperImg(labelList(SelectNotYet(idx))) + 1;
    exemplarIndex(k) = SelectNotYet(idx);
    
    k = k + 1;
end
