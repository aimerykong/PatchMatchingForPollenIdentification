function [exemplarIndex] = DiscriminativeExemplarSelection_LazyGreedy(SimilarityGraph, LocationAffinityGraph, FullConnectGraph, labelList, K, lambdaL, lambdaD, lambdaE, lambdaB)
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

%% configuration
fprintf('seek for %d most desirable exemplars...\n', K);

nClass = length(unique(labelList)); % number of classes

SelectedFlag = zeros( 1, size(SimilarityGraph, 2) );
CountNumRFperImg = zeros(1, nClass);
exemplarIndex = zeros(1, K);

nClass = length(unique(labelList));

storeIteration = zeros(size(SimilarityGraph,1), K);
flag = true;

%% initialize gain list
gainList = zeros(1,size(SimilarityGraph,1));
for i = 1:length(gainList)
    %% term1
    tmpRows = SimilarityGraph(:, i);
    tmpSimilarity = max(tmpRows,[],2);
    term1 = sum(tmpSimilarity);
    Hcur = term1; % first term    
    %% term2
    tmpRows = LocationAffinityGraph(:, i);
    tmpSimilarity = max(tmpRows,[],2);
    term2 = sum(tmpSimilarity);    
    Hcur = Hcur + lambdaL*term2; % second term
    %% term3
    tmpColumns = SimilarityGraph(:, i);
    classLabels = repmat(labelList(:), 1, size(tmpColumns,2));
    classLabels(classLabels==0) = 0;
    count = zeros(nClass, size(tmpColumns,2) );
    for c = 1:nClass
        count(c,:) = sum(classLabels==c,1);
    end
    term3 = max(count,[],1);
    term3 = sum(term3);
    term3 = term3/nClass - 1;
    Hcur = Hcur + lambdaD*term3;
    %% term4
    term4 = sum(log(1+1));
    Hcur = Hcur + lambdaE*term4;        
    %% term5    
    term5 = 0;
    tmpColumns = FullConnectGraph(:, i);  
    [~, simCluster] = max(tmpColumns,[],2);
    for k = 1:length(i)
        tmp = length(find(simCluster==k));
        if tmp == 0
            tmp = 0;
        else
            tmp = tmp/size(simCluster,1);
            tmp = -tmp*log(tmp);
        end
        term5 = term5 + tmp;
    end
    term5 = term5  - length(i);
    Hcur = Hcur + lambdaB*term5; 
    %%
    gainList(i) = Hcur;
end

%% maximize the submodular function in lazy greedy way
k=0;
preSelectID = -1;
curF = 0;
while flag  
    if k>=K
        break;
    end
    SelectNotYet = find(SelectedFlag ~= 1);
    [Delta_e, ei] = max(gainList(SelectNotYet));
    if ei == preSelectID
        delta = Delta_e;
    else
        
        %% new F value if adding this point to the exemplar list
        % first term -- representation power
        SelectedFlagNew = SelectedFlag;
        SelectedFlagNew(SelectNotYet(ei)) = 1;
        SelectAlreadyNew = find(SelectedFlagNew == 1);
        if isempty(SelectAlreadyNew)
            term1 = 0;
        else
            tmpRows = SimilarityGraph(:, SelectAlreadyNew);
            tmpSimilarity = max(tmpRows,[],2);
            term1 = sum(tmpSimilarity);
        end
        newF = term1; % first term
        
        % second term -- location coverage power
        SelectAlreadyNew = find(SelectedFlagNew == 1);
        if isempty(SelectAlreadyNew)
            term2 = 0;
        else
            tmpRows = LocationAffinityGraph(:, SelectAlreadyNew);
            tmpSimilarity = max(tmpRows,[],2);
            term2 = sum(tmpSimilarity);
        end
        newF = newF + lambdaL*term2; % first term
        
        % third term -- discriminative term
        tmpColumns = SimilarityGraph(:, SelectAlreadyNew);
        classLabels = repmat(labelList(:), 1, size(tmpColumns,2));
        classLabels(classLabels==0) = 0;
        count = zeros(nClass, size(tmpColumns,2) );
        for c = 1:nClass
            count(c,:) = sum(classLabels==c,1);
        end
        term3 = max(count,[],1);
        term3 = sum(term3);
        term3 = term3/nClass - length(SelectAlreadyNew) ;
        newF = newF + lambdaD*term3;
        
        % plus equal-size term (fourth one)
        term4 = sum(log(1+CountNumRFperImg));
        newF = newF + lambdaE*term4;
        
        % plus balance cluster size term (fifth term)
        term5 = 0;
        tmpColumns = FullConnectGraph(:, SelectAlreadyNew);
        [~, simCluster] = max(tmpColumns,[],2);
        for i = 1:length(SelectAlreadyNew)
            tmp = length(find(simCluster==i));
            if tmp == 0
                tmp = 0;
            else
                tmp = tmp/size(simCluster,1);
                tmp = -tmp*log(tmp);
            end
            term5 = term5 + tmp;
        end
        term5 = term5  - length(SelectAlreadyNew);
        newF = newF + lambdaB*term5;
        
        delta = newF-curF; %[!!!! perform gain !!!!!]
        gainList(SelectNotYet(ei)) = delta;
    end
    
    %% checking and updating
    if delta >= max(gainList(SelectNotYet))
        if delta <= 0
            break;
        else            
            k = k + 1;
            SelectedFlag(SelectNotYet(ei)) = 1;
            gainList(SelectNotYet(ei)) = 0;
            exemplarIndex(k) = SelectNotYet(ei);  
            %% current F value
            % first term -- representation power
            SelectAlready = find(SelectedFlag == 1);
            if isempty(SelectAlready)
                term1 = 0;
            else
                tmpRows = SimilarityGraph(:, SelectAlready);
                tmpSimilarity = max(tmpRows,[],2);
                term1 = sum(tmpSimilarity);
            end
            curF = term1; % first term
            
            % second term -- location coverage power
            SelectAlready = find(SelectedFlag == 1);
            if isempty(SelectAlready)
                term2 = 0;
            else
                tmpRows = LocationAffinityGraph(:, SelectAlready);
                tmpSimilarity = max(tmpRows,[],2);
                term2 = sum(tmpSimilarity);
            end
            curF = curF + lambdaL*term2; % first term
            
            % third term -- discriminative term
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
            curF = curF + lambdaD*term3;
            
            % plus equal-size term (fourth one)
            term4 = sum(log(1+CountNumRFperImg));
            curF = curF + lambdaE*term4;
            
            % plus balance cluster size term (fifth term)
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
            curF = curF + lambdaB*term5;
        end
    end
end
