% This code demonstrates the proposed objective function for discovering
% the most informative receptive fields in images.
% Specifically, a synthetic dataset is generated, including three cluters
% to stand for three images, and the most correlated image regions or the
% common foreground objects are represented by the points lying in the
% intersection of the three clusters.
%
%
% Readers can directly run main_syntheticData.m to see the results, and the
% intermediate results/figures are automatically stored in 'figures' folder.
%
%
% For details, readers are suggested to refer to the following report:
%       Shu Kong, "Collaborative Receptive Field Learning", arXiv, 2014
%
% 
% The code is writen by
%           Shu Kong (Aimery)
%           aimerykong@gmail.com
%           http://www.ics.uci.edu/~skong2/
%           Dec. 2013, release version available on Jan. 13, 2014
%              updated: Nov. 16, 2015
%              updated: May 2, 2016

clear;
close all;
clc;

%% prelude
%
% some questions to consider:
%   1. sparse graph or non-sparse graph
%

%sigma = .2; % sigma controls the Gaussian kernel that convert pair-wise distance to similarity

K = 21;  % select K desired receptive fields

lambdaL = 0.6;
lambdaD = 1; % discriminative term (2) -- make the selected points more "class-pure" 
lambdaE = 20; % equal-size term (3) -- a very large lambda1 means that the points are being selected from class to class iteratively
lambdaB = 2; % balance term (4) -- make the size of clusters represented by selected point similar

ChoiceGraphBuilding = 'overall_knn'; % which way to build the graph, overall_knn or individual_knn or nothing
knn = 7; % knn-graph to measure RF similarity
eta = .4; % soft-threshold the similarity graph to remove some "dissimilar" pairs

distBound = 0.30;

pointNum = 140;

%% data generation and visualization colors
switch2Display = true;   % set this switch as 'true' to display the results

%strOptPoint = {'rd', 'md', 'b^', 'c^', 'ks', 'ys', 'gx'};
strOptPoint = {'rd', 'b^', 'ks', 'gx'};
strOptExemplar = {'ro', 'bo', 'ko'};
gtColor =  'gp';
OrgPoint = [0.1; 0.1];

X = cell(1, 3);
X{1} = sqrt(0.01) * randn(2, pointNum);
X{1}(1, :) = 1.1*X{1}(1, :);
theta = -50 / 180 * pi;
rotMat = [cos(theta), -sin(theta); sin(theta), cos(theta) ];
X{1} = rotMat*X{1};
X{1} = repmat( [0; 0], 1, size(X{1}, 2) ) + X{1};

X{2} = sqrt(0.01) * randn(2, pointNum);
X{2}(1, :) = 1.5*X{2}(1, :);
theta = 60 / 180 * pi;
rotMat = [cos(theta), -sin(theta); sin(theta), cos(theta) ];
X{2} = rotMat*X{2};
X{2} = repmat( [-0.05; -0.15], 1, size(X{2}, 2) ) + X{2};

X{3} = sqrt(0.01) * randn(2, pointNum);
X{3}(1, :) = 2.7*X{3}(1, :);
X{3}(2, :) = 2.0*X{3}(2, :);
theta = 30 / 180 * pi;
rotMat = [cos(theta), -sin(theta); sin(theta), cos(theta) ];
X{3} = rotMat*X{3};
X{3} = repmat( [0.1; 0.1], 1, size(X{3}, 2) ) + X{3};

%% visualize the original dataset
%%{
tempPoint = [];
ttt = 1;
if switch2Display
    figure; hold on; grid on;
    for c = 1:length(X)
        a = X{c} - repmat(OrgPoint, 1, size(X{c},2));
        a = sqrt(sum(a.^2, 1) );
        a = find(a > distBound);        
        tt(ttt) = plot( X{c}(1, :), X{c}(2, :), strOptPoint{c}  );
        ttt = ttt+1;
        
        tempPoint = [tempPoint  a(1)];
        
%         a = setdiff(1:size(X{c},2), a);
%         tt(ttt) = plot( X{c}(1, a), X{c}(2, a), strOptPoint{c}, 'linewidth', 4, 'MarkerSize', 4 );
%         ttt= ttt+1;
        %plot( X{c}(1, a), X{c}(2, a), gtColor  );
        
        tempPoint = [tempPoint  a(1)];
    end
    titleName = 'syntheticData';
    title( titleName );
    legend( [tt(1), tt(2), tt(3) ],...
        'class-1 point','class-2 point','class-3 point', 'Location', 'SouthEastoutside');
%     legend( [tt(1), tt(2), tt(3), tt(4), tt(5), tt(6), ],...
%         'class-1 background point','class-1 foreground point','class-2 background point','class-2 foreground point','class-3 foreground point','class-3 background point',...
%         'Location', 'SouthEast');
    hold off;
%     saveas(gcf, ['./figures/' titleName '.eps'], 'psc2');
%     saveas(gcf,['./figures/' titleName '.fig']);
end
%}

%% fetch the data, extracting RF candidates for all images from specific category
numImg = length(X);
%pointNum = pointNum + pointNumGT;
fprintf('fetch the data (receptive fields candidates)...\n');
imgIdx = zeros( 1, numImg*pointNum );
for i = 1:numImg
    imgIdx(1+pointNum*(i-1):i*pointNum ) = i;
end

%% create the similarity-graph (to maximize)
%{
% build a full graph that measures pairwise contrast/distance of all images
fprintf('build pairwise-distance graph...\n');
Xmat = cell2mat(X);
SimilarityGraph = Xmat'*Xmat; % maximize reward graph -- contrast
% clear X A;

Xsqr = diag(SimilarityGraph);
Xsqr = Xsqr(:);
SimilarityGraph = repmat(Xsqr, 1, size(SimilarityGraph,2)) ...
    +repmat(Xsqr', size(SimilarityGraph,2), 1) ...
    - 2*SimilarityGraph; % Euclidean distance, maximize sum of pairwise distance

SimilarityGraph_full= SimilarityGraph;

switch ChoiceGraphBuilding
    case 'overall_knn',
        for i = 1:size(SimilarityGraph, 2)
            costList = SimilarityGraph(i, :);
            [costListSorted, costIdx]= sort(costList, 'ascend'); % smaller distance
            SimilarityGraph(i, :) = 0;
            SimilarityGraph(i, costIdx(1:knn*1)) = costListSorted(1:knn*1);
        end
    case 'individual_knn',
        for i = 1:size(SimilarityGraph, 2)
            idx = imgIdx(i);
            tmp = SimilarityGraph(i, :);
            SimilarityGraph(i, :) = 0;
            
            for j = 1:imgIdx(end)
                %if j == idx
                    %continue;
%                     a = find(imgIdx==idx);
%                     costList = tmp(a);
%                     [costList costIdx]= sort(costList, 'ascend');
%                     SimilarityGraph( i, a(costIdx(1: round(knn ))) ) = costList( 1:round(knn ) );
                %else
                    a = find(imgIdx==j);
                    costList = tmp(a);
                    [costList costIdx]= sort(costList, 'ascend');
                    SimilarityGraph(i, a(costIdx(1:knn))) = costList(1:knn);
                %end
            end
            tmpVal = SimilarityGraph(i, :);
            tmpIdx = find(tmpVal~=0);
            tmpVal = tmpVal(tmpIdx);
            SimilarityGraph(i, find(SimilarityGraph(i,:) < eta*mean(tmpVal) ) ) = 0;
            
        end
    case 'nothing',
           SimilarityGraph = SimilarityGraph;
end

SimilarityGraph = (SimilarityGraph + SimilarityGraph') / 2;
%SimilarityGraph = SimilarityGraph ./ max(SimilarityGraph(:));
%SimilarityGraph = SimilarityGraph_backup; % also works!

% normalize the similarity graph
fprintf('normalize the similarity graph...\n');
%SimilarityGraph = SimilarityGraph ./ max(SimilarityGraph(:));
sigma = 2* mean(SimilarityGraph(:)); % adaptively selecting sigma

% convert it to the similarity graph
fprintf('convert into similarity graph...\n');
a = find(SimilarityGraph~=0);
SimilarityGraph(a) = exp(-SimilarityGraph(a) / sigma);


sigma2 = 2* mean(SimilarityGraph_full(:)); % adaptively selecting sigma
a = find(SimilarityGraph_full~=0);
SimilarityGraph_full(a) = exp(-SimilarityGraph_full(a) / sigma2);

%SimilarityGraph = sparse(SimilarityGraph);
%}

%%
Xmat = cell2mat(X);
[SimilarityGraph, SimilarityGraph_full]= graphConstruction(Xmat, imgIdx, ChoiceGraphBuilding, knn);
RFpositive = DiscriminativeExemplarSelection(SimilarityGraph, SimilarityGraph_full, SimilarityGraph_full, imgIdx, K, lambdaL, lambdaD, lambdaE, lambdaB);

%% visualize the selected RFs
% RFpositive; % find(SelectedFlag);
if switch2Display
    figure; hold on; grid on;
    ttt = 1;
    tt = [];
    for c = 1:3
        tt(ttt) = plot( X{c}(1, :), X{c}(2, :), strOptPoint{c}  );
        ttt = ttt+1;
    end
    IdxImg = imgIdx(RFpositive); % image index
    for i = 1:length(RFpositive)
        plot( Xmat(1, RFpositive(i)), Xmat(2, RFpositive(i)), strOptExemplar{IdxImg(i)}, 'linewidth', 3, 'MarkerSize', 11);        
    end
    
    a = imgIdx(RFpositive);
    b = find(a==1);
    b = b(1);
    b = RFpositive(b);
    tt(ttt) = plot( X{1}(1, b), X{1}(2, b), strOptExemplar{1}, 'linewidth', 3, 'MarkerSize', 11  );
    ttt = ttt+1;
    
    b = find(a==2);
    b = b(1);
    b = RFpositive(b) - pointNum;
    tt(ttt) = plot( X{2}(1, b), X{2}(2, b), strOptExemplar{2}, 'linewidth', 3, 'MarkerSize', 11  );
    ttt = ttt+1;
    
    b = find(a==3);
    b = b(1);
    b = RFpositive(b) - 2*pointNum;
    tt(ttt) = plot( X{3}(1, b), X{3}(2, b), strOptExemplar{3}, 'linewidth', 3, 'MarkerSize', 11  );
    ttt = ttt+1;
    
    titleName = 'selected exemplars';
    hold off; title(titleName);
    legend( [tt(1), tt(4), tt(2), tt(5), tt(3), tt(6), ],...
        'class-1 point', 'class-1 selected point', 'class-2 point', 'class-2 selected point', 'class-3 point','class-3 selected point',...
        'Location', 'SouthEastoutside');
    hold off;
%     saveas(gcf, ['./figures/' titleName '.eps'], 'psc2');
%     saveas(gcf,['./figures/' titleName '.fig'])
%     saveas(gcf,['./figures/' titleName '.jpg']);
end

%% visualize the marginal gain at each iteration
%{
if switch2Display
    az = 30;
    el = 15;
    for iter = 1:size(storeIteration, 2)
        ttt = 1;
        tt = [];
        curFigure = figure(3); hold on; grid on; 
        titleName = ['marginalGainIteration-' num2str(iter)];
        title( titleName );
        for c =1:3
            tt(ttt) = plot3( X{c}(1, :), X{c}(2, :), storeIteration(1+pointNum*(c-1):pointNum*c, iter), strOptPoint{c}  );
            ttt = ttt+1;
        end
        for i = 1:iter
           tt(ttt) = plot3( Xmat(1, RFpositive(i)), Xmat(2, RFpositive(i)), storeIteration( RFpositive(i), iter ), ...
                strOptExemplar{IdxImg(i)}, 'linewidth', 3, 'MarkerSize', 11);
            ttt = ttt+1;
        end
        view(az, el);
%         saveas(gcf, ['./figures/' titleName '.eps'], 'psc2');
%         saveas(gcf,['./figures/' titleName '.fig']);
        saveas(gcf,['./figures/' titleName '.jpg']);
       	close(curFigure)
    end
end
%}


