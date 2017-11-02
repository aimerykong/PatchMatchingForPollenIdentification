function [graph, graph_full]= graphConstruction(Xmat, label, ChoiceGraphBuilding, knn)
%
%
%
%
%
% 
% Copyright (c) 2015 Shu Kong <skong2@uci.edu>
%           http://www.ics.uci.edu/~skong2/
%              updated: Nov. 17, 2015

if nargin<3
    ChoiceGraphBuilding = 'overall_knn';
end

if nargin<4
    knn = Inf;
end

%% create the similarity-graph (to maximize)
% build a full graph that measures pairwise contrast/distance of all images
fprintf('build pairwise-distance graph...\n');

graph = Xmat'*Xmat; % maximize reward graph -- contrast
% clear X A;

Xsqr = diag(graph);
Xsqr = Xsqr(:);
graph = repmat(Xsqr, 1, size(graph,2)) ...
    +repmat(Xsqr', size(graph,2), 1) ...
    - 2*graph; % Euclidean distance, maximize sum of pairwise distance

graph_full = graph;

switch ChoiceGraphBuilding
    case 'overall_knn',
        for i = 1:size(graph, 2)
            if knn ~= Inf
                costList = graph(i, :);
                [costListSorted, costIdx]= sort(costList, 'ascend'); % smaller distance
                graph(i, :) = 0;
                graph(i, costIdx(1:knn*1)) = costListSorted(1:knn*1);
            end
        end
    case 'individual_knn',
        for i = 1:size(graph, 2)
            idx = label(i);
            tmp = graph(i, :);
            graph(i, :) = 0;
            
            for j = 1:label(end)
                %if j == idx
                    %continue;
%                     a = find(label==idx);
%                     costList = tmp(a);
%                     [costList costIdx]= sort(costList, 'ascend');
%                     graph( i, a(costIdx(1: round(knn ))) ) = costList( 1:round(knn ) );
                %else
                    a = find(label==j);
                    costList = tmp(a);
                    [costList costIdx]= sort(costList, 'ascend');
                    graph(i, a(costIdx(1:knn))) = costList(1:knn);
                %end
            end
            tmpVal = graph(i, :);
            tmpIdx = find(tmpVal~=0);
            tmpVal = tmpVal(tmpIdx);
            graph(i, find(graph(i,:) < eta*mean(tmpVal) ) ) = 0;
            
        end
    case 'nothing',
           graph = graph;
end

graph = (graph + graph') / 2;
%graph = graph ./ max(graph(:));

% normalize the similarity graph
fprintf('normalize the similarity graph...\n');
%graph = graph ./ max(graph(:));
sigma = 2* mean(graph(:)); % adaptively selecting sigma

% convert it to the similarity graph
fprintf('convert into similarity graph...\n');
a = find(graph~=0);
graph(a) = exp(-graph(a) / sigma);


sigma2 = 2* mean(graph_full(:)); % adaptively selecting sigma
a = find(graph_full~=0);
graph_full(a) = exp(-graph_full(a) / sigma2);


