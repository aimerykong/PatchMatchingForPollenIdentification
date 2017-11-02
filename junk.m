clc;

scXTrain = zeros( size(scFeatTrain{1},1), length(scFeatTrain));
YTrain = zeros(numel(pollenName), length(scFeatTrain));
locTrain = zeros(numel(pollenName), length(scFeatTrain));
lcXTrain = zeros(size(locFeatTrain{1},1), length(scFeatTrain));
 
alpha = 1;


for i = 1:length(scFeatTrain)
    A = scFeatTrain{i};
    a = scFeatTrainLoc{i};
    b = scFeatTrainValidIndex{i};
    a = sqrt(sum(a(:,scFeatTrainValidIndex{i}).^2,1));
    thresh = mean(a);
% %     thresh = 0.1;
    inner = find(a<thresh); 
%     if isempty(inner) 
%         inner = 1:length(a); 
%     end
%     outer = find(a>=thresh);  
%     if isempty(outer) 
%         outer = 1:length(a); 
%     end
%     
% %     inner = find(a>-1);
% %     outer = find(a>-1);    
%     scXTrain(:,i) = [ alpha*mean(A(:,inner),2); (1-alpha)*mean(A(:,outer),2) ];
    scXTrain(:,i) = mean(A(:,inner),2);
%     scXTrain(:,i) = mean(scFeatTrain{i},2);
    %scXTrain(:,i) = max(scFeatTrain{i},[],2);
%     scXTrain(:,i) = [mean(scFeatTrain{i},2); max(scFeatTrain{i},[],2)];
    
    lcXTrain(:,i) = locFeatTrain{i};
    YTrain(trainClassLabel(i),i) = 1;
end
% scXTrain = scXTrain ./ repmat( sqrt(sum(scXTrain.^2,1)), size(scXTrain,1), 1 );

scXVal = zeros( size(scFeatVal{1},1), length(scFeatVal));
YVal = zeros(numel(pollenName), length(scFeatVal));
lcXval = zeros(size(locFeatVal{1},1), length(scFeatVal));
for i = 1:length(scFeatVal)    
    A = scFeatVal{i};
    a = sqrt(sum(scFeatTestLoc{i}.^2,1));
    thresh = mean(a);
% %     thresh = 0.1;
    inner = find(a<thresh); 
%     if isempty(inner) 
%         inner = 1:length(a); 
%     end
%     outer = find(a>=thresh);  
%     if isempty(outer) 
%         outer = 1:length(a); 
%     end
% %     inner = find(a>-1);
% %     outer = find(a>-1); 
%     scXVal(:,i) = [ alpha*mean(A(:,inner),2); (1-alpha)*mean(A(:,outer),2) ];
    scXVal(:,i) = mean(A(:,inner),2);
%     scXVal(:,i) = mean(scFeatVal{i},2);
%     scXVal(:,i) = max(scFeatVal{i},[],2);
%     scXVal(:,i) = [mean(scFeatVal{i},2); max(scFeatVal{i},[],2)];

    lcXval(:,i) = locFeatVal{i};
    YVal(testDataClassLabel(i),i) = 1;
end
% scXVal = scXVal ./ repmat( sqrt(sum(scXVal.^2,1)), size(scXVal,1), 1 );

scXTrain = [scXTrain; 0.92*lcXTrain]; % location statistics for shape
scXVal = [scXVal; 0.92*lcXval]; %

% scXTrain = scXTrain ./ repmat( sqrt(sum(scXTrain.^2,1)), size(scXTrain,1), 1 );
% scXVal = scXVal ./ repmat( sqrt(sum(scXVal.^2,1)), size(scXVal,1), 1 );

%% svm on the sparse codes
addpath(genpath('../toolbox/libsvm-3.20/matlab'));

ccList = [];
accList = [];
for cc = 30:1:60
    model = libsvm_svmtrain(trainClassLabel', scXTrain', sprintf('-s 0 -c %.2f -t 0', cc)); % linear
    [predLabelVal, accuracy, dec_values] = libsvm_svmpredict(testDataClassLabel(:), scXVal', model); % test the training data
    predLabelVal = predLabelVal';
    acc = mean(testDataClassLabel(:) == predLabelVal(:));
    accList(end+1) = acc;
    ccList(end+1) = cc;
    fprintf('\non validation set\n\taccuracy=%.4f (lambda1=%.4f, lambdaLoc=%.4f, T=%d)\n\n', acc, lambda1, lambdaLoc, T);
end

[a,b] =max(accList);
fprintf('\nbest acc: %.4f (c=%d)\n', a,ccList(b));

