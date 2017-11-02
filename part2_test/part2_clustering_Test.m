clear
close all
clc

%% parameters
K = 2; % k-mediod clustering
numTrain = 65; % the number of training data in each class
reSZ = [40, 40];
% dirDataset = './mixedSamples_2Dshape_binaryMask'; % using binary mask
dirDataset = './mixedSamples_2Dshape_thumbnail'; % using thumbnail images

%% fetch data of modern pollen grains
dirDataset_orgImg = './dataset';
categoryList = dir( strcat(dirDataset,'/mixed*'));

dataMatTest = [];
imList = {};
imListTest = {};
labelTest = [];

fprintf('fetch data...\n');
for categID = 1:numel(categoryList)
    fprintf('\t%s...\n', categoryList(categID).name);
	imList{end+1} = dir( fullfile(dirDataset, categoryList(categID).name, '*.bmp') );
    
    for imID = 1:numel(imList{end})       
        im = imread( fullfile(dirDataset, categoryList(categID).name, imList{end}(imID).name) );        
        im = double(im)/255;
        
        labelTest(end+1) = categID;
        imListTest{end+1} = fullfile(dirDataset, categoryList(categID).name, imList{end}(imID).name);
        dataMatTest = [dataMatTest, im(:)];
    end
end
fprintf('\tdone!\n');
dataMatTest = reshape(dataMatTest, [reSZ, size(dataMatTest,2)] );
%imageNameList = [imList{1};imList{2};imList{2}];
fprintf('done\n');

%% calculate the affinity matrix between every pair of images
shapeType = strfind(dirDataset,'_');
shapeType = dirDataset(shapeType(end)+1:end);
affinityMatName = ['affinityMatrix_' shapeType '.mat'];

fprintf('calculate/load the affinity matrix between every pair of images...\n');
if ~exist( affinityMatName, 'file')    
    affMat = zeros(size(dataMatTest,3));
    bestThyI2J = zeros(size(dataMatTest,3));
    for i = 1:size(dataMatTest,3)
        fprintf('\timage-%d is being matched to all others...\n', i);
        for j = i+1:size(dataMatTest,3)
            [mindist, thy, im3new]= minDistRot2D(dataMatTest(:,:,i), dataMatTest(:,:,j));
            bestThyI2J(j,i) = thy; % rotate image-j by thy to best match image-i
            affMat(i,j) = mindist;
        end
    end
    fprintf('\ndone!\n');    
    affMat = affMat + affMat';    
    save( affinityMatName, 'affMat', 'dataMatTrain', 'dataMatVal', 'imList', ...
        'imListTrain', 'imListVal', 'labelTrain', 'labelVal', 'bestThyI2J');
else
    load(affinityMatName);
end
%% clustering
shapeType = strfind(dirDataset,'_');
shapeType = dirDataset(shapeType(end)+1:end);
mediodsMatName = ['mediods_' shapeType '_K' num2str(K) '.mat'];

if exist( mediodsMatName, 'file')
    load( mediodsMatName );
    
    C = reshape(C3D, [size(C3D,1)*size(C3D,2), size(C3D,3)]);
    Cdisp = showdict(C, [reSZ(1), reSZ(2)], ceil(sqrt(K)), ceil(sqrt(K)));    
%     figure; imshow(Cdisp); title('2D shape centroids');
    
    figure;
    rNUM = ceil(sqrt(K));
    cNUM = ceil(K/rNUM);
    for i = 1:K
        subplot(rNUM, cNUM, i);
        imshow(ImgMedoid{i});
        title(['2D shape centroids' num2str(i)]);
    end
else
    
    %% k-medoid on whole set
    [medoidIdx, cidx] = kmedioids(affMat, K);
    
    C3D = dataMatTest(:,:,cidx);
    C = reshape(C3D, [size(C3D,1)*size(C3D,2), size(C3D,3)]);
    Cdisp = showdict(C, [size(im,1), size(im,2)], ceil(sqrt(K)), ceil(sqrt(K)));
    
%     figure; imshow(Cdisp); title('2D shape centroids');
    countK = zeros(K, 1);
    fprintf('\nstatistics for training set\n');
    for i = 1:K
        countK(i) = length(find(medoidIdx==i));
    end
    disp(countK);
    
    %% first pass to get the medoids of original images
    fprintf('\nfirst pass to get the medoids of original images...');
    
    ImgMedoid = cell(K,1);
    for i = 1:K
        idx = cidx(i);
        [junk, NAME, EXT] = fileparts(imListTest{idx});
        im = imread(  fullfile(dirDataset_orgImg , categoryList( labelTest(idx) ).name, strcat(NAME,'.jpg')) );
        ImgMedoid{i} = im;
    end
    
    % visualize the medoids of pollen images
    figure;
    rNUM = ceil(sqrt(K));
    cNUM = ceil(K/rNUM);
    for i = 1:K
        subplot(rNUM, cNUM, i);
        imshow(ImgMedoid{i});
        title(['2D shape centroids' num2str(i)]);
    end
    save(mediodsMatName, 'ImgMedoid', 'cidx', 'C3D', 'C');
end
%% get clustering result of fossil pollen grains and save
fprintf('\nsaving for visualization...\n');
dirDataSave = ['./mixedSamples_' shapeType '_visualize'];
dirDataSaveCanonicalShape = ['./mixedSamples_' shapeType '_canonicalShape'];

%{
% training data 
trainMedoidIdx = zeros(1,numel(imListTest));
for i = 1:numel(imListTest)
    fprintf('image-%d/%d...\n', i,numel(imListTest));
    [junk, NAME, ~] = fileparts(imListTest{i});
    im = imread(  fullfile(dirDataset_orgImg, categoryList( labelTest(i) ).name, strcat(NAME,'.jpg')) );
    
    disList = -1*ones(1,K);
    bestThyI2J = zeros(1,K);
    for k = 1:K
        [mindist, thy, im3new]= minDistRot2D(C3D(:,:,k), dataMatTest(:,:,i));
        bestThyI2J(k) = thy; % rotate image-j by thy to best match image-i
        disList(k) = mindist;
    end
    [valjunk, idjunk] = min(disList);
    trainMedoidIdx(i) = idjunk;
    im = makeImSquareByPadding(double(im)/255);
    im = imrotate(im, thy*180/pi);
    
    if ~isdir( fullfile(dirDataSave, categoryList(labelTest(i)).name) )
        mkdir(fullfile(dirDataSave, categoryList(labelTest(i)).name));
    end
    if ~isdir( fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(trainMedoidIdx(i)))) )
        mkdir(fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(trainMedoidIdx(i)))));
    end
    
    [~, curImgName, curExt] = fileparts(imListTest{i});
    imwrite( im, fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(trainMedoidIdx(i))), [curImgName, curExt]) );
    
    if ~isdir( fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name) )
        mkdir(fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name) );
    end
    
    imwrite( im, fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name, ...
        strcat(curImgName, '_K', num2str(trainMedoidIdx(i)), '.jpg')   ) );
end
%}

% test set
testMedoidIdx = zeros(1,numel(imListTest));
for i = 1:numel(imListTest)
    fprintf('image-%d/%d...\n', i,numel(imListTest));
    [junk, NAME, ~] = fileparts(imListTest{i});
    im = imread(  fullfile(dirDataset_orgImg, categoryList( labelTest(i) ).name, strcat(NAME,'.jpg')) );
    
    disList = -1*ones(1,K);
    bestThyI2J = zeros(1,K);
    for k = 1:K
        [mindist, thy, im3new]= minDistRot2D(C3D(:,:,k), dataMatTest(:,:,i));
        bestThyI2J(k) = thy; % rotate image-j by thy to best match image-i
        disList(k) = mindist;
    end
    [valjunk, idjunk] = min(disList);
    testMedoidIdx(i) = idjunk;
    im = makeImSquareByPadding(double(im)/255);
    im = imrotate(im, thy*180/pi);
    
    if ~isdir( fullfile(dirDataSave, categoryList(labelTest(i)).name) )
        mkdir(fullfile(dirDataSave, categoryList(labelTest(i)).name));
    end
    if ~isdir( fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(testMedoidIdx(i)))) )
        mkdir(fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(testMedoidIdx(i)))));
    end
    
    [~, curImgName, curExt] = fileparts(imListTest{i});
    imwrite( im, fullfile(dirDataSave, categoryList(labelTest(i)).name, strcat('cluster_', num2str(testMedoidIdx(i))), [curImgName, curExt]) );
    
    
    if ~isdir( fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name) )
        mkdir(fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name) );
    end
    
%     [PATHSTR,NAME,EXT] = fileparts(imageNameList(i).name);
    imwrite( im, fullfile(dirDataSaveCanonicalShape, categoryList(labelTest(i)).name, ...
        strcat(curImgName, '_K', num2str(testMedoidIdx(i)), '.jpg')   ) );
end



