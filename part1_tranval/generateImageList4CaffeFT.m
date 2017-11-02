clear
clc

databaseDIR = './database';

speciesName = {'critchfieldii fossil', 'glauca fossil', 'mariana fossil'};

label = [];
imgPath = {};

trainNames = {};
testNames = {};
trainLabels = [];
testLabels = [];
numTrain = 65;

%% split into training/testing sets
trIdx = 1;
tsIdx = 1;
for sid = 1:length(speciesName)
    imList = dir( fullfile(databaseDIR, speciesName{sid}, '*.jpg') );
    itmp = 0;
    for i = 1:numTrain
        itmp = itmp + 1;
        trainNames{end+1} = fullfile(databaseDIR, speciesName{sid}, imList(itmp).name);
        trainLabels(end+1) = sid;
    end
    for i = 1+numTrain:length(imList)
        itmp = itmp + 1;
        testNames{end+1} = fullfile(databaseDIR, speciesName{sid}, imList(itmp).name);
        testLabels(end+1) = sid;
    end    
end

%% shuffle
randIdx = randperm(length(trainNames));
trainNames = trainNames(randIdx);
trainLabels = trainLabels(randIdx);

%% output to file
fn = 'imgListTrain4Caffe.txt';
fn = fopen(fn, 'w');
for i = 1:length(trainNames)
    fprintf(fn, '%s %d\n', trainNames{i}, trainLabels(i));
end
fclose(fn);

fn = 'imgListTest4Caffe.txt';
fn = fopen(fn, 'w');
for i = 1:length(testNames)
    fprintf(fn, '%s %d\n', testNames{i}, testLabels(i));
end
fclose(fn);





