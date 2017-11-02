% shuffle data for later use -- split into training&validation sets
%
%
clear 
close all
clc;

%% 
dirDataset = './database';
categoryList = dir( strcat(dirDataset,'/* fossil'));
randIdxMat = cell(1,length(categoryList));
for categID = 1:numel(categoryList)
    imList = dir( fullfile(dirDataset, categoryList(categID).name , '*.jpg') );
    
    randTMP = randperm(length(imList));
    randIdxMat{categID} = randTMP;
    
    for imID = 1:numel(imList)                
        fprintf('%s\n', imList(imID).name);
        srcName = fullfile(dirDataset, categoryList(categID).name, imList(imID).name);
        dstName = fullfile(dirDataset, categoryList(categID).name, [ num2str(randTMP(imID)) '_' imList(imID).name]);
        movefile( srcName, dstName );
    end
end

save('part0.mat', 'randIdxMat');