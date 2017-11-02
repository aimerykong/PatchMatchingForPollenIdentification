clear 
close all
clc;

addpath('../toolbox/exportFig');

%% get the hist as a side oracle over a set of training data
im = {};
histRecord = 0;

% critch examples
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/15_CropLBS-4-2__K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/22_CropPC-21_Slide20Pos5(GOOD)_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/24_ZCropLBS-1_Slide16Pos3(OK)_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/44_CropLBS-4-2_Slide12Pos12(GOOD)_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/83_ZCropLBS-1_Slide16Pos5(GOOD)_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/84_CropPC-21_Slide20Pos1(OK)_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/92_CropPC-21_Slide20Pos4(GOOD)_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/99_ZCropLBS-1_Slide11Pos1(OK)_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/105_CropPC-19_Slide23Pos23(OK)_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/critchfieldii fossil/107_CropPC-21_Slide5Pos1(OK)_K2.jpg');

% glauca examples
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/205_Nelson Lake NE954 Pos 56_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/276_Nelson Lake NE906 Pos 90_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/278_Nelson Lake NE1010 Pos 149_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/295_Nelson Lake NE906 Pos 2_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/306_Nelson Lake NE778 Pos 319_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/308_Nelson Lake NE954 Pos 244_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/363_Nelson Lake NE978 Pos 172 cropped_K2.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/394_Nelson Lake NE978 Pos 173 cropped_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/396_Nelson Lake NE978 Pos 74 cropped_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/glauca fossil/405_Nelson Lake NE954 Pos 227_K1.jpg');

% mariana examples
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/228_Nelson Lake NE778 Pos 82_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/230_Nelson Lake NE 1010 Pos 36 cropped_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/268_Nelson Lake NE1010 Pos 69_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/273_Nelson Lake NE1010 Pos 76_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/279_Nelson Lake NE954 Pos 212_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/292_Nelson Lake NE830 Pos 67_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/298_Nelson Lake NE 1010 Pos 51 cropped_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/307_Nelson Lake NE830 Pos 75_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/332_Nelson Lake NE1010 Pos 97_K1.jpg');
im{end+1}=imread('/home/skong/PollenProject/subproject1_fossilOnlyForThreeClassification/database_binaryMask_canonicalShape/mariana fossil/348_Nelson Lake NE906 Pos 147_K1.jpg');

for i = 1:length(im)
    A = double(im{i});
    idx = find(A>0);
    b = hist( A(idx), 100);
    histRecord = histRecord + b;
end
histRecord = histRecord ./ length(im);

%%
% copy conanical shape folder to this one, and preprocess images w.r.t
% global contrast
datasetDir = 'database_binaryMask_canonicalShape_globalContrastNorm';
classPath = dir(datasetDir);
for i = 3:length(classPath)
    imList = dir( fullfile(datasetDir, classPath(i).name, '*.jpg') );
    fprintf('%s...\n', classPath(i).name);
    for imgID = 1:length(imList)
        filename = imList(imgID).name;
        im = imread( fullfile(datasetDir, classPath(i).name, filename) );
        im2 = im;
        idx = find(im2>0);
        im2(idx) = histeq(im2(idx), histRecord);
        
        imwrite(im2,  fullfile(datasetDir, classPath(i).name, filename) );
        
%         subplot(1,2,1);
%         imshow(im); title('original image');
%         subplot(1,2,2);
%         imshow(im2); title('histeq');
%         fprintf('.');
    end
end




