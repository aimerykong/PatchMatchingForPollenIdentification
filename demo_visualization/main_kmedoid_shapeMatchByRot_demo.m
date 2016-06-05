clear
close all
clc

%% fetch data
dirDataset = './database_2Dshape';
dirDataset_orgImg = './database';
categoryList = dir(dirDataset);
categoryList = categoryList(3:end);
categoryList = categoryList([1,3]);

im2NameList = {'P glauca B1500 Pos 5.tif_Files', 'P glauca B1503 Pos 104', 'P glauca B1500 Pos 202.tif_Files', 'P glauca B1503 Pos 120'};
im2ID = 2;

im1 = imread( fullfile(dirDataset, 'glauca - modern', 'P glauca B1500 Pos 193.tif_Files.bmp') );
im2 = imread( fullfile(dirDataset, 'glauca - modern', strcat(im2NameList{im2ID}, '.bmp') ) );

im1Org = imread(fullfile(dirDataset_orgImg, 'glauca - modern', 'P glauca B1500 Pos 193.tif_Files.jpg'));
im2Org = imread(fullfile(dirDataset_orgImg, 'glauca - modern', strcat(im2NameList{im2ID}, '.jpg')));

im1 = mat2gray(im1);
im2 = mat2gray(im2);

%% rotate image 2 on fixed angles to match the target image 1
diagLength = ceil(sqrt(sum(size(im1).^2)))+4;

im1new = zeros(diagLength);
cIm1New = round(size(im1new)/2);
im1new( 1+cIm1New(1)-floor(size(im1,1)/2):cIm1New(1)+ceil(size(im1,1)/2), 1+cIm1New(2)-floor(size(im1,2)/2):cIm1New(2)+ceil(size(im1,2)/2)) = im1;

rotList = [0:10:359] *pi/180 ; % pre-defined rotation bins

%% function for calculating the distance
%[mindist, thy, im3new]= minDistRot2D(im1, im2, rotList);

%%{
[idx1, idx2] = ind2sub( size(im2), 1:numel(im2));
X = [idx1; idx2];
distList = zeros(2*length(rotList),1);
count = 1;
for rotId = 1:length(rotList)
    thy = rotList(rotId);
    %{
    RotMat = [  cos(thy)   -sin(thy) ; sin(thy)   cos(thy) ];
    
    Y = RotMat*X;
    
    Z = round(bsxfun(@minus, Y, min(Y,[],2))+1);
    im3 = zeros( max(Z(1,:)), max(Z(2,:)) );
    z = sub2ind(size(im3), Z(1,:), Z(2,:));
    im3(z) = im2(:);               
    %im3 = fillObjectHoles(im3);
    %}
    im3 = imrotate(im2, thy*180/pi, 'nearest');
    
    im3new = zeros(diagLength);
    cIm3New = round(size(im3new)/2);
    im3new( 1+cIm3New(1)-floor(size(im3,1)/2):cIm3New(1)+ceil(size(im3,1)/2), 1+cIm3New(2)-floor(size(im3,2)/2):cIm3New(2)+ceil(size(im3,2)/2)) = im3;
        
    distList(count) = norm(im1new-im3new, 'fro');
    count = count + 1;
end

%% flip left-right and match again
im2flip = fliplr(im2);
for rotId = 1:length(rotList)
    thy = rotList(rotId);
    %{
    RotMat = [  cos(thy)   -sin(thy) ; sin(thy)   cos(thy) ];
    
    Y = RotMat*X;
    
    Z = round(bsxfun(@minus, Y, min(Y,[],2))+1);
    im3 = zeros( max(Z(1,:)), max(Z(2,:)) );
    z = sub2ind(size(im3), Z(1,:), Z(2,:));
    im3(z) = im2flip(:);
       
    im3 = fillObjectHoles(im3);
    %}
    im3 = imrotate(im2flip, thy*180/pi, 'nearest');
    
    im3new = zeros(diagLength);
    cIm3New = round(size(im3new)/2);
%     im3new( 1+cIm3New(1)-floor(size(im4,1)/2):cIm3New(1)+ceil(size(im4,1)/2), 1+cIm3New(2)-floor(size(im4,2)/2):cIm3New(2)+ceil(size(im4,2)/2)) = im4;
    im3new( 1+cIm3New(1)-floor(size(im3,1)/2):cIm3New(1)+ceil(size(im3,1)/2), 1+cIm3New(2)-floor(size(im3,2)/2):cIm3New(2)+ceil(size(im3,2)/2)) = im3;
    distList(count) = norm(im1new-im3new, 'fro');
    count = count + 1;
end

%% choose the smallest distance among all possible rotations
[minVal, minIdx] = min(distList);
flipFlag = 0;
if minIdx > length(distList)/2
    minIdx = minIdx - length(distList)/2;
    flipFlag = 1;
end

thy = rotList(minIdx);
%{
RotMat = [  cos(thy)   -sin(thy) ; sin(thy)   cos(thy) ];
Y = RotMat*X;

Z = round(bsxfun(@minus, Y, min(Y,[],2))+1);
im3 = zeros( max(Z(1,:)), max(Z(2,:)) );
z = sub2ind(size(im3), Z(1,:), Z(2,:));
%}
if flipFlag
    
    im3 = imrotate(im2flip, thy*180/pi, 'nearest');
%     im3(z) = im2flip(:); im3 = fillObjectHoles(im3);
else
    im3 = imrotate(im2, thy*180/pi, 'nearest');
%     im3(z) = im2(:); im3 = fillObjectHoles(im3);
end

im3new = zeros(diagLength);
cIm3New = round(size(im3new)/2);
im3new( 1+cIm3New(1)-floor(size(im3,1)/2):cIm3New(1)+ceil(size(im3,1)/2), 1+cIm3New(2)-floor(size(im3,2)/2):cIm3New(2)+ceil(size(im3,2)/2)) = im3;
%}

%% visualize the result -- the shape matching
figure(1);
subplot(2,3,1);
imshow(im1); title('im1 - the target to match');
subplot(2,3,2);
imshow(im2); title('im2 - the novel image');
subplot(2,3,3);
imshow(im2flip); title('fliped im2 (left-right)');

subplot(2,3,4);
imshow(im1new); title('im1 - the target with padding');
subplot(2,3,5);
imshow(im3); title('im3 - best rotation on im2');
subplot(2,3,6);
imshow(im3new); title('im3 - with padding');

fprintf('dist(im1,im2) without rotation: %.4f\n', norm(im1-im2,'fro'));
[minVal] = min(distList(1:length(distList)/2));
fprintf('dist(im1,rot(im2)) with rotation: %.4f\n', minVal);
[minVal] = min(distList(length(distList)/2+1:end));
fprintf('dist(im1,rot(flip(im2))) with flip and rotation: %.4f\n', minVal);

%% visualize the result -- the comparison between rotated testing image and the target image
figure(2);

im1OrgSquare = zeros(max([ceil(sqrt(sum(size(im1Org).^2))) ceil(sqrt(sum(size(im2Org).^2)))]));
im2OrgSquare = im1OrgSquare;

cIm1Org = round(size(im1OrgSquare)/2);
cIm2Org = round(size(im2OrgSquare)/2);
im2OrgRot = 0*im2OrgSquare;

im1OrgSquare( 1+cIm1Org(1)-floor(size(im1Org,1)/2):cIm1Org(1)+ceil(size(im1Org,1)/2), 1+cIm1Org(2)-floor(size(im1Org,2)/2):cIm1Org(2)+ceil(size(im1Org,2)/2)) = im1Org;
im2OrgSquare( 1+cIm2Org(1)-floor(size(im2Org,1)/2):cIm2Org(1)+ceil(size(im2Org,1)/2), 1+cIm2Org(2)-floor(size(im2Org,2)/2):cIm2Org(2)+ceil(size(im2Org,2)/2)) = im2Org;

if flipFlag
    im2Org = fliplr(im2Org);
end

cim2Org = round(size(im2Org)/2);
[idx1, idx2] = ind2sub( size(im2Org), 1:numel(im2Org));
X = [idx1; idx2];

thy = rotList(minIdx);
A = imrotate(im2Org, thy*180/pi, 'bicubic');
im2OrgRot( 1+cIm2Org(1)-floor(size(A,1)/2):cIm2Org(1)+ceil(size(A,1)/2), 1+cIm2Org(2)-floor(size(A,2)/2):cIm2Org(2)+ceil(size(A,2)/2)) = A;

RotMat = [  cos(thy)   -sin(thy) ; sin(thy)   cos(thy) ];
%{
Y = RotMat*bsxfun(@minus, X, cim2Org(:));
Y = bsxfun(@plus, Y, round([size(im2OrgSquare,1);size(im2OrgSquare,2)]/2));

Z = round(Y+1);
z = sub2ind(size(im2OrgRot), Z(1,:), Z(2,:));
im2OrgRot(z) = im2Org(:);
%}

subplot(1,3,1);
imshow(uint8(im2OrgSquare)); title('image 2');
subplot(1,3,2);
imshow(uint8(im1OrgSquare)); title('image 1');
subplot(1,3,3);
imshow(uint8(im2OrgRot)); title('rotated image 2');

%% visualize the result -- randomly draw lines to see how well they are matching
if flipFlag
    imVis = cat(2, fliplr(im2OrgSquare/255), im1OrgSquare/255, im2OrgRot/255);
else
    imVis = cat(2, im2OrgSquare/255, im1OrgSquare/255, im2OrgRot/255);
end

imVis(:, size(im2OrgSquare,2):size(im2OrgSquare,2)+1) = 1;
imVis(:, 2*size(im2OrgSquare,2):2*size(im2OrgSquare,2)+1) = 1;

figure(3);
imshow(imVis); title('randomly matching patches');

patchNum = 10;
patchSize = 52;
thresh = 0.15; % generate meaningful patches randomly
[xList,yList] = genMeaningfulPatches(im1OrgSquare, patchNum, patchSize, thresh);
X = [yList';xList'];

Y = RotMat\bsxfun(@minus, X, cIm1Org(:)); % rotate on the image center
%Y = bsxfun(@plus, Y, round([size(im2OrgSquare,1);size(im2OrgSquare,2)]/2));
Y = round(bsxfun(@plus, Y, cIm1Org(:)));


x = xList+size(im1OrgSquare,1);
y = yList;
for i = 1:patchNum
    % draw boxes on image 1
    line([x(i)-patchSize/2, x(i)+patchSize/2-1, x(i)+patchSize/2-1, x(i)-patchSize/2, x(i)-patchSize/2],  [y(i)-patchSize/2, y(i)-patchSize/2, y(i)+patchSize/2-1, y(i)+patchSize/2-1, y(i)-patchSize/2], 'linewidth', 2, 'color', 'm' )
    
    % draw lines and boxes on rotated image 2, a mirror from image1 to
    % rotated image 2
    line([x(i)+size(im1OrgSquare,1)-patchSize/2, x(i)+patchSize/2-1+size(im1OrgSquare,1), x(i)+patchSize/2-1+size(im1OrgSquare,1), x(i)+size(im1OrgSquare,1)-patchSize/2, x(i)+size(im1OrgSquare,1)-patchSize/2],  [y(i)-patchSize/2, y(i)-patchSize/2, y(i)+patchSize/2-1, y(i)+patchSize/2-1, y(i)-patchSize/2], 'linewidth', 2, 'color', 'c' )    
    line( [size(im1OrgSquare,1)+xList(i), 2*size(im1OrgSquare,1)+xList(i)], [yList(i), yList(i)], 'linewidth', 2, 'color', 'r');        
    
    % draw lines and boxes on the original image 2
    line( [size(im1OrgSquare,1)+xList(i), Y(2,i)], [yList(i), Y(1,i)], 'linewidth', 2, 'color', 'b');    
    line([Y(2,i)-patchSize/2, Y(2,i)-1+patchSize/2, Y(2,i)-1+patchSize/2, Y(2,i)-patchSize/2, Y(2,i)-patchSize/2],  [Y(1,i)-patchSize/2, Y(1,i)-patchSize/2, Y(1,i)-1+patchSize/2, Y(1,i)-1+patchSize/2, Y(1,i)-patchSize/2], 'linewidth', 2, 'color', 'c' )    
end





