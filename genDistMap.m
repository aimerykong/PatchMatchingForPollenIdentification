function W = genDistMap(dictLoc, imSize, feaSize, locationScaleFactor)
%
%
%
%
%

% imSize = curImgMat.imSize;
% feaSize = curImgMat.feaSize(1:2);
% dictLocBackup = dictLoc;

dictLoc = dictLoc/locationScaleFactor;

W = zeros(feaSize(1), feaSize(2), size(dictLoc,2));
yy = 1:feaSize(1);
xx = 1:feaSize(2);
yy = repmat( yy(:), [1, feaSize(2)] );
xx = repmat( xx, [feaSize(1), 1] );
patchLoc = yy;
patchLoc = repmat(patchLoc, [1,1,2]);
patchLoc(:,:,2) = xx;


patchLoc = patchLoc-1;
patchLoc(:,:,1) = patchLoc(:,:,1) ./ feaSize(1) .* imSize(1);
patchLoc(:,:,2) = patchLoc(:,:,2) ./ feaSize(2) .* imSize(2);
patchLoc = patchLoc + 1;

patchLoc(:,:,1) = patchLoc(:,:,1) - imSize(1)/2;
patchLoc(:,:,2) = patchLoc(:,:,2) - imSize(2)/2;

for i = 1:size(dictLoc,2)
    Wtmp = patchLoc;
    Wtmp(:,:,1) = Wtmp(:,:,1) - dictLoc(1,i);
    Wtmp(:,:,2) = Wtmp(:,:,2) - dictLoc(2,i);
    Wtmp = sum(abs(Wtmp),3);
    Wtmp = Wtmp.^0.5;
%     Wtmp = sqrt( sum((Wtmp.^2),3) );
    W(:,:,i) = Wtmp;
end
A = max(W,[],3);
W = bsxfun(@rdivide, W, A);
