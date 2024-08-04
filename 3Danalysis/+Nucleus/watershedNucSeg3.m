function [fim, t1] = watershedNucSeg3(im, mask, metaDataDS, filledNuc)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
nZslices = metaDataDS.analysisInfo.nZslices;
zPixUM = metaDataDS.analysisInfo.zPixUM;
nucThres = metaDataDS.analysisInfo.nucThres;
clearBorder = metaDataDS.analysisInfo.clearBorder;
minNucVol = metaDataDS.analysisInfo.minNucVol;
% mask2 = activecontour(imcomplement(im), mask, 10, 'edge', 'SmoothFactor', 1, 'ContractionBias', 1);
imScale = rescale(im);
imFilt = imgaussfilt3(imadjustn(imScale));
t1 = adaptthresh(imFilt, nucThres, 'ForegroundPolarity', 'bright', 'Statistic', 'gaussian'); 
if filledNuc == 0
    bw1 = imcomplement(imbinarize(imFilt, t1));    
elseif filledNuc == 1
%     bw1 = imbinarize(imFilt, t1);
    bw1 = mask;
end

bw1 = imfill(bw1, 'holes');
bw1 = bwareaopen(bw1, ceil(minNucVol/3));
D = -bwdist(~bw1);
D(bw1)= -inf;
DL = watershed(D);
bgm = DL == 0;

bw4 = zeros(size(bw1));
for i=1:size(bw1, 3)
    bw4(:,:,i) = any((bw1(:,:,i)==1 & bgm(:,:,i)==0),3);
end
bw5 = imfill(bw4, 'holes');
bw5 = bwareaopen(bw5, ceil(minNucVol/3));
if clearBorder
    bw5 = imclearborder(bw5);
end
fim = bw5;
end