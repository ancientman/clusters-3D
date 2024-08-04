function [fim] = segmentFilledNuc3(im, metaDataDS, filledNuc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elementSize = metaDataDS.analysisInfo.elementSize;
nucThres = metaDataDS.analysisInfo.nucThres;
clearBorder = metaDataDS.analysisInfo.clearBorder;
minNucVol = metaDataDS.analysisInfo.minNucVol;
imScale = rescale(im);
imScale = imadjustn(imScale);
imFilt1 = medfilt3(imScale);
imFilt2 = imgaussfilt3(imFilt1);
se = strel('cuboid',elementSize*[1 1 1]);
bw1 = imerode(imFilt2, se);
bw2 = imreconstruct(bw1, imScale);
bw3 = imdilate(bw2, se);
bw4 = imreconstruct(imcomplement(bw3),imcomplement(bw2));
bw4 = imcomplement(bw4);
% bw5 = imbinarize(imgaussfilt3(bw4, 7)); %use 2 for normal 7 for low
bw5 = imbinarize(imgaussfilt3(bw4, 7)); 
bw6 = imclose(bw5,se);
bw6 = imerode(bw6,se);

[bw7, ~] = Nucleus.watershedNucSeg3(imFilt2, bw6, metaDataDS, filledNuc);
% bw7 = bw6;

if clearBorder
    bw7 = imclearborder(bw7);
end

% windowSize = 5; 
% expFactor = 1.15;
% kernel = ones(windowSize, windowSize, windowSize) / (windowSize^expFactor);
% imBlur2 = convn(bw7, kernel, 'same');
% bw8 = imBlur2 > nucThres;
% bw8 = imfill(bw8, 'holes');
% bw8 = bwareaopen(bw8, ceil(minNucVol));
fim = bw7;
end

% elementSize = metaDataDS.analysisInfo.elementSize;
% nucThres = metaDataDS.analysisInfo.nucThres;
% clearBorder = metaDataDS.analysisInfo.clearBorder;
% imScale = rescale(im);
% imScale = imadjustn(imScale);
% imFilt1 = medfilt3(imScale);
% imFilt2 = imgaussfilt3(imFilt1);
% se = strel('cuboid',elementSize*[1 1 1]);
% bw1 = imerode(imFilt2, se);
% bw2 = imreconstruct(bw1, imScale);
% bw3 = imdilate(bw2, se);
% bw4 = imreconstruct(imcomplement(bw3),imcomplement(bw2));
% bw4 = imcomplement(bw4);
% bw5 = imbinarize(imgaussfilt3(bw4, 7)); %use 2 for normal 7 for low
% bw6 = imclose(bw5,se);
% bw6 = imerode(bw6,se);
% % windowSize = 11; 
% % expFactor = 1.15;
% % kernel = ones(windowSize, windowSize, windowSize) / (windowSize^expFactor);
% % imBlur2 = convn(bw6, kernel, 'same');
% % bw7 = imBlur2 > nucThres;
% [bw7, ~] = Nucleus.watershedNucSeg3(imFilt2, bw6, metaDataDS, filledNuc);
% if clearBorder
%     bw7 = imclearborder(bw7);
% end
% fim = bw7;
% end