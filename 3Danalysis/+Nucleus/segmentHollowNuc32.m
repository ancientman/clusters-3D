function [fim] = segmentHollowNuc32(im, metaDataDS, filledNuc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a crude segmentation technique that produces smoothass boundaries to
% segmented nuclei. 
% This corresponds to option "watershed" of nucleusDetectMethod
% This function segments the nuclei specifically by watershed. 
% This can be used for most applications.
% It first blurs the input image then performs watershed operation. 
% Then uses a convolution of the blurred watershed labels to create the
% final segmented image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucThres = metaDataDS.analysisInfo.nucThres;
elementSize = metaDataDS.analysisInfo.elementSize;
minNucVol = metaDataDS.analysisInfo.minNucVol;
clearBorder = metaDataDS.analysisInfo.clearBorder;
imScale = rescale(im);
imFilt1 = medfilt3(imScale);
imFilt2 = imgaussfilt3(imFilt1, 7);
imBin1 = imbinarize(imFilt2);
bw1 = ~imBin1;
bw1 = imfill(bw1, 'holes');
bw1 = bwareaopen(bw1, ceil(minNucVol/3));

t1 = adaptthresh(imFilt2, nucThres, 'ForegroundPolarity', 'bright', 'Statistic', 'gaussian'); 
imFilt3 = imgaussfilt3(imadjustn(imScale));
bw2 = imcomplement(imbinarize(imFilt2, t1));    

windowSize = 15; 
expFactor = 1.5;
kernel = ones(windowSize, windowSize, windowSize) / (windowSize^expFactor);
bw2Blur = convn(bw2, kernel, 'same');

bw2 = rescale(bw2Blur) > 0.2; %%%%%% change this is necessary (original 0.1)
bw2 = imfill(bw2, 'holes');
bw2 = bwareaopen(bw2, ceil(minNucVol/3));

D = -bwdist(~bw2);
D(bw2)= -inf;
DL = watershed(D);
bgm = DL == 0;
bw4 = zeros(size(bw2));
for i=1:size(bw2, 3)
%     bw4(:,:,i) = any((bw1(:,:,i)==1 & bgm(:,:,i)==0),3);
    bw4(:,:,i) = any((bw2(:,:,i)==1 & bgm(:,:,i)==0),3);
end
bw5 = imfill(bw4, 'holes');
bw5 = bwareaopen(bw5, ceil(minNucVol/3));
if clearBorder
    bw5 = imclearborder(bw5);
end
bw5 = imfill(bw5, 'holes');
se = strel('cuboid', elementSize*[1, 1, 2]);
Io = imopen(bw5, se);
If = imfill(Io, 'holes');
nucBW = bwareaopen(If, ceil(minNucVol/3));    
for i = 1:size(nucBW, 3)
    nucBW(:,:,i) = imfill(nucBW(:,:,i), 'holes');
end
fim = nucBW;
end
