function [fim] = segmentHollowNuc3(im, metaDataDS, filledNuc)
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
imScale = rescale(im);
imScale = imadjustn(imScale);
imFilt1 = medfilt3(imScale);
imFilt2 = imgaussfilt3(imFilt1);
imBin0 = imextendedmax(imFilt2, 0.7);
imBin1 = imbinarize(imFilt2, 0.4);
imBin2 = imbinarize(imFilt2, 0.5);
imBin3 = imbinarize(imFilt2, 0.6);
imBin = imBin0+imBin1+imBin2+imBin3;
bw1 = ~imBin;
bw1 = imfill(bw1, 'holes');
windowSize = 5; 
expFactor = 1.5;
kernel = ones(windowSize, windowSize, windowSize) / (windowSize^expFactor);
imBlur2 = convn(bw1, kernel, 'same');
bw2 = imBlur2 > nucThres;
[bw3, ~] = Nucleus.watershedNucSeg3(imFilt2, bw2, metaDataDS, filledNuc);
bw3 = imfill(bw3, 'holes');
se = strel('cuboid', elementSize*[1, 1, 2]);
Io = imopen(bw3, se);
If = imfill(Io, 'holes');
nucBW = bwareaopen(If, ceil(minNucVol/3));    
fim = nucBW;
end
