function [fim] = segmentHollowNucThin3(im, metaDataDS, filledNuc)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otsuFactor = 0.2;
imBlur1 = imgaussfilt3(medfilt3(im), 3);
% imBin1 = edge3(imBlur1,'approxcanny',otsuFactor, 3);
% imBin1 = imfill(imBin1, 'holes'); % Spots
imMorph1 = imtophat(imBlur1, strel('sphere', 11));
imBin1 = imbinarize(rescale(imMorph1), otsuFactor); % Spots
% im2 = im;
% im2(imBin1) = 0;
imBlur2 = imBlur1;
imBlur2(imBin1) = 0;
imBin2 = imbinarize(1-rescale(imBlur2));
imBin2 = bwareaopen(imBin2, minNucVol);
imBin2 = imfill(imBin2, 'holes');
[bw3, ~] = Nucleus.watershedNucSeg3(imBlur2, imBin2, metaDataDS, filledNuc); % Nuc
windowWidth = 21;
kernel = ones(windowWidth) / windowWidth^2;
% boundValue = max(imBlur2, [], 'all', 'omitnan');
boundValue = ceil(mean(maxk(max(max(imBlur2, [], 3, 'omitnan'), [], 2),50)));
imBlur2(1:windowWidth,:,:) = boundValue; 
imBlur2(end-windowWidth:end,:,:) = boundValue; 
imBlur2(:,1:windowWidth,:) = boundValue; 
imBlur2(:,end-windowWidth:end,:) = boundValue; 

if size(imBlur2, 3)>2
    windowWidth = 1;
    imBlur2(:,:,1:windowWidth) = boundValue; 
    imBlur2(:,:,end-windowWidth:end) = boundValue;
end
imBlur3 = conv2(double(imBlur2), kernel, 'same');
imBlur3 = 1-rescale(imBlur3);
imBin2 = imbinarize((imBlur3), 0.75);
imBin2 = imfill(imBin2, 'holes'); 
imBin2 = bwareaopen(imBin2, minNucVol); % Nuc
imBin2(1,:,:) = 0; imBin2(end,:,:) = 0; imBin2(:,1,:) = 0; imBin2(:,end,:) = 0; 
if size(imBin2, 3)>2
    imBin2(:,:,1) = 0; imBin2(:,:,end) = 0;
end
% imBin2 = imfill(imBin2, 'holes');

D = -bwdist(~imBin2);
D(imBin2)= -inf;
DL = watershed(D);
bgm = DL == 0;

windowSize = 21; 
expFactor = 1.5;
kernel = ones(windowSize, windowSize, windowSize) / (windowSize^expFactor);
imBlur3 = convn(imBin2, kernel, 'same');
imBin3 = imBlur3 > nucThres;
imBin3 = imfill(imBin3, 'holes'); 
imBin3 = bwareaopen(imBin3, minNucVol); % Nuc
imBin3(bgm) = 0;
imBin3 = bwareaopen(imBin3, minNucVol); % Nuc
fim = imBin3;
end
