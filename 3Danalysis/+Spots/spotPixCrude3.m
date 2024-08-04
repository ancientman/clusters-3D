function [bwSpot, nucLabelSpot] = spotNucFilled(im, nucMask, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called if "deviationofGradient" is used
%   Use tophat with extreme caution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imScale = rescale(im);
CC = bwconncomp(nucMask);
nNuc = CC.NumObjects;
labelMat = labelMatrix(CC);
spotFilter = metaDataDS.analysisInfo.spotFilter;
el = metaDataDS.analysisInfo.elementSize;
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
meanNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc, 1);
switch spotFilter
    case 1% gaussian
        thresholdFraction = 0.9;
        imBlur = imgaussfilt(im);
    case 2% bilateral
        thresholdFraction = 0.9;
        imBlur = imbilatfilt(im);
    case 3% tophat
        thresholdFraction = 0.9;
        se = strel('disk', el, 8);
        imBlur = imtophat(imScale,se);
%         imBlur = imadjust(imBlur);
    case 4% median
        warning('No spot detect type defined, using median');
        thresholdFraction = 0.9;
        imBlur = medfilt2(im);
    otherwise
        warning('No spot detect type defined, using median');
        thresholdFraction = 0.9;
        imBlur = im;
end
imAdd = imgradient(imScale, 'sobel') + medfilt2((imBlur));
for i = 1:nNuc
    meanNuc(i) = mean2(imAdd(labelMat==i));
    thresNuc(i) = prctile(imAdd(labelMat==i & imAdd>meanNuc(i)),...
        (thresholdFraction*100),'all');    
    imSpotTemp = imAdd>thresNuc(i);    
    imSpotTemp(labelMat ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end
bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end


