function [fim, L] = watershedSpot(im, nucMask, mask, metaDataDS)
% mask = bwareaopen(mask,minSpotSize);
CC = bwconncomp(nucMask);
labelMat = labelmatrix(CC);
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
xyPsfUM = metaDataDS.analysisInfo.xyPsfUM;
spotFilter = metaDataDS.analysisInfo.spotFilter;
el = metaDataDS.analysisInfo.elementSize;
minSpotSizeUM = metaDataDS.analysisInfo.minSpotSize/1000;
if minSpotSizeUM==0 || minSpotSizeUM<xyPsfUM% Convert to pixels here
    minSpotSizeUM = ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
else
    minSpotSizeUM = ceil((pi*(minSpotSizeUM/2)^2)/(xPixUM*yPixUM));
end
im = rescale(im);
im(mask==0) = 0;
switch spotFilter
    case 1% gaussian (recomended for global spots)
        imBlur = imgaussfilt(im);        
    case 2% bilateral (recomended for global spots)
        imBlur = imbilatfilt(im);
    case 3% tophat (NOT recomended for global spots)
        se = strel('disk', el, 8);
        imBlur = imtophat(im,se); 
        imBlur = imadjust(imBlur);        
    case 4% median (NOT recomended for global spots)
        imBlur = medfilt2(im);        
    otherwise% none
        imBlur = im;
end
bw = imBlur>0;
bw = bwareaopen(bw,minSpotSizeUM);
D = bwdist(~bw);
D = -D;
D(~bw) = -Inf;
L = watershed(D);
bw2 = mask;
bw2(L==0) = 0;
bw3 = bw2;
bw3 = bwareaopen(bw2, minSpotSizeUM);
bw3 = imclearborder(bw3);
spotLabel = labelMat.*uint8(bw3);
fim = bw3;
L = spotLabel;
end