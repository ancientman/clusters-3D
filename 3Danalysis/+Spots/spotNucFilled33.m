function [bwSpot, nucLabelSpot] = spotNucFilled33(im, nucLabel, metaDataDS,  removeBorderSpots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Deviation from mean
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nNuc = max(nucLabel,[],'all');
elementSize = metaDataDS.analysisInfo.elementSize;
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
thresNuc = zeros(nNuc, 1);

% spotFilter = metaDataDS.analysisInfo.spotFilter;
meanNuc = zeros(nNuc,1);

thresholdFraction = 0.8;

%% Technique based on deviation from the mean
for i = 1:nNuc
    meanNuc(i) = mean2(im(nucLabel==i));
    thresNuc(i) = 2.*meanNuc(i);
    imSpotTemp = im>thresNuc(i);    
    imSpotTemp(nucLabel ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end

%% Watershed the spots
% imSpotW = watershedSpot(im, imSpot, metaDataDS);
imSpotW = imclearborder(imSpotW);
imNucLabelSpotW = imSpotW.*nucLabel;

if  removeBorderSpots == 0
    bwSpot = imSpot;
    nucLabelSpot = imNucLabelSpot;
else
    bwSpot = imSpotW;
    nucLabelSpot = imNucLabelSpotW;
end
end

function [fim] = watershedSpot(im, mask, metaDataDS)
el = metaDataDS.analysisInfo.elementSize;
% im = rescale(im);
im(mask==0) = 0;

imBlur = imgaussfilt(im);        
se = strel('disk', el, 8);
imBlur = imtophat(im,se); 

bw = imBlur>0;
D = bwdist(~bw);
D = -D;
D(~bw) = -Inf;
L = watershed(D);
bw2 = mask;
bw2(L==0) = 0;
bw3 = bw2;
bw3 = imclearborder(bw3);
fim = bw3;
end