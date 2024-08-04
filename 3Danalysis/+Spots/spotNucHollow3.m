function [bwSpot, nucLabelSpot] = spotNucHollow3(im, nucLabel, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called if "deviationofGradient" is used
%   Use tophat with extreme caution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imUseType = metaDataDS.analysisInfo.imUseType;
howManySigmas = 3;
thresholdFraction = 0.95; %% Default value is 0.90
labels = sort(nonzeros(HelperFunctions.count_unique(nucLabel)));
nNuc = length(labels);
meanNuc = zeros(nNuc,1);
sigmaNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc, 1);
minMRNAspotVol = metaDataDS.analysisInfo.minMRNAspotVol;
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
for i = 1:nNuc
    meanNuc(i) = mean2(im(nucLabel==labels(i)));
    sigmaNuc(i) = std2(im(nucLabel==labels(i)));
    thresNuc(i) = prctile(im(nucLabel==labels(i) & im>(howManySigmas*sigmaNuc(i)+meanNuc(i))),...
        (thresholdFraction*100),'all');
    imSpotTemp = im>thresNuc(i);
    imSpotTemp(nucLabel~=labels(i)) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
end
bw1 = imbinarize(imSpot, 0.7);
bw2 = bwareaopen(bw1, minMRNAspotVol/3);
% bw2 = imclearborder(bw2);


% Assuming transcription spots are twice as large and appears over twice as many frames
% if imUseType==3 % use "sharp"
%     bw2 = bwareaopen(bw1, minMRNAspotVol/3);
% else
%     bw2 = bwareaopen(bw1, minMRNAspotVol);
% end

bwSpot = bw2;
nucLabelSpot = bw2.*imNucLabelSpot;
end
