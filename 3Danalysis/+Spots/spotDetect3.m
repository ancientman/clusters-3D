function [spotLabel, spotPropStruct] = spotDetect3(im, nucLabel, channelFilled, metaDataDS,  removeBorderSpots)
minSpotVol = metaDataDS.analysisInfo.minSpotVol;
imUseType = metaDataDS.analysisInfo.imUseType;
% colorChannels = metaDataDS.imagingInfo.channelCount;
if imUseType==1 % use "raw"
    imUse = im;
elseif imUseType==2 % use "smooth"
    imUse = Preprocess.smoothRaw3(im, metaDataDS);
elseif imUseType==3 % use "sharp"
    imUse = Preprocess.sharpRaw3(im, metaDataDS);
else % use "raw"
    imUse = im;
end 

if channelFilled==1
    [bwSpot, ~] = Spots.spotNucFilled3(imUse, nucLabel, metaDataDS); % based on deviations from base
%     [bwSpot, ~] = Spots.spotNucFilled32(imUse, nucLabel, metaDataDS,  removeBorderSpots); % based on deviations from base
%     [bwSpot, ~] = Spots.spotNucFilled33(imUse, nucLabel, metaDataDS,  removeBorderSpots); % based on deviations from base
else
%     [bwSpot, ~] = Spots.spotNucHollow3(imUse, nucLabel, metaDataDS); % based on deviations from base
    [bwSpot, ~] = Spots.spotNucHollow32(imUse, nucLabel, metaDataDS); % based on deviations from base
%  [bwSpot, ~] = Spots.spotNucHollow33(imUse(:,:,:,t), nucLabel(:,:,:,t), metaDataDS); % for bottleneck
    if removeBorderSpots ~=0
        bwSpot = imclearborder(bwSpot);
    end
end
 maskSpot = bwSpot;
% if colorChannels==1
%     bigSpotMask = bwareaopen(bwSpot, 10*minSpotVol);
%     [watershedSpotMask] = Spots.watershedSpot3(bigSpotMask, metaDataDS);
%     maskSpot = xor(xor(bwSpot, bigSpotMask), watershedSpotMask);
% else
%     maskSpot = bwSpot;
% end
maskSpot = imfill(maskSpot, 'holes');
if channelFilled== 0
    if imUseType==3
        maskSpot = bwareaopen(maskSpot,minSpotVol/3); % smaller mask for deconvolved
    else
        maskSpot = bwareaopen(maskSpot,minSpotVol);
    end
end
% labeledSpots = int(nucLabel).*int(maskSpot);
labeledSpots = nucLabel;
labeledSpots(maskSpot==0) = 0;
[spotPropStruct] = Spots.spotProp3(nucLabel, im, labeledSpots, metaDataDS);
spotLabel = labeledSpots;
end



