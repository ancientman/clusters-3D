function [spotLabel, spotPropStruct] = spotDetect3New(im, nucLabel, channelFilled, metaDataDS,  removeBorderSpots)
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

spotLabel = zeros(size(im));
spotPropStruct = struct([]);
timePoints = size(im, 4); 

for t = 1:timePoints
    if channelFilled==1
        [bwSpot, ~] = Spots.spotNucFilled3(imUse(:,:,:,t), nucLabel(:,:,:,t), metaDataDS); % based on deviations from base
%         [bwSpot, ~] = Spots.spotNucFilled32(imUse(:,:,:,t), nucLabel(:,:,:,t), metaDataDS,  removeBorderSpots); % based on deviations from base
    else
    %     [bwSpot, ~] = Spots.spotNucHollow3(imUse, nucLabel, metaDataDS); % based on deviations from base
        [bwSpot, ~] = Spots.spotNucHollow32(imUse(:,:,:,t), nucLabel(:,:,:,t), metaDataDS); % based on deviations from base
%         [bwSpot, ~] = Spots.spotNucHollow33(imUse(:,:,:,t), nucLabel(:,:,:,t), metaDataDS); % for bottleneck
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
    if channelFilled == 0
        if imUseType == 3
            maskSpot = bwareaopen(maskSpot,minSpotVol/3); % smaller mask for deconvolved
        else
            maskSpot = bwareaopen(maskSpot,minSpotVol);
        end
        maskSpot = imclearborder(maskSpot); %%%
        labeledSpots = maskSpot;
    elseif channelFilled == 1
        if imUseType == 3
            
            maskSpot = bwareaopen(maskSpot,minSpotVol/3); % smaller mask for deconvolved
        else
            maskSpot = bwareaopen(maskSpot,minSpotVol);
        end
        maskSpot = imclearborder(maskSpot);
        labeledSpots = nucLabel(:,:,:,t).*maskSpot;
    end
    
    [spotPropStruct{t}] = Spots.spotProp3(nucLabel(:,:,:,t), im(:,:,:,t), labeledSpots, metaDataDS);
    spotLabel(:,:,:,t) = labeledSpots;
end

end
