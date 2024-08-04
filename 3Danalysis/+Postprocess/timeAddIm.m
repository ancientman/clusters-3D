function [addNucLabelProp, addSpotLabelProp] = timeAddIm(im, nucLabelStack, spotLabelStack, metaDataDS, channelFilled, removeBorderSpots)
timePoints = size(im, 4);
zSlices = size(im, 3);
nNuc = length(nonzeros(HelperFunctions.count_unique(nucLabelStack(:,:,:,:))));
nSpot = length(nonzeros(HelperFunctions.count_unique(spotLabelStack(:,:,:,:))));
bgMean = double(zeros(timePoints, nNuc));
bgStd = double(zeros(timePoints, nNuc));
imSpotAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
imThresAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
imNucAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
imThresMaskAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
timeAddINucLabel = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
addNucLabelProp = struct([]);
addSpotLabelProp = struct([]);

if  nNuc ==0 && nSpot ~= 0
    for i = 1:nSpot
        for t=1:timePoints
            imTemp = rescale(im(:,:,:,t));
            imTemp(spotLabelStack(:,:,:,t)~=i) = 0;
            imTemp = rescale(medfilt3(imTemp)); % Normalize
            imSpotAdd = imThresAdd + double(imTemp);
            imThresAdd = imSpotAdd;
        end
    end
else
    for i = 1:nNuc
        imJustNuc = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
        for t=1:timePoints
            imTemp = im(:,:,:,t);
            imTemp = (medfilt3(imTemp)); 
            imTemp(nucLabelStack(:,:,:,t)~=i) = 0;
            for z = 1:zSlices
                imTemp(:,:,z) = rescale(imTemp(:,:,z));
            end
            imTempNucNoSpots = imTemp;  
            imTempSpots = imTemp;
            imTempSpots(spotLabelStack(:,:,:,t)~=i) = 0; % spot pixels from a nucleus
            imTempNucNoSpots(spotLabelStack(:,:,:,t)==i) = 0; % background pixels from a nucleus
            bgMean(t, i) = mean(nonzeros(imTempNucNoSpots),'all', 'omitnan');
            bgStd(t, i) = std(nonzeros(imTempNucNoSpots),0,'all', 'omitnan');
            imTempThres = imTemp - (bgMean(t, i) + 2.0*bgStd(t, i));
            imTempThres(imTempThres<0) = 0;
            imTempThres(isnan(imTempThres)) = 0;
            imThresAdd = imThresAdd + double(imTempThres);
            imSpotAdd = imSpotAdd + rescale(imTempSpots);
            imJustNuc = imJustNuc + imTemp;
            imNucAdd = imNucAdd + imTemp;
        end
        timeAddINucLabel(imJustNuc~=0) = i;
        imThresMaskAdd = imThresMaskAdd+imbinarize(imThresAdd);
    end
end

timeAddSpot = imSpotAdd; %rescale(imSpotAdd);
timeAddThres = imThresAdd; %rescale(imThresAdd);
timeAddNuc = imNucAdd; %rescale(imNucAdd);

[addNucLabelProp] = Nucleus.nucProp3(bwconncomp(timeAddINucLabel), timeAddNuc, metaDataDS);       
[spotLabel, addSpotLabelProp] = Spots.spotDetect3New(timeAddThres, timeAddINucLabel, channelFilled, metaDataDS,  removeBorderSpots);

% sliceViewer((imStackNew + bwperim(imbinarize(squeeze(nucLabelStack(:,:,:,ceil(zSlices/2))))/5)));

end