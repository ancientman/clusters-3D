function [nucAddLabelProp, addSpotLabelProp, thresAddLabelProp, thresAdd] = timeRegAddIm(im, nucLabelStack, nucProp, spotLabelStack, metaDataDS, channelFilled, removeBorderSpots)
if channelFilled == 1
    howManySDs = 3;
elseif channelFilled == 0
    howManySDs = 5;
end
timePoints = size(im, 4);
zSlices = size(im, 3);
nNuc = max(nucLabelStack(:,:,:,1), [], 'all');
nSpot = length(nonzeros(HelperFunctions.count_unique(spotLabelStack(:,:,:,:))));
bgMean = double(zeros(timePoints, nNuc));
bgStd = double(zeros(timePoints, nNuc));
spotAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
thresAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
nucAdd = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
nucAddLabel = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
spotAddLabel = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
thresAddLabel =  double(zeros(size(im, 1),size(im, 2), size(im, 3)));
addSpotLabelProp = struct([]);
nucAddLabelProp = struct([]);
thresAddLabelProp = struct([]);

if  nNuc ==0 && nSpot ~= 0 % No nuc filter, no alignment needed
    for i = 1:nSpot
        for t=1:timePoints
            imTemp = rescale(im(:,:,:,t));
            imTemp(spotLabelStack(:,:,:,t)~=i) = 0;
            imTemp = rescale(medfilt3(imTemp)); % Normalize
            spotAdd = thresAdd + double(imTemp);
            thresAdd = spotAdd;
        end
    end
else % Align nuc for fast interval data only
    if channelFilled == 1
        for i = 1:nNuc        
            tempNuc = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            thresAddTemp = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            spotAddTemp = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            nucAddTemp = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            nucAddNoRegTemp = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            for t=1:timePoints
                imTemp = im(:,:,:,t);
                imTemp = (medfilt3(imTemp)); 
                imTemp(nucLabelStack(:,:,:,t)~=i) = 0;
                for z = 1:zSlices
                    imTemp(:,:,z) = rescale(imTemp(:,:,z));
                end
                tempNucNoSpots = imTemp;
                tempSpots = imTemp;
                tempSpots(spotLabelStack(:,:,:,t)~=i) = 0; % spot pixels from a nucleus
                tempNucNoSpots(spotLabelStack(:,:,:,t)==i) = 0; % background pixels from a nucleus
                bgMean(t, i) = mean(nonzeros(tempNucNoSpots),'all', 'omitnan');
                bgStd(t, i) = std(nonzeros(tempNucNoSpots),0,'all', 'omitnan');
                tempThres = imTemp - (bgMean(t, i) + howManySDs*bgStd(t, i));
                tempThres(tempThres<0) = 0;
                tempThres(isnan(tempThres)) = 0;
                tempThres = rescale(tempThres);
%                 if t ==1
%                     nucCentTemp1 = nucProp{t}(i).center;
%                 elseif t>1
%                     nucCentTemp = nucProp{t}(i).center;
%                     centShift = round(nucCentTemp1 - nucCentTemp);
%                     if metaDataDS.analysisInfo.analyze3D == 1
%                         tempNuc = imtranslate(imTemp, centShift);
%                         tempSpots = imtranslate(tempSpots, centShift);
%                         tempThres = imtranslate(tempThres, centShift);
%                     elseif metaDataDS.analysisInfo.analyze3D == 0
%                         tempNuc = imtranslate(imTemp, centShift(1:2));
%                         tempSpots = imtranslate(tempSpots, centShift(1:2));
%                         tempThres = imtranslate(tempThres, centShift(1:2));
%                     end
%                 end
                thresAddTemp = thresAddTemp + double(tempThres);
                spotAddTemp = rescale(spotAddTemp + rescale(tempSpots));
                nucAddTemp = rescale(nucAddTemp + rescale(tempNuc)); 
                nucAddNoRegTemp = rescale(nucAddNoRegTemp + imTemp);
            end
            spotAdd = spotAdd + spotAddTemp;
            thresAdd = thresAdd + thresAddTemp;
            nucAdd = nucAdd + nucAddTemp;
            nucAddLabel(nucAddTemp~=0) = i;
    %         timeAddINucLabel(imNucAddNoRegTemp~=0) = i;        
            thresAddLabel(thresAddTemp~=0) = i;
            spotAddLabel(spotAddTemp~=0) = i;
        end
    elseif channelFilled==0
        %%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:nNuc        
            spotAddTemp = double(zeros(size(im, 1),size(im, 2), size(im, 3)));
            for t=1:timePoints
                tempSpots = im(:,:,:,t);
                tempSpots = (medfilt3(tempSpots)); 
                tempSpots(spotLabelStack(:,:,:,t)~=i) = 0;
                tempSpots = rescale(tempSpots);
                if t ==1
                    nucCentTemp1 = nucProp{t}(i).center;
                elseif t>1
                    nucCentTemp = nucProp{t}(i).center;
                    centShift = round(nucCentTemp1 - nucCentTemp);
                    if metaDataDS.analysisInfo.analyze3D == 1
                        tempSpots = imtranslate(tempSpots, centShift);
                    elseif metaDataDS.analysisInfo.analyze3D == 0
                        tempSpots = imtranslate(tempSpots, centShift(1:2));
                    end
                end
                spotAddTemp = rescale(spotAddTemp + rescale(tempSpots));
            end
            spotAdd = spotAdd + spotAddTemp;    
            spotAddLabel(spotAddTemp~=0) = i;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
    end    
end


if channelFilled == 0
    [addSpotLabelProp] = Spots.spotProp3(nucLabelStack, spotAdd, spotAddLabel, metaDataDS);     
    thresAdd = spotAdd;
elseif channelFilled == 1
    timeAddThres = thresAdd; 
    timeAddNuc = nucAdd; 
    [nucAddLabelProp] = Nucleus.nucProp3(bwconncomp(nucAddLabel), timeAddNuc, metaDataDS);       
    [addSpotLabel, addSpotLabelProp] = Spots.spotDetect3New(nucAdd, nucAddLabel, channelFilled, metaDataDS,  removeBorderSpots);
    [thresAddLabel, thresAddLabelProp] = Spots.spotDetect3New(timeAddThres, nucAddLabel, channelFilled, metaDataDS,  removeBorderSpots);
end
end