function [seriesMetaDataDS] = InitializeOtherParams(seriesMetaDataDS)
xPixUM = seriesMetaDataDS.analysisInfo.xPixUM;
yPixUM = seriesMetaDataDS.analysisInfo.yPixUM;
zPixUM = seriesMetaDataDS.analysisInfo.zPixUM;
startZ = seriesMetaDataDS.analysisInfo.startZ;
endZ = seriesMetaDataDS.analysisInfo.endZ;
Zslices = seriesMetaDataDS.imagingInfo.stackSizeZ;
bicoidChannel = seriesMetaDataDS.imagingInfo.bicoidChannel;
mrnaSpotSizeFactor = seriesMetaDataDS.analysisInfo.mrnaSpotSizeFactor;
xyPsfUM = seriesMetaDataDS.imagingInfo.XYpsf/1000;
zPsfUM = seriesMetaDataDS.imagingInfo.Zpsf/1000;
minSpotSizeUM = (seriesMetaDataDS.analysisInfo.minSpotSize/1000);

if startZ==0 && endZ==0
    startZ = 1;
    endZ = Zslices;
    nZslices = Zslices;
elseif startZ==0 && endZ>0
    startZ = 1;
    if endZ<Zslices
        nZslices = endZ;
    else
    endZ = Zslices;
    nZslices = Zslices;
    end
elseif startZ>0 && endZ>0
    if startZ<Zslices && startZ<endZ && endZ<=Zslices
        nZslices = endZ-startZ+1;
    else
        error('check starting and or ending z slice');
    end
end
seriesMetaDataDS.imagingInfo.startZ = startZ;
seriesMetaDataDS.imagingInfo.endZ = endZ;
seriesMetaDataDS.analysisInfo.nZslices = nZslices; 

nucleusFeatureSize = seriesMetaDataDS.analysisInfo.nucleusFeatureSize;% already in microns
minNucAreaPix = ceil(3*(nucleusFeatureSize^2)/(xPixUM*yPixUM)); %all in microns
seriesMetaDataDS.analysisInfo.minNucSize = minNucAreaPix;%%%%%%%% Assignment in pixels here

if (2*nucleusFeatureSize)<=(nZslices*zPixUM)
    minNucVolPix = ceil(4*nucleusFeatureSize^3/(xPixUM*yPixUM*zPixUM)); %all in microns
elseif (2*nucleusFeatureSize)>(nZslices*zPixUM)
    minNucVolPix = ceil(minNucAreaPix*(2*nucleusFeatureSize/zPixUM));
end

seriesMetaDataDS.analysisInfo.minNucVol = ceil(minNucVolPix);
seriesMetaDataDS.analysisInfo.xyPsfUM = xyPsfUM;
seriesMetaDataDS.analysisInfo.zPsfUM = zPsfUM;
seriesMetaDataDS.analysisInfo.minSpotSizeUM = minSpotSizeUM;

if minSpotSizeUM==0
    minSpotAreaPix = 0;
    minMRNAspotAreaPix = mrnaSpotSizeFactor*ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
elseif minSpotSizeUM>0 && minSpotSizeUM<xyPsfUM% Convert to pixels here
%     minSpotAreaPix = ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
    minSpotAreaPix = ceil((pi*(minSpotSizeUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
else
    minSpotAreaPix = ceil((pi*(minSpotSizeUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
end
seriesMetaDataDS.analysisInfo.minSpotArea = ceil(minSpotAreaPix); %%%%%%%% assignment in pixels
seriesMetaDataDS.analysisInfo.minMRNAspotArea = ceil(minMRNAspotAreaPix); %%%%%%%% assignment in pixels

if zPixUM<zPsfUM
    minSpotVolPix = ceil(minSpotAreaPix)*ceil(zPsfUM/zPixUM); %all in microns
    minSpotVolPix = ceil(minSpotAreaPix); % Only one plane spot for 2 colors
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
else
    minSpotVolPix = 2*ceil(minSpotAreaPix);
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
end
seriesMetaDataDS.analysisInfo.minSpotVol = ceil(minSpotVolPix); %%%%%%%% assignment in pixels
seriesMetaDataDS.analysisInfo.minMRNAspotVol = ceil(minMRNAspotVolPix); %%%%%%%% assignment in pixels
end