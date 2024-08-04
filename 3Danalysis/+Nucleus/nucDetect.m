function [fim] = nucDetect(im, filledNuc, metaDataDS)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
xyPsfUM = metaDataDS.analysisInfo.xyPsfUM;
zPsfUM = metaDataDS.analysisInfo.zPsfUM;
imUseType = metaDataDS.analysisInfo.imUseType;
elementSize = metaDataDS.analysisInfo.elementSize;
clearBorder = metaDataDS.analysisInfo.clearBorder;
nucleusDetectMode = metaDataDS.analysisInfo.nucleusDetectMode;
nucThres = metaDataDS.analysisInfo.nucThres;
nucleusFeatureSize = metaDataDS.analysisInfo.nucleusFeatureSize;% already in microns
centroidShiftConstraint = metaDataDS.analysisInfo.centroidShiftConstraint; %permissible shift between centroids in consecutive time points
smoothFilterType = metaDataDS.analysisInfo.nucSmoothFilter;
smoothingParam = metaDataDS.analysisInfo.smoothingParam;
deconvolutionType = metaDataDS.analysisInfo.deconvolutionType;
minNucSize = ceil(nucleusFeatureSize^2/(3*xPixUM*yPixUM)); %all in microns

if imUseType==1 % use "raw"
    imUse = im;
elseif imUseType==2 % use "smooth"
    imUse = Preprocess.smoothRawImage(im, smoothFilterType, smoothingParam);
    if smoothFilterType == 7
        warning ('using abandoned filter type: DoG');
    end
elseif imUseType==3 % use "sharp"
    sharpParam = ceil(xyPsfUM/(2.0*xPixUM));
    imUse = Preprocess.sharpRawImage(im, deconvolutionType, sharpParam);
else % use "raw"
    imUse = im;
end 

if nucleusDetectMode==1 
    if filledNuc==1
        nucMaskFilled = Nucleus.segmentFilledNuc(imUse, elementSize, minNucSize, nucThres, clearBorder);
        nucBW = nucMaskFilled;
    elseif filledNuc==0
        nucMaskHollow = Nucleus.segmentHollowNuc(imUse, elementSize, minNucSize, nucThres, clearBorder);
        nucBW = nucMaskHollow;
    end
end
% [nucBW, thres] = Nucleus.noisyNuclearMask(imUse, nucBW, nucThres, minNucSize, clearBorder);
% nucBW = bwareaopen(nucBW, minNucSize);

se = strel('disk', elementSize);
Io = imopen(nucBW, se);
If = imfill(Io, 'holes');
if filledNuc==0
    nucBW = bwareaopen(If, minNucSize/3);    
    fim = nucBW;
else
    fim = If;
end
end