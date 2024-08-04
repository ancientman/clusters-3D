function [fim, nucProp] = nucDetect3(im, filledNuc, metaDataDS)
imUse = im;

if filledNuc==1
    nucMaskFilled = Nucleus.segmentFilledNuc3(imUse, metaDataDS, filledNuc);
%     nucMaskFilled = Nucleus.segmentFilledNuc32(imUse, metaDataDS, filledNuc);
    nucBW = nucMaskFilled;
elseif filledNuc==0
    if metaDataDS.analysisInfo.nZslices>3
        nucMaskHollow = Nucleus.segmentHollowNuc3(imUse, metaDataDS, filledNuc);
    else
    nucMaskHollow = Nucleus.segmentHollowNuc32(imUse, metaDataDS, filledNuc);
%     nucMaskHollow = Nucleus.segmentHollowNucThin3(imUse, metaDataDS, filledNuc);
    end
    nucBW = nucMaskHollow;
end
fim = nucBW;
[nucProp] = Nucleus.nucProp3(bwconncomp(nucBW), im, metaDataDS);
end