function [adjustFilledNucLabelMat, hollowLabelMat, nucFilledProp] = adjustNucLabeler3(filledMask, hollowHull, im, primeHollowLabelMat, metaDataDS, t)
centroidShiftConstraint = metaDataDS.analysisInfo.centroidShiftConstraint;
minNucVol = metaDataDS.analysisInfo.minNucVol;
hollowCC = bwconncomp(hollowHull);
if t==1
    hollowLabelMat = labelmatrix(hollowCC);
    [hollowLabels] = HelperFunctions.count_unique(nonzeros(hollowLabelMat));
    adjustFilledNucMask = filledMask.*hollowHull;
    adjustFilledNucMask = bwareaopen(adjustFilledNucMask, ceil(minNucVol/3));
    adjustFilledNucLabelMat = uint8(adjustFilledNucMask).*uint8(hollowLabelMat);
    [filledLabels] = HelperFunctions.count_unique(nonzeros(adjustFilledNucLabelMat));
    orphanHollowLabels = setxor(hollowLabels, filledLabels);

    if~isempty(orphanHollowLabels)
        for j=1:size(orphanHollowLabels,1)
            hollowHull(hollowCC.PixelIdxList{orphanHollowLabels(j)}) = 0;
        end
    end
    hollowCC = bwconncomp(hollowHull);
    hollowLabelMat = labelmatrix(hollowCC);
    adjustFilledNucLabelMat = uint8(adjustFilledNucMask).*uint8(hollowLabelMat);
else    
    
    [primeHollowLabels] = HelperFunctions.count_unique(nonzeros(primeHollowLabelMat));
    totalPrimeHollowNuc = max(primeHollowLabels);
    primeHollowStat = regionprops3(primeHollowLabelMat, 'Centroid'); 
    primeHollowNucCenter = primeHollowStat.Centroid;
    hollowStat = regionprops3(hollowCC, 'Centroid', 'VoxelIdxList');
    hollowNucCenter = hollowStat.Centroid;
    hollowNucVoxel = hollowStat.VoxelIdxList;
    hollowLabelMat = zeros(size(hollowHull));
    for i=1:totalPrimeHollowNuc
        centerDist = pdist2(primeHollowNucCenter(i,:), hollowNucCenter, 'euclidean');
        [minDist,minInd] = min(centerDist);
        if minDist<centroidShiftConstraint
            hollowLabelMat(hollowNucVoxel{minInd}) = i;
        end
    end
    hollowHull(~hollowLabelMat) = 0;
    [hollowLabels] = HelperFunctions.count_unique(nonzeros(hollowLabelMat));
    adjustFilledNucMask = filledMask.*hollowHull;
    adjustFilledNucMask = bwareaopen(adjustFilledNucMask, ceil(minNucVol/3));
    adjustFilledNucLabelMat = uint8(adjustFilledNucMask).*uint8(hollowLabelMat);
    [filledLabels] = nonzeros(HelperFunctions.count_unique(adjustFilledNucLabelMat));
    orphanHollowLabels = setxor(hollowLabels, filledLabels);
    if ~isempty(orphanHollowLabels)
        for i=1:length(orphanHollowLabels)
            hollowHull(hollowCC.PixelIdxList{orphanHollowLabels(i)}) = 0;
            hollowLabelMat(hollowCC.PixelIdxList{orphanHollowLabels(i)}) = 0;
        end
    end    
end
[nucFilledProp] = Nucleus.nucProp3(bwconncomp(adjustFilledNucMask), im, metaDataDS);
end


