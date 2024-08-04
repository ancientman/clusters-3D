function [spotPropNL] = spotProp(nucLabel, spotMask)
sCC = bwconncomp(spotMask);
spotProp = regionprops(sCC, 'Area', 'Centroid', 'PixelIdxList' );
spotCentroid = cat(1, spotProp.Centroid);
spotArea = cat(1, spotProp.Area);
spotPixIdx = cat(1, spotProp.PixelIdxList);
nucCentroid = struct2cell(regionprops(nucLabel, 'Centroid'));
nucPixIdx = struct2cell(regionprops(nucLabel, 'PixelIdxList'));
nNuc = length(nucPixIdx);

spotIdxNL = cellfun(@(x) intersect(x,spotPixIdx), nucPixIdx, 'UniformOutput',false);
spotAreaNL = cellfun(@(x) length(x), spotIdxNL, 'UniformOutput',false);
% [spotSubNLy, ~] = cellfun(@(x) ind2sub(sCC.ImageSize, x), spotIdxNL, 'UniformOutput',false);
% [~, spotSubNLx] = cellfun(@(x) ind2sub(sCC.ImageSize, x), spotIdxNL, 'UniformOutput',false);
spotCentroidNL = cell(1, nNuc);
tempMask = zeros(sCC.ImageSize);
for i = 1:nNuc
    tempMask(spotIdxNL{i}) = 1;
    prop = regionprops(bwconncomp(tempMask), 'Centroid');
    spotCentroidNL{i} = cat(1, prop.Centroid);
    tempMask = zeros(sCC.ImageSize);
end
spotPropNL.spotIdxNL = spotIdxNL;
spotPropNL.spotAreaNL = spotAreaNL;
spotPropNL.spotCentroidNL = spotCentroidNL;
end