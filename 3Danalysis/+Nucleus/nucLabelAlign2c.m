function [outLabelStack, nucLabelProp] = nucLabelAlign2c(c2BinStack, c2Stack, c1LabelStack, c1Prop, metaDataDS)

end

function [nucProp] = nucProp(nucCC)
nNuc = nucCC.NumObjects;
nucProp = struct([]);
for i=1:nNuc
    s = regionprops3(nucCC, "Centroid", "VoxelIdxList", "VoxelList");
    nucProp(i).center = s.Centroid;
    nucProp(i).voxIdx = s.VoxelIdxList;
    nucProp(i).voxList = s.VoxelList;
end

end