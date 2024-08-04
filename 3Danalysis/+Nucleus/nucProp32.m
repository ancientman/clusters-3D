function [nucProp3] = nucProp32(nucLabelMat, im, metaDataDS)
%%%%%%%%%%%%%%%%%%%
% has issues: check
%%%%%%%%%%%%%%%%%%%
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
voxVolume = xPixUM*yPixUM*zPixUM;
nNuc = max(nucLabelMat,[], 'all');
CC = cell(1, nNuc);
nucProp3 = struct([]);
for i=1:nNuc
    CC{i} = bwconncomp(nucLabelMat==i);
    s = regionprops3(CC{i}, im, "Centroid","PrincipalAxisLength",...
        "Volume", "EquivDiameter",  "BoundingBox", "VoxelIdxList", "VoxelList", "VoxelValues", "WeightedCentroid");
    nucProp3(i).center = s.Centroid;
    nucProp3(i).paLen = s.PrincipalAxisLength;
    nucProp3(i).bb = s.BoundingBox;
    nucProp3(i).vol = s.Volume; % in voxels
    nucProp3(i).dia = s.EquivDiameter; % in voxels
    nucProp3(i).volUM = voxVolume*s.Volume; % in um^3
    nucProp3(i).voxIdx = s.VoxelIdxList;
    nucProp3(i).voxList = s.VoxelList;
    nucProp3(i).voxVal = s.VoxelValues; % Value of the voxels in the region
    nucProp3(i).voxValCenter= s.WeightedCentroid; %Center of the region based on location and intensity value
    nucProp3(i).meanVal = mean(nucProp3(i).voxVal{:}); %mean of all nuclear values
end
end