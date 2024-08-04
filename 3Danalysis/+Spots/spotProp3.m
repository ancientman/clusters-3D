function [spotProp3] = spotProp3(nucLabel, im, labeledSpotMat, metaDataDS)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
voxVolume = xPixUM*yPixUM*zPixUM;
nNuc = max(nucLabel,[],'all');

% nucLabels = nonzeros(unique(nucLabel));
% spotLabels = nonzeros(unique(labeledSpotMat));
% nNuc = max(nucLabels);

CC = cell(1, nNuc);
spotProp3 = struct([]);
for i=1:nNuc
    if ismember(i, labeledSpotMat)
        CC{i} = bwconncomp(labeledSpotMat==i);
        s = regionprops3(CC{i}, im, "Centroid","PrincipalAxisLength", ...
            "Volume", 'EquivDiameter',  "BoundingBox", "VoxelIdxList", 'VoxelList', "VoxelValues", "WeightedCentroid");
        spotProp3(i).center = s.Centroid;
        spotProp3(i).paLen = s.PrincipalAxisLength;
        spotProp3(i).bb = s.BoundingBox;
        spotProp3(i).vol = s.Volume; % in voxels
        spotProp3(i).dia = s.EquivDiameter; % in voxels
        spotProp3(i).volUM = voxVolume*s.Volume; % in um^3
        spotProp3(i).voxIdx = s.VoxelIdxList;   
        spotProp3(i).voxList = s.VoxelList;
        spotProp3(i).voxVal = s.VoxelValues; % Value of the voxels in the region
        spotProp3(i).voxValCenter= s.WeightedCentroid; %Center of the region based on location and intensity value
    else
        spotProp3(i).center = [];
        spotProp3(i).paLen = [];
        spotProp3(i).bb = [];
        spotProp3(i).vol = [];
        spotProp3(i).dia = [];
        spotProp3(i).volUM = [];
        spotProp3(i).voxIdx = [];
        spotProp3(i).voxList = [];
        spotProp3(i).voxVal = [];
        spotProp3(i).voxValCenter= [];
    end
end
end
