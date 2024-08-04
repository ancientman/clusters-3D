function [spotLabel, spotProp] = spotNoNuc3(im, channelFilled, mask, metaDataDS)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
voxVolume = xPixUM*yPixUM*zPixUM;
minSpotVol = metaDataDS.analysisInfo.minSpotVol;
imUseType = metaDataDS.analysisInfo.imUseType;
clearBorder = metaDataDS.analysisInfo.clearBorder;
thresholdFraction = 0.8; %%%%%%%%%%%%%%%%%%%%
spotProp = struct([]);

% colorChannels = metaDataDS.imagingInfo.channelCount;
if imUseType==1 % use "raw"
    imUse = im;
elseif imUseType==2 % use "smooth"
    imUse = Preprocess.smoothRaw3(im, metaDataDS);
elseif imUseType==3 % use "sharp"
    imUse = Preprocess.sharpRaw3(im, metaDataDS);
else % use "raw"
    imUse = im;
end 

if channelFilled == 1
    [spotVal, totalSpots] = Spots.adapThres3(im, mask);
    spotLabel = imbinarize(spotVal);
else
    imFilt = medfilt3(rescale(imUse));

    H1 = fspecial3('gaussian',[13, 13, size(imFilt, 3) ], 13);
    H2 = fspecial3('gaussian',[13, 13, size(imFilt, 3) ], 9);
    dog = H2 - H1;
    imDogFilt = rescale(convn(imFilt,dog,'same'));
    spotBin = imbinarize(imDogFilt, thresholdFraction);
    
    clearBorder = 1; %%%%%%%%%%%%%%%%%%
    if clearBorder==1
        spotBin = imclearborder(spotBin);
    end
    spotBin = bwareaopen(spotBin, ceil(minSpotVol/2));
    
    spotLabel = bwlabeln(spotBin);
end

CC = bwconncomp(spotLabel);
s = regionprops3(CC, im, "Centroid","PrincipalAxisLength", ...
    "Volume", 'EquivDiameter',  "BoundingBox", "VoxelIdxList", 'VoxelList', "VoxelValues", "WeightedCentroid");
spotProp(1).center = s.Centroid;
spotProp(1).paLen = s.PrincipalAxisLength;
spotProp(1).bb = s.BoundingBox;
spotProp(1).vol = s.Volume; % in voxels
spotProp(1).dia = s.EquivDiameter; % in voxels
spotProp(1).volUM = voxVolume*s.Volume; % in um^3
spotProp(1).voxIdx = s.VoxelIdxList;   
spotProp(1).voxList = s.VoxelList;
spotProp(1).voxVal = s.VoxelValues; % Value of the voxels in the region
spotProp(1).voxValCenter= s.WeightedCentroid; %Center of the region based on location and intensity value
end