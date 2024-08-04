function [nucBinStackProject, nucStackProp] = nucBinStacker(nucBinStackProject, nucBin, lastZ, metaDataDS)
nucBinStackProject = or(nucBinStackProject, nucBin);
nucStackProp.centroid = 0;
nucStackProp.area = 0;
if lastZ
	minNucSize = metaDataDS.analysisInfo.minNucSize;
    nucBinStackProject = imfill(nucBinStackProject, 'holes');
    clearBorder = metaDataDS.analysisInfo.clearBorder;
    if clearBorder
        nucBinStackProject = imclearborder(nucBinStackProject);
        nucBinStackProject = bwareaopen(nucBinStackProject, minNucSize);
    else
        nucBinStackProjectBorder = nucBinStackProject - imclearborder(nucBinStackProject);
        nucBinStackProjectIn = bwareaopen(imclearborder(nucBinStackProject), minNucSize);
        nucBinStackProject = or(nucBinStackProjectBorder, nucBinStackProjectIn);
    end
    stat = regionprops(bwconncomp(nucBinStackProject), 'Area', 'Centroid');
    nucStackProp.nCC = bwconncomp(nucBinStackProject);
    nucStackProp.centroid = cat(1, stat.Centroid);
    nucStackProp.area = cat(1, stat.Area);
end
end