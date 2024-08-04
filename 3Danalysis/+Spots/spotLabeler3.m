function [maskSpot, labelSpot] = spotLabeler3(nucLabel, bwSpot, minSpotVol)
CC = bwconncomp(nucBinCC);
labelMat = labelmatrix(CC);
maskSpot = bwareaopen(bwSpot,minSpotVol);
labelSpot = labelMat.*maskSpot;
end
