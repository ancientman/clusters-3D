function [combineNucBin, combineNucProp] = nucBinStackCombiner(c1NucBin, c2NucBin, lastZ, metaDataDS)
combineNucBin = and(c1NucBin, c2NucBin);
stat = regionprops(bwconncomp(combineNucBin), 'Area', 'Centroid');
combineNucProp.centroid = cat(1, stat.Centroid);
combineNucProp.area = cat(1, stat.Area);
end
