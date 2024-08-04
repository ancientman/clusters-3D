function [adjustFilledNucMask, adjustFilledNucLabel] = adjustFilledNucSeg(filledNucMask, hollowNucMask)
hollowNucConvHull = bwconvhull(hollowNucMask,'objects');
hollowCC = bwconncomp(hollowNucConvHull);
hollowNucLabel = labelmatrix(hollowCC);
adjustFilledNucMask = filledNucMask.*hollowNucConvHull;
adjustFilledNucLabel = uint8(filledNucMask).*uint8(hollowNucLabel);
end