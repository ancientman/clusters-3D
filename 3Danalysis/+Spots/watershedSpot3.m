function [fim] = watershedSpot3(bw, metaDataDS)
elementSize = metaDataDS.analysisInfo.elementSize;
bw = imtophat(bw, strel('cube', 11));
D = bwdist(~bw);
D = -D;
D(~bw) = -Inf;
L = watershed(D);
bw2 = bw;
bw2(L==0) = 0;
fim = bw2;
end