function [hotSpotStruct] = globalHotSpotFinder(im, channelFilled, nucPropStruct, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bw = Nucleus.nucDetect(im, thres, nucleusDetectMode, nucEl, minNucSize, clearBorder, 0, 1);
labelMat = bwlabel(bw);
nucLabelStats = regionprops(labelMat, 'Area', 'Centroid');
[newGlobalNucLabelMat] = Spots.globalHotSpotLabelCorrect(labelMat, nucLabelStats, nucPropStruct);
labelMat = newGlobalNucLabelMat; % reassign the labels after checking if consistent with the raw frames
nucLabel = nonzeros(unique(labelMat, 'sorted'));
nucLabelStats = regionprops(labelMat, 'Area', 'Centroid');
spotDetect = 4; % Mode #4 is dedicated for hotspot detection
spotFilter = 1; % Gaussian (use only gaussian or bilateral)
% minSpotSize = 2*minSpotSize; % Double the size of individual spots
% Segment global hotspots
[bwHotSpot, hotSpotLabelMatNL] = Spots.spotDetect(imScale, labelMat, spotEl, spotDetect, spotFilter, minSpotSize);
% sort nuclei with detected spots
uniqSpotLabelNL = nonzeros(unique(hotSpotLabelMatNL, 'sorted'));
noSpotNuc = setdiff(nucLabel, uniqSpotLabelNL);
if ~isempty(noSpotNuc)
    uniqSpotLabelNL = HelperFunctions.vectorInsertAfter(uniqSpotLabelNL, 0, noSpotNuc);
    uniqSpotLabelNL = sort(uniqSpotLabelNL);
end
% compute the properties of individual spots, and asign them to cell array,
% with each array element representing one nucleus
spotCentroidUniq = cell(1, length(uniqSpotLabelNL));
spotAreaUniq = cell(1, length(uniqSpotLabelNL));
spotPixListUniq = cell(1, length(uniqSpotLabelNL));
spotIdxListUniq = cell(1, length(uniqSpotLabelNL));
spotNLtemp = hotSpotLabelMatNL;
for i=1:length(uniqSpotLabelNL)
    spotNLtemp(spotNLtemp~=uniqSpotLabelNL(i)) = 0;
    spotNLtemp(spotNLtemp==uniqSpotLabelNL(i)) = 1;
    CC = bwconncomp(spotNLtemp);
    stats = regionprops(CC, 'Centroid', 'Area', 'PixelList', 'PixelIdxList');
    spotCentroidUniq{i} = cat(1,stats.Centroid);
    spotAreaUniq{i} = cat(1,stats.Area);
    spotPixListUniq{i} = cat(1, stats.PixelList);
    spotIdxListUniq{i} = cat(1, stats.PixelIdxList);
    spotNLtemp = hotSpotLabelMatNL;
end
hotspotAreaNL = regionprops(hotSpotLabelMatNL,'Area');
hotSpotStruct.globalNucLabel = labelMat;
hotSpotStruct.globalNucAreaNL = cat(1, nucLabelStats.Area);
hotSpotStruct.globalNucCentroidNL = cat(1, nucLabelStats.Centroid);
hotSpotStruct.hotspotBW = bwHotSpot;
hotSpotStruct.hotspotNL = hotSpotLabelMatNL;
hotSpotStruct.hotspotAreaNL = cell2mat(struct2cell(hotspotAreaNL));
hotSpotStruct.hotspotCentroidUniq = spotCentroidUniq;
hotSpotStruct.hotspotAreaUniq = spotAreaUniq;
hotSpotStruct.hotspotPixListUniq = spotPixListUniq;
hotSpotStruct.hotspotIdxListUniq = spotIdxListUniq;
hotSpotStruct.noSpotNuc = noSpotNuc;
% plotHotSpotBoundaries(im, bwHotSpot);
end

function plotHotSpotBoundaries(im, bw)
[B,L] = bwboundaries(bw,'noholes');
figure('Color', [1 1 1]);
imshow(im,[]);
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
end
hold off;
end