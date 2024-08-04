function [c2SpotLabelStack, c2SpotProp] = c2SpotLabelAssign(c2SpotBinStack, c2Stack, c1NucLabelStack, c1NucProp, metaDataDS)
c2SpotLabelStack = zeros(size(c2SpotBinStack));
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;

maxSpotDistUM = 3.5; % max ms2 spot distance from nuc center

c2SpotProp = struct([]);

for t = 1:size(c2SpotBinStack, 4)
    if any(c2SpotBinStack, 'all')
        c2CC = bwconncomp(c2SpotBinStack(:,:,:,t));
        c2Lab = zeros(size(c2SpotBinStack(:,:,:,t)));
        s = regionprops3(c2CC, "Centroid", "VoxelIdxList");
        c2SpotCent = s.Centroid;
        nNuc = length(c1NucProp{t});
        for i = 1:nNuc
            if i == 1
                nucCent = c1NucProp{t}(i).center;
            else
                nucCent = vertcat(nucCent, c1NucProp{t}(i).center);
            end
        end
        if ~isempty(nucCent)
            nucCentUM = [xPixUM*nucCent(:, 1), yPixUM*nucCent(:, 2), zPixUM*nucCent(:, 3)];
        else
            nucCentUM = [NaN, NaN, NaN];
        end
        if ~isempty(c2SpotCent)
            c2SpotCentUM = [xPixUM*c2SpotCent(:, 1), yPixUM*c2SpotCent(:, 2), zPixUM*c2SpotCent(:, 3)];
        else
            c2SpotCentUM = [NaN, NaN, NaN];
        end
        centDist = pdist2(c2SpotCentUM, nucCentUM);  
        centDist(centDist>maxSpotDistUM) = NaN;
        idx = find(~isnan(centDist));
        spotNoTemp = zeros(size(idx));
        nucNoTemp = zeros(size(idx));
        for i = 1:length(idx)
            [spotNoTemp(i), nucNoTemp(i)] = ind2sub([c2CC.NumObjects, nNuc],idx(i));
        end
        
        % Find same spots multi nuclei
        [~, w] = unique(spotNoTemp, 'stable' );
        dupInd = setdiff( 1:numel(spotNoTemp), w );
        dupSpot = spotNoTemp(dupInd);
        % Find multiple spots in the same nucleus
        [~, w] = unique(nucNoTemp, 'stable' );
        dupInd = setdiff( 1:numel(nucNoTemp), w );
        dupNucSpot = spotNoTemp(dupInd);
        
        allRejectedSpots = [dupSpot; dupNucSpot];
        allRejectedSpots = unique(allRejectedSpots);
        acceptedSpots = setdiff(spotNoTemp, allRejectedSpots);
        
        for i = 1:c2CC.NumObjects
            if ismember(i, acceptedSpots)
                c2Lab(s.VoxelIdxList{(i)}) = nucNoTemp(spotNoTemp ==i);
            end
        end
        c2SpotLabelStack(:,:,:,t) = c2Lab;        
    end    
%         centDist = pdist2(nucCentUM, c2SpotCentUM);     
%         [minVal, idx] = min(centDist,[], 2);         
%         tempIdx = idx;
%         [~, ia, ~] = unique(idx, 'stable');
%         repEl = idx(setdiff(1:nNuc, ia));
%         for i = 1:size(nucCent, 1)
%             if ismember(tempIdx(i), repEl)
%                 tempIdx(i) = 0;
%             end
%             if minVal(i) < maxSpotDistUM && tempIdx(i)~=0
%                 c2Lab(s.VoxelIdxList{tempIdx(i)}) = i;
%             end
%         end
%         c2SpotLabelStack(:,:,:,t) = c2Lab;        
%     end    
    [c2SpotProp{t}] = Spots.spotProp3(c1NucLabelStack(:,:,:,t), c2Stack(:,:,:,t), c2SpotLabelStack(:,:,:,t), metaDataDS);
    
%     for i = 1:length(c2SpotProp{t})
%         if i == 1
%             spotCent = c2SpotProp{t}(i).center;
%         else
%             spotCent = vertcat(spotCent, c2SpotProp{t}(i).center);
%         end
%     end
%     viewLabels(nucCent, spotCent, c1NucLabelStack(:,:,:,t), c2SpotProp{t}, c2Lab);
end
end

function viewLabels(nucCent, spotCent, c1NucLabelStack, c2SpotProp, c2Lab)
nucSpot = zeros(size(c2Lab));
c2LabView = zeros([size(c2Lab, 1), size(c2Lab, 2), 3, size(c2Lab, 3)]);
nucSpotPerimView = zeros([size(c2Lab, 1), size(c2Lab, 2), 3, size(c2Lab, 3)]);
spotLabels = nonzeros(unique(c2Lab));
if ~isempty(spotLabels)
    for i = 1:size(c2Lab, 3)  
        text = vertcat(cellfun(@(x) num2str(x), num2cell(spotLabels), 'un', 0), cellfun(@(x) num2str(x), num2cell(1:length(nucCent))', 'un', 0));
        textPos = [spotCent(:, 1:2); nucCent(:, 1:2)];
        boxCol = vertcat(repmat({'red'}, size(spotCent, 1), 1), repmat({'green'}, size(nucCent, 1), 1));

        nuc = c1NucLabelStack(:,:,i);
        spot = c2Lab(:,:,i);
        nn = label2rgb(imbinarize(nuc),'jet','k');
        ss = label2rgb(imbinarize(spot),'jet','k');
        nucSpot(:,:,i) = (bwperim(nuc) | bwperim(spot));
        c2LabView(:,:,:,i) =  insertText(label2rgb(nuc), round(nucCent(:, 1:2)), cellfun(@(x) num2str(x), num2cell(1:length(nucCent)), 'un', 0),'AnchorPoint','LeftBottom','FontSize',20);
        nucSpotPerimView(:,:,:,i) = insertText(imfuse(nn, ss, 'ColorChannels', 'red-cyan'), textPos, text,'AnchorPoint','LeftBottom','FontSize',20, 'TextColor', boxCol);
    %     nucText =  insertText(label2rgb(nuc), round(nucCent(:, 1:2)), cellfun(@(x) num2str(x), num2cell(1:length(nucCent)), 'un', 0),'AnchorPoint','LeftBottom','FontSize',20);
    %     spotText =  insertText(label2rgb(spot), round(spotCent(:, 1:2)), cellfun(@(x) num2str(x), num2cell(spotLabels), 'un', 0),'AnchorPoint','LeftBottom','FontSize',20);
    end
end
end

















