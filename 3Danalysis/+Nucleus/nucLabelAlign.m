function [outLabelStack, nucLabelProp] = nucLabelAlign(imBinStack, im, metaDataDS)

defaultLabStack = zeros(size(imBinStack));
outLabelStack = zeros(size(imBinStack));
shiftLim = metaDataDS.analysisInfo.centroidShiftConstraint;
nTimePoints = size(imBinStack, 4);
centList = cell(nTimePoints, 1);
voxList = cell(nTimePoints, 1);
labIdxList = cell(nTimePoints, 1);
acceptedSpots = cell(nTimePoints, 1);
idRTemp = cell(nTimePoints, 1);
idCTemp = cell(nTimePoints, 1);
nucLabelProp = struct([]);


% get the centers and indices of the connected components
for t = 1:nTimePoints
    labStackTemp = labelmatrix(bwconncomp(imBinStack(:,:,:,t)));
    nLabels = max(labStackTemp, [], 'all');
    centTemp = zeros(nLabels, 3);
    defaultLabStack(:,:,:,t) = labStackTemp;
    for i = 1:nLabels
        CC = labStackTemp;
        CC(labStackTemp ~= i) = 0;
        CC(labStackTemp == i) = 1;
        s = regionprops3(CC,  "Centroid", "VoxelIdxList");
        voxListTemp = s.VoxelIdxList;
        voxList{t}{i} = voxListTemp;
        centTemp(i,:) = s.Centroid;
        centList{t} = round(centTemp);
    end    
    % Distance matrix of the centers in each frame with the previous
    if t==1
        centDist = pdist2(centList{1}, centList{t});        
    else
        centDist = pdist2(centList{t-1}, centList{t});
    end
    % Only keep the distances meeting limit criterion
    centDist(centDist>shiftLim) = NaN;
    idx = find(~isnan(centDist)); % indices of all permissible centers
    idRTemp{t} = zeros(size(idx)); % Earlier time (rows of cent dist)
    idCTemp{t} = zeros(size(idx)); % Current time (cols of cent dist)
    for i = 1:length(idx)
        [idRTemp{t}(i), idCTemp{t}(i)] = ind2sub(size(centDist),idx(i));
    end
    % same {t} nucleus mapping to multiple (1) nuclei
    [~, w] = unique(idRTemp{t}, 'stable' );
    dupInd = setdiff( 1:numel(idRTemp{t}), w );
%     dupNuc1 = nuc1Temp{t}(dupInd);
     idCTemp{t}(dupInd) = 0; % set to zero
    % same (1) nucleus mapping to multiple {t} nuclei
    [~, w] = unique(idCTemp{t}, 'stable' );
    dupInd = setdiff( 1:numel(idCTemp{t}), w );
    idCTemp{t}(dupInd) = 0; % set to zero
% %     dupNucT = nuc1Temp{t}(dupInd);
%     allRejectedSpots = [dupNuc1; dupNucT];
%     allRejectedSpots = unique(allRejectedSpots);
%     acceptedSpots{t} = setdiff(nuc1Temp{t}, allRejectedSpots);
end

% for t = 1:nTimePoints
%     if t == 1
%         commonAccepted = nonzeros(acceptedSpots{t});
%     else
%         commonAccepted = nonzeros(intersect(commonAccepted, acceptedSpots{t}));
%     end
% end

idNuc = cell(nTimePoints, 1);

for t = 1:nTimePoints
    
    for i = 1:length(idCTemp{t})
        if t < 3
%             if ismember(nuc1Temp{t}(i), commonAccepted) && ismember(nucTTemp{t}(i), commonAccepted)
               idNuc{t}(idCTemp{t}(i)) = idRTemp{t}(i);
%             else
%                 nucIDs{t}(nucTTemp{t}(i)) = 0;
%             end
        else
%             if ismember(nucTTemp{t-1}(i), commonAccepted) 
try 
    idNuc{t}(idCTemp{t}(i)) = idNuc{t-1}(idRTemp{t}(i));
catch
    aaa = t
end
%             else
%                 nucIDs{t}(nucTTemp{t}(i)) = 0;
%             end
        end
    end 
end

for t = 1:nTimePoints
    if t == 1
        commonAccepted = nonzeros(idNuc{t});
    else
        commonAccepted = nonzeros(intersect(idNuc{t}, idNuc{t-1}));
    end
end

for t = 1:nTimePoints
    [~,diffIdx] = setdiff(idNuc{t}, commonAccepted, 'stable');
    idNuc{t}(diffIdx) = 0;
end

if max(idNuc{1}) ~= length(nonzeros(idNuc{1}))        
    tempOldSort = nonzeros(sort(idNuc{t}));
    tempNewSort = 1:length(tempOldSort);        
else
    tempOldSort = nonzeros(sort(idNuc{t}));  
    tempNewSort = tempOldSort;
end

nucIDf = idNuc; %cell(nTimePoints, 1);
for t = 1:nTimePoints
    outLabTemp = zeros(size(labStackTemp));
   for i = 1:length(idNuc{t})
       tempIdx = find(tempOldSort==idNuc{t}(i));
       if ~isempty(tempIdx)
        nucIDf{t}(i) = tempNewSort(tempIdx);
       end
       outLabTemp(voxList{t}{(i)}{1}) = nucIDf{t}(i);
   end
   outLabelStack(:,:,:,t) = outLabTemp;
   [nucLabelProp{t}] = Nucleus.nucProp32(outLabelStack(:,:,:,t), im(:,:,:,t), metaDataDS);     
end


% for t = 1:nTimePoints
%     outLabTemp = zeros(size(labStackTemp));
%     for i = 1:length(nucTTemp{t})
%         if ismember(i, commonAccepted)
%             outLabTemp(voxList{t}{nucTTemp{t}(i)}{1}) = nuc1Temp{t}(i);
%         end
%     end 
%     outLabelStack(:,:,:,t) = outLabTemp;
%     [nucLabelProp{t}] = nucProp3(outLabelStack(:,:,:,t), im(:,:,:,t), metaDataDS);     
% %      [nucLabelProp{t}] = nucProp3(defaultLabStack(:,:,:,t), im(:,:,:,t), metaDataDS);     
% end

%%%%%%%%%%%% Visualize labels %%%%%%%%%%%%%%%%%%%%

% viewLabels(nucLabelProp, outLabelStack);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     centDist = pdist2(centList{1}, centList{t});
%     [minDist, minIdx] = min(centDist,[], 2);
%     minIdx(minDist>shiftLim) = 0;
%     labIdxList{t} = minIdx;
%     
%     for i = 1:length(labIdxList{t})
%         if labIdxList{t}(i) ~=0
%             outLabTemp(voxList{t}{labIdxList{t}(i)}{1}) = i;
%         end
%     end    
%     
%     outLabelStack(:,:,:,t) = outLabTemp;
%     if t<3
%         commonLabel = intersect(labIdxList{1}, labIdxList{t});
%     else
%         commonLabel = intersect(labIdxList{t}, commonLabel);
%     end   
% 
% for t = 1:nTimePoints
%     outTemp = outLabelStack(:,:,:,t);
%     delPoints = setdiff(nonzeros(unique(outTemp)), commonLabel, 'stable'); 
%     for i = 1:size(delPoints)
%          outTemp(outTemp==delPoints(i)) = 0;
%     end  
%     labels = nonzeros(unique(outTemp));    
%     Align the labels from 1 to n
%     for i = 1:length(labels)
%         outTemp(outTemp==labels(i)) = i;        
%     end
%     outLabelStack(:,:,:,t) = 0;
%     outLabelStack(:,:,:,t) = outTemp;
%     [nucLabelProp{t}] = nucProp3(outLabelStack(:,:,:,t), im(:,:,:,t), metaDataDS);       
% end

end

function [nucProp3] = nucProp3(nucLabelMat, im, metaDataDS)
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

function viewLabels(nucLabelProp, outLabelStack)
timePoints = length(nucLabelProp);
zSlices = size(outLabelStack, 3);
nucCentXY = cell(1, timePoints);
nucCentZ = cell(1, timePoints);
imagePlane = cell(1, timePoints);
displayLabImage = zeros([size(outLabelStack, [1, 2]), timePoints]);

for t = 1:timePoints
    for i = 1:length(nucLabelProp{t})
        if i == 1
            nucCent = nucLabelProp{t}(i).center;
        else
            nucCent = vertcat(nucCent, nucLabelProp{t}(i).center);
        end
    end
    nucCentXY{t} = round(nucCent(:,1:2));
    nucCentZ{t} = round(nucCent(:,3));
    zPlane = mle(nucCentZ{t});
    zPlane = round(zPlane(1));
    if zPlane<1
        zPlane = 1;
    elseif zPlane>zSlices
        zPlane = zSlices;
    end
    imagePlane{t} = zPlane;
    displayLabImage(:,:,t) = outLabelStack(:,:,imagePlane{t},t);
    nucText =  insertText(label2rgb(displayLabImage(:,:,t)), round(nucCentXY{t}), cellfun(@(x) num2str(x), num2cell(1:length(nucCent)), 'un', 0),'AnchorPoint','LeftBottom','FontSize',20);
    displayLabImage(:,:,t) = rgb2gray(nucText);
end
% figure; sliceViewer(displayLabImage);
end