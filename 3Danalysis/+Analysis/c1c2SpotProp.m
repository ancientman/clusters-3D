function [C1C2SpotStruct] = c1c2SpotProp(c1NucStruct, c1SpotStruct, c2SpotStruct, metaDataDS, t)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
if (length(c1SpotStruct)~=length(c2SpotStruct))
    error('unequal spot labels at t = %d\n', t);
end
nNuc = length(c1SpotStruct);
C1C2SpotDistUM = cell(1, nNuc);
C1C2SpotVolUM = cell(1, nNuc);
C1C2SpotVal = cell(1, nNuc);
C1C2NucVal = cell(1, nNuc);
for i=1:nNuc
    if ~isempty(c2SpotStruct(i).center) && ~isempty(c1SpotStruct(i).center)
        if length(c2SpotStruct(i).vol)==1
            c1SpotCentUM = [xPixUM*c1SpotStruct(i).voxValCenter(:,1), yPixUM*c1SpotStruct(i).voxValCenter(:,2), zPixUM*c1SpotStruct(i).voxValCenter(:,3)];
            c2SpotCentUM = [xPixUM*c2SpotStruct(i).voxValCenter(:,1), yPixUM*c2SpotStruct(i).voxValCenter(:,2), zPixUM*c2SpotStruct(i).voxValCenter(:,3)];
            C1C2SpotDistUM{i} = pdist2(c1SpotCentUM, c2SpotCentUM);
            C1C2SpotVolUM{i} = c1SpotStruct(i).volUM;
            if ~isempty(c1SpotStruct(i).voxVal)
                C1C2SpotVal{i} = cellfun(@mean, c1SpotStruct(i).voxVal);
            end
            C1C2NucVal{i} = mean(c1NucStruct(i).voxVal{:});            
        end
    end
end
C1C2SpotStruct.C1C2SpotDistUM = C1C2SpotDistUM;
C1C2SpotStruct.C1C2SpotVolUM = C1C2SpotVolUM;
C1C2SpotStruct.C1C2SpotVal = C1C2SpotVal;
C1C2SpotStruct.C1C2NucVal = C1C2NucVal;
end