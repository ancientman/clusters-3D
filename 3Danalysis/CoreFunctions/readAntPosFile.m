function [emDim] = readAntPosFile(path)
reader = bfGetReader(path);
omeGlobalMeta = reader.getGlobalMetadata();
seriesLength = str2double(omeGlobalMeta.get(strcat('Information|Image|SizeS #1')));	
if seriesLength ~=2
    seriesLength = str2double(omeGlobalMeta.get(strcat('Information|Image|SizeS')));	
    if seriesLength ~=2
        msg = ("not a two position file\n");
        error(msg)
    end
end
positionList = cell(seriesLength, 1);
for i = 1:seriesLength    
    positionList{i}.xPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|X #', int2str(i))));
    positionList{i}.yPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Y #', int2str(i))));
    positionList{i}.zPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Z #', int2str(i))));       
end
emAnt = [positionList{1}.xPosition, positionList{1}.yPosition];
emPos = [positionList{2}.xPosition, positionList{2}.yPosition];
emLength = pdist([emAnt; emPos]);
emAngle = atan((emAnt(2)-emPos(2))/(emAnt(1)-emPos(1)));
emDim.ant = emAnt;
emDim.pos = emPos;
emDim.len = emLength;
emDim.angle = emAngle;
% save([resultFolder,'/emDim.mat'],'emDim');
end