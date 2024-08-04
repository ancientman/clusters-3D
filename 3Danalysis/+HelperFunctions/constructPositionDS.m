function constructPositionDS()
seriesFolder = '\\tigress-cifs.princeton.edu\fileset-gregor\projects\Rahul\Airyscan\bcdXp2wt\20201011_bcdXp2wt\';
seriesFile = '20201011_bcdXp2wt_nc14_em2_1_Airyscan Processing.czi';
antPosFile = '20201011_bcdXp2wt_antpos_em1_1.czi';
[filepath,fileName,ext] = cellfun(@fileparts, append({[seriesFolder], [seriesFile]}), 'un', 0);
fileName = erase(fileName, '_Airyscan Processing');
fileName = append({'DS_'}, fileName);

resultFolder = strcat('E:\Data\airyscan\Test\2c_Pos',filesep,fileName);
emDim = readAntPosPosition(strcat(seriesFolder, antPosFile)); 
reader = bfGetReader(strcat(seriesFolder, seriesFile));
omeGlobalMeta = reader.getGlobalMetadata();

seriesLength = str2double(omeGlobalMeta.get(strcat('Information|Image|SizeS #1')));	
seriesLength(isnan(seriesLength)) = 1;
positionList = cell(seriesLength, 1);

for i = 1:seriesLength
    positionList{i}.xPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|X #', int2str(i))));
    positionList{i}.yPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Y #', int2str(i))));
    positionList{i}.zPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Z #', int2str(i))));     
    positionList{i}.distFromAnt = pdist([emDim.ant; positionList{i}.xPosition,  positionList{i}.yPosition]);
    positionList{i}.angleWithAnt = atan((positionList{i}.yPosition-emDim.ant(2))/(positionList{i}.xPosition-emDim.ant(1)));
    positionList{i}.angleWithAxis =positionList{i}.angleWithAnt - emDim.angle;
    positionList{i}.distAlongAP = positionList{i}.distFromAnt*cos(positionList{i}.angleWithAxis);
    positionList{i}.fracDistEL = positionList{i}.distAlongAP/emDim.len;
end
if ~exist(resultFolder{2}, 'file')
    mkdir(resultFolder{2});
end
save([resultFolder{2},'/positionDS.mat'],'positionList');
end

function [emDim] = readAntPosPosition(path)
reader = bfGetReader(path);
omeGlobalMeta = reader.getGlobalMetadata();
seriesLength = str2double(omeGlobalMeta.get(strcat('Information|Image|SizeS #1')));	
if seriesLength ~=2
    msg = ("not a two position file\n");
    error(msg)
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
end