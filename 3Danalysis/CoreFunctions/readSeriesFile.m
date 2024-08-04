function [stackList, positionList, emDim] = readSeriesFile(seriesMetaDataDS)
antPosFilePath = seriesMetaDataDS.expInfo.antPosFilePath;
seriesFilePath = seriesMetaDataDS.expInfo.rawDataFilePath;
CZIdata = HelperFunctions.loadCZI(seriesFilePath);
emDim = readAntPosFile(antPosFilePath); 
reader = bfGetReader(seriesFilePath);
omeGlobalMeta = reader.getGlobalMetadata();
positionList = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
positionStack = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
stackList = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
bufferTop = seriesMetaDataDS.analysisInfo.bufferTop;% turns planes from Z = 1:bufferTop to zeros
bufferBot = seriesMetaDataDS.analysisInfo.bufferBot; % turns planes from Z = end-bufferBot:end to zeros
for i = 1:seriesMetaDataDS.imagingInfo.seriesLength
    positionStack{i} = CZIdata{i, 1}(:, 1);
    stackList{i} = permute(reshape(cell2mat(positionStack{i})', [seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.imagingInfo.stackSizeZ]), [2, 1, 3]);
%     stackList{i} = permute(reshape(cell2mat(positionStack{i})', [seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.imagingInfo.stackSizeZ, seriesMetaDataDS.imagingInfo.channelCount]), [2, 1, 3]);

    if bufferTop
        stackList{i}(:,:,1:bufferTop) = 0;
    end
    if bufferBot
        stackList{i}(:,:,end-bufferBot:end) = 0;
    end   

    if seriesMetaDataDS.imagingInfo.seriesLength<10
        positionList{i}.xPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|X #', num2str(i, '%01d'))));
        positionList{i}.yPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Y #', num2str(i, '%01d'))));
        positionList{i}.zPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Z #', num2str(i, '%01d')))); 
    else
        positionList{i}.xPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|X #', num2str(i, '%02d'))));
        positionList{i}.yPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Y #', num2str(i, '%02d'))));
        positionList{i}.zPosition = str2double(omeGlobalMeta.get(strcat('Information|Image|S|Scene|Position|Z #', num2str(i, '%02d'))));     
    end
    
    positionList{i}.distFromAnt = pdist([emDim.ant; positionList{i}.xPosition,  positionList{i}.yPosition]);
    positionList{i}.angleWithAnt = atan((positionList{i}.yPosition-emDim.ant(2))/(positionList{i}.xPosition-emDim.ant(1)));
    positionList{i}.angleWithAxis =positionList{i}.angleWithAnt - emDim.angle;
    positionList{i}.distAlongAP = positionList{i}.distFromAnt*cos(positionList{i}.angleWithAxis);
    positionList{i}.fracDistEL = positionList{i}.distAlongAP/emDim.len;
end

% save('E:\Analysis\1c1P\3D\6x\DS\DS_20201223_bcd6x_nc14_em1_1_Airyscan Processing/positionListDS.mat','positionList');
% save('E:\Analysis\1c1P\3D\6x\DS\DS_20201223_bcd6x_nc14_em1_1_Airyscan Processing/embryoDimDS.mat','emDim');
end