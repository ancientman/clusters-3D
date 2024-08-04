function [mcp_ch, bcd_ch, positionList, emDim] = readSeriesFile2(seriesMetaDataDS)
antPosFilePath = seriesMetaDataDS.expInfo.antPosFilePath;
seriesFilePath = seriesMetaDataDS.expInfo.rawDataFilePath;
CZIdata = HelperFunctions.loadCZI(seriesFilePath);
emDim = readAntPosFile(antPosFilePath); 
reader = bfGetReader(seriesFilePath);
omeGlobalMeta = reader.getGlobalMetadata();
positionList = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
positionStack = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
bcd_ch = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
mcp_ch = cell(seriesMetaDataDS.imagingInfo.seriesLength, 1);
bufferTop = seriesMetaDataDS.analysisInfo.bufferTop;% turns planes from Z = 1:bufferTop to zeros
bufferBot = seriesMetaDataDS.analysisInfo.bufferBot; % turns planes from Z = end-bufferBot:end to zeros
for i = 1:seriesMetaDataDS.imagingInfo.seriesLength
    if seriesMetaDataDS.imagingInfo.channelCount == 2
        imTemp = CZIdata{i, 1}(:, 1);
        z = 1;
        for j = 1:size(imTemp, 1)
            if seriesMetaDataDS.imagingInfo.bicoidChannel == 1
                if rem(j, 2) == 0
                   mcp_ch{i}(:,:,z) = imTemp{j, 1};
                    z = z+1;
                elseif rem(j, 2) == 1
                    bcd_ch{i}(:,:,z) = imTemp{j, 1};
                else
                        warning('work on the channel assignment!')
                end
            elseif seriesMetaDataDS.imagingInfo.bicoidChannel == 2
                if rem(j, 2) == 0
                   bcd_ch{i}(:,:,z) = imTemp{j, 1};
                    z = z+1;
                elseif rem(j, 2) == 1
                    mcp_ch{i}(:,:,z) = imTemp{j, 1};
                else
                        warning('work on the channel assignment!')
                end
            end                

        end

        if bufferTop
            mcp_ch{i}(:,:,1:bufferTop) = 0;
            bcd_ch{i}(:,:,1:bufferTop) = 0;
        end
        if bufferBot
            mcp_ch{i}(:,:,1:bufferBot) = 0;
            bcd_ch{i}(:,:,1:bufferBot) = 0;
        end 
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
    
    if ~isempty(nonzeros(positionList{i}.xPosition)) && ~isempty(nonzeros(positionList{i}.yPosition))
        positionList{i}.distFromAnt = pdist([emDim.ant; positionList{i}.xPosition,  positionList{i}.yPosition]);
        positionList{i}.angleWithAnt = atan((positionList{i}.yPosition-emDim.ant(2))/(positionList{i}.xPosition-emDim.ant(1)));
        positionList{i}.angleWithAxis =positionList{i}.angleWithAnt - emDim.angle;
        positionList{i}.distAlongAP = positionList{i}.distFromAnt*cos(positionList{i}.angleWithAxis);
        positionList{i}.fracDistEL = positionList{i}.distAlongAP/emDim.len;
    else
        positionList{i}.xPosition = [];
        positionList{i}.yPosition = [];
        positionList{i}.zPosition = [];
        positionList{i}.distFromAnt = [];
        positionList{i}.angleWithAnt = [];
        positionList{i}.angleWithAxis = [];
        positionList{i}.distAlongAP = [];
        positionList{i}.fracDistEL = [];
    end
end

% save('E:\Analysis\1c1P\3D\6x\DS\DS_20201223_bcd6x_nc14_em1_1_Airyscan Processing/positionListDS.mat','positionList');
% save('E:\Analysis\1c1P\3D\6x\DS\DS_20201223_bcd6x_nc14_em1_1_Airyscan Processing/embryoDimDS.mat','emDim');
end