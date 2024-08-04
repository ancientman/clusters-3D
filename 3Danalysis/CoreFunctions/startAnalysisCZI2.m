function startAnalysisCZI2(mcp_ch, bcd_ch, seriesMetaDataDS, positionList, emDim)
removeBorderSpots = 1;
resultFolder = seriesMetaDataDS.expInfo.procImagesFolderName;
timePoints = seriesMetaDataDS.analysisInfo.totalTimePoints;

bcdStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);
mcpStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);
bcdNucBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);
bcdNucBinHullStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);
mcpNucBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);
mcpNucBinHullStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, seriesMetaDataDS.analysisInfo.nZslices, timePoints]);

bcdSpotProp = struct([]);
bcdNucProp = struct([]);
mcpSpotProp = struct([]);
mcpNucProp = struct([]);
bcdMcpSpotProp = struct([]);
bcdMcpIntensityProp = struct([]);
bcdNucCentIntensityProp = struct([]);

removeNucMask = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX]);

for p = 1:length(mcp_ch) % scenes
    fprintf('Processing series #%d\n', p);
    for t = 1:timePoints  
        bcdStack(:,:,:,t) = bcd_ch{p}(:,:,seriesMetaDataDS.imagingInfo.startZ:seriesMetaDataDS.imagingInfo.endZ);
        channelFilled = 1; [bcdNucBinStack(:,:,:,t), ~] = Nucleus.nucDetect3(bcdStack(:,:,:,t), channelFilled, seriesMetaDataDS);
        %   Remove undesired nuclei by drawing polygons around them
        if seriesMetaDataDS.analysisInfo.clearBorder == 0
            [imNucRemoved, removeNucMask] = Preprocess.removeNuc3(bcdNucBinStack(:,:,:,t), t, removeNucMask);
        else            
            imNucRemoved = imclearborder(bcdNucBinStack(:,:,:,t));
        end
        bcdNucBinStack(:,:,:,t) = bwareaopen(imNucRemoved, ceil(seriesMetaDataDS.analysisInfo.minNucVol/3));
        channelFilled = 1; bcdNucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(bcdNucBinStack(:,:,:,t), channelFilled, seriesMetaDataDS); 
        
        mcpStack(:,:,:,t) = mcp_ch{p}(:,:,seriesMetaDataDS.imagingInfo.startZ:seriesMetaDataDS.imagingInfo.endZ);
    end

    [bcdNucLabelStack, bcdNucProp{p}] = Nucleus.nucLabelAlign(bcdNucBinStack, bcdStack, seriesMetaDataDS);
    [mcpNucLabelStack, mcpNucProp{p}] = Nucleus.nucLabelAlign(bcdNucBinStack, mcpStack, seriesMetaDataDS);

    channelFilled = 1; [bcdSpotLabelStack, bcdSpotProp{p}] = Spots.spotDetect3New(bcdStack, ...
        bcdNucLabelStack, channelFilled, seriesMetaDataDS, removeBorderSpots);
    
    channelFilled = 0; [mcpSpotBinStack, ~] = Spots.spotDetect3New(mcpStack, ...
        mcpNucLabelStack, channelFilled, seriesMetaDataDS, removeBorderSpots);
    [mcpSpotLabelStack, mcpSpotProp{p}] = Spots.c2SpotLabelAssign(mcpSpotBinStack, ...
        mcpStack, bcdNucLabelStack, bcdNucProp{p}, seriesMetaDataDS);

    for t = 1:timePoints 
        [bcdMcpSpotProp{p}{t}] = Analysis.c1c2SpotProp(bcdNucProp{p}{t}, bcdSpotProp{p}{t}, mcpSpotProp{p}{t}, seriesMetaDataDS, t);
    end   

    bcdMcpIntensityProp{p}  = Postprocess.c1c2Intensity3D(bcdStack, mcpSpotLabelStack, mcpSpotProp{p}, bcdNucLabelStack, bcdSpotLabelStack, seriesMetaDataDS);
    bcdNucCentIntensityProp{p}  = Postprocess.c1nucCentIntensity3D(bcdStack, bcdNucLabelStack, bcdNucProp{p}, seriesMetaDataDS);
end

save([resultFolder,'/c1SpotPropDS.mat'],'bcdSpotProp', '-v7.3');
save([resultFolder,'/c1NucPropDS.mat'],'bcdNucProp', '-v7.3');

save([resultFolder,'/c2SpotPropDS.mat'],'mcpSpotProp', '-v7.3');
save([resultFolder,'/c2NucPropDS.mat'],'mcpNucProp', '-v7.3');

save([resultFolder,'/c1c2SpotPropDS.mat'],'bcdMcpSpotProp');
save([resultFolder,'/c1c2IntensityPropDS.mat'],'bcdMcpIntensityProp')
save([resultFolder,'/c1nucCentIntensityPropDS.mat'],'bcdNucCentIntensityProp')

save([resultFolder,'/positionListDS.mat'],'positionList', '-v7.3');
save([resultFolder,'/embryoDimDS.mat'],'emDim', '-v7.3');

save([resultFolder,'/seriesMetaDataDS.mat'],'seriesMetaDataDS', '-v7.3');

Analysis.spotPropCZI2Cplot1Lite(resultFolder)
% Analysis.nucPropCZI2Cplot1(resultFolder)
end