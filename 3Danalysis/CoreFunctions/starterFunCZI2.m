function starterFunCZI2(analysisFolder)
[seriesMetaDataDS] =  readInputParams(analysisFolder);% constructs datastruct from expinfo.txt with user input
[seriesMetaDataDS] = readCZIParams(seriesMetaDataDS);
[seriesMetaDataDS] = InitializeOtherParams(seriesMetaDataDS);
save(append(analysisFolder,filesep,'seriesMetaDataDS.mat'),'seriesMetaDataDS','-V7.3')
[mcp_ch, bcd_ch, positionList, emDim] = readSeriesFile2(seriesMetaDataDS);
createCZIDSFolder(seriesMetaDataDS);
% startAnalysisCZI(stackList, seriesMetaDataDS, positionList, emDim);% main analysis function
startAnalysisCZI2(mcp_ch, bcd_ch, seriesMetaDataDS, positionList, emDim);% main analysis function
end