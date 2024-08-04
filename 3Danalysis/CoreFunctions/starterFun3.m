function starterFun3(analysisFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the starting function to run
% Create experiment information data structure
% run analysis file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metaDataDS = makeMetaDataDS(analysisFolder);% constructs datastruct from expinfo.txt with user input
metaDataDS = makeMetaDataDSNew(analysisFolder);% constructs datastruct from expinfo.txt with user input
save(append(analysisFolder,filesep,'metaDataDS.mat'),'metaDataDS','-V7.3')
[metaDataDS] = createOutputFolder(metaDataDS);% creates a copy of the input image with split frames

% if metaDataDS.analysisInfo.analyze3D==0
%     startAnalysis(metaDataDS);% main analysis function
% else
%     startAnalysis3(metaDataDS);% main analysis function with 3d
%     startAnalysis3New(metaDataDS);% main analysis function with 3d time matching
%     startAnalysis3New2(metaDataDS);% main analysis function with 3d (alt)
    startAnalysis3New3(metaDataDS);% updated from startAnalysis3New for color channel swap
% end
end