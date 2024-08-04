function startAnalysisCZI(stackList, seriesMetaDataDS, positionList, emDim)
colorChannels = seriesMetaDataDS.imagingInfo.channelCount;
resultFolder = seriesMetaDataDS.expInfo.procImagesFolderName;
timePoints = seriesMetaDataDS.analysisInfo.totalTimePoints;
zSlices = seriesMetaDataDS.analysisInfo.nZslices; 
bicoidChannel = seriesMetaDataDS.imagingInfo.bicoidChannel;
% c1ImagesFolderName = seriesMetaDataDS.expInfo.c1ImagesFolderName;
c1Stack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
c1NucBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
c1NucLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
c1SpotBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
c1SpotLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
c1SpotProp = struct([]);
c1NucProp = struct([]);

removeBorderSpots = 1;

% if ~exist(append(c1ImagesFolderName, filesep, 'zStack'), 'file')
% 	mkdir(append(c1ImagesFolderName, filesep, 'zStack'));
% end

if colorChannels>1
%     c2ImagesFolderName = seriesMetaDataDS.expInfo.c2ImagesFolderName;
%     if ~exist(append(c2ImagesFolderName, filesep, 'zStack'), 'file')
%         mkdir(append(c2ImagesFolderName, filesep, 'zStack'));
%     end
    c2Stack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c2NucBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c2NucBinHullStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c2SpotBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c2SpotLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c1AdjNucLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);    
    c2NucLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c2SpotProp = struct([]);
    c2NucProp = struct([]);
    primeC2NucLabelStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
end

if colorChannels>2
    c3ImagesFolderName = seriesMetaDataDS.expInfo.c3ImagesFolderName;
    if ~exist(append(c3ImagesFolderName, filesep, 'zStack'), 'file')
        mkdir(append(c3ImagesFolderName, filesep, 'zStack'));
    end
    c3Stack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c3nucBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
    c3spotBinStack = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX, zSlices, timePoints]);
end

removeNucMask = zeros([seriesMetaDataDS.imagingInfo.stackSizeY, seriesMetaDataDS.imagingInfo.stackSizeX]);
c1c2SpotProp = struct([]);

saveTiffOptions.append = 1; %options for saveTiff function
saveTiffOptions.message = 0; %options for saveTiff function
for p = 1:length(stackList)
    fprintf('Processing series #%d\n', p);
    for t = 1:timePoints  
        for c=1:colorChannels
            if c==bicoidChannel
                channelFilled = 1; % TF channel
            else
                channelFilled = 0; % TF channel
            end
            cId = append('c',num2str(c));
            tId = append('t',sprintf('%04d',t));
            filePrefix = {cId, tId};
            filePrefix = strjoin(filePrefix,'_');      

            if c==1
                c1Stack(:,:,:,t) = stackList{p}(:,:,seriesMetaDataDS.imagingInfo.startZ:seriesMetaDataDS.imagingInfo.endZ);
                [c1NucBinStack(:,:,:,t), ~] = Nucleus.nucDetect3(c1Stack(:,:,:,t), channelFilled, seriesMetaDataDS);

                %% Remove undesired nuclei by drawing polygons around them
                if seriesMetaDataDS.analysisInfo.clearBorder == 0
                    [imNucRemoved, removeNucMask] = Preprocess.removeNuc3(c1NucBinStack(:,:,:,t), t, removeNucMask);
                else
                    
                    imNucRemoved = imclearborder(c1NucBinStack(:,:,:,t));
                end
                c1NucBinStack(:,:,:,t) = bwareaopen(imNucRemoved, ceil(seriesMetaDataDS.analysisInfo.minNucVol/3));
                [c1NucProp{p}{t}] = Nucleus.nucProp3(bwconncomp(c1NucBinStack(:,:,:,t)), c1Stack(:,:,:,t), seriesMetaDataDS);

                if colorChannels==1
                    c1NucLabelStack(:,:,:,t) = labelmatrix(bwconncomp(c1NucBinStack(:,:,:,t)));
                    [c1SpotLabelStack(:,:,:,t), c1SpotProp{p}{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
                        c1NucLabelStack(:,:,:,t), channelFilled, seriesMetaDataDS, removeBorderSpots);
                end

%             elseif c==2
%                 c2Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
%                 channelFilled = 0;
%                 [c2NucBinStack(:,:,:,t), c2NucProp{t}] = Nucleus.nucDetect3(c2Stack(:,:,:,t), channelFilled, seriesMetaDataDS); 
% 
%                 c2NucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(c2NucBinStack(:,:,:,t), channelFilled, seriesMetaDataDS);            
%                 [c1AdjNucLabelStack(:,:,:,t), c2NucLabelStack(:,:,:,t), c1NucProp{t}] = Nucleus.adjustNucLabeler3(c1NucBinStack(:,:,:,t), ...
%                     c2NucBinHullStack(:,:,:,t), primeC2NucLabelStack, seriesMetaDataDS, t);
%                 if t==1
%                     primeC2NucLabelStack = c2NucLabelStack(:,:,:,1);
%                 end
%                     channelFilled = 1;
%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
%                         c1AdjNucLabelStack(:,:,:,t), channelFilled, seriesMetaDataDS);
%                     channelFilled = 0;
%                     [c2SpotLabelStack(:,:,:,t), c2SpotProp{t}] = Spots.spotDetect3(c2Stack(:,:,:,t), ...
%                             c2NucLabelStack(:,:,:,t), channelFilled, seriesMetaDataDS);
%                     [c1c2SpotProp{t}] = Analysis.calcSpotDist(c1SpotProp{t}, c2SpotProp{t}, seriesMetaDataDS, t);
% 
%             elseif c==3
    %             c3Stack(:,:,:,t) = Helperfunctions.loadTiff(append(c3ImagesFolderName, filesep, filePrefix, '.tif'));
    %             c3NucBinStack(:,:,:,t) = Nucleus.nucDetect3(c3Stack(:,:,:,t), channelFilled, metaDataDS);
            end
        end     

%         HelperFunctions.saveTiff(single(c1SpotLabelStack(:,:,:,t)), ...
%             append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Spot_', tId, '.tif'), ...
%             saveTiffOptions);

        if colorChannels==1
    %         HelperFunctions.saveTiff(single(c1NucBinStack(:,:,:,t)), ...
    %             append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Nuc_', tId, '.tif'), ...
    %             saveTiffOptions);
        elseif colorChannels>1
    %         HelperFunctions.saveTiff(single(c1AdjNucLabelStack(:,:,:,t)), ...
    %             append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Nuc_', tId, '.tif'), ...
    %             saveTiffOptions);
    %         HelperFunctions.saveTiff(single(c2NucBinStack(:,:,:,t)), ...
    %             append(c2ImagesFolderName, filesep, 'zStack', filesep, 'c2Nuc_', tId, '.tif'), ...
    %             saveTiffOptions);      
    
%             HelperFunctions.saveTiff(single(c2SpotLabelStack(:,:,:,t)), ...
%                 append(c2ImagesFolderName, filesep, 'zStack', filesep, 'c2Spot_', tId, '.tif'), ...
%                 saveTiffOptions);
        end
    end
%     figure; kk = bwperim(c1NucBinStack(:,:,ceil(zSlices/2)) | bwperim(c1SpotLabelStack(:,:,ceil(zSlices/2))));  imshow(imfuse(rescale(c1Stack(:,:,ceil(zSlices/2))), kk));
sliceViewer(bwperim(c1NucBinStack) | bwperim(imbinarize(c1SpotLabelStack)));
end
save([resultFolder,'/c1SpotPropDS.mat'],'c1SpotProp', '-v7.3');
save([resultFolder,'/c1NucPropDS.mat'],'c1NucProp', '-v7.3');
save([resultFolder,'/positionListDS.mat'],'positionList', '-v7.3');
save([resultFolder,'/embryoDimDS.mat'],'emDim', '-v7.3');
% if colorChannels>1
%     save([resultFolder,'/c2NucPropDS.mat'],'c2NucProp');
%     save([resultFolder,'/c2SpotPropDS.mat'],'c2SpotProp');
%     save([resultFolder,'/c1c2SpotPropDS.mat'],'c1c2SpotProp');
% end
save([resultFolder,'/seriesMetaDataDS.mat'],'seriesMetaDataDS', '-v7.3');
% PlottingFunctions.C1SpotPropCZIplot1(resultFolder)
PlottingFunctions.C1SpotPropCZIplot1Lite2(resultFolder)
PlottingFunctions.C1NucPropCZIplot1(resultFolder)
end