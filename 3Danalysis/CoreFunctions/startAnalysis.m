function startAnalysis(metaDataDS)
colorChannels = metaDataDS.imagingInfo.colorChannels;
imUseType = metaDataDS.analysisInfo.imUseType;
rawImagesFolderName  = metaDataDS.expInfo.rawImagesFolderName;
resultFolderName = metaDataDS.expInfo.procImagesFolderName;
timePoints = metaDataDS.analysisInfo.totalTimePoints;
zSlices = metaDataDS.analysisInfo.nZslices; 
bicoidChannel = metaDataDS.imagingInfo.bicoidChannel;
totalFrames = colorChannels*zSlices*timePoints;
c1ImagesFolderName = metaDataDS.expInfo.c1ImagesFolderName;
c1Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucBinStackProject = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X]);
c1NucPropZ = struct([]);
c1SpotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1SpotPropNL = struct([]);

c1SpotProp3 = struct([]);
if ~exist(append(c1ImagesFolderName, filesep, 'zStack'), 'file')
	mkdir(append(c1ImagesFolderName, filesep, 'zStack'));
end

if colorChannels>1
    c2ImagesFolderName = metaDataDS.expInfo.c2ImagesFolderName;
    if ~exist(append(c2ImagesFolderName, filesep, 'zStack'), 'file')
        mkdir(append(c2ImagesFolderName, filesep, 'zStack'));
    end
    c2Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2NucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2NucBinStackProject = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X]);
    c2NucPropZ = struct([]);
    c2SpotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2SpotPropNL = struct([]);
    c12NucBinOverlayStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, 3, timePoints]);
    c12SpotBinOverlayStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, 3, timePoints]);
    c1AdjNucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    
    c2SpotProp3 = struct([]);
end

if colorChannels>2
    c3ImagesFolderName = metaDataDS.expInfo.c3ImagesFolderName;
    if ~exist(append(c3ImagesFolderName, filesep, 'zStack'), 'file')
        mkdir(append(c3ImagesFolderName, filesep, 'zStack'));
    end
    c3Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c3nucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c3spotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
end

saveTiffOptions.append = 1; %options for saveTiff function
saveTiffOptions.message = 0; %options for saveTiff function

for t = 1:timePoints    
    for z = 1:zSlices
        for c=1:colorChannels
            if c==bicoidChannel
                channelFilled = 1;
            else
                channelFilled = 0;
            end
            cId = append('c',num2str(c));
            tId = append('t',sprintf('%04d',t));
            zId = append('z',sprintf('%04d',z));
            filePrefix = {cId, tId, zId};
            filePrefix = strjoin(filePrefix,'_');            
            if c==1
                c1Stack(:,:,z,t) = imread(append(c1ImagesFolderName, filesep, filePrefix, '.tif'));                
                if imUseType==1 % use "raw"
                    imUse = c1Stack(:,:,z,t);
                elseif imUseType==2 % use "smooth"
                    imUse = Preprocess.smoothRawImage(c1Stack(:,:,z,t), metaDataDS);
                elseif imUseType==3 % use "sharp"
                    imUse = Preprocess.sharpRawImage(c1Stack(:,:,z,t), metaDataDS);
                else % use "raw"
                    imUse = c1Stack(:,:,z,t);
                end
                c1Stack(:,:,z,t) = imUse;
                c1NucBinStack(:,:,z,t) = Nucleus.nucDetect(c1Stack(:,:,z,t), channelFilled, metaDataDS);
                if colorChannels==1
                    c1NucLabel = labelmatrix(bwconncomp(c1NucBinStack(:,:,z,t)));
                    [c1SpotBinStack(:,:,z,t)] = Spots.spotDetect(c1Stack(:,:,z,t), ...
                        c1NucBinStack(:,:,z,t), c1NucLabel, channelFilled, metaDataDS);    
                end
            elseif c==2
                c2Stack(:,:,z,t) = imread(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
                if imUseType==1 % use "raw"
                    imUse = c2Stack(:,:,z,t);
                elseif imUseType==2 % use "smooth"
                    imUse = Preprocess.smoothRawImage(c2Stack(:,:,z,t), metaDataDS);
                elseif imUseType==3 % use "sharp"
                    imUse = Preprocess.sharpRawImage(c2Stack(:,:,z,t), metaDataDS);
                else % use "raw"
                    imUse = c2Stack(:,:,z,t);
                end
                c2Stack(:,:,z,t) = imUse;    
                c2NucBinStack(:,:,z,t) = Nucleus.nucDetect(c2Stack(:,:,z,t), channelFilled, metaDataDS); 
                c2NucLabelTemp = labelmatrix(bwconncomp(bwconvhull(c2NucBinStack(:,:,z,t), 'objects')));%%%%%%%%%%
                [c2SpotBinStack(:,:,z,t)] = Spots.spotDetect(c2Stack(:,:,z,t), c2NucBinStack(:,:,z,t), c2NucLabelTemp, channelFilled, metaDataDS);
                [c1AdjNucBinStack(:,:,z,t), c1NucLabel] = Nucleus.adjustFilledNucSeg(c1NucBinStack(:,:,z,t), c2NucBinStack(:,:,z,t));
                channelFilled = 1;
                [c1SpotBinStack(:,:,z,t)] = Spots.spotDetect(c1Stack(:,:,z,t), c1AdjNucBinStack(:,:,z,t), c1NucLabel, channelFilled, metaDataDS);
%             c12NucBinOverlayStack(:,:,z,:,t) = (imfuse(c1AdjNucBinStack(:,:,z,t), (c2NucBinStack(:,:,z,t))));
%             c12SpotBinOverlayStack(:,:,z,:,t) = (imfuse(c1SpotBinStack(:,:,z,t), (c2SpotBinStack(:,:,z,t))));
            elseif c==3
                c3Stack(:,:,z,t) = imread(append(c3ImagesFolderName, filesep, filePrefix, '.tif'));
                c3NucBinStack(:,:,z,t) = Nucleus.nucDetect(c3Stack(:,:,z,t), channelFilled, metaDataDS);
            end
        end        
        if colorChannels>1                    
        end 
        if t==1
            if z == zSlices
                lastZ = 1;
            else
                lastZ = 0;
            end
            if colorChannels==1
                [c1NucBinStackProject, c1NucStackProp] = Nucleus.nucBinStacker(c1NucBinStackProject, c1NucBinStack(:,:,z,t), lastZ, metaDataDS);
            elseif colorChannels>1
                [c1NucBinStackProject, c1NucStackProp] = Nucleus.nucBinStacker(c1NucBinStackProject, c1NucBinStack(:,:,z,t), lastZ, metaDataDS);
                [c2NucBinStackProject, c2NucStackProp] = Nucleus.nucBinStacker(c2NucBinStackProject, c2NucBinStack(:,:,z,t), lastZ, metaDataDS);
                if lastZ ==1
                    [combineNucBin, combineNucProp] = Nucleus.nucBinStackCombiner(c1NucBinStackProject, c2NucBinStackProject, lastZ, metaDataDS);
                end
            end
        end
    end
    [c1SpotBinStack(:,:,:,t), c1SpotProp3{t}] = Spots.spotPruner(c1SpotBinStack(:,:,:,t), metaDataDS);
    HelperFunctions.saveTiff(single(c1SpotBinStack(:,:,:,t)), ...
        append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Spot_', tId, '.tif'), ...
        saveTiffOptions);
    if colorChannels==1
        HelperFunctions.saveTiff(single(c1NucBinStack(:,:,:,t)), ...
            append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Nuc_', tId, '.tif'), ...
            saveTiffOptions);
    elseif colorChannels>1
        HelperFunctions.saveTiff(single(c1AdjNucBinStack(:,:,:,t)), ...
            append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Nuc_', tId, '.tif'), ...
            saveTiffOptions);
        HelperFunctions.saveTiff(single(c2NucBinStack(:,:,:,t)), ...
            append(c2ImagesFolderName, filesep, 'zStack', filesep, 'c2Nuc_', tId, '.tif'), ...
            saveTiffOptions);
        [c2SpotBinStack(:,:,:,t), c2SpotProp3{t}] = Spots.spotPruner(c2SpotBinStack(:,:,:,t), metaDataDS);        
        HelperFunctions.saveTiff(single(c2SpotBinStack(:,:,:,t)), ...
            append(c2ImagesFolderName, filesep, 'zStack', filesep, 'c2Spot_', tId, '.tif'), ...
            saveTiffOptions);
    end
end
end