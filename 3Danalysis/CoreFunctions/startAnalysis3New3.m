function startAnalysis3New3(metaDataDS)

%% %%%%%%%%%%%% Some global options %%%%%%%%%%%%%%
skipNuc = 0;%% skip nucleus segmentation: 1 for skip; 0 for no skip
skipChannel = 0; %% skip a particular channel: 1 for yes; 0 for no
trackNuc = 0; %% Track nucleus: 1 for yes; 0 for no
removeBorderSpots = 0; %% remove border spots
saveTiffOptions.append = 1; %options for saveTiff function
saveTiffOptions.message = 0; %options for saveTiff function

%% %%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
colorChannels = metaDataDS.imagingInfo.colorChannels;
resultFolder = metaDataDS.expInfo.procImagesFolderName;
timePoints = metaDataDS.analysisInfo.totalTimePoints;
zSlices = metaDataDS.analysisInfo.nZslices; 
bicoidChannel = metaDataDS.imagingInfo.bicoidChannel;

%% %%%%%% Datastructure Initialization for color channel #1 %%%%%%%%%%%%%
c1ImagesFolderName = metaDataDS.expInfo.c1ImagesFolderName;
if ~exist(append(c1ImagesFolderName, filesep, 'zStack'), 'file')
	mkdir(append(c1ImagesFolderName, filesep, 'zStack'));
end
c1Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1SpotLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
removeNucMask = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X]);
c1SpotProp = struct([]);
c1NucProp = struct([]);

%% %%%%%%%%%% Datastructure for color channel #2 %%%%%%%%%%%%%%%
if colorChannels>1
    c2ImagesFolderName = metaDataDS.expInfo.c2ImagesFolderName;
    if ~exist(append(c2ImagesFolderName, filesep, 'zStack'), 'file')
        mkdir(append(c2ImagesFolderName, filesep, 'zStack'));
    end
    c2Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2NucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c1NucBinHullStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2NucBinHullStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c1NucHullLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2SpotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2SpotLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c3SpotLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c1AdjNucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);    
    c2NucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c2SpotProp = struct([]);
    c2NucProp = struct([]);
    c1c2SpotProp = struct([]);
    c1c2IntensityProp = struct([]);
    radValBin = struct([]);
    corrCoeff = struct([]);
    peakLoc = struct([]);
    primeC2NucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices]);   
end

%% %%%%%%%%%% Datastructure for color channel #3 (placeholder) %%%%%%%%%%%%%%%
if colorChannels>2
    c3ImagesFolderName = metaDataDS.expInfo.c3ImagesFolderName;
    if ~exist(append(c3ImagesFolderName, filesep, 'zStack'), 'file')
        mkdir(append(c3ImagesFolderName, filesep, 'zStack'));
    end
    c3Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c3nucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
    c3spotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
end

%% %%%%%%%%%%%% Main time loop starts %%%%%%%%%%%%%% %%
for t = 1:timePoints  
    for c=1:colorChannels % total color channels
        cId = append('c',num2str(c));
        tId = append('t',sprintf('%04d',t));
        filePrefix = {cId, tId};
        filePrefix = strjoin(filePrefix,'_'); 
        
        if skipNuc == 1 % If nuclear segmentation is skipped
            if c==1 
                 c1Temp = HelperFunctions.loadTiff(append(c1ImagesFolderName, filesep, filePrefix, '.tif'));
                 
                if t==1
                     [imNucRemoved, maskRemoved] = Preprocess.removeNuc3(c1Temp, t, removeNucMask);
                     maskRemoved = repmat(~maskRemoved, 1, 1, size(c1Stack, 3));
                end
                
                [c1NucBinTemp, ~] = Nucleus.nucDetect3(c1Temp, channelFilled, metaDataDS);            
                c1NucBinTemp(~maskRemoved) = 0;
                c1NucBinStack(:,:,:,t) = bwareaopen(c1NucBinTemp, ceil(metaDataDS.analysisInfo.minNucVol));
                c1NucLabelStack(:,:,:,t) = labelmatrix(bwconncomp(c1NucBinStack(:,:,:,t))); 
                c1Temp(~maskRemoved) = 0;
                c1Stack(:,:,:,t) = c1Temp;     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotNoNuc3(c1Stack(:,:,:,t), channelFilled, maskRemoved, metaDataDS);
            elseif c==2 
                c2Temp = HelperFunctions.loadTiff(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
                c2Temp(~maskRemoved) = 0;
                c2Stack(:,:,:,t)  = c2Temp;
                [c2SpotLabelStack(:,:,:,t), c2SpotProp{t}] = Spots.spotNoNuc3(c2Stack(:,:,:,t), channelFilled, maskRemoved, metaDataDS);                     
            elseif c==3 
                c3Temp = HelperFunctions.loadTiff(append(c3ImagesFolderName, filesep, filePrefix, '.tif'));
                c3Temp(~maskRemoved) = 0;
                c3Stack(:,:,:,t)  = c3Temp;
                [c3SpotLabelStack(:,:,:,t), c2SpotProp{t}] = Spots.spotNoNuc3(c3Stack(:,:,:,t), channelFilled, maskRemoved, metaDataDS);     
            end
        %% %%%%%%%%% If nuclear segmentation carried out %%%%%%%%%%%%%%%%%    
        elseif skipNuc == 0  
            if c==1 && skipChannel ~= 1
                if bicoidChannel == 1
                    cId = append('c',num2str(1));
                    tId = append('t',sprintf('%04d',t));
                    filePrefix = {cId, tId};
                    filePrefix = strjoin(filePrefix,'_'); 
                    c1Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c1ImagesFolderName, filesep, filePrefix, '.tif'));
                    channelFilled = 1;
                elseif bicoidChannel == 2                    
                    cId = append('c',num2str(2));
                    tId = append('t',sprintf('%04d',t));
                    filePrefix = {cId, tId};
                    filePrefix = strjoin(filePrefix,'_'); 
                    c1Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
                    channelFilled = 1;
                end
                [c1NucBinStack(:,:,:,t), ~] = Nucleus.nucDetect3(c1Stack(:,:,:,t), channelFilled, metaDataDS);            
                
                if channelFilled == 1    
                    % Remove undesired nuclei by drawing polygons around them
                    [imNucRemoved, removeNucMask] = Preprocess.removeNuc3(c1NucBinStack(:,:,:,t), t, removeNucMask);
                    c1NucBinStack(:,:,:,t) = bwareaopen(imNucRemoved, ceil(metaDataDS.analysisInfo.minNucVol/3));
                end
                
                if colorChannels>1 && channelFilled == 1    
                    c1NucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(c1NucBinStack(:,:,:,t), channelFilled, metaDataDS);    
                elseif colorChannels == 1
                    continue
                else
%                     error('some error in channel assignment, bugger');
                end

            elseif c==2 && skipChannel ~= 2
                if bicoidChannel == 1                    
                    cId = append('c',num2str(2));
                    tId = append('t',sprintf('%04d',t));
                    filePrefix = {cId, tId};
                    filePrefix = strjoin(filePrefix,'_'); 
                    c2Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
                    channelFilled = 0;   
                elseif bicoidChannel == 2                    
                    cId = append('c',num2str(1));
                    tId = append('t',sprintf('%04d',t));
                    filePrefix = {cId, tId};
                    filePrefix = strjoin(filePrefix,'_'); 
                    c2Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c1ImagesFolderName, filesep, filePrefix, '.tif'));
                    channelFilled = 0;
                end
                [c2NucBinStack(:,:,:,t), ~] = Nucleus.nucDetect3(c2Stack(:,:,:,t), channelFilled, metaDataDS);                   
                c2NucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(c2NucBinStack(:,:,:,t), channelFilled, metaDataDS);            
                
            elseif c==3 && skipChannel ~= 3
            end
        end
    end   
    fprintf('time iteration = %d of %d\n',t,timePoints);
end
%   Main time loop ends
%_______________________________________________________________________________________________

for c=1:colorChannels % total color channels    
    if c==bicoidChannel % index for the filled channel
        channelFilled = 1; % Current color is Bicoid/histone channel
    else
        channelFilled = 0; % current color is TF channel
    end

    if c == 1 %% always the bicoid channel
        channelFilled = 1;
        [c1NucLabelStack, c1NucProp] = Nucleus.nucLabelAlign(c1NucBinStack, c1Stack, metaDataDS);
        [c1SpotLabelStack, c1SpotProp] = Spots.spotDetect3New(c1Stack, ...
            c1NucLabelStack, channelFilled, metaDataDS, removeBorderSpots);        

    elseif c == 2 %% always the mcp channel
        channelFilled = 0;
        [c2SpotBinStack, ~] = Spots.spotDetect3New(c2Stack, ...
            c1NucLabelStack, channelFilled, metaDataDS, removeBorderSpots);
        [c2SpotLabelStack, c2SpotProp] = Spots.c2SpotLabelAssign(c2SpotBinStack, c2Stack, c1NucLabelStack, c1NucProp, metaDataDS);
        for t = 1:timePoints 
            [c1c2SpotProp{t}] = Spots.c1c2SpotProp(c1NucProp{t}, c1SpotProp{t}, c2SpotProp{t}, metaDataDS, t);
        end        
    end
end

if colorChannels>1 && xor(bicoidChannel==1, bicoidChannel==2) 
    if skipNuc == 0
%         c1c2IntensityProp  = Postprocess.c1IntensityC2(c1Stack, c2SpotLabelStack, c2SpotProp, c1NucLabelStack, c1SpotLabelStack);
        c1c2IntensityProp  = Postprocess.c1c2Intensity3D(c1Stack, c2SpotLabelStack, c2SpotProp, c1NucLabelStack, c1SpotLabelStack, metaDataDS);
        c1nucCentIntensityProp  = Postprocess.c1nucCentIntensity3D(c1Stack, c1NucLabelStack, c1NucProp, metaDataDS);
%         spotDistTime = Postprocess.distTime(c1c2SpotProp);
    else 
%             c1IntensityFromC2 = c1IntensityC2(c1Stack, c2Stack, c2SpotProp, c1NucLabelStack);
    end
end  
%_______________________________________________________________________________________________

save([resultFolder,'/c1SpotPropDS.mat'],'c1SpotProp');
save([resultFolder,'/c1NucPropDS.mat'],'c1NucProp');
if colorChannels>1
    save([resultFolder,'/c2NucPropDS.mat'],'c2NucProp');
    save([resultFolder,'/c2SpotPropDS.mat'],'c2SpotProp');
    save([resultFolder,'/c1c2SpotPropDS.mat'],'c1c2SpotProp');
    save([resultFolder,'/c1c2IntensityPropDS.mat'],'c1c2IntensityProp')
    save([resultFolder,'/c1nucCentIntensityPropDS.mat'],'c1nucCentIntensityProp')
end
save([resultFolder,'/metaDataDS.mat'],'metaDataDS');
end



function spotTrack(spotProp)

distMax = 100;
timePoints = length(spotProp);
centLoc = cell(timePoints, 1);
nNuc = length(spotProp{1});
spotCentNuc = cell(1, nNuc);
ind = cell(1, nNuc);
meanCentLoc =  cell(1, nNuc);
relCentLoc =  cell(1, nNuc);
for i = 1:nNuc
    for t = 1:timePoints
        if ~isempty(spotProp{t}(i).center)
            centLoc{t} = spotProp{t}(i).center;
        end
    end
    spotCentNuc{i} = centLoc;
    spotCentNucArr = vertcat(spotCentNuc{i}{:});
    ind{i} = all(spotCentNucArr, 2);
    meanCentLoc{i} = mean(spotCentNucArr(ind{i},:));
    relCentLoc{i} = spotCentNucArr(ind{i},:) - repmat(meanCentLoc{i}, length(nonzeros(ind{i})),1);
    relCentLoc{i} = [0.042, 0.042, 0.2].*relCentLoc{i};
    p1 = plot3(relCentLoc{i}(:, 1), relCentLoc{i}(:, 2), relCentLoc{i}(:, 3));
end
% meanCentLoc = mean(centLoc(ind,:));
% relCentLoc = (centLoc(ind,:)) - repmat(mean(centLoc(ind,:)), length(nonzeros(ind)),1);
% relCentLoc = [0.042, 0.042, 0.2].*relCentLoc;
% p1 = plot3(relCentLoc(:, 1), relCentLoc(:, 2), relCentLoc(:, 3));
p1.LineWidth = 1.5;
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
xlabel('x {\mu}m', 'Rotation', 20);
ylabel('y {\mu}m', 'Rotation', -30);
zlabel('z {\mu}m');
grid on;
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 12;
hold on;
end




