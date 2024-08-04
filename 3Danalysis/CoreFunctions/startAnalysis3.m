function startAnalysis3(metaDataDS)
colorChannels = metaDataDS.imagingInfo.colorChannels;
resultFolder = metaDataDS.expInfo.procImagesFolderName;
timePoints = metaDataDS.analysisInfo.totalTimePoints;
zSlices = metaDataDS.analysisInfo.nZslices; 
bicoidChannel = metaDataDS.imagingInfo.bicoidChannel;
c1ImagesFolderName = metaDataDS.expInfo.c1ImagesFolderName;
c1Stack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1ThresStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1NucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1SpotBinStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1SpotLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices, timePoints]);
c1SpotProp = struct([]);
c1NucProp = struct([]);
c1AddNucProp = struct([]);
c1AddSpotProp = struct([]);

skipNuc = 0;%% skip nucleus segmentation: 1 for skip; 0 for no skip
skipChannel = 0; %% skip a particular channel: 1 for yes; 0 for no
trackNuc = 0; %% Track nucleus: 1 for yes; 0 for no
removeBorderSpots = 0; %% remove border spots

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
    primeC2NucLabelStack = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X, zSlices]);   
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

removeNucMask = zeros([metaDataDS.imagingInfo.Y, metaDataDS.imagingInfo.X]);
c1c2SpotProp = struct([]);

saveTiffOptions.append = 1; %options for saveTiff function
saveTiffOptions.message = 0; %options for saveTiff function

for t = 1:timePoints  
    for c=1:colorChannels % total color channels
        if c==bicoidChannel % index for the filled channel
            channelFilled = 1; % Current color is Bicoid/histone channel
        else
            channelFilled = 0; % current color is TF channel
        end
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
            
        elseif skipNuc == 0  % If nuclear segmentation carried out
            if c==1 && skipChannel ~= 1
                c1Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c1ImagesFolderName, filesep, filePrefix, '.tif'));
                [c1NucBinStack(:,:,:,t), ~] = Nucleus.nucDetect3(c1Stack(:,:,:,t), channelFilled, metaDataDS);            
                % Remove undesired nuclei by drawing polygons around them
                [imNucRemoved, removeNucMask] = Preprocess.removeNuc3(c1NucBinStack(:,:,:,t), t, removeNucMask);
                c1NucBinStack(:,:,:,t) = bwareaopen(imNucRemoved, ceil(metaDataDS.analysisInfo.minNucVol/3));
                
                if colorChannels==1 && channelFilled == 1                    
                    c1NucLabelStack(:,:,:,t) = labelmatrix(bwconncomp(c1NucBinStack(:,:,:,t))); 
                    [c1NucProp{t}] = Nucleus.nucProp3(bwconncomp(c1NucBinStack(:,:,:,t)), c1Stack(:,:,:,t), metaDataDS); 
%                     [c1ThresStack(:,:,:,t), ~] = Spots.adapThres3( c1Stack(:,:,:,t), c1NucLabelStack(:,:,:,t));       
%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1ThresStack(:,:,:,t), ...
%                         c1NucLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);                       
                        
                    [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
                        c1NucLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);
                    
                elseif colorChannels>1 && channelFilled == 1             
                    [c1NucProp{t}] = Nucleus.nucProp3(bwconncomp(c1NucBinStack(:,:,:,t)), c1Stack(:,:,:,t), metaDataDS);           
                    c1NucLabelStack(:,:,:,t) = labelmatrix(bwconncomp(c1NucBinStack(:,:,:,t)));
                    c1NucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(c1NucBinStack(:,:,:,t), channelFilled, metaDataDS);    
                    c1NucHullLabelStack(:,:,:,t) = labelmatrix(bwconncomp(c1NucBinHullStack(:,:,:,t)));
%                     [c1ThresStack(:,:,:,t), ~] = Spots.adapThres3( c1Stack(:,:,:,t), c1NucLabelStack(:,:,:,t));%         


%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1ThresStack(:,:,:,t), ...
%                         c1NucLabelStack(:,:,:,t), channelFilled, metaDataDS);                    
%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
%                         c1NucLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);
                else
                    error('some error in channel assignment, bugger');
                end

            elseif c==2 && skipChannel ~= 2
                c2Stack(:,:,:,t) = HelperFunctions.loadTiff(append(c2ImagesFolderName, filesep, filePrefix, '.tif'));
                channelFilled = 0;
                [c2NucBinStack(:,:,:,t), c2NucProp{t}] = Nucleus.nucDetect3(c2Stack(:,:,:,t), channelFilled, metaDataDS);             
                c2NucBinHullStack(:,:,:,t) = Nucleus.nucHollowHull3(c1NucBinStack(:,:,:,t), channelFilled, metaDataDS);            
                [c1AdjNucLabelStack(:,:,:,t), c2NucLabelStack(:,:,:,t), c1NucProp{t}] = Nucleus.adjustNucLabeler3(c1NucBinStack(:,:,:,t), ...
                    c2NucBinHullStack(:,:,:,t), c1Stack(:,:,:,t), primeC2NucLabelStack, metaDataDS, t);
                if t==1
                    primeC2NucLabelStack = c2NucLabelStack(:,:,:,1);
                else
%                      [c1AdjNucLabelStack(:,:,:,t), c2NucLabelStack(:,:,:,t), c1NucProp{t}, c2NucProp{t}] = Nucleus.timeAdjustNucLabeler3(c1NucBinStack(:,:,:,t), ...
%                     c2NucBinHullStack(:,:,:,t), primeC2NucLabelStack, metaDataDS, t);
                    
                end
    %             c1AdjNucLabelStack(:,:,:,t) = c1NucBinStack(:,:,:,t); %%%%%%%%%%%%%%%%%%%%%
                channelFilled = 1;
%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
%                         c1AdjNucLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);

                    [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1Stack(:,:,:,t), ...
                        c1NucHullLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);
%                     [c1SpotLabelStack(:,:,:,t), c1SpotProp{t}] = Spots.spotDetect3(c1ThresStack(:,:,:,t), ...
%                         c1NucHullLabelStack(:,:,:,t), channelFilled, metaDataDS, removeBorderSpots);
                    channelFilled = 0;
%                     [c2SpotLabelStack(:,:,:,t), c2SpotProp{t}] = Spots.spotDetect3(c2Stack(:,:,:,t), ...
%                             c2NucLabelStack(:,:,:,t), channelFilled, metaDataDS,  removeBorderSpots);                    
                    [c2SpotLabelStack(:,:,:,t), c2SpotProp{t}] = Spots.spotDetect3(c2Stack(:,:,:,t), ...
                            c1NucHullLabelStack(:,:,:,t), channelFilled, metaDataDS,  removeBorderSpots);
%                     [c1c2SpotProp{t}] = Analysis.c1c2SpotProp(c1NucProp{t}, c1SpotProp{t}, c2SpotProp{t}, metaDataDS, t);

            elseif c==3
    %             c3Stack(:,:,:,t) = Helperfunctions.loadTiff(append(c3ImagesFolderName, filesep, filePrefix, '.tif'));
    %             c3NucBinStack(:,:,:,t) = Nucleus.nucDetect3(c3Stack(:,:,:,t), channelFilled, metaDataDS);
            end
        end
    end
    
    HelperFunctions.saveTiff(single(c1SpotLabelStack(:,:,:,t)), ...
        append(c1ImagesFolderName, filesep, 'zStack', filesep, 'c1Spot_', tId, '.tif'), ...
        saveTiffOptions);
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
        HelperFunctions.saveTiff(single(c2SpotLabelStack(:,:,:,t)), ...
            append(c2ImagesFolderName, filesep, 'zStack', filesep, 'c2Spot_', tId, '.tif'), ...
            saveTiffOptions);
    end
    fprintf('time iteration = %d of %d\n',t,timePoints);
end
%   Main time loop ends
%_______________________________________________________________________________________________

%% Visualization
%_______________________________________________________________________________________________
% imshow(zeros(size(c1Stack(:,:,:,1))));
% figure('Color', 'w')
% imshow(zeros(size(c1Stack(:,:,2))));
% imshow(bwperim(c1NucBinHullStack(:,:,2)));
% hold on;
% % imshow(rescale(imbinarize(c1SpotLabelStack(:,:,2)).*c1Stack(:,:,2)),[]);
% % imshow(rescale(imbinarize(c1SpotLabelStack(:,:,2)).*c1Stack(:,:,2))+ 0.25.*bwperim(c1NucBinHullStack(:,:,2)),[]);
% hold on;
% distC1C2 = [];
% volLim = 20;
% for t = 1:length(c1SpotProp)
%     hold on;
%     for i = 2:length(c1SpotProp{t})
%         hold on;
%         if ~isempty(c1SpotProp{t}(i).voxValCenter)
%             center1X = c1SpotProp{t}(i).voxValCenter(:, 1) - (c1NucProp{t}(i).voxValCenter(1) - c1NucProp{t}(i).voxValCenter(1));
%             center1Y = c1SpotProp{t}(i).voxValCenter(:, 2) - (c1NucProp{t}(i).voxValCenter(2) - c1NucProp{1}(i).voxValCenter(2));
%             center1X(c1SpotProp{t}(i).vol<volLim) = [];
%             center1Y(c1SpotProp{t}(i).vol<volLim) = [];
%             plot(center1X, center1Y, 'g*');
%             hold on;
%         end
%         if colorChannels>1 && ~isempty(c2SpotProp{t}(i).voxValCenter)
%             center2X = c2SpotProp{t}(i).voxValCenter(:, 1) - (c1NucProp{t}(i).voxValCenter(1) - c1NucProp{t}(i).voxValCenter(1));
%             center2Y = c2SpotProp{t}(i).voxValCenter(:, 2) - (c1NucProp{t}(i).voxValCenter(2) - c1NucProp{1}(i).voxValCenter(2));
%             plot(center2X, center2Y, 'mv');
%             hold on;
%             spotC1 = [center1X, center1Y];
%             spotC2 = [center2X, center2Y];
%             distC1C2 = vertcat(distC1C2, min(pdist2(spotC2, spotC1))); %%%%%% check
%         end
%     end
% end
% hold on;
% B = bwboundaries(c1NucBinStack(:,:,:,1),'noholes');
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'y--', 'LineWidth', 1.5)
% end

%% Limits for single nucleus
%%---------------------------------------------------
% xMin = min(boundary(:,2));
% xMax = max(boundary(:,2));
% yMin = min(boundary(:,1));
% yMax = max(boundary(:,1));
% 
% xDiff = xMax - xMin;
% yDiff = yMax - yMin;
% 
% len = ceil(max(xDiff, yDiff)/10)*10;
% 
% xMid = round((xMin+xMax)/2);
% xMin = max(xMid - ceil(len/2), 1);
% xMax = min(xMid + ceil(len/2), size(c1Stack(:,:,:,1), 2));
% 
% yMid = round((yMin+yMax)/2);
% yMin = max(yMid - ceil(len/2), 1);
% yMax = min(yMid + ceil(len/2), size(c1Stack(:,:,:,1), 1));
% 
% xlim([xMin, xMax])
% ylim([yMin, yMax])
% %%---------------------------------------------------
% 
% x0=10;
% y0=10;
% width=200;
% height=200;
% set(gcf,'position',[x0,y0,width,height]);
% set(gca, 'color', 'k');
% set(gca,'visible','off');
% hold off;
% 
% figure('Color', 'w')
% p1 = plot(1.5.*(1:length(distC1C2)), distC1C2.*(0.042), 'k');
% p1.LineWidth = 1.0;
% set(gca, 'YScale', 'log')
% grid on
% ylabel('Distance (\mu m)');
% xlabel('Time (s)');

%_______________________________________________________________________________________________

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

% [c1NucLabelStack, c1NucProp] = Nucleus.nucLabelAlign(c1NucBinStack, c1Stack, metaDataDS);
% [c1SpotLabelStack, c1SpotProp] = Spots.spotDetect3New(c1Stack, ...
%     c1NucLabelStack, channelFilled, metaDataDS, removeBorderSpots);

%   Watch
%_______________________________________________________________________________________________
% figure; sliceViewer(imbinarize(squeeze((c2SpotLabelStack))));
% figure; sliceViewer(imbinarize(squeeze((c1SpotLabelStack))));

%   Create time add image
%_______________________________________________________________________________________________

% for c=1:colorChannels % total color channels
%         if c==bicoidChannel % index for the filled channel
%             channelFilled = 1; % Current color is Bicoid/histone channel
%             
%             [c1AddNucProp, c1AddSpotProp] = Postprocess.timeAddIm(c1Stack, c1NucLabelStack, c1SpotLabelStack, metaDataDS, channelFilled, removeBorderSpots);
%         else
%             channelFilled = 0; % current color is TF channel
%         end
% end
% 
% if colorChannels ==1
% %     c1SpotLabelStack = imclearborder(c1SpotLabelStack);
%     [c1AddSpot, c1AddThres, c1AddNuc, c1AddINucLabel, c1AddNucProp] = Postprocess.timeAddIm(c1Stack, c1NucLabelStack, c1SpotLabelStack, metaDataDS, channelFilled, removeBorderSpots);
%     [c1AddSpotProp] = Postprocess.globalHotSpotFinder(c1AddThres, channelFilled, c1AddNucProp, metaDataDS);
%     if bicoidChannel ~= 1
% %         spotTrack(c1SpotProp);
%     end
% elseif colorChannels>1
% %     c1SpotLabelStack = imclearborder(c1SpotLabelStack);
%     [c1AddSpot, c1AddThres, c1AddNuc, c1AddINucLabel] = Postprocess.timeAddIm(c1Stack, c1AdjNucLabelStack, c1SpotLabelStack);
%     if bicoidChannel == 1
% %         spotTrackFilled(c1SpotProp);
%     end
%     if bicoidChannel ~= 1
% %         spotTrack(c1SpotProp); 
%     end
%     [c2AddSpot, c2AddThres, c2AddNuc, c2AddINucLabel] = Postprocess.timeAddIm(c2Stack, c1AdjNucLabelStack, c2SpotLabelStack);
%     if bicoidChannel ~= 2
% %         spotTrack(c2SpotProp);
%     end
%     
%     if xor(bicoidChannel==1, bicoidChannel==2) 
%         if skipNuc == 0
%             Postprocess.c1IntensityC2(c1Stack, c2SpotLabelStack, c2SpotProp, c1NucLabelStack, c1SpotLabelStack);
%             spotDistTime = distTime(c1c2SpotProp);
%         else 
% %             c1IntensityFromC2 = c1IntensityC2(c1Stack, c2Stack, c2SpotProp, c1NucLabelStack);
%         end
%     end
%     
% elseif colorChannels>2
%     [c1AddSpot, c1AddThres, c1AddNuc, c1AddINucLabel] = Postprocess.timeAddIm(c1Stack, c1AdjNucLabelStack, c1SpotLabelStack);
%     if bicoidChannel ~= 1
% %         spotTrack(c1SpotProp);
%     end
%     [c2AddSpot, c2AddThres, c2AddNuc, c2AddINucLabel] = Postprocess.timeAddIm(c2Stack, c1AdjNucLabelStack, c2SpotLabelStack);
%     if bicoidChannel ~= 2
% %         spotTrack(c2SpotProp);
%     end
%     [c3AddSpot, c3AddThres, c3AddNuc, c3AddINucLabel] = Postprocess.timeAddIm(c3Stack, c1AdjNucLabelStack, c3SpotLabelStack);
%     if bicoidChannel ~= 3
% %         spotTrack(c3SpotProp);
%     end
% end
%_______________________________________________________________________________________________

%   Visualize
% figure; imshow(imoverlay(c1TimeStack, bwperim(imbinarize(c2TimeStack))));
%_______________________________________________________________________________________________

save([resultFolder,'/c1SpotPropDS.mat'],'c1SpotProp');
save([resultFolder,'/c1NucPropDS.mat'],'c1NucProp');
save([resultFolder,'/c1NucCentIntensityPropDS.mat'],'c1nucCentIntensityProp');
if colorChannels>1
    save([resultFolder,'/c2NucPropDS.mat'],'c2NucProp');
    save([resultFolder,'/c2SpotPropDS.mat'],'c2SpotProp');
    save([resultFolder,'/c1c2SpotPropDS.mat'],'c1c2SpotProp');
    save([resultFolder,'/c1c2IntensityPropDS.mat'],'c1c2IntensityProp');
end
save([resultFolder,'/metaDataDS.mat'],'metaDataDS');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% needs a alot of work %%%%%%%%%%%%%
function spotTrackFilled(spotProp)
distMax = 100;
volCutOff = 0.005; % in UM
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


function spotDist = distTime(c1c2SpotProp)
timePoints = length(c1c2SpotProp);
nNuc = length(c1c2SpotProp{1}.C1C2SpotVolUM);
minDist = cell(1, nNuc);
for i = 1:nNuc
    minDist{i} = zeros(timePoints, 1);
    for t = 1:timePoints    
%     volTemp = c1c2SpotProp{t}.C1C2SpotVolUM{i};
%     [volTemp, volIdx] = sort(volTemp);
    distTemp = c1c2SpotProp{t}.C1C2SpotDistUM{i};
%     distTemp = distTemp(volIdx);
    if ~isempty(distTemp)
        minDist{i}(t) = min(distTemp);
    end
    end
figure('Color', 'w');
p1 = plot (minDist{i});
p1.Color = [0.1 0.1 0.1];
p1.LineWidth = 1.5;
p1.LineStyle = '-';
legend('off');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
grid ('on');
xlabel('Time (s)');
ylabel('Dist. from mRNA center {\mu}m');
title('Distance of closest bicoid from mRNA');
hold on;
x = 0:1:timePoints;
c=mean(minDist{i});
const = @(x)(c).*x.^(0);
p2 = plot(x, const(x));
p2.Color = [0.7 0.2 0.2];
p2.LineWidth = 1.5;
p2.LineStyle = '-.';
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end
spotDist = minDist;
end

function c1IntensityC2(c1Stack, c2SpotLabelStack, c2SpotProp, c1NucLabelStack, c1SpotLabelStack)
bgSub = 1;
maxcorrPlane = 1;
distMax = 25;
im = c1Stack;
timePoints = size(im, 4);
nNuc = length(nonzeros(HelperFunctions.count_unique(c1NucLabelStack(:,:,:,:))));

bgMean = double(zeros(timePoints, nNuc));
bgStd = double(zeros(timePoints, nNuc));

timePoints = size(c1Stack, 4);
nNuc = max(c1NucLabelStack, [], 'all');
if nNuc == 0
    nNuc = 1;
end
spotCentNuc = cell(1, nNuc);
peakLoc =  cell(1, nNuc);
corrCoeff = cell(1, nNuc);
mockIm = zeros(size(c1Stack(:,:,1,1)));
for i = 1:nNuc
    spotCentNuc{i} = zeros(timePoints, 3);
end

for t=1:timePoints
    for i = 1:size(c2SpotProp{t}, 2)
        if ~isempty(c2SpotProp{t}(i).center)
            spotCentNuc{i}(t,:) = round(c2SpotProp{t}(i).center);
        end
    end
end

radValBin= cell(1, nNuc);
for i = 1:nNuc
    radValBin{i} = zeros(distMax, timePoints);
    peakLoc{i} = zeros(timePoints, 1);
    corrCoeff{i} = zeros(timePoints, 1);
    iter = 0;
    flag = 0;
    for t = 1:timePoints
        if spotCentNuc{i}(t,3) ~=0
            zTemp = spotCentNuc{i}(t,3);
            labelMatTemp = c1NucLabelStack(:,:,zTemp,t);
            labelMatTemp(labelMatTemp~=i) = 0;
            c1Temp = c1Stack(:,:,zTemp,t).*labelMatTemp;           
            % For subtraction
            %--------------------------------------------------------------
            if bgSub == 1
                c1Temp2 = c1Temp;          
                c1Temp2(c1SpotLabelStack(:,:,zTemp,t)==i) = 0;       
                bgMean(t, i) = mean(nonzeros(c1Temp2),'all', 'omitnan');
                bgStd(t, i) = std(nonzeros(c1Temp2),0,'all', 'omitnan');
                imTempNew = c1Temp - (bgMean(t, i) + 2*bgStd(t, i));
                imTempNew(imTempNew<0) = 0;
                imTempNew(isnan(imTempNew)) = 0;     
                c1Temp = rescale(imTempNew);
            end
            %--------------------------------------------------------------
            spotCentTemp = spotCentNuc{i}(t,1:2);
            
            xEdge = (spotCentTemp(1) - distMax) : (spotCentTemp(1) + distMax);
            yEdge = (spotCentTemp(2) - distMax) : (spotCentTemp(2) + distMax);
            [xNuc,yNuc] = meshgrid(xEdge, yEdge);
            if  ~isempty(setdiff(xNuc,1:size(mockIm, 1))) || ~isempty(setdiff(yNuc,1:size(mockIm, 2)))
                continue
            else            
                mockIm(sub2ind(size(mockIm), yNuc, xNuc)) = 1;
                c1Temp(~mockIm) = 0;
                c1Crop = imcrop(c1Temp, [(spotCentTemp(1) - distMax), (spotCentTemp(2) - distMax), 2*distMax, 2*distMax]);
                c1Crop = medfilt2(c1Crop);
                if bgSub == 1
                    c1Crop = rescale(c1Crop); % For subtraction             
                else
                    c1Crop = c1Crop/c1Crop(distMax+1, distMax+1); % For no subtraction
                end               
%                 corrMat = corr2D(c1Crop);
                
                [valRadBinAvg, valRadBinSem] = cart2PolAvg(distMax, c1Crop);
%                 radValBin{i} =  horzcat(radValBin{i},valRadBinAvg);
                radValBin{i}(:,t) =  valRadBinAvg;
                iter = iter + 1;
            end
        end
        if any(radValBin{i}(:,t))
            if flag == 0
                c1CropPrime = c1Crop;
            end
            corrCoeff{i}(t) = corr2(c1CropPrime, c1Crop);
%             [pks,loc] = findpeaks(radValBin{i}(:,t));
            [maxVal,loc] = max(radValBin{i}(:,t));
            peak1Loc = min(loc);
            if peak1Loc == 0
                peak1Loc = 1;
            end
            peakLoc{i}(t) = peak1Loc;
            flag = flag+1;
        end
    end

     %  Visualizations%%%%%%%%%%%%%%
%     linePlot1(radValBin{i});
%     linePlot2(0.042*peakLoc{i});
%     linePlot2(corrCoeff{i});
end
end


function [valRadBinAvg, valRadBinSem] = cart2PolAvg(rMax, imCrop)
xVals = ones(1, 2*rMax+1)'*(-rMax:rMax);
yVals = (-rMax:rMax)'*ones(1, 2*rMax+1);
[theta, rho, valTemp] = cart2pol(xVals,yVals,  imCrop);  % convert x, y to polar coordinates
aR = reshape(rho,1, (2*rMax+1)^2);
aValTemp = reshape(valTemp,1, (2*rMax+1)^2);
[rrTemp,ind] = sort(aR);
corrRadVal = aValTemp(ind);
% rad= 0:floor(max(rrTemp));
rad= 0:rMax;
[n, bin] = histcounts(rrTemp, rad);
valRadBinAvg = zeros(length(bin)-1, 1);
valRadBinSem = zeros(length(bin)-1, 1);
for j = 1:length(bin)-1
    valRadBinAvg(j) = sum(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sum(rrTemp>bin(j) & rrTemp<=bin(j+1));
    valRadBinSem(j) = std(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sqrt(sum(rrTemp>bin(j) & rrTemp<=bin(j+1)));                        
end
end


function linePlot1(radValBin)
figure('color', 'w');
for i=1:size(radValBin, 2)
    p1= plot(0.042*(1:size(radValBin, 1)), radValBin(:,i));
    p1.LineWidth = 1.5;
    p1.Color = [0.7 0.7 0.7];
    p1.LineStyle = '-';    
    hold on;
end

nonZeroRows = any(radValBin, 1);
radValBinAvg = mean(radValBin(:,nonZeroRows), 2);
radValBinErr = std(radValBin(:,nonZeroRows), 0, 2);

p2= plot(0.042*(1:size(radValBin, 1)), radValBinAvg);
p2.LineWidth = 2;
p2.Color = [0.8 0.3 0.3];
p2.LineStyle = '-';
ylabel('Relative bcd intensity (a. u.)');
xlabel('Apparent radius ({\mu}m)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;

plotErr(radValBinAvg, radValBinErr);
end

function plotErr(val, err)
figure('color', 'w');
pixScale = 0.042; 
patchTop = val+err;
patchTop = reshape(patchTop,1,[]);
patchBot = val-err;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch=[pixScale*(1:length(val)),fliplr(pixScale*(1:length(val)))];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = [0.3, 0.7, 0.6];
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.6;
hold on;
pl = plot(pixScale*(1:length(val)), val);
pl.LineWidth = 1.5;
pl.Color = [0.4 0.4 0.4];
pl.LineStyle = '--';
ylabel('Relative bcd intensity (a. u.)');
xlabel('Apparent radius (r)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function linePlot2(radValBin)
figure('color', 'w');

p1= plot(1.1*(1:size(radValBin, 1)), radValBin);
p1.LineWidth = 1.5;
p1.Color = [0.3 0.3 0.3];
p1.LineStyle = '-';   
ylabel('Correlation coefficient');
% ylabel('Brightest Bicoid from hb ({\mu}m');
xlabel('Time (s)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function corrMat = corr2D(cropIm)
cropZeroPad = 100;
centCropPad = zeros(size(cropIm, 1) + 2*cropZeroPad,  size(cropIm, 2) + 2*cropZeroPad);
centCropPad(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = cropIm;
mask = zeros(size(centCropPad));
mask(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = 1;
nParticles = sum(cropIm, 'all');  % number of particles within mask
nPixMask = sum(mask, 'all');      % area of mask
NP = real(fftshift(ifft2(abs(fft2(mask)).^2))); % Normalization for correct boundary conditions
G1 = nPixMask^2/nParticles^2*real(fftshift(ifft2(abs(fft2(centCropPad)).^2)))./NP; % 2D G(r) with proper normalization
corrMat = imcrop(G1, [floor(size(G1, 2)/2+1)-cropZeroPad, floor(size(G1, 1)/2+1)-cropZeroPad, 2*cropZeroPad, 2*cropZeroPad]);  %only return valid part of G
[valAvg, valSem] = cart2PolCorr(cropZeroPad, corrMat);
corrCentRadAvg = valAvg;
corrCentRadSem = valSem;
plotErr(corrCentRadAvg, corrCentRadSem);
end

function [valRadBinAvg, valRadBinSem] = cart2PolCorr(rMax, corrCut)
xVals = ones(1, 2*rMax+1)'*(-rMax:rMax);
yVals = (-rMax:rMax)'*ones(1, 2*rMax+1);
[theta, rho, valTemp] = cart2pol(xVals,yVals,  corrCut);  % convert x, y to polar coordinates
aR = reshape(rho,1, (2*rMax+1)^2);
aValTemp = reshape(valTemp,1, (2*rMax+1)^2);
[rrTemp,ind] = sort(aR);
corrRadVal = aValTemp(ind);
rad= 0:floor(max(rrTemp));
[n, bin] = histcounts(rrTemp, rad);
valRadBinAvg = zeros(length(bin)-1, 1);
valRadBinSem = zeros(length(bin)-1, 1);
for j = 1:length(bin)-1
    valRadBinAvg(j) = sum(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sum(rrTemp>bin(j) & rrTemp<=bin(j+1));
    valRadBinSem(j) = std(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sqrt(sum(rrTemp>bin(j) & rrTemp<=bin(j+1)));                        
end
end
