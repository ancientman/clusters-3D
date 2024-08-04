function allSpotPropPlot2TimeNew(txtFilePath, geneName)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   uses two colour data to calculate distance between the 
%   transciption spot and the bicoid spots.
%   also calculates the diameter of the bicoid spots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtFilePathChar = convertStringsToChars(txtFilePath);

minVolFactor = 3*3*3;
minVol = minVolFactor*0.043*0.043*0.2;

minBoundingBoxZ = 1; % earlier 3
minBoundingBoxXY = 1*1;

if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
nucProp = cell(1, totalSubDirs);
c1c2SpotStruct = cell(1, totalSubDirs);
c2SpotStruct = cell(1, totalSubDirs);
timeFrames = cell(1, totalSubDirs);
nucVolUM = cell(1, totalSubDirs);
nucVal = cell(1, totalSubDirs);
metaDataStruct = cell(1, totalSubDirs);
c1SpotProp = cell(1, totalSubDirs);
c2SpotProp = cell(1, totalSubDirs);
cellC1C2SpotDistUM = cell(1, totalSubDirs);
cellC1C2SpotVolUM = cell(1, totalSubDirs);
c1SpotsPerNuc = cell(1, totalSubDirs);
c1SpotsDist = cell(1, totalSubDirs);
sortC1C2spotDistUM = cell(1, totalSubDirs);
sortDistIdx = cell(1, totalSubDirs);
distSortC1C2SpotVolUM = cell(1, totalSubDirs);
distSortC1SpotVal = cell(1, totalSubDirs);
close1DistNucTime = cell(1, totalSubDirs);
close1VolNucTime = cell(1, totalSubDirs);
close2VolNucTime = cell(1, totalSubDirs);
close1ValNucTime = cell(1, totalSubDirs);
close2ValNucTime = cell(1, totalSubDirs);
c2SpotVolNucTime = cell(1, totalSubDirs);
close1DistEm = cell(1, totalSubDirs);
close1VolEm = cell(1, totalSubDirs);
close2VolEm = cell(1, totalSubDirs);
close1ValEm = cell(1, totalSubDirs);       
close2ValEm = cell(1, totalSubDirs);
c2SpotVolEm = cell(1, totalSubDirs);
c2SpotMolEm = cell(1, totalSubDirs);


allNucVal = [];
allC1C2SpotDistUM = [];
allC1C2SpotVolUM = [];
allCIC2meanInterC1Dist = [];
allC1C2Close1SpotDistUM = [];
allC1C2Close2SpotDistUM = [];
allC1C2Close1SpotVolUM = [];
allC1C2Close2SpotVolUM = [];
allC1C2Close1SpotVal = [];
allC1C2Close2SpotVal = [];
allC2SpotVal = [];
allC2SpotVol = [];
c1c2Close1SpotDistUM = cell(1, totalSubDirs);
c1c2Close2SpotDistUM = cell(1, totalSubDirs);
c1c2Close1SpotVolUM = cell(1, totalSubDirs);
c1c2Close2SpotVolUM = cell(1, totalSubDirs);
c1c2Close1SpotVal = cell(1, totalSubDirs);
c1c2Close2SpotVal = cell(1, totalSubDirs);
c1SpotVal = cell(1, totalSubDirs);
c1SpotVol = cell(1, totalSubDirs);
c2SpotVal = cell(1, totalSubDirs);
c2SpotVol = cell(1, totalSubDirs);
c1SpotValNucSub = cell(1, totalSubDirs);
c1c2SpotDist = cell(1, totalSubDirs);
c1c2SpotDistSort = cell(1, totalSubDirs);
c1c2SpotValSort = cell(1, totalSubDirs);
c1c2SpotVolSort = cell(1, totalSubDirs);
numSpot = cell(1, totalSubDirs);
c2SpotMol = cell(1, totalSubDirs);
c2SpotMolNucTime = cell(1, totalSubDirs);

customFun1 = @(idx, d) d(idx);
customFun2 = @(d) d(1);
customFun3 = @(d) d(2);
init = 1;
for i=init:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        meta = load(append(dataFiles.(fileID), filesep, 'metaDataDS.mat'));
        scaling = [meta.metaDataDS.analysisInfo.yPixUM, meta.metaDataDS.analysisInfo.xPixUM, meta.metaDataDS.analysisInfo.zPixUM];
        nucProp{i} = load(append(dataFiles.(fileID), filesep, 'c1NucPropDS.mat'));
        c1SpotProp{i} = load(append(dataFiles.(fileID), filesep, 'c1SpotPropDS.mat'));
        c2SpotProp{i} = load(append(dataFiles.(fileID), filesep, 'c2SpotPropDS.mat'));

        timeFrames{i} = length(nucProp{i}.c1NucProp); % check   
        c1SpotVal{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1SpotVol{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c2SpotVal{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c2SpotVol{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1SpotValNucSub{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1c2SpotDist{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1c2SpotDistSort{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1c2SpotValSort{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        c1c2SpotVolSort{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        numSpot{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));

        nucVolUM{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        nucVal{i} = cell(1, length(c1SpotProp{i}.c1SpotProp{1}));
        
        for j = 1:length(c1SpotProp{i}.c1SpotProp{1}) % total nuclei
            nucVal{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            nucVolUM{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1SpotVal{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1SpotVol{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c2SpotVal{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c2SpotVol{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1SpotValNucSub{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1c2SpotDist{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1c2SpotDistSort{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1c2SpotValSort{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            c1c2SpotVolSort{i}{j} = cell(1, length(nucProp{i}.c1NucProp));
            numSpot{i}{j} = cell(1, length(nucProp{i}.c1NucProp));

            for t = 1:timeFrames{i} % total time
                if ~isempty(nucProp{i}.c1NucProp{t})
                    if numel(c2SpotProp{i}.c2SpotProp{t}(j).vol) == 1
                        nucVolUM{i}{j}{t} = nucProp{i}.c1NucProp{t}(j).volUM;  
%                         nucVal{i}{t}{j} = nucProp{i}.c1NucProp{t}(j).meanVal;  
                        nucVal{i}{j}{t} = mean(nucProp{i}.c1NucProp{t}(j).voxVal{1});  
                        c2SpotVal{i}{j}{t} = mean(c2SpotProp{i}.c2SpotProp{t}(j).voxVal{1});
                        c2SpotVol{i}{j}{t} = c2SpotProp{i}.c2SpotProp{t}(j).volUM;
                        c1ValTemp = [];
                        c1VolTemp = [];
                        c1ValSubTemp = [];
                        c1DistTemp = [];
                        for k = 1:length(c1SpotProp{i}.c1SpotProp{t}(j).voxVal) % total spots (unfiltered)
%                         if c1SpotProp{i}.c1SpotProp{t}(j).bb(k,4)*c1SpotProp{i}.c1SpotProp{t}(j).bb(k,5) >= 36 && c1SpotProp{i}.c1SpotProp{t}(j).bb(k,6)>=2
%                             if c1SpotProp{i}.c1SpotProp{t}(j).volUM(k) >= minVol
%                             if c1SpotProp{i}.c1SpotProp{t}(j).bb(k,6)>=minBoundingBoxZ
                            if c1SpotProp{i}.c1SpotProp{t}(j).bb(k,4)*c1SpotProp{i}.c1SpotProp{t}(j).bb(k,5)>=minBoundingBoxXY
%                                 c1ValTemp = vertcat(c1ValTemp, mean(c1SpotProp{i}.c1SpotProp{t}(j).voxVal{k}));
%                                 c1VolTemp = vertcat(c1ValTemp, c1SpotProp{i}.c1SpotProp{t}(j).volUM{k});
%                                 c1ValSubTemp = vertcat(c1ValSubTemp, vertcat(c1SpotValNucSub{i}{t}{j}, mean(c1SpotProp{i}.c1SpotProp{t}(j).voxVal{k}) - nucVal{i}{t}{j}));
%                                 c1DistTemp = vertcat(c1DistTemp, pdist2(scaling.*c1SpotProp{i}.c1SpotProp{t}(j).center(k,:), scaling.*c2SpotProp{i}.c2SpotProp{t}(j).center));
                                c1SpotVal{i}{j}{t}(k,:) = mean(c1SpotProp{i}.c1SpotProp{t}(j).voxVal{k});
                                c1SpotVol{i}{j}{t}(k,:) = c1SpotProp{i}.c1SpotProp{t}(j).volUM(k);
                                c1SpotValNucSub{i}{j}{t}(k,:) = mean(c1SpotProp{i}.c1SpotProp{t}(j).voxVal{k}) - nucVal{i}{j}{t};
                                c1c2SpotDist{i}{j}{t}(k,:) = pdist2(scaling.*c1SpotProp{i}.c1SpotProp{t}(j).center(k,:), scaling.*c2SpotProp{i}.c2SpotProp{t}(j).center);                                                            
                            else
%                                 c1ValTemp = vertcat(c1ValTemp, NaN);
%                                 c1VolTemp = vertcat(c1VolTemp, NaN);
                                c1c2SpotDist{i}{j}{t}(k,:) = NaN;
                                c1SpotVal{i}{j}{t}(k,:) = NaN;
                                c1SpotVol{i}{j}{t}(k,:) = NaN;
                                c1SpotValNucSub{i}{j}{t}(k,:) = NaN;
                            end
                        end         
                    else
                        nucVolUM{i}{j}{t} = [];
                        nucVal{i}{j}{t} = [];
                        c2SpotVal{i}{j}{t} = [];
                        c2SpotVol{i}{j}{t} = [];
                        c1SpotVal{i}{j}{t} = [];
                        c1SpotVol{i}{j}{t} = [];
                        c1SpotValNucSub{i}{j}{t} = [];
                        c1c2SpotDist{i}{j}{t} = [];                                                            
                    end

                    tempDist = c1c2SpotDist{i}{j}{t};
                    tempDist = tempDist(~isnan(tempDist));                    
                    [c1c2SpotDistSort{i}{j}{t}, sortIdx] = sort(tempDist);
                    numSpot{i}{j}{t} = numel(tempDist);
                    tempVal = c1SpotValNucSub{i}{j}{t};
                    tempVal = tempVal(~isnan(tempVal));
                    tempVal = tempVal(sortIdx);

                    if isempty(tempVal)
                        nucVal{i}{j}{t} = [];
                        nucVolUM{i}{j}{t} = [];
                    end
                    c1c2SpotValSort{i}{j}{t} = tempVal;
                    tempVol = c1SpotVol{i}{j}{t};
                    tempVol = tempVol(~isnan(tempVol));
                    tempVol = tempVol(sortIdx);
                    c1c2SpotVolSort{i}{j}{t} = tempVol;


                end
            end
        end
    end
end

TFSpotProp.geneName = geneName; %" example: p2 all weak";
fprintf("\ngene name = %s\n\n", TFSpotProp.geneName);
TFSpotProp.nucVal = nucVal;
TFSpotProp.nucVol = nucVolUM; %in um
TFSpotProp.numSpot = numSpot;
TFSpotProp.distSortUM = c1c2SpotDistSort;
TFSpotProp.c1ValSort = c1c2SpotValSort; % nuc subbed val
TFSpotProp.c1VolSort = c1c2SpotVolSort;
TFSpotProp.c2Val = c2SpotVal; 
TFSpotProp.c2Vol = c2SpotVol; % in um

resultFolder = 'G:\Dropbox (Princeton)\bcd_ss_paper\data\aa_val\geneCombineSpot2';
resultFolder = 'G:\Tyrone_analysis\newtest\test';
fileName = append('\', TFSpotProp.geneName, '_DS.mat');
save(append(resultFolder, fileName),'TFSpotProp');
end

function plotDouble(xArr, y1, y2Mean, y2Dev)

fig = figure('Color', 'w');
leftColor = [0.9 0.3 0.3];
rightColor = [0.3 0.7 0.5];

yyaxis left;
pp = plot(xArr, y1);
pp.LineWidth = 1.5;
pp.Color = leftColor;
pp.LineStyle = '-';

ylabel('MS2 Intensity');
set(gca,'YTick',[]);
set(gca,'ycolor',leftColor) ;
hold on;

yyaxis right;
patchTop = y2Mean+y2Dev;
patchTop = reshape(patchTop,1,[]);
patchBot = y2Mean-y2Dev;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch = [xArr,fliplr(xArr)];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = rightColor;
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.4;
hold on;
pl = plot(xArr, y2Mean);
pl.LineWidth = 1.5;
pl.Color = rightColor;
pl.LineStyle = '-';

ylabel('Cluster dist. {\mu}m');
% set(gca,'YTick',[]);
set(gca,'ycolor',rightColor) ;
hold on;
xlabel('Time (s)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=500;
plotHeight=200;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(fig,'defaultAxesColorOrder',[leftColor; rightColor]);
hold off;
end