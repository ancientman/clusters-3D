function allSpotPropPlot1Time(txtFilePath, geneName)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   uses two colour data to calculate distance between the 
%   transciption spot and the bicoid spots.
%   also calculates the diameter of the bicoid spots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtFilePathChar = convertStringsToChars(txtFilePath);

minVolFactor = 4*4*2;
minVol = minVolFactor*0.042*0.042*0.2;

if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
nucStruct = cell(1, totalSubDirs);
c1c2SpotStruct = cell(1, totalSubDirs);
c2SpotStruct = cell(1, totalSubDirs);
timeFrames = cell(1, totalSubDirs);
nucVol = cell(1, totalSubDirs);
nucVal = cell(1, totalSubDirs);
metaDataStruct = cell(1, totalSubDirs);
c1SpotPropStruct = cell(1, totalSubDirs);
c2SpotPropStruct = cell(1, totalSubDirs);
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
c2SpotVol = cell(1, totalSubDirs);
c2SpotVal = cell(1, totalSubDirs);
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
%         meta = load(append(dataFiles.(fileID), filesep, 'seriesMetaDataDS.mat'));
        minVolFactor = (meta.metaDataDS.imagingInfo.XYpsf)^2*meta.metaDataDS.imagingInfo.Zpsf/(1000^3*(meta.metaDataDS.analysisInfo.xPixUM)^2*meta.metaDataDS.analysisInfo.zPixUM);
        minVol = minVolFactor*meta.metaDataDS.analysisInfo.xPixUM*meta.metaDataDS.analysisInfo.yPixUM*meta.metaDataDS.analysisInfo.zPixUM;
        nucStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1NucPropDS.mat'));
        c1SpotPropStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1SpotPropDS.mat'));
        timeFrames{i} = length(c1SpotPropStruct{i}.c1SpotProp);
        c1SpotVal{i} = cell(1, timeFrames{i});
        
        for t = 1:timeFrames{i}
            if ~isempty(nucStruct{i}.c1NucProp{t})
                nucVol{i}{t} = vertcat(nucStruct{i}.c1NucProp{t}.volUM);   
                nucVal{i}{t} = vertcat(nucStruct{i}.c1NucProp{t}.meanVal);   
                for n  =1:length(nucVol{i}{t})
                    spotsPerNuc{i}(t,n) = length(c1SpotPropStruct{i}.c1SpotProp{t}(n).volUM);
                    spotsAvgVol{i}(t,n) = mean(c1SpotPropStruct{i}.c1SpotProp{t}(n).volUM);
                    spotsStdVol{i}(t,n) = std(c1SpotPropStruct{i}.c1SpotProp{t}(n).volUM);
                    spotsAvgVal{i}(t,n) = mean(cellfun(@mean, c1SpotPropStruct{i}.c1SpotProp{t}(n).voxVal));
                    spotsStdVal{i}(t,n) = std(cellfun(@std, c1SpotPropStruct{i}.c1SpotProp{t}(n).voxVal));
                    totValTemp = cellfun(@mean, c1SpotPropStruct{i}.c1SpotProp{t}(n).voxVal).*(c1SpotPropStruct{i}.c1SpotProp{t}(n).volUM);
                    spotsAvgTotVal{i}(t,n) = mean(totValTemp);
                    spotsStdTotVal{i}(t,n) = std(totValTemp);
                end
            end
        end
        nucVol{i} = horzcat(nucVol{i}{:})';
        nucVal{i} = horzcat(nucVal{i}{:})';
        spotsPerNucNorm{i} = (spotsPerNuc{i} - spotsPerNuc{i}(1,:))./spotsPerNuc{i}(1,:);
        spotsAvgVolNorm{i} = (spotsAvgVol{i} - spotsAvgVol{i}(1,:))./spotsAvgVol{i}(1,:);
        spotsAvgValNorm{i} = (spotsAvgVal{i} - spotsAvgVal{i}(1,:))./spotsAvgVal{i}(1,:);
    end
end

nucValNorm = horzcat(nucVal{:});
nucValNorm = (nucValNorm - nucValNorm(1,:))./nucValNorm(1,:);
spotsPerNucNormMean = mean(horzcat(spotsPerNucNorm{:}), 2);
spotsPerNucNormSem = std(horzcat(spotsPerNucNorm{:}), 0, 2)./sqrt(size(horzcat(spotsPerNucNorm{:}), 2));
spotsAvgVolNormMean = mean(horzcat(spotsAvgVolNorm{:}), 2);
spotsAvgVolNormSem = std(horzcat(spotsAvgVolNorm{:}), 0, 2)./sqrt(size(horzcat(spotsPerNucNorm{:}), 2));
spotsAvgValNormMean = mean(horzcat(spotsAvgValNorm{:}), 2);
spotsAvgValNormSem = std(horzcat(spotsAvgValNorm{:}), 0, 2)./sqrt(size(horzcat(spotsPerNucNorm{:}), 2));


% spotsAvgValNormAdj = horzcat(spotsAvgVal{:})./horzcat(nucVal{:});
% spotsAvgValNormAdjMean = mean(spotsAvgValNormAdj, 2);
% spotsAvgValNormAdjSem = std(spotsAvgValNormAdj, 0, 2)./sqrt(size(spotsAvgValNormAdj, 2));


spotsAvgValNormAdj = horzcat(spotsAvgVal{:})./horzcat(nucVal{:});
spotsAvgValNormAdj = (spotsAvgValNormAdj - spotsAvgValNormAdj(1,:))./(spotsAvgValNormAdj(1,:));
spotsAvgValNormAdjMean = mean(spotsAvgValNormAdj, 2);
spotsAvgValNormAdjSem = std(spotsAvgValNormAdj, 0, 2)./sqrt(size(spotsAvgValNormAdj, 2));

tArr = 27.6.*(1:20);
figure('color', 'w');
p1 = plotter1(tArr, spotsPerNucNormMean, spotsPerNucNormSem);
hold on;
p2 = plotter1(tArr, spotsAvgValNormAdjMean, spotsAvgValNormAdjSem);
hold on;
p3 = plotter1(tArr, spotsAvgVolNormMean, spotsAvgVolNormSem);

legTex{1} = "Cluster density";
legTex{2} = "Adjusted cluster intensity";
legTex{3} = "Cluster volume";
leg =legend([p1, p2, p3], legTex{:});
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off')
end
         
function p1 = plotter1(xArr, yMean, yDev)
p1 = errorbar(xArr, yMean, yDev);
p1.Marker = 'o';
p1.LineWidth = 1;
% pl.Color = color;
p1.LineStyle = '--';
p1.CapSize = 0;
hold on;
% %-----------------------------------------
ax = gca;
ylabel('Cluster Paramteres');
xlabel('Time (s)');
% xlim([0.17, 0.6]);
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=500;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
% %---------------------------------------
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