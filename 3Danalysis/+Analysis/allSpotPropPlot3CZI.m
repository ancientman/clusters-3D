function allSpotPropPlot3CZI(txtFilePath, geneName)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   uses two colour CZI data to calculate distance between the 
%   transciption spot and the bicoid spots.
%   also calculates the diameter of the bicoid spots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtFilePathChar = convertStringsToChars(txtFilePath);

if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
nucStruct = cell(1, totSubDirs);
bcdFitDS = cell(1, totSubDirs);
c1SpotStruct = cell(1, totSubDirs);
c2SpotStruct = cell(1, totSubDirs);
timeFrames = cell(1, totSubDirs);
nucVol = cell(1, totSubDirs);
nucVal = cell(1, totSubDirs);
metaDataStruct = cell(1, totSubDirs);
c1SpotStruct = cell(1, totSubDirs);
c2SpotStruct = cell(1, totSubDirs);
cellC1C2SpotDistUM = cell(1, totSubDirs);
cellC1C2SpotVolUM = cell(1, totSubDirs);
c1SpotsPerNuc = cell(1, totSubDirs);
c1SpotsDist = cell(1, totSubDirs);
sortC1C2spotDistUM = cell(1, totSubDirs);
sortDistIdx = cell(1, totSubDirs);
distSortC1C2SpotVolUM = cell(1, totSubDirs);
distSortC1SpotVal = cell(1, totSubDirs);
close1DistNucTime = cell(1, totSubDirs);
close1VolNucTime = cell(1, totSubDirs);
close2VolNucTime = cell(1, totSubDirs);
close1ValNucTime = cell(1, totSubDirs);
close2ValNucTime = cell(1, totSubDirs);
c2SpotVolNucTime = cell(1, totSubDirs);
close1DistEm = cell(1, totSubDirs);
close1VolEm = cell(1, totSubDirs);
close2VolEm = cell(1, totSubDirs);
close1ValEm = cell(1, totSubDirs);       
close2ValEm = cell(1, totSubDirs);
c2SpotVolEm = cell(1, totSubDirs);
c2SpotMolEm = cell(1, totSubDirs);

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
c1c2Close1SpotDistUM = cell(1, totSubDirs);
c1c2Close2SpotDistUM = cell(1, totSubDirs);
c1c2Close1SpotVolUM = cell(1, totSubDirs);
c1c2Close2SpotVolUM = cell(1, totSubDirs);
c1c2Close1SpotVal = cell(1, totSubDirs);
c1c2Close2SpotVal = cell(1, totSubDirs);
c1SpotVal = cell(1, totSubDirs);
c2SpotVol = cell(1, totSubDirs);
c2SpotVal = cell(1, totSubDirs);
c2SpotMol = cell(1, totSubDirs);
c2SpotMolNucTime = cell(1, totSubDirs);

customFun1 = @(idx, d) d(idx);
customFun2 = @(d) d(1);
customFun3 = @(d) d(2);
init = 1;
% ----------------------------------------------

minDistAll = [];
minDistValSubAll = [];
nucValAll = [];

for i=init:totSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})        
        bcdFitDS{i} = load(append(dataFiles.(fileID), filesep, 'spotFitPropDS.mat'));
        if ~isempty(bcdFitDS{i}.fitDS)
            bcdFitPerNucDS{i} = load(append(dataFiles.(fileID), filesep, 'bcdSpotPerNuc.mat'));
            totNuc = length(bcdFitPerNucDS{i}.spotPropPerNuc);
    
            for j = 1:totNuc
                minDistAll = vertcat(minDistAll, bcdFitPerNucDS{i}.spotPropPerNuc{j}.nucMcpSpotMinDist);
                minDistValSubAll = vertcat(minDistValSubAll, bcdFitPerNucDS{i}.spotPropPerNuc{j}.nucMcpSpotMinDistValSub);
                nucValAll = vertcat(nucValAll, bcdFitPerNucDS{i}.spotPropPerNuc{j}.nucVal);
            end
        end
    end
end

TFSpotProp.geneName = geneName; %" example: p2 all weak";
fprintf("\ngene name = %s\n\n", TFSpotProp.geneName);
% TFSpotProp.minFitDistUM = minFitDistAll;
TFSpotProp.minDistUM = minDistAll;
TFSpotProp.minDistValSub = minDistValSubAll;
TFSpotProp.nucVal = nucValAll;

resultFolder = 'G:\Tyrone_analysis\test\CZI2c';
fileName = append('\', TFSpotProp.geneName, '_DS.mat');
save(append(resultFolder, fileName),'TFSpotProp');
end


% 
% 
% 
%         totUseNuc = length(vertcat(bcdFitDS{i}.fitDS{1}.nuc{:}));
%         bcdMcpFitDistNuc = cell(1, totNuc);
%         bcdMcpDistNuc = cell(1, totNuc);
%         minFitDist{i} = cell(1, totNuc);
%         bcdTotValNuc = cell(1, totNuc);
%         minTotVal{i} = cell(1, totNuc);
%         for j = 1:totNuc
%             if(~isempty(bcdFitDS{i}.fitDS{1}.nuc{j}))
%                 totSpot = length(bcdFitDS{i}.fitDS{1}.nuc{j}.spotFitProp);
%                 for k = 1:totSpot
%                     if ~isempty(bcdFitDS{i}.fitDS{1}.nuc{j}.spotFitProp{k})
%                         bcdMcpFitDistNuc{j} = vertcat(bcdMcpFitDistNuc{j}, bcdFitDS{i}.fitDS{1}.nuc{j}.spotFitProp{k}.bcdMcpFitDist);
%                         bcdTotValNuc{j} = vertcat(bcdTotValNuc{j}, bcdFitDS{i}.fitDS{1}.nuc{j}.spotFitProp{k}.spot2DTotVal);
%                         
%                         if c1SpotStruct{i}.bcdSpotProp{1}{1}(j).bb(k, 4)*c1SpotStruct{i}.bcdSpotProp{1}{1}(j).bb(k, 5) >= bcdSpotXYpixelLim || ...
%                             c1SpotStruct{i}.bcdSpotProp{1}{1}(j).bb(k, 6)>= bcdSpotZpixelLim % only spots above the limit   
% %                             bcdMcpDistNuc{j} = vertcat(bcdMcpDistNuc{j}, c1SpotStruct{i}.bcdSpotProp.bcdSpotProp{i}{1}(j).bcdMcpDist(k));
%                             bcdMcpDistNuc{j} = vertcat(bcdMcpDistNuc{j}, c1c2SpotStruct{i}.bcdMcpSpotProp{1}{1}.C1C2SpotDistUM{j}(k));
%                        
%                         end
%                     end
%                 end
% 
%                 [minFitDist{i}{j}, minFitIdx] = min(bcdMcpFitDistNuc{j});
%                 minFitDistAll = vertcat(minFitDistAll, minFitDist{i}{j});                    
%                 minTotVal{i}{j} = bcdTotValNuc{j}(minFitIdx);
%                 minTotValAll = vertcat(minTotValAll, minTotVal{i}{j});
% 
%                 
% 
%                 [minDist{i}{j}, minIdx] = min(bcdMcpDistNuc{j});
%                 minDistAll = vertcat(minDistAll, minDist{i}{j});
%             end
%         end
%     end
% end


% TFSpotProp.geneName = geneName; %" example: p2 all weak";
% fprintf("\ngene name = %s\n\n", TFSpotProp.geneName);
% TFSpotProp.minFitDistUM = minFitDistAll;
% TFSpotProp.minDistUM = minDistAll;
% TFSpotProp.minTotVal = minTotValAll;
% 
% resultFolder = 'G:\Tyrone_analysis\test\CZI2c';
% fileName = append('\', TFSpotProp.geneName, '_DS.mat');
% save(append(resultFolder, fileName),'TFSpotProp');



% xArr = 93*(1:length(c2SpotMolNucTime{1}(5,:)));
% plotDouble(xArr, c2SpotMolNucTime{1}(5,:), close1DistNucTime{1}(5,:), xArr.*0)
% plotDouble(xArr, c2SpotMolNucTime{1}(1,:), close1DiaNucTime{1}(1,:), xArr.*0)

% 
% allNucVal = nonzeros(allNucVal(~isnan(allNucVal)));
% allC1C2SpotDistUM = nonzeros(allC1C2SpotDistUM(~isnan(allC1C2SpotDistUM)));
% allC1C2SpotVolUM = nonzeros(allC1C2SpotVolUM(~isnan(allC1C2SpotVolUM)));
% allCIC2meanInterC1Dist = nonzeros(allCIC2meanInterC1Dist(~isnan(allCIC2meanInterC1Dist)));
% allC1C2Close1SpotDistUM = nonzeros(allC1C2Close1SpotDistUM(~isnan(allC1C2Close1SpotDistUM)));
% allC1C2Close2SpotDistUM = nonzeros(allC1C2Close2SpotDistUM(~isnan(allC1C2Close2SpotDistUM)));
% allC1C2Close1SpotVolUM = nonzeros(allC1C2Close1SpotVolUM(~isnan(allC1C2Close1SpotVolUM)));
% allC1C2Close2SpotVolUM = nonzeros(allC1C2Close2SpotVolUM(~isnan(allC1C2Close2SpotVolUM)));
% allC1C2SpotDiaUM = power((2*allC1C2SpotVolUM), (1/3));
% allC1C2Close1SpotDiaUM = power(2*(allC1C2Close1SpotVolUM), (1/3));
% allC1C2Close2SpotDiaUM = power(2*(allC1C2Close2SpotVolUM), (1/3));
% allC1C2Close1SpotVal = nonzeros(allC1C2Close1SpotVal(~isnan(allC1C2Close1SpotVal)));
% allC1C2Close2SpotVal = nonzeros(allC1C2Close2SpotVal(~isnan(allC1C2Close2SpotVal)));
% allC2SpotVal = nonzeros(allC2SpotVal(~isnan(allC2SpotVal)));
% allC2SpotVol = nonzeros(allC2SpotVol(~isnan(allC2SpotVol)));
% 
% TFSpotProp.geneName = geneName; %" example: p2 all weak";
% fprintf("\ngene name = %s\n\n", TFSpotProp.geneName);
% TFSpotProp.nucVal = allNucVal;
% TFSpotProp.close1DistUM = allC1C2Close1SpotDistUM;
% TFSpotProp.close2DistUM = allC1C2Close2SpotDistUM;
% TFSpotProp.close1DiaUM = allC1C2Close1SpotDiaUM;
% TFSpotProp.close2DiaUM = allC1C2Close2SpotDiaUM;
% TFSpotProp.close1Val = allC1C2Close1SpotVal;
% TFSpotProp.close2Val = allC1C2Close2SpotVal;
% TFSpotProp.c2Val = allC2SpotVal;
% TFSpotProp.c2Vol = allC2SpotVol;
% 
% resultFolder = 'C:\xData\Rawdata_RM\test\test_spot';
% fileName = append('\', TFSpotProp.geneName, '_DS.mat');
% save(append(resultFolder, fileName),'TFSpotProp');
% 
% plotColor = [0 1 1];
% plotHist(allC1C2Close1SpotDistUM, geneName, plotColor);
% 
% % Analysis.barPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1);
% % Analysis.barPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2);
% 
% % Analysis.boxPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1, 'light');
% % Analysis.boxPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2, 'light');
% 
% % Analysis.scatterBoxPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1, 'light');
% % Analysis.scatterBoxPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2, 'light');
% end

function plotHist(data, type, plotColor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram plotter and fitter
% >>>>>>>>>Inputs<<<<<<<<<<
% data
% type: string for plot title
% plotColor: values in rgb, eg [1 0 1] for plotting.
%......................................................................................................%
% The initialization parameters will depend on the data passed. Make sure
% to change all before running the program.
% theme: choice dark or light themes
% binType: manually determine bin spacing etc. Check before running
% fitType: Gaussian mixture or Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binType = 2; 
fitType = 'gaussian'; % oisson or gaussian or halfnormal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color', 'w');

if binType == 0 % fix bin numbers
    [xMin, xMax] = bounds(data);
    BE = linspace(xMin, xMax, 10);
    hIn = histogram(data, BE, 'Normalization', 'probability');
elseif binType == 1 % fix bin position
    BE = 0:0.1:3; 
    hIn = histogram(data, BE, 'Normalization', 'probability');
elseif binType == 2  % auto
    hIn = histogram(data, 'Normalization', 'probability');
end

%% Plot 1
axes1 = gca;
set(axes1,'FontSize',12,'LineWidth',1.5,'XColor',...
[0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
set(hIn,'Parent',axes1, 'LineWidth',1.5,...
    'EdgeColor',plotColor, 'FaceColor',plotColor./2);

hold on;

%% For gaussian mixture fit
if strcmp(fitType, 'gaussian')
    trials = 1; % change max trails here
    pdIn = cell(1,trials);
    AIC = zeros(1,trials);
    options = statset('MaxIter',500);
    hold on;
    for k = 1:trials
        pdIn{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
        AIC(k)= pdIn{k}.AIC;
    end
    [minAIC,numComponents] = min(AIC);
    bestModel = pdIn{numComponents};
    pdIn = bestModel;
    % pdIn = fitgmdist(data,3);
    pdfIn = pdf(pdIn,data);
    pdfIn = pdfIn*sum(hIn.Values * hIn.BinWidth); %pdf times area under the histogram curve
    [data, idxIn] = sort(data); %sort raw data, det indices
    pdfIn = pdfIn(idxIn); %sort y per those indices
    meanIn = sort(pdIn.mu);
    pIn = plot(data,pdfIn,'-', 'linewidth', 2);
    set(pIn, 'Parent',axes1, 'LineWidth',2,'LineStyle','-.', 'Color',plotColor);
    
    hold on;

    %% For poisson fit
elseif strcmp(fitType, 'poisson')
    pdIn = cell(1,2);
    pdIn = fitdist(data,'Poisson');
    if binType == 0
        p = ceil(xMin):ceil(xMax);
        pdfIn = poisspdf(p, pdIn.lambda);
        pdfIn = pdfIn*sum(hIn.Values * hIn.BinWidth);
    elseif binType == 1
        p = BE;
        increment = (1/(BE(2) - BE(1)));
        pdfIn = poisspdf((p.*increment), increment*pdIn.lambda);
    end   
%     pIn = scatter(p,pdfIn,20, plotColor, 'filled', 'o', 'MarkerEdgeColor',[0 0 0]);
    pIn = plot(BE, pdfIn, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor',plotColor);
    hold on;
elseif strcmp(fitType, 'halfnormal')
    pdIn = fitdist(data,'hn');
    pdfIn = pdf(pdIn, BE);
    pdfIn = pdfIn*sum(hIn.Values * hIn.BinWidth);
    pIn = plot(BE,pdfIn,'-', 'linewidth', 2, 'Color', plotColor);
end

% % global plot properties
xLabel = '{distance from mRNA hotspot {\mu}m}';
yLabel = 'p';
ylabel(yLabel);
xlabel(xLabel);
title (type);
box(axes1,'on');
gca;
grid off;
if binType==0
    xlim([xMin, xMax]);
else
    xlim([0, 3]);
end
hold off;

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