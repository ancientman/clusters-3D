function combine2CSpotPlotter5Mod(folderPath)
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x)load(append(folderPath, filesep, x)), fileNames, 'un', 0);

distLim =  [0.45 0.45];
% distLim = [0.67 0.67 0.67];
distLim = [1 1];
nucValUB = 400;
nucValLB = 100;

minVolFactor = 4*4*3;
% minVolFactor = 2*2*1;
minVol = minVolFactor*0.043*0.043*0.2;

nucValCutOff = [0 1000]; % Lower, Upper

totalKMeans = 10;
dataCutOff = 40;


colorStruct{6} = [210, 59, 227; 0 0 0];
colorStruct{5} = [138, 62, 8; 0 0 0];
colorStruct{4} = [213,62,79; 0, 0, 0];
colorStruct{3} = [153,112,171; 0, 0, 0];
colorStruct{2} = [102,194,165; 0, 0, 0];
colorStruct{1} = [50,136,189; 0, 0, 0];

%------- Manual Reordering of data Sequence --------%

groupOrder = [2 1];

% groupOrder = [1 2 3];
groupOrder = 1:length(struct2C); % default
%------------------------------------------------------------------%

geneName = strings(1, length(struct2C));
% geneNameStruct = {'Strong', 'Weak'};
groupColor = zeros(length(struct2C), 3);

c1ValAll = cell(1, length(struct2C));
c1ValMean = cell(1, length(struct2C));
c1VolAll = cell(1, length(struct2C));
c1volUnfiltAll = cell(1, length(struct2C));
distAll = cell(1, length(struct2C));
nucValRepAll = cell(1, length(struct2C));
close1DistAll = cell(1, length(struct2C));
close1ValAll = cell(1, length(struct2C));
close1VolAll = cell(1, length(struct2C));
c2ValRepAll = cell(1, length(struct2C));
numSpotsAll = cell(1, length(struct2C));

c1ValClose1 = cell(1, length(struct2C));
c1VolClose1 = cell(1, length(struct2C));
c1DistClose1 = cell(1, length(struct2C));

for i = 1:length(struct2C)    
    j = groupOrder(i);
    geneName(i) = struct2C{j}.TFSpotProp.geneName;   
    groupColor(i,:) = colorStruct{j}(1,:);
%     geneName(i) = geneNameStruct{j};
    c1Vol{i} = vertcat(struct2C{j}.TFSpotProp.c1VolSort{:}); % split to embryos (in um)
    nucVal{i} = vertcat(struct2C{j}.TFSpotProp.nucVal{:}); % split to embryos    
    c2Val{i} = vertcat(struct2C{j}.TFSpotProp.c2Val{:}); % split to embryos    
    c1Val{i} = vertcat(struct2C{j}.TFSpotProp.c1ValSort{:}); % split to embryos (nuc values are subbed)
    c1VolMean{i} = cell(1, length(c1Vol{i}));
    c1ValMean{i} = cell(1, length(c1Val{i}));
    c1TotValMean{i} = cell(1, length(c1Val{i}));
    c1Dist{i} = vertcat(struct2C{j}.TFSpotProp.distSortUM{:}); % split to embryos
    for m=1:length(c1Val{i}) % total nuclei in each embryo
%         nucVal{i}{m} = cellfun(@(x,y) x(y), nucVal{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0);        
        c1volUnfiltAll{i} = vertcat(c1volUnfiltAll{i}, vertcat(c1Vol{i}{m}{:}));
        c1Vol{i}{m} = cellfun(@(x,y) x(y), c1Vol{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0);    
        c1VolMean{i}{m} = cellfun(@(x) mean(x), c1Vol{i}{m}, 'un', 0);    
        c1VolAll{i} = vertcat(c1VolAll{i}, vertcat(c1Vol{i}{m}{:}));
        c1VolClose1{i} = vertcat(c1VolClose1{i}, cellfun(@(x) x(1), c1Vol{i}{m}(cellfun(@(x) ~isempty(x), c1Vol{i}{m})))');
        
        c1Val{i}{m} = cellfun(@(x,y) x(y), c1Val{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0);        

        c1ValMean{i}{m} = cellfun(@(x) mean(x), c1Val{i}{m}, 'un', 0);        
        c1ValAll{i} = vertcat(c1ValAll{i}, vertcat(c1Val{i}{m}{:}));
        c1ValClose1{i} = vertcat(c1ValClose1{i}, cellfun(@(x) x(1), c1Val{i}{m}(cellfun(@(x) ~isempty(x), c1Val{i}{m})))');          
        
        c1TotValMean{i}{m} = cellfun(@mean, cellfun(@(x,y) x.*y, c1Val{i}{m}, c1Vol{i}{m}, 'un', 0), 'un', 0);

        c1Dist{i}{m} = cellfun(@(x,y) x(y), c1Dist{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0); 
        distAll{i} = vertcat(distAll{i}, vertcat(c1Dist{i}{m}{:}));
        c1DistClose1{i} = vertcat(c1DistClose1{i}, cellfun(@(x) x(1), c1Dist{i}{m}(cellfun(@(x) ~isempty(x), c1Dist{i}{m})))');

        numSpotsAll{i} = vertcat(numSpotsAll{i}, nonzeros(cellfun(@length, c1Val{i}{m}')));
        
        if ~isempty(nonzeros(cellfun(@length, c1Val{i}{m}')))
            nucRepTemp = arrayfun(@(x,y) repelem(x, y), vertcat(nucVal{i}{m}{:}), nonzeros(cellfun(@length, nucVal{i}{m})'), 'un', 0);
            nucRepTemp = cellfun(@transpose, nucRepTemp, 'un', 0);
            nucRepTemp = vertcat(nucRepTemp{:});
        else
            nucRepTemp = [];
        end
        nucValRepAll{i} = vertcat(nucValRepAll{i}, nucRepTemp);
    end

    nucVal{i} = horzcat(nucVal{i}{:});
    nucVal{i} = vertcat(nucVal{i}{:});
    c2Val{i} = horzcat(c2Val{i}{:});
    c2Val{i} = vertcat(c2Val{i}{:});
    c1ValMean{i} = horzcat(c1ValMean{i}{:});
    c1ValMean{i} = vertcat(c1ValMean{i}{:});
    c1VolMean{i} = horzcat(c1VolMean{i}{:});
    c1VolMean{i} = vertcat(c1VolMean{i}{:});
    c1TotValMean{i} = horzcat(c1TotValMean{i}{:});
    c1TotValMean{i} = vertcat(c1TotValMean{i}{:});

    nucValTrim{i} = nucVal{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c2ValTrim{i} = c2Val{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1ValMeanTrim{i} = c1ValMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1VolMeanTrim{i} = c1VolMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1TotValMeanTrim{i} = c1TotValMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
       
    coupledDist{i} = c1DistClose1{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    coupledVal{i} = c1ValClose1{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    coupledVol{i} = c1VolClose1{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    coupledNucVal{i} = nucVal{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1ValMeanTrimCoupled{i} = c1ValMean{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1VolMeanTrimCoupled{i} = c1VolMean{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1TotValMeanTrimCoupled{i} = c1TotValMean{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));

    uncoupledDist{i} = c1DistClose1{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    uncoupledVal{i} = c1ValClose1{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    uncoupledVol{i} = c1VolClose1{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    uncoupledNucVal{i} = nucVal{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1ValMeanTrimUncoupled{i} = c1ValMean{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1VolMeanTrimUncoupled{i} = c1VolMean{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1TotValMeanTrimUncoupled{i} = c1TotValMean{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));

    nucValRepAllTrim{i} = nucValRepAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));
    c1ValAllTrim{i} = c1ValAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));
    c1VolAllTrim{i} = c1VolAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));

    coupledC2Val{i} = c2Val{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    uncoupledC2Val{i} = c2Val{i}(c1DistClose1{i}>distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
end

nucValBin = cell(1, length(struct2C));
coupledNucValBin = cell(1, length(struct2C));
uncoupledNucValBin = cell(1, length(struct2C));
coupledDistBin = cell(1, length(struct2C));
coupledValBin = cell(1, length(struct2C));
coupledVolBin = cell(1, length(struct2C));
coupledTotValBin = cell(1, length(struct2C));
coupledDiaBin = cell(1, length(struct2C));

coupledAllValBin = cell(1, length(struct2C));
coupledAllVolBin = cell(1, length(struct2C));
coupledAllDiaBin = cell(1, length(struct2C));
coupledAllTotValBin = cell(1, length(struct2C));

uncoupledDistBin = cell(1, length(struct2C));
uncoupledValBin = cell(1, length(struct2C));
uncoupledVolBin = cell(1, length(struct2C));
coupledC2ValBin = cell(1, length(struct2C));
uncoupledC2ValBin = cell(1, length(struct2C));

f1 = figure('Color', 'w');
f2 = figure('Color', 'w');
hold on;

for i = 1:length(struct2C)   
    nucValBin{i} = vertcat(coupledNucVal{i}, uncoupledNucVal{i});    
    [nucValBinSort, allIdx, groupLengthSort] = binInxFun(vertcat(coupledNucVal{i}, uncoupledNucVal{i}), totalKMeans, dataCutOff);
    for j = 1:length(nucValBinSort)-1
        coupledNucValBin{i}{j} = coupledNucVal{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        uncoupledNucValBin{i}{j} = uncoupledNucVal{i}(uncoupledNucVal{i}>nucValBinSort(j) & uncoupledNucVal{i}<=nucValBinSort(j+1));

%         coupledAllC1ValBin{i}{j} = 


        coupledDistBin{i}{j} = coupledDist{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        coupledValBin{i}{j} = coupledVal{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        coupledVolBin{i}{j} = coupledVol{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        coupledDiaBin{i}{j} = (6/pi.*(coupledVolBin{i}{j})).^(1/3);
        coupledTotValBin{i}{j} = coupledValBin{i}{j}.*coupledVolBin{i}{j};

        coupledAllValBin{i}{j} = c1ValMeanTrimCoupled{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        coupledAllVolBin{i}{j} = c1VolMeanTrimCoupled{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        coupledAllDiaBin{i}{j} = (6/pi.*(coupledAllVolBin{i}{j})).^(1/3);
        coupledAllTotValBin{i}{j} = c1TotValMeanTrimCoupled{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));

        uncoupledDistBin{i}{j} = uncoupledDist{i}(uncoupledNucVal{i}>nucValBinSort(j) & uncoupledNucVal{i}<=nucValBinSort(j+1));
        uncoupledValBin{i}{j} = uncoupledVal{i}(uncoupledNucVal{i}>nucValBinSort(j) & uncoupledNucVal{i}<=nucValBinSort(j+1));
        uncoupledVolBin{i}{j} = uncoupledVol{i}(uncoupledNucVal{i}>nucValBinSort(j) & uncoupledNucVal{i}<=nucValBinSort(j+1));
        coupledC2ValBin{i}{j} = coupledC2Val{i}(coupledNucVal{i}>nucValBinSort(j) & coupledNucVal{i}<=nucValBinSort(j+1));
        uncoupledC2ValBin{i}{j} = uncoupledC2Val{i}(uncoupledNucVal{i}>nucValBinSort(j) & uncoupledNucVal{i}<=nucValBinSort(j+1));
    end
    
    coupledNucValMean{i} = cellfun(@mean, coupledNucValBin{i});    
    coupledDistMean{i} = cellfun(@mean, coupledDistBin{i});
    coupledValMean{i} = cellfun(@mean, coupledValBin{i});
    coupledVolMean{i} = cellfun(@mean, coupledVolBin{i});
    coupledDiaMean{i} = cellfun(@mean, coupledDiaBin{i});
    coupledC2ValMean{i} = cellfun(@mean, coupledC2ValBin{i});       
    
    coupledNucValStd{i} = cellfun(@std, coupledNucValBin{i});
    coupledDistStd{i} = cellfun(@std, coupledDistBin{i});
    coupledValStd{i} = cellfun(@std, coupledValBin{i});
    coupledVolStd{i} = cellfun(@std, coupledVolBin{i});
    coupledDiaStd{i} = cellfun(@std, coupledDiaBin{i});
    coupledC2ValStd{i} = cellfun(@std, coupledC2ValBin{i});
    
    markerType(1) = 'o';
    markerType(2) = 'v';



[pl1{1}, pl2{1}] = errPlot32(coupledNucValBin{i}, coupledTotValBin{i}, dataCutOff, colorStruct{1}(1,:), markerType(i), f1, f2);
hold on;
[pl1{2}, pl2{2}] = errPlot32(coupledNucValBin{i}, coupledValBin{i}, dataCutOff, colorStruct{2}(1,:), markerType(i), f1, f2);
hold on;
[pl1{4}, pl2{4}] = errPlot32(coupledNucValBin{i}, coupledDiaBin{i}, dataCutOff, colorStruct{4}(1,:), markerType(i), f1, f2);
hold on;
[pl1{5}, pl2{5}] = errPlot32(coupledNucValBin{i}, coupledC2ValBin{i}, dataCutOff, colorStruct{6}(1,:), markerType(i), f1, f2);

% [pl1{1}, pl2{1}] = errPlot32(coupledNucValBin{i}, coupledAllTotValBin{i}, dataCutOff, colorStruct{1}(1,:), markerType(i),f1, f2);
% hold on;
% [pl1{2}, pl2{2}] = errPlot32(coupledNucValBin{i}, coupledAllValBin{i}, dataCutOff, colorStruct{2}(1,:), markerType(i),f1, f2);
% hold on;
% [pl1{4}, pl2{4}] = errPlot32(coupledNucValBin{i}, coupledAllDiaBin{i}, dataCutOff, colorStruct{4}(1,:), markerType(i),f1, f2);
% hold on;
% [pl1{5}, pl2{5}] = errPlot32(coupledNucValBin{i}, coupledC2ValBin{i}, dataCutOff, colorStruct{6}(1,:), markerType(i),f1, f2);

legTex = {
    append("Cluster total intensity")...
    append("Cluster mean intensity"),...
    append("Cluster diameter")...
    append("MS2 mean intensity")...
    };

set (0, "CurrentFigure", f1)
leg =legend(horzcat(pl1{:}), legTex{:});
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
ylabel("{\sigma} (c)/c")
set(gca, 'YScale', 'log')
xlim([0 1])

set (0, "CurrentFigure", f2)
leg =legend(horzcat(pl2{:}), legTex{:});
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
ylabel("{\sigma}/ {\mu}")
xlim([0 1])
end

% title([append(geneName(1), markerType(1)); append(geneName(2), markerType(2))]);

end

function [valBinSort, idx, groupLenSort] = binInxFun(data, nBins, minLen)
[groupLen, valBin, idx] = histcounts(data, nBins);
valBin = (valBin(1:end-1) + diff(valBin) / 2)';
[valBinSort, sortIdx] = sort(valBin);

tempIdx = idx;
for i = 1:nBins
    idx(tempIdx==sortIdx(i)) = (i);    
end
groupLenSort = groupLen';
groupLenSort(sortIdx) = groupLenSort;
valBinSort(groupLenSort<minLen) = [];
end

function [pl1, pl2] = errPlot32(xCell, yCell, countLimit, color, markerType, f1, f2)
color = color(1,:)./255;
nBoot = 300;

xCell = xCell(~cellfun(@isempty, yCell));
yCell = yCell(~cellfun(@isempty, yCell));

yCell = cellfun(@(x) x(~isnan(x)), yCell, 'un', 0);
xCell = cellfun(@(x) x(~isnan(x)), xCell, 'un', 0);

maskAP = cellfun(@(x) (length(x)>countLimit), yCell);

yCellMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yCell(maskAP), 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yCellMean);
yCellStd = cellfun(@(x) bootstrp(nBoot,@std,x), yCell(maskAP), 'un', 0); % std of bins
yStdMean = cellfun(@mean, yCellStd);
yStdStd = cellfun(@std, yCellStd);
xCellMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xCell(maskAP), 'un', 0);
xMeanMean = cellfun(@mean, xCellMean);

f = fit(xMeanMean', yMeanMean', 'poly1');
% pf = plot(f);
% pf(1).Color = [color, 0.5];
% pf.LineStyle = '--';
% pf.LineWidth = 1;
% hold on;

coeff = coeffvalues(f);
slope = coeff(1);
% estX = (yMeanMean)./slope(1);
estXErr = (yStdMean)./abs(slope(1));
estXErrNorm = estXErr./xMeanMean;

% %-----------------------------------------
set(0, "CurrentFigure", f1)
pl1 = plot(((xMeanMean - min(xMeanMean))./min(xMeanMean))', estXErrNorm);
pl1.Marker = markerType;
pl1.LineWidth = 1;
pl1.Color = color;
pl1.LineStyle = '--';
hold on;
ax = gca;
ylabel('\sigma C');
xlabel('Nuclear concentration (c/c_{min})');

ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=400;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
% %-----------------------------------------

% %-----------------------------------------
set(0, "CurrentFigure", f2)
pl2 = errorbar(((xMeanMean - min(xMeanMean))./min(xMeanMean))', yStdMean./yMeanMean, yStdStd./yMeanMean);
pl2.Marker = markerType;
pl2.LineWidth = 1;
pl2.Color = color;
pl2.LineStyle = '--';
pl2.CapSize = 0; 
hold on;
ax = gca;
ylabel('\sigma C');
xlabel('Nuclear concentration (c/c_{min})');

ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=400;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end



function [pl1, pl2] = errPlot33(xCell, yCell, countLimit, color, markerType, f1, f2)
color = color(1,:)./255;
nBoot = 300;

% xCell = xCell(~cellfun(@isempty, yCell));
% yCell = yCell(~cellfun(@isempty, yCell));
% 
% yCell = cellfun(@(x) x(~isnan(x)), yCell, 'un', 0);
% xCell = cellfun(@(x) x(~isnan(x)), xCell, 'un', 0);
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xCell = xCell(~cellfun(@isempty, yCell));
yCell = yCell(~cellfun(@isempty, yCell));
lenCell = cellfun(@length, yCell);
yCellMin = mean(yCell{min(find(lenCell>countLimit))}, 'omitnan');
yCell = cellfun(@(x) (x-yCellMin)./yCellMin, yCell, 'un', 0);
% yCell = cellfun(@(x) log(x), yCell, 'un', 0);

xCellMin = mean(xCell{min(find(lenCell>countLimit))}, 'omitnan');
xCell = cellfun(@(x) (x-xCellMin)./xCellMin, xCell, 'un', 0);
% xCell = cellfun(@(x) log(x), xCell, 'un', 0);

% anovaFun(xCell, yCell);
% boxDataCell = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], yCell, xCell, 'un', 0); % adding labels to the cells
% boxData = vertcat(boxDataCell{:});


meanY = cellfun(@mean, yCell, 'un', 0);

% errY = cellfun(@(x) std(x)/sqrt(length(x)), yCell, 'un', 0);
errY = cellfun(@(x) std(x), yCell, 'un', 0);
%~~~~Don't change the order of these lines~~~~
errY(cellfun(@(x) any(isnan(x)), meanY)) = [];
xCellNan = xCell;
xCellNan(cellfun(@(x) any(isnan(x)), meanY)) = [];
meanX = cellfun(@mean, xCellNan, 'un', 0);
meanY(cellfun(@(x) any(isnan(x)), meanY)) = [];

xMeanMean = vertcat(meanX{:});
yMeanMean = vertcat(meanY{:});
yStdMean = vertcat(errY{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = fit(xMeanMean', yMeanMean', 'poly1');
f = fit(xMeanMean, yMeanMean, 'poly1');
% pf = plot(f);
% pf(1).Color = [color, 0.5];
% pf.LineStyle = '--';
% pf.LineWidth = 1;
% hold on;

coeff = coeffvalues(f);
slope = coeff(1);
% estX = (yMeanMean)./slope(1);
estXErr = (yStdMean)./abs(slope(1));
estXErrNorm = estXErr./xMeanMean;

% %-----------------------------------------
set(0, "CurrentFigure", f1)
pl1 = plot(((xMeanMean - min(xMeanMean))./min(xMeanMean))', estXErrNorm);
pl1.Marker = markerType;
pl1.LineWidth = 1;
pl1.Color = color;
pl1.LineStyle = '--';
hold on;
ax = gca;
ylabel('\sigma C');
xlabel('Nuclear concentration (c/c_{min})');

ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=400;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
% %-----------------------------------------

% %-----------------------------------------

maskAP = cellfun(@(x) (length(x)>countLimit), yCell);

yCellMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yCell(maskAP), 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yCellMean);
yCellStd = cellfun(@(x) bootstrp(nBoot,@std,x), yCell(maskAP), 'un', 0); % std of bins
yStdMean = cellfun(@mean, yCellStd);
yStdStd = cellfun(@std, yCellStd);
xCellMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xCell(maskAP), 'un', 0);
xMeanMean = cellfun(@mean, xCellMean);

set(0, "CurrentFigure", f2)
pl2 = errorbar(((xMeanMean - min(xMeanMean))./min(xMeanMean)), yStdMean./yMeanMean, yStdStd./yMeanMean);
pl2.Marker = markerType;
pl2.LineWidth = 1;
pl2.Color = color;
pl2.LineStyle = '--';
pl2.CapSize = 0; 
hold on;
ax = gca;
ylabel('\sigma C');
xlabel('Nuclear concentration (c/c_{min})');

ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=400;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end

