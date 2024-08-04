function emC1Plotter1(folderPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Works for single color multiposition data generated from CZI file format.
%   Can combine data from 2xa, 6x etc.
%   Plots various properties based on position bins 
%   Input#1: folder that contains all the combinedDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------
% Anterior nucleus =  30 molecules/um^3                     
% Anterior mean intensity = 150 (a.u.)                            
% Intensity to molecule conversion 30/150 = 1/5
% multiply all intensities by 1/5 = 0.2
%-------------------------------------------------------------------
factor = 0.2; % default = 0.2 from thomas (explanation above)
bins = 5; % use 13
minLen = 30;
nBoot = 100;

nucMolsAt20pc = 60000; %28450; %14225; % bcd 2x %  for 6x 
% nucMolsAt20pc = 100000; % bcd 6x

posCoord = [0.3000 0.2000 0.6000 0.8000];

posCutOff = [0.15, 0.7];
valCutOff = [0, 600];


cd(folderPath);
files=dir('combinedEmDSNew*');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 

colorStruct{1} = [185, 0, 91; 0, 0, 0];
colorStruct{1} = [10, 10, 10; 0, 0, 0];
colorStruct{2} = [59, 68, 246;0, 0, 0];

colorStruct{5} = [138, 62, 8; 0 0 0];
colorStruct{4} = [213,62,79; 0, 0, 0];
colorStruct{3} = [153,112,171; 0, 0, 0];
colorStruct{2} = [102,194,165; 0, 0, 0];
colorStruct{1} = [50,136,189; 0, 0, 0];

color = colorStruct(1:length(struct2C));

% names = {'Bicoid 2 Copy', 'Bicoid 6 copy'};
names = cellfun(@(x) x.emSpotProp.name, struct2C, 'un', 0);

nucPos = cellfun(@(x) x.emSpotProp.position, struct2C, 'un', 0);
nucPos = cellfun(@(x) vertcat(x{:}), nucPos{1}, 'un', 0); % each cell is an embryo
nucPosAll = vertcat(nucPos{:});

nucPosRep = cellfun(@(x) x.emSpotProp.positionRep, struct2C, 'un', 0);
nucPosRep = cellfun(@(x) vertcat(x{:}), nucPosRep{1}, 'un', 0); % each cell is an embryo
nucPosRepAll = vertcat(nucPosRep{:});

nucVal = cellfun(@(x) x.emSpotProp.nucValMod, struct2C, 'un', 0); % each cell is an embryo
nucVal = nucVal{1};
% nucVal = cellfun(@(x) vertcat(x{:}), nucVal, 'un', 0); % use if nuc val is used not nuc val mod
nucValAll = vertcat(nucVal{:}); % use if nuc val is used not nuc val mod
% nucValAll = horzcat(nucVal{:});
% nucValAll = horzcat(nucValAll{:});
% nucValAll = nucValAll';

nucValRep = cellfun(@(x) x.emSpotProp.nucValRep, struct2C, 'un', 0);
nucValRep = cellfun(@(x) vertcat(x{:}), nucValRep{1}, 'un', 0); % each cell is an embryo
nucValRepAll = vertcat(nucValRep{:});

%%%%%%%%%%%%%%%% Nuc val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNucVal = [165,42,42];
fNucVal = figure('Color', 'w');
set(0, "CurrentFigure", fNucVal)
hold on;
% for i=2:length(nucVal)
%     [hNucVal{i}, lambDaNucVal(i,:), ~] = scatterFitPlot(log(nucVal{i}), nucPos{i}, posCutOff, colorNucVal, 'lineNoFit');
% end
benPlot(nucPosAll, log(nucValAll))
% [N,c]  = hist3([nucPosAll, log(nucValAll)], 'Nbins', [30, 30]);
% N(N==0) = NaN;
% [x1, x2] = meshgrid(sort(c{1}),sort(c{2}));
% scatter(x1(:), x2(:), 20, N(:), 'filled');
% % colorbar
% colormap sky
hold on;
[hNucValAll, ~, ~] = scatterFitPlot(log(nucValAll), nucPosAll, posCutOff, 	colorNucVal, 'errNoFit');
hold on;
[hNucValAll, lambdaNucValAll, ~] = scatterFitPlot(log(nucValAll), nucPosAll, posCutOff,	colorNucVal, 'justFit');
xlim([0.15 0.65]);
ylim([3 6]);
ylim([1 7]);
ylabel('log(I_{nuc})');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaNucValAll(1), lambdaNucValAll(2));
% T = text(0.2, 5.5, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit total value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
colorTotVal = [51, 133, 141];
spotTotVal = cellfun(@(x) x.emSpotProp.spotTotValFilt, struct2C, 'un', 0);
for i=1:length(spotTotVal{1})
    spotTotVal{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotTotVal{1}{i})';
end
spotTotVal = spotTotVal{1};
spotTotValAll = vertcat(spotTotVal{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotTotValRep = cellfun(@(x) x.emSpotProp.spotTotValFilt, struct2C, 'un', 0);
spotTotValRep = spotTotValRep{1};
spotTotValRep = cellfun(@(x) vertcat(x{:}), spotTotValRep, 'un', 0);
spotTotValRepAll = vertcat(spotTotValRep{:});
% nucValEdges = [50 100 150 200 250];
% [xBin, yBin] = binMaker(nucValRepAll, spotTotValRepAll, nucValEdges, minLen);
% kernelPlotPosBin2(arrayfun(@(x) x, nucValEdges(1:end-1), 'un', 0), yBin, "I_{clust.}", '2x', colorTotVal);
% xlim([0 100])
% xlabel('I_{tot.}');
% ylabel('I_{nuc}');
% box on

% [xBin, yBin] = binMaker(nucPosAll, spotTotValAll, bins, minLen);
% kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "I_{clust.}", '2x', colorTotVal);


vTotVal = figure('Color', 'w');
set(0, "CurrentFigure", vTotVal)
ax = gca;
ax.Position = posCoord;
benPlot(nucValAll, spotTotValAll);
hold on;

% % cMapNucTotVal = colo
colorTotVal = [0 0 0];
[vhFrac, ~, ~] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'errNoFit');
hold on;
[~, ~, rSqSpotTotValAll, slopeSpotTotValAll] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'justFit');
xlim([50 200]);
ylabel('Cluster Total Value');
xlabel('I_{nuc}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pSpotTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", pSpotTotVal)
set(0, "CurrentFigure", fNucVal)
ax = gca;
ax.Position = posCoord;
hold on;
% for i=2:length(spotTotVal)
%     [phSpotTotVal{i}, lambdaSpotTotVal(i,:)] = scatterFitPlot(log(spotTotVal{i}), nucPos{i}, posCutOff, colorTotVal, 'lineNoFit');
% end
hold on;

% [N,c]  = hist3([nucPosAll, log(spotTotValAll)], 'Nbins', [30, 30]);
% N(N==0) = NaN;
% [x1, x2] = meshgrid(sort(c{1}),sort(c{2}));
% scatter(x1(:), x2(:), 20, N(:), 'filled');
% % colorbar
% colormap sky
% hold on;

benPlot(nucPosAll, log(spotTotValAll));
hold on;

[phSpotTotValAll, ~, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff,[0 0 0], 'errNoFit');
hold on;
[~, lambdaSpotTotValAll, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff, [0 0 0], 'justFit');
xlim([0.15 0.65]);
ylim([1 5])
ylim([1 7])
ylabel('log(I_{tot.})');
xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotTotValAll(1), lambdaSpotTotValAll(2));
% T = text(0.2, 2, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');


%%%%%%%%%%%%%%%% Nuc sum val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p
% NucTotVal = figure('Color', 'w');
nucTotVal = cellfun(@(x) x.emSpotProp.nucTotVal, struct2C, 'un', 0); % each cell is an embryo
nucTotVal = nucTotVal{1}';
nucTotVal = cellfun(@(x) cell2mat(x'), nucTotVal, 'un', 0);
nucTotValAll = vertcat(nucTotVal{:});
% 
% set(0, "CurrentFigure", pNucTotVal)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(nucTotVal)
%     [phNucTotVal{i}, lambdaNucTotVal(i,:), ~] = scatterFitPlot(log(nucTotVal{i}), nucPos{i}, posCutOff, [70 70 70], 'lineNoFit');
% end
% % hold on;
% [phNucTotValAll, ~, ~] = scatterFitPlot(log(nucTotValAll), nucPosAll, posCutOff, [70 70 70], 'errNoFit');
% hold on;
% [~, lambdaNucTotValAll, ~] = scatterFitPlot(log(nucTotValAll), nucPosAll, posCutOff, [70 70 70], 'justFit');
% xlim([0.15 0.65]);
% ylim([12 20])
% ylabel('log(\Sigma I_{nuc.})');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaNucTotValAll(1), lambdaNucTotValAll(2));
% T = text(0.2, 20, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%
% vNucTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", vNucTotVal)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(nucTotVal)
%     [vhNucTotVal{i}, ~, rSqNucTotVal(i,:)] = scatterFitPlot(nucTotVal{i}, nucVal{i}, valCutOff, color1(i,1:3), 'fitOff');
% end
% hold on;
% [vhNucTotValAll, ~, rSqNucTotValAll] = scatterFitPlot(nucTotValAll, nucValAll, valCutOff, [70 70 70], 'fitOn');
% xlim([50 200]);
% ylabel('\Sigma I_{nuc.}');
% xlabel('I_{nuc}');

%%%%%%%%%%%%%%%% Spot sum val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotTotVal = cellfun(@(x) x.emSpotProp.spotRawTotValSum, struct2C, 'un', 0); % each cell is an embryo
spotTotVal = spotTotVal{1}';
spotTotVal = cellfun(@(x) cell2mat(x'), spotTotVal, 'un', 0);
spotTotValAll = vertcat(spotTotVal{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSpotTotVal = figure('Color', 'w');
set(0, "CurrentFigure", pSpotTotVal)
ax = gca;
ax.Position = posCoord;
hold on;
for i=2:length(spotTotVal)
    [phSpotTotVal{i}, lambdaSpotTotVal(i,:)] = scatterFitPlot(log(spotTotVal{i}), nucPos{i}, posCutOff, [0, 112, 187], 'lineNoFit');
end
hold on;
[phSpotTotValAll, ~, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff, [0, 112, 187], 'errNoFit');
hold on;
[~, lambdaSpotTotValAll, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff, [0, 112, 187], 'justFit');
xlim([0.15 0.65]);
ylim([12 16])
ylabel('log(\Sigma I_{clust.})');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotTotValAll(1), lambdaSpotTotValAll(2));
% T = text(0.2, 15, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vSpotTotVal = figure('Color', 'w');
set(0, "CurrentFigure", vSpotTotVal)
ax = gca;
ax.Position = posCoord;
% cMapNucTotVal = colormap(summer(length(nucTotVal)));
hold on;
for i=2:length(spotTotVal)
    [vhSpotTotVal{i}, ~, rSqSpotTotVal(i,:)] = scatterFitPlot(spotTotVal{i}, nucVal{i}, valCutOff, [0, 112, 187], 'lineNoFit');
end
hold on;
[vhSpotTotValAll, ~, ~] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, [0, 112, 187], 'errNoFit');
hold on;
[~, ~, rSqSpotTotValAll] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, [0, 112, 187], 'justFit');
xlim([50 200]);
ylabel('\Sigma I_{clust.}');
xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotTotValAll);
% T = text(75, 0.5*10^6, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%% Bg val sum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bgValSum = cellfun(@(x) x.emSpotProp.bgRawValSum, struct2C, 'un', 0); % each cell is an embryo
bgValSum = bgValSum{1}';
bgValSum = cellfun(@(x) cell2mat(x'), bgValSum, 'un', 0);
bgValSumAll = vertcat(bgValSum{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pBgTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", pBgTotVal)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(bgValSum)
%     [hBgValSum{i}, lambdaBgValSum(i,:)] = scatterFitPlot(log(bgValSum{i}), nucPos{i}, posCutOff, [222, 49, 99], 'lineNoFit');
% end
% hold on;
% [phBgValSumAll, ~, ~] = scatterFitPlot(log(bgValSumAll), nucPosAll, posCutOff, [222, 49, 99], 'errNoFit');
% hold on;
% [~, lambdaBgValSumAll, ~] = scatterFitPlot(log(bgValSumAll), nucPosAll, posCutOff, [222, 49, 99], 'justFit');
% xlim([0.15 0.65]);
% ylim([12 20])
% ylabel('log(\Sigma I_{non-clust.})');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaBgValSumAll(1), lambdaBgValSumAll(2));
% T = text(0.2, 20, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vBgTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", vBgTotVal)
% ax = gca;
% ax.Position = posCoord;
% % cMapNucTotVal = colormap(summer(length(nucTotVal)));
% hold on;
% for i=2:length(bgValSum)
%     [hBgValSum{i}, ~, rSqBgTotVal(i,:)] = scatterFitPlot(bgValSum{i}, nucVal{i}, valCutOff, [222, 49, 99], 'lineNoFit');
% end
% hold on;
% [vhBgValSumAll, ~, ~] = scatterFitPlot(bgValSumAll, nucValAll, valCutOff, [222, 49, 99], 'errNoFit');
% hold on;
% [~, ~, rSqBgValSumAll] = scatterFitPlot(bgValSumAll, nucValAll, valCutOff, [222, 49, 99], 'justFit');
% xlim([50 200]);
% ylabel('\Sigma I_{non-clust.}');
% xlabel('I_{nuc}');
% 
% str=sprintf('R^2 = %1.2f', rSqBgValSumAll);
% T = text(75, 2*10^7, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%% clustered fraction intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawSpotFracInt = cellfun(@(x) x.emSpotProp.rawSpotFracInt, struct2C, 'un', 0); % each cell is an embryo
rawSpotFracInt = rawSpotFracInt{1}';
rawSpotFracInt = cellfun(@(x) cell2mat(x'), rawSpotFracInt, 'un', 0);
rawSpotFracIntAll = vertcat(rawSpotFracInt{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pRawSpotFracInt = figure('Color', 'w');
set(0, "CurrentFigure", pRawSpotFracInt)
ax = gca;
ax.Position = posCoord;
hold on;
for i=2:length(rawSpotFracInt)
    [hRawSpotFracInt{i}, lambdaRawSpotFracInt(i,:)] = scatterFitPlot(log(rawSpotFracInt{i}), nucPos{i}, posCutOff, [0, 112, 187], 'lineNoFit');
end
hold on;
[phRawSpotFracInt, ~, ~] = scatterFitPlot(log(rawSpotFracIntAll), nucPosAll, posCutOff, [0, 112, 187], 'errNoFit');
hold on;
[~, lambdaRawSpotFracIntAll, ~] = scatterFitPlot(log(rawSpotFracIntAll), nucPosAll, posCutOff, [0, 112, 187], 'justFit');
xlim([0.15 0.65]);
ylim([12 20])
ylabel('log(I_{clust.})');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaRawSpotFracIntAll(1), lambdaRawSpotFracIntAll(2));
% T = text(0.2, 20, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vRawSpotFracIntAll = figure('Color', 'w');
set(0, "CurrentFigure", vRawSpotFracIntAll)
ax = gca;
ax.Position = posCoord;
% cMapNucTotVal = colormap(summer(length(nucTotVal)));
hold on;
% for i=2:length(rawSpotFracInt)
%     [hRawSpotFracIntAll{i}, ~, rSqRawSpotFracIntAll(i,:)] = scatterFitPlot(rawSpotFracInt{i}, nucVal{i}, valCutOff, [0, 112, 187], 'lineNoFit');
% end
% hold on;
% [vhEnrichRatioAll, ~, ~] = scatterFitPlot(rawSpotFracIntAll, nucValAll, valCutOff, [70 70 70], 'errNoFit');
% hold on;
[~, ~, rSqRawSpotFracIntAll, slopeRawSpotFracIntAll] = scatterFitPlot(rawSpotFracIntAll, nucValAll, valCutOff, [70 70 70], 'justFit');
xlim([50 200]);
ylabel('I_{clust.}');
xlabel('I_{nuc}');

% str=sprintf('R^2 = %1.2f', rSqRawSpotFracIntAll);
% T = text(75, 2*10^7, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%% enrichRatio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% enrichRatio = cellfun(@(x) x.emSpotProp.enrichRatio, struct2C, 'un', 0); % each cell is an embryo
% enrichRatio = enrichRatio{1}';
% enrichRatio = cellfun(@(x) cell2mat(x'), enrichRatio, 'un', 0);
% enrichRatioAll = vertcat(enrichRatio{:});
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pEnrichRatio = figure('Color', 'w');
% set(0, "CurrentFigure", pEnrichRatio)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(enrichRatio)
%     [hEnrichRatio{i}, lambdaEnrichRatio(i,:)] = scatterFitPlot(log(enrichRatio{i}), nucPos{i}, posCutOff, [222, 49, 99], 'lineNoFit');
% end
% hold on;
% [phEnrichRatio, ~, ~] = scatterFitPlot(log(enrichRatioAll), nucPosAll, posCutOff, [222, 49, 99], 'errNoFit');
% hold on;
% [~, lambdaEnrichRatioAll, ~] = scatterFitPlot(log(enrichRatioAll), nucPosAll, posCutOff, [222, 49, 99], 'justFit');
% xlim([0.15 0.65]);
% ylim([12 20])
% ylabel('log(enrich)');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaEnrichRatioAll(1), lambdaEnrichRatioAll(2));
% T = text(0.2, 20, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vEnrichRatioAll = figure('Color', 'w');
% set(0, "CurrentFigure", vEnrichRatioAll)
% ax = gca;
% ax.Position = posCoord;
% % cMapNucTotVal = colormap(summer(length(nucTotVal)));
% hold on;
% for i=2:length(enrichRatio)
%     [hEnrichRatioAll{i}, ~, rSqEnrichRatioAll(i,:)] = scatterFitPlot(enrichRatio{i}, nucVal{i}, valCutOff, [222, 49, 99], 'lineNoFit');
% end
% hold on;
% [vhEnrichRatioAll, ~, ~] = scatterFitPlot(enrichRatioAll, nucValAll, valCutOff, [222, 49, 99], 'errNoFit');
% hold on;
% [~, ~, rSqEnrichRatioAll] = scatterFitPlot(enrichRatioAll, nucValAll, valCutOff, [222, 49, 99], 'justFit');
% xlim([50 200]);
% ylabel('Cluster enrichment factor');
% xlabel('I_{nuc}');
% 
% str=sprintf('R^2 = %1.2f', rSqEnrichRatioAll);
% T = text(75, 2*10^7, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%% cluster fraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorFrac = [87, 55, 43];
pFrac = figure('Color', 'w');
frac = cellfun(@(x, y) x./y, spotTotVal, nucTotVal, 'un', 0); % each cell is an embryo
fracAll = vertcat(frac{:});

set(0, "CurrentFigure", pFrac)
ax = gca;
ax.Position = posCoord;
cMapNucTotVal = colormap(summer(length(nucTotVal)));
hold on;
for i=2:length(frac)
    [phFrac{i}, lambdaFrac(i,:)] = scatterFitPlot((frac{i}), nucPos{i}, posCutOff, colorFrac, 'lineNoFit');
end
hold on;
[phFrac, ~, ~] = scatterFitPlot((fracAll), nucPosAll, posCutOff, colorFrac, 'errNoFit');
xlim([0.15 0.65]);
ylabel('Cluster fraction');
xlabel('x/L');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vFrac = figure('Color', 'w');
set(0, "CurrentFigure", vFrac)
ax = gca;
ax.Position = posCoord;
% % cMapNucTotVal = colormap(summer(length(nucTotVal)));
% hold on;
% for i=2:length(frac)
%     scatterFitPlot(frac{i}, nucVal{i}, valCutOff, colorFrac, 'lineNoFit');
% end
% hold on;
colorFrac = [0 0 0];
[vhFrac, ~, ~] = scatterFitPlot(fracAll, nucValAll, valCutOff, colorFrac, 'errNoFit');
hold on;
[vhFrac, ~, ~] = scatterFitPlot(fracAll, nucValAll, valCutOff, colorFrac, 'justFit');
xlim([50 200]);
ylim([0.025 0.05])
ylabel('Cluster fraction');
xlabel('I_{nuc}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%~Error Plots~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('color', 'w');
f2 = figure('color', 'w');
hold on;
plErrNuc = binErrPlot(nucValAll, nucPosAll, posCutOff, [165,42,42], bins, minLen, nBoot, f1, f2, 'pos', lambdaNucValAll);
hold on;
plErrRawSpot = binErrPlot(rawSpotFracIntAll, nucPosAll, posCutOff, [0, 112, 187], bins, minLen, nBoot, f1, f2, 'pos', lambdaRawSpotFracIntAll);
% hold on;
% binErrPlot(bgValSumAll, nucPosAll, posCutOff, [222, 49, 99], bins, minLen, nBoot, f1, f2, 'pos', lambdaBgValSumAll);

hh = [plErrNuc; plErrRawSpot];

set(0, "CurrentFigure", f1)
legend(hh, 'Nuclear intensity', 'Clustered fraction intensity')
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.4])
ylabel('\sigma / \mu')
xlabel('x/L');
set(0, "CurrentFigure", f2)
xlim([0.15 0.65])
ylim([0 0.1])
ylabel('\sigma_{x/L}')
xlabel('x/L');

f1 = figure('color', 'w');
f2 = figure('color', 'w');
hold on;
binErrPlot(nucValAll, nucValAll, valCutOff, [165,42,42], bins, minLen, nBoot, f1, f2, 'conc', lambdaNucValAll);
hold on;
binErrPlot(rawSpotFracIntAll, nucValAll, valCutOff, [0, 112, 187], bins, minLen, nBoot, f1, f2, 'conc', lambdaRawSpotFracIntAll);
% hold on;
% binErrPlot(bgValSumAll, nucValAll, valCutOff, [222, 49, 99], bins, minLen, nBoot, f1, f2, 'conc', lambdaBgValSumAll);

set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([50 200])
ylim([0 0.4])
ylabel('\sigma / \mu')
xlabel('I_{nuc}');
set(0, "CurrentFigure", f2)
xlim([50 200])
ylim([0 0.1])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');

%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cluster density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorDen = [44,87,117];
spotPerNuc = cellfun(@(x) x.emSpotProp.spotCount, struct2C, 'un', 0);
spotPerNuc = cellfun(@(x) vertcat(x{:}), spotPerNuc{1}, 'un', 0); % each cell is an embryo
spotPerNucAll = vertcat(spotPerNuc{:});


nucValEdges = [50 100 150 200 250];
[xBin, yBin] = binMaker(nucValAll, spotPerNucAll, nucValEdges, minLen);
kernelPlotPosBin2(arrayfun(@(x) x, nucValEdges(1:end-1), 'un', 0), yBin, "#", '2x', colorDen);
xlim([20 100])

pSpotPerNuc = figure('Color', 'w');
benPlot(nucValAll, spotPerNucAll);
hold on;
[~, ~, rSqSpotPerNucAll, ~] = scatterFitPlot((spotPerNucAll), nucValAll, valCutOff, colorDen, 'justFit');
set(0, "CurrentFigure", pSpotPerNuc)
hold on;
% for i=2:length(spotPerNuc)
%     [phSpotPerNuc{i}, ~] = scatterFitPlot((spotPerNuc{i}), nucPos{i}, posCutOff, colorDen, 'lineNoFit');
% end
% hold on;
[phSpotPerNucAll, ~, ~] = scatterFitPlot(spotPerNucAll, nucPosAll, posCutOff, colorDen, 'errNoFit');
xlim([0.15 0.65]);
ylim([20 100])
ylabel('Cluster count');
xlabel('x/L');

[~, lambdaSpotPerNucAll, ~] = scatterFitPlot(log(spotPerNucAll), nucPosAll, posCutOff, colorDen, 'justFit');
f1 = figure ('color', 'w');
f2 = figure ('color', 'w');
binErrPlot(spotPerNucAll, nucPosAll, posCutOff, colorDen, bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotPerNucAll);
set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.4])
ylabel('\sigma / \mu')
xlabel('I_{nuc}');
set(0, "CurrentFigure", f2)
xlim([50 200])
ylim([0 0.1])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vSpotPerNuc = figure('Color', 'w');
set(0, "CurrentFigure", vSpotPerNuc)
hold on;
for i=2:length(spotPerNuc)
    [phSpotPerNuc{i}, ~] = scatterFitPlot((spotPerNuc{i}), nucVal{i}, valCutOff, colorDen, 'lineNoFit');
end
hold on;
[vhSpotPerNucAll, ~, ~] = scatterFitPlot(spotPerNucAll, nucValAll, valCutOff, colorDen, 'errNoFit');
xlim([50 200]);
ylim([20 100])
ylabel('Cluster count');
xlabel('I_{nuc}');
hold on;
[~, ~,rSqSpotPerNucAll, slopeSpotPerNucAll] = scatterFitPlot((spotPerNucAll), nucValAll, valCutOff, colorDen, 'justFit');
f1 = figure ('color', 'w');
f2 = figure ('color', 'w');
binErrPlot(spotPerNucAll, nucValAll, valCutOff, colorDen, bins, minLen, nBoot, f1, f2, 'val', slopeSpotPerNucAll);
set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.4])
ylabel('\sigma / \mu')
xlabel('I_{nuc}');
set(0, "CurrentFigure", f2)
xlim([50 200])
ylim([0 0.1])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% raw dia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% spotRawDia = cellfun(@(x) x.emSpotProp.spotRawDiaFilt, struct2C, 'un', 0);
% for i=1:length(spotRawDia{1})
%     spotRawDia{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotRawDia{1}{i})';
% end
% spotRawDia = spotRawDia{1};
% spotRawDiaAll = vertcat(spotRawDia{:});
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [xBin, yBin] = binMaker(nucPosAll, spotRawDiaAll, bins, minLen);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "dia ({\mu}m)", '2x', [84,39,136]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pSpotRawDia = figure('Color', 'w');
% set(0, "CurrentFigure", pSpotRawDia)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotRawDia)
%     [phSpotRawDia{i}, lambdaSpotRawDia(i,:)] = scatterFitPlot(log(spotRawDia{i}), nucPos{i}, posCutOff, [84,39,136], 'lineNoFit');
% end
% hold on;
% [phSpotRawDiaAll, ~, ~] = scatterFitPlot(log(spotRawDiaAll), nucPosAll, posCutOff, [84,39,136], 'errNoFit');
% hold on;
% [~, lambdaSpotRawDiaAll, ~] = scatterFitPlot(log(spotRawDiaAll), nucPosAll, posCutOff, [84,39,136], 'justFit');
% xlim([0.15 0.65]);
% ylim([-1.5 -1])
% ylabel('log(\Phi_{clust.})');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotRawDiaAll(1), lambdaSpotRawDiaAll(2));
% T = text(0.2, -1.4, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vSpotRawDia = figure('Color', 'w');
% set(0, "CurrentFigure", vSpotRawDia)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotRawDia)
%     [hSpotRawDia{i}, ~, rSqSpotRawDia(i,:)] = scatterFitPlot(spotRawDia{i}, nucVal{i}, valCutOff, [84,39,136], 'lineNoFit');
% end
% hold on;
% [vhSpotRawDiaAll, ~, ~] = scatterFitPlot(spotRawDiaAll, nucValAll, valCutOff, [84,39,136], 'errNoFit');
% hold on;
% [~, ~, rSqSpotRawDiaAll] = scatterFitPlot(spotRawDiaAll, nucValAll, valCutOff, [84,39,136], 'justFit');
% xlim([50 200]);
% ylim([0.25 0.35])
% ylabel('\Sigma \Phi_{clust.}');
% xlabel('I_{nuc}');
% 
% str=sprintf('R^2 = %1.2f', rSqSpotRawDiaAll);
% T = text(75, 0.27, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mean int %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spotMeanInt = cellfun(@(x) x.emSpotProp.spotMeanIntFilt, struct2C, 'un', 0);
% for i=1:length(spotMeanInt{1})
%     spotMeanInt{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotMeanInt{1}{i})';
% end
% spotMeanInt = spotMeanInt{1};
% spotMeanIntAll = vertcat(spotMeanInt{:});
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [xBin, yBin] = binMaker(nucPosAll, spotMeanIntAll, bins, minLen);
% kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "I_{clust.}", '2x', [216,179,101]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pSpotMeanInt = figure('Color', 'w');
% set(0, "CurrentFigure", pSpotMeanInt)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotMeanInt)
%     [phSpotMeanInt{i}, lambdaSpotMeanInt(i,:)] = scatterFitPlot(log(spotMeanInt{i}), nucPos{i}, posCutOff, [216,179,101], 'lineNoFit');
% end
% hold on;
% [phSpotMeanIntAll, ~, ~] = scatterFitPlot(log(spotMeanIntAll), nucPosAll, posCutOff, [216,179,101], 'errNoFit');
% hold on;
% [~, lambdaSpotMeanIntAll, ~] = scatterFitPlot(log(spotMeanIntAll), nucPosAll, posCutOff, [216,179,101], 'justFit');
% xlim([0.15 0.65]);
% ylim([3 6])
% ylabel('log(I_{clust.})');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotMeanIntAll(1), lambdaSpotMeanIntAll(2));
% T = text(0.2, 3.5, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% vSpotMeanInt = figure('Color', 'w');
% set(0, "CurrentFigure", vSpotMeanInt)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotMeanInt)
%     [hSpotMeanInt{i}, ~, rSqSpotMeanInt(i,:)] = scatterFitPlot(spotMeanInt{i}, nucVal{i}, valCutOff, [216,179,101], 'lineNoFit');
% end
% hold on;
% [vhSpotMeanIntAll, ~, ~] = scatterFitPlot(spotMeanIntAll, nucValAll, valCutOff, [216,179,101], 'errNoFit');
% hold on;
% [~, ~, rSqSpotMeanIntAll] = scatterFitPlot(spotMeanIntAll, nucValAll, valCutOff, [216,179,101], 'justFit');
% xlim([50 200]);
% ylim([0 400])
% ylabel('\Sigma \Phi_{clust.}');
% xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotMeanIntAll);
% T = text(75, 100, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% raw total value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% spotRawTotVal = cellfun(@(x) x.emSpotProp.spotRawTotVal, struct2C, 'un', 0);
% for i=1:length(spotRawTotVal{1})
%     spotRawTotVal{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotRawTotVal{1}{i})';
% end
% spotRawTotVal = spotRawTotVal{1};
% spotRawTotValAll = vertcat(spotRawTotVal{:});
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [xBin, yBin] = binMaker(nucPosAll, spotRawTotValAll, bins, minLen);
% kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "I_{clust.}", '2x', [90,180,172]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pSpotRawTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", pSpotRawTotVal)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotRawTotVal)
%     [phSpotRawTotVal{i}, lambdaSpotRawTotVal(i,:)] = scatterFitPlot(log(spotRawTotVal{i}), nucPos{i}, posCutOff, [90,180,172], 'lineNoFit');
% end
% hold on;
% [phSpotRawTotValAll, ~, ~] = scatterFitPlot(log(spotRawTotValAll), nucPosAll, posCutOff, [90,180,172], 'errNoFit');
% hold on;
% [~, lambdaSpotRawTotValAll, ~] = scatterFitPlot(log(spotRawTotValAll), nucPosAll, posCutOff, [90,180,172], 'justFit');
% xlim([0.15 0.65]);
% ylim([7 11])
% ylabel('log(I_{tot_{clust.}})');
% xlabel('x/L');
% 
% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotRawTotValAll(1), lambdaSpotRawTotValAll(2));
% T = text(0.2, 8, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% vSpotRawTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", vSpotRawTotVal)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotRawTotVal)
%     [hSpotRawTotVal{i}, ~, rSqSpotRawTotVal(i,:)] = scatterFitPlot(spotRawTotVal{i}, nucVal{i}, valCutOff, [90,180,172], 'lineNoFit');
% end
% hold on;
% [vhSpotRawTotValAll, ~, ~] = scatterFitPlot(spotRawTotValAll, nucValAll, valCutOff, [90,180,172], 'errNoFit');
% hold on;
% [~, ~, rSqSpotRawTotValAll] = scatterFitPlot(spotRawTotValAll, nucValAll, valCutOff, [90,180,172], 'justFit');
% xlim([50 200]);
% ylim([0 30000])
% ylabel('I_{tot_{clust.}}');
% xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotRawTotValAll);
% T = text(75, 7000, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit dia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorDia = [163, 146, 162];

spotFitDia = cellfun(@(x) x.emSpotProp.spotFitDia, struct2C, 'un', 0);
for i=1:length(spotFitDia{1})
    spotFitDiaMean{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitDia{1}{i})';
    spotFitDiaStd{i} = cellfun(@(x) std(x, 'omitnan'), spotFitDia{1}{i})';
    spotFitDiaVar{i} = cellfun(@(x) var(x, 'omitnan'), spotFitDia{1}{i})';
end
spotFitDia = spotFitDiaMean;
spotFitDiaAll = vertcat(spotFitDia{:});
spotFitDiaStd = spotFitDiaStd;
spotFitDiaStdAll = vertcat(spotFitDiaStd{:});
spotFitDiaVar = spotFitDiaVar;
spotFitDiaVarAll = vertcat(spotFitDiaVar{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spotFitDiaRep = cellfun(@(x) x.emSpotProp.spotFitDiaFilt, struct2C, 'un', 0);
spotFitDiaRep = spotFitDiaRep{1};
spotFitDiaRep = cellfun(@(x) vertcat(x{:}), spotFitDiaRep, 'un', 0);
spotFitDiaRepAll = vertcat(spotFitDiaRep{:});
[xBin, yBin] = binMaker(nucPosAll, spotFitDiaAll, bins, minLen);
kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "\D (\mu m)", '2x', colorDia);
nucValEdges = [0 50 100 150 200 250];
[xBin, yBin] = binMaker(nucValAll, spotFitDiaVarAll, nucValEdges, minLen);
% [xBin, yBin] = binMaker(nucValRepAll, spotFitDiaRepAll, nucValEdges, minLen);
kernelPlotPosBin2(arrayfun(@(x) x, nucValEdges(1:end-1), 'un', 0), yBin, "D", '2x', colorDia);
xlim([0 0.8])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSpotFitDia = figure('Color', 'w');
colorDia     = [0 0 0 ];
set(0, "CurrentFigure", pSpotFitDia)
ax = gca;
ax.Position = posCoord;
benPlot(nucPosAll, log(spotFitDiaAll));
hold on;
% for i=2:length(spotFitDia)
%     [phSpotFitDia{i}, lambdaSpotFitDia(i,:)] = scatterFitPlot(log(spotFitDia{i}), nucPos{i}, posCutOff, colorDia, 'lineNoFit');
% end
hold on;
[phSpotFitDiaAll, ~, ~] = scatterFitPlot(log(spotFitDiaAll), nucPosAll, posCutOff, colorDia, 'errNoFit');
hold on;
[~, lambdaSpotFitDiaAll, ~] = scatterFitPlot(log(spotFitDiaAll), nucPosAll, posCutOff, colorDia, 'justFit');
xlim([0.15 0.65]);
ylim([-1.5 -1])
ylabel('\Phi');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotFitDiaAll(1), lambdaSpotFitDiaAll(2));
% T = text(0.2, -1.1, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vSpotFitDia = figure('Color', 'w');
set(0, "CurrentFigure", vSpotFitDia)
ax = gca;
ax.Position = posCoord;

hold on;
% for i=2:length(spotFitDia)
%     [hSpotFitDia{i}, ~, rSqSpotFitDia(i,:)] = scatterFitPlot(spotFitDia{i}, nucVal{i}, valCutOff, colorDia, 'lineNoFit');
% end

benPlot(nucValAll, spotFitDiaAll)
hold on;
[vhSpotFitDiaAll, ~, ~] = scatterFitPlot(spotFitDiaAll, nucValAll, valCutOff, [0 0 0], 'errNoFit');
hold on;
[~, ~, rSqSpotFitDiaAll, slopeSpotFitDiaAll] = scatterFitPlot(spotFitDiaAll, nucValAll, valCutOff, [0 0 0], 'justFit');
xlim([0 200]);
ylim([-inf inf])
ylabel('cluster Diameter');
xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotFitDiaAll);
% T = text(75, 0.35, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit Val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
colorVal = [143, 150, 80]; 
spotVal = cellfun(@(x) x.emSpotProp.spotValFilt, struct2C, 'un', 0);
for i=1:length(spotVal{1})
    spotValStd{1}{i} = cellfun(@(x) std(x, 'omitnan'), spotVal{1}{i})';
    spotVal{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotVal{1}{i})';
end
spotVal = spotVal{1};
spotValAll = vertcat(spotVal{:});
spotValStd = spotValStd{1};
spotValStdAll = vertcat(spotValStd{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotValRep = cellfun(@(x) x.emSpotProp.spotValFilt, struct2C, 'un', 0);
spotValRep = spotValRep{1};
spotValRep = cellfun(@(x) vertcat(x{:}), spotValRep, 'un', 0);
spotValRepAll = vertcat(spotValRep{:});
nucValEdges = [0 50 100 150 200 250];
[xBin, yBin] = binMaker(nucValRepAll, spotValRepAll, nucValEdges, minLen);
kernelPlotPosBin2(arrayfun(@(x) x, nucValEdges(1:end-1), 'un', 0), yBin, "I_{clust.}", '2x', [0 0 0]);
xlim([0 700])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSpotVal = figure('Color', 'w');
set(0, "CurrentFigure", pSpotVal)
% set(0, "CurrentFigure", fNucVal)
ax = gca;
ax.Position = posCoord;
benPlot(nucPosAll, log(spotValAll));
hold on;
% for i=2:length(spotVal)
%     [phSpotVal{i}, lambdaSpotVal(i,:)] = scatterFitPlot(log(spotVal{i}), nucPos{i}, posCutOff, colorVal, 'lineNoFit');
% end
% hold on;
[phSpotValAll, ~, ~] = scatterFitPlot(log(spotValAll), nucPosAll, posCutOff, colorVal, 'errNoFit');
hold on;
[~, lambdaSpotValAll, ~] = scatterFitPlot(log(spotValAll), nucPosAll, posCutOff, colorVal, 'justFit');
xlim([0.15 0.65]);
ylim([3 7])
ylim([1 7])
ylabel('log(I_{r})');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotValRepAll(1), lambdaSpotValRepAll(2));
% T = text(0.2, 3.7, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vSpotVal = figure('Color', 'w');
set(0, "CurrentFigure", vSpotVal)
ax = gca;
ax.Position = posCoord;
benPlot(nucValAll, spotValAll);
hold on;
% for i=2:length(spotVal)
%     [hSpotVal{i}, ~, rSqSpotVal(i,:)] = scatterFitPlot(spotVal{i}, nucVal{i}, valCutOff, colorVal, 'lineNoFit');
% end
% hold on;
[vhSpotValAll, ~, ~] = scatterFitPlot(spotValAll, nucValAll, valCutOff, [0 0 0], 'errNoFit');
hold on;
[~, ~, rSqSpotValAll, slopeSpotValAll] = scatterFitPlot(spotValAll, nucValAll, valCutOff, [0 0 0], 'justFit');
xlim([50 200]);
ylim([0 600])
ylabel('I_{clust.}');
xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotValAll);
% T = text(75, 50, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dia vs val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNow = [255,179,186];
colorNow = [0 0 0];
spotValCutOff = [0.00001, 1000];
spotDiaCutOff = [0.0001 1];% 
fsh = figure('color', 'w');
% sh = histScatPlot(spotValAll, spotFitDiaAll, colorNow);
% set(fsh, 'CurrentAxes', sh(1))
benPlot(spotValAll, spotFitDiaAll)
hold on;
hold on;
[vhSpotValAll, ~, ~] = scatterFitPlot(spotFitDiaAll, spotValAll, spotValCutOff, colorNow, 'errNoFit');
hold on;
[~, ~, rSq, slope] = scatterFitPlot(spotFitDiaAll, spotValAll, spotValCutOff, colorNow, 'justFit');
xlim([0 500])
ylim([0.2 0.6])
ylabel('\Phi (\mu m)');
xlabel('I_{p}');
% str=sprintf('slope = %1.2f', slope);
% T = text(75, 0.35, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%~ fit_bg~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorBg = [224, 187, 228];
spotBg = cellfun(@(x) x.emSpotProp.spotBgFilt, struct2C, 'un', 0); % each cell is an embryo

for i=1:length(spotBg{1})
    spotBg{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotBg{1}{i})';
end
spotBg = spotBg{1};
spotBgAll = vertcat(spotBg{:});
[xBin, yBin] = binMaker(nucPosAll, spotBgAll, bins, minLen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotBgRep = cellfun(@(x) x.emSpotProp.spotBgFilt, struct2C, 'un', 0);
spotBgRep = spotBgRep{1};
spotBgRep = cellfun(@(x) vertcat(x{:}), spotBgRep, 'un', 0);
spotBgRepAll = vertcat(spotBgRep{:});
[xBin, yBin] = binMaker(nucPosRepAll, spotBgRepAll, bins, minLen);
kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(cellfun(@mean, xBin), 2), 'un', 0), yBin, "I_{clust.}", '2x', colorBg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pBg = figure('Color', 'w');
set(0, "CurrentFigure", pBg)
ax = gca;
ax.Position = posCoord;
benPlot(nucPosAll, log(spotBgAll));
hold on;
for i=2:length(spotBg)
    [hSpotBg{i}, lambdaSpotBg(i,:)] = scatterFitPlot(log(spotBg{i}), nucPos{i}, posCutOff, colorBg, 'lineNoFit');
end
hold on;
[phspotBg, ~, ~] = scatterFitPlot(log(spotBgAll), nucPosAll, posCutOff, colorBg, 'errNoFit');
hold on;
[~, lambdaSpotBgAll, ~] = scatterFitPlot(log(spotBgAll), nucPosAll, posCutOff, colorBg, 'justFit');
xlim([0.15 0.65]);
ylim([3 6])
ylabel('log(I_{bg})');
xlabel('x/L');

% str=sprintf('\x03bb = %1.2f \x00B1 %1.2f', lambdaSpotBgAll(1), lambdaSpotBgAll(2));
% T = text(0.2, 3.7, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vSpotBg = figure('Color', 'w');
set(0, "CurrentFigure", vSpotBg)
ax = gca;
ax.Position = posCoord;
benPlot(nucValAll, (spotBgAll));
hold on;
% for i=2:length(spotBg)
%     [hSpotBg{i}, ~, rSqSpotBg(i,:)] = scatterFitPlot(spotBg{i}, nucVal{i}, valCutOff, colorBg, 'lineNoFit');
% end
% hold on;
[vhSpotBgAll, ~, ~] = scatterFitPlot(spotBgAll, nucValAll, valCutOff, [0 0 0], 'errNoFit');
hold on;
[~, ~, rSqSpotBgAll, slopeSpotBgAll] = scatterFitPlot(spotBgAll, nucValAll, valCutOff, [0 0 0], 'justFit');
xlim([50 200]);
ylim([0 250]);
ylabel('I_{bg}');
xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotBgAll);
% T = text(75, 50, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Peak vs bg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color', 'w');
histogram(spotValAll./ spotBgAll, 'DisplayStyle','stairs','Normalization','probability');
xlim([0 10]);
ylim([0 0.2])
xlabel('I_{peak}/I_{bg}');
ylabel('Probability');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

fsh = figure('color', 'w');
colorNow = [151, 102, 102];
hold on;
sh = histScatPlot(spotBgAll, spotValAll, colorNow);
set(fsh, 'CurrentAxes', sh(1))
hold on;
% for i=2:length(spotBg)
%     [~, ~] = scatterFitPlot(spotValRep{i}, spotBgRep{i}, [0 inf], colorNow, 'lineNoFit');
% end
% hold on;
[~, ~, ~] = scatterFitPlot(spotValRepAll, spotBgRepAll, [0 inf], colorNow, 'errNoFit');
hold on;
[~, ~, rSqSpotValBg, slopeSpotValBg] = scatterFitPlot(spotValAll, spotBgAll, [0 inf], colorNow, 'justFit');
xlim([0 500]);
ylim([0 500])
xlabel('I_{bg}');
ylabel('I_{peak}');
% str=sprintf('slope = %1.2f \x00B1 %1.2f', slopeSpotValBg(1), slopeSpotValBg(2));
% T = text(75, 50, str);
% set(T, 'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ('color', 'w');
colorNow = [80 80 80];
[~, ~, ~] = scatterFitPlot(spotValAll./spotBgAll, nucValAll, [0 inf], colorNow, 'errNoFit');
hold on;
[~, ~, rSqSpotValBg, slopeSpotValBg] = scatterFitPlot(spotValAll./spotBgAll, nucValAll, [0 inf], colorNow, 'justFit');
xlim([0 inf]);
ylim([0 inf])
xlabel('I_{nuc}');
ylabel('I_{peak}/ I_{bg}');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vSpotRawTotVal = figure('Color', 'w');
set(0, "CurrentFigure", vSpotRawTotVal)
ax = gca;
ax.Position = posCoord;
hold on;
for i=2:length(spotTotVal)
    [hSpotTotVal{i}, ~, rSqSpotTotVal(i,:)] = scatterFitPlot(spotTotVal{i}, nucVal{i}, valCutOff, colorTotVal, 'lineNoFit');
end
hold on;
[vhSpotTotValAll, ~, ~] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'errNoFit');
hold on;
[~, ~, rSqSpotTotValAll, slopeSpotTotValAll] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'justFit');
xlim([50 200]);
ylim([0 50])
ylabel('I_{tot.}');
xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotTotValAll);
% T = text(75, 10, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');


%%%%%%%%%%%%%%%%%%%%%%% dia vs tot val %%%%%%%%%%%%%%%%%%%%%
% colorNow = [0 0 0];
% spotTotValCutOff = [0.00001, 1000];
% spotDiaCutOff = [0.0001 1];% 
% fsh = figure('color', 'w');
% % sh = histScatPlot(spotValAll, spotFitDiaAll, colorNow);
% % set(fsh, 'CurrentAxes', sh(1))
% % hold on;
% % [vhSpotValAll, ~, ~] = scatterFitPlot(spotFitDiaAll, spotTotValAll, spotValCutOff, colorNow, 'errNoFit');
% % hold on;
% [~, ~, rSq, slope] = scatterFitPlot(spotFitDiaAll, spotTotValAll, spotValCutOff, colorNow, 'justFit');
% hold on;
% xlim([0 inf])
% ylim([0.2 0.6])
% ylabel('\Phi (\mu m)');
% xlabel('I_{tot}');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%~ Fit Error Plots ~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('color', 'w');
f2 = figure('color', 'w');
% hold on;
% binErrPlot(nucValAll, nucPosAll, posCutOff, colorNucVal, bins, minLen, nBoot, f1, f2, 'pos', lambdaNucValAll);
% hold on;
% binErrPlot(spotFitDiaAll, nucPosAll, posCutOff, colorDia, bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotFitDiaAll);
% hold on;
% binErrPlot(spotTotValAll, nucPosAll, posCutOff, colorTotVal, bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotTotValAll);
hold on;
% binErrPlot(spotBgAll, nucPosAll, posCutOff, [0 0 0], bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotBgAll);
% hold on;
% binErrPlot(spotValAll, nucPosAll, posCutOff, colorVal, bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotValAll);

set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.5])
ylabel('\sigma / \mu')
xlabel('x/L');
set(0, "CurrentFigure", f2)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.15])
ylabel('\sigma_{x/L}')
xlabel('x/L');

f1 = figure('color', 'w');
f2 = figure('color', 'w');
hold on;
binErrPlot(spotFitDiaAll, nucValAll, valCutOff, colorDia, bins, minLen, nBoot, f1, f2, 'conc', slopeSpotFitDiaAll);
hold on;
% binErrPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, bins, minLen, nBoot, f1, f2, 'conc', slopeSpotTotValAll);
% hold on;
binErrPlot(spotValAll, nucValAll, valCutOff, [0 0 0], bins, minLen, nBoot, f1, f2, 'conc', slopeSpotValAll);
% hold on;
% binErrPlot(spotBgAll, nucValAll, valCutOff, [0 0 0], bins, minLen, nBoot, f1, f2, 'conc', slopeSpotBgAll);

set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([50 200])
ylim([0 0.5])
ylabel('\sigma / \mu')
xlabel('I_{nuc}');

set(0, "CurrentFigure", f2)
ax = gca;
ax.Position = posCoord;
xlim([50 200])
% set(gca, 'YScale', 'log')
ylim([0 1])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');


% %~~~~~~~~~~~~~~~~~~~~~integration time to read by cluster precisely~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % ................................






aEn = 0.34*10*10^(-3);


D_fast = 10; %unit: mu m^2/s
c_nuc_at0_low = 60; %unit: molecules per mu m^2
c_nuc_at0_hi = 573; % shelby (60,000 molecules per nuc at anterior)
I_nuc_at0 = 522.2; % for bcd 2x calculated from intercept)
D_nuc = 7.5; % abu arish
nuc_vol = 65; % um^3 (2.5 um nuc dia)

cNucAt0 = 840; %value from shelby

concFactor = cNucAt0*exp(-0.2)/mean(nucVal{2});% 60 molecules per um m^3 at mean of position bin #2
nucConc = cellfun(@(x) x.*concFactor, nucVal, 'un', 0);
nucConcAll = vertcat(nucConc{:});
% timeTo10Clust = cellfun(@(x,y)(300./(D_nuc.*x.*y)), spotFitDia, nucConc, 'un', 0 );
% timeTo10ClustAll = vertcat(timeTo10Clust{:});

%%%%%%%%%%%%%%%%%%%
minLen = 40;
[~, allIdx, ~] = binInxFun(nucValAll, 10, minLen);
for i = 1:length(unique(allIdx))
    nTemp = nucValAll(allIdx==i);   
    pTemp = nucPosAll(allIdx==i);  
    cTemp =nucConcAll(allIdx==i);
    dTemp =spotFitDiaAll(allIdx==i);
    tTemp = spotTotValAll(allIdx==i);
    sTemp = spotPerNucAll(allIdx==i);
    fTemp = fracAll(allIdx==i);
    iTemp = spotValAll(allIdx==i);
    bTemp = spotBgAll(allIdx==i);
    if length(nTemp)>=minLen
        nucConcBin{i} = cTemp;
        diaBin{i} = dTemp;        
        nucValBin{i} = nTemp;
        nucPosBin{i} = pTemp;
        totValBin{i} = tTemp;
        spotPerNucBin{i} = sTemp;
        fracBin{i} = fTemp;
        valBin{i} = iTemp;
        bgBin{i} = bTemp;
    end
end

nucPosMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), nucPosBin, 'un', 0)); 
nucPosStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), nucPosBin, 'un', 0)); 

nucValMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), nucValBin, 'un', 0)); 
nucValStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), nucValBin, 'un', 0)); 

valMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), valBin, 'un', 0)); 
valStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), valBin, 'un', 0)); 

bgMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), bgBin, 'un', 0)); 
bgStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), bgBin, 'un', 0));

diaMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), diaBin, 'un', 0)); 
diaStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), diaBin, 'un', 0)); 

spotPerNucMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), spotPerNucBin, 'un', 0)); 
spotPerNucStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), spotPerNucBin, 'un', 0)); 

timeGainMean = diaMean.*valMean./nucValMean./(6*aEn);
timeGainStd = timeGainMean.*sqrt((diaStd./diaMean).^2+(valStd./valMean).^2+(nucValStd./nucValMean).^2);


% timeGainMean = cellfun(@(x,y,z) x.*y./z./aEn, diaMean, valMean,  nucValMean, 'un', 0);
% timeGainStd = cellfun(@(f, a, da, b, db, c, dc) f.*sqrt((da./a)^2+(db./d)^2+(dc./c)^2), ...
%     timeGainMean, nucDiaMean, nucDiaStd, nucValMean, nucValStd, valMean, valStd, 'un', 0);

figure('color', 'w')
pe = errorbar(nucPosMean, timeGainMean, timeGainStd, timeGainStd, nucPosStd, nucPosStd);
pe.Marker = 'square';
pe.MarkerSize = 2;
pe.MarkerEdgeColor = [0 0 0];
pe.LineWidth = 1.5;
pe.Color = [0 0 0];
pe.LineStyle = 'none';
pe.CapSize = 0;
xlabel('x/L')
ylabel('g_T')
xlim([0.15 0.65]);

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  


%%%%%%%%%%% Method1: calculate the number of molecules in a cluster %%%%%%%%%%%%%%%
nucVol = (4*pi/3).*(2.5)^3; %um^3
nucConAt0 = 38000; % molecules of Bcd per nucleus.
lambda = lambdaNucValAll(1);
molPerUMmean = (nucConAt0*exp(-((nucPosMean(end)-0.1)/lambda))).*nucValMean./(nucVol.*nucValMean(end)); % molecules per um^3
molPerUMstd = molPerUMmean.*nucValStd./nucValMean; % molecules per um^3

volMean = (4*pi/3).*(diaMean).^3./8;
volStd = volMean.*3.*(diaStd)./diaMean;

molsPerClusterMean = molPerUMmean.*volMean.*valMean./nucValMean;
molsPerClusterStd = molsPerClusterMean.*((molPerUMstd/molPerUMmean).^2 + (volStd./volMean).^2 + (valStd./valMean).^2 + (bgStd./bgMean).^2).^0.5;
 

figure('color', 'w')
pe = errorbar(nucPosMean, molsPerClusterMean, molsPerClusterStd, molsPerClusterStd, nucPosStd, nucPosStd);
% pe = errorbar(nucValMean, molsPerClusterMean, molsPerClusterStd, molsPerClusterStd, nucValStd, nucValStd);
pe.Marker = 'square';
pe.MarkerSize = 2;
pe.MarkerEdgeColor = [0 0 0];
pe.LineWidth = 1.5;
pe.Color = [0 0 0];
pe.LineStyle = 'none';
pe.CapSize = 0;
xlabel('I_{nuc}')
ylabel('N_clust')

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  


fracMean = (molsPerClusterMean.*spotPerNucMean)./((nucConAt0*exp(-((nucPosMean(end)-0.1)/lambda))).*(nucValMean./nucValMean(end)));
fracStd = fracMean.*((molsPerClusterStd./molsPerClusterMean).^2 + (spotPerNucStd./spotPerNucMean).^2).^0.5;

figure('color', 'w')
pe = errorbar(nucPosMean, fracMean, fracStd, fracStd, nucPosStd, nucPosStd);
% pe = errorbar(nucValMean, fracMean, fracStd, fracStd, nucValStd, nucValStd);
pe.Marker = 'square';
pe.MarkerSize = 2;
pe.MarkerEdgeColor = [0 0 0];
pe.LineWidth = 1.5;
pe.Color = [0 0 0];
pe.LineStyle = 'none';
pe.CapSize = 0;
xlabel('I_{nuc}')
ylabel('Frac_{clust}')
xlim([50 200])
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  
%%%%%%%%%%% Method2: calculate the number of molecules in a cluster %%%%%%%%%%%%%%%


nucConcMean = (cNucAt0 /I_nuc_at0).*(nucValMean);
nucConcStd =  (cNucAt0 /I_nuc_at0).*(nucValStd);

% nucConcMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), nucConcBin, 'un', 0)); 
% nucConcStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), nucConcBin, 'un', 0)); 

fracMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x),fracBin, 'un', 0)); 
fracStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), fracBin, 'un', 0)); 

spotPerNucMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), spotPerNucBin, 'un', 0)); 
spotPerNucStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), spotPerNucBin, 'un', 0)); 

molsPerClusterMean = nucConcMean.*nuc_vol.*fracMean./spotPerNucMean;
molsPerClusterStd = molsPerClusterMean.*((nucConcStd/nucConcMean).^2 + (fracStd./fracMean).^2 + (spotPerNucStd./spotPerNucMean).^2).^0.5;

figure('color', 'w')
pe = errorbar(nucPosMean, molsPerClusterMean, molsPerClusterStd, molsPerClusterStd, nucPosStd, nucPosStd);
pe = errorbar(nucValMean, molsPerClusterMean, molsPerClusterStd, molsPerClusterStd, nucValStd, nucValStd);
pe.Marker = 'square';
pe.MarkerSize = 2;
pe.MarkerEdgeColor = [0 0 0];
pe.LineWidth = 1.5;
pe.Color = [0 0 0];
pe.LineStyle = 'none';
pe.CapSize = 0;
xlabel('x/L}')
ylabel('N_clust')
xlim([0.15 0.65]);

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotTotMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), totValBin, 'un', 0)); 
spotTotStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), totValBin, 'un', 0));

timeToErrMean = 6./(D_nuc.*diaMean.*nucConcMean).*(spotTotStd./spotTotMean).^(-2);
timeToErrStd = timeToErrMean.*((diaStd./diaMean).^2 + (nucConcStd./nucConcMean).^2).^(0.5);

timeTo10ClustMean = (300./(D_nuc.*nucConcMean.*diaMean));
timeTo10ClustStd = timeTo10ClustMean.*((diaStd./diaMean).^2 + (nucConcStd./nucConcMean).^2).^(0.5);

figure('color', 'w')
p10 = errorbar(nucPosMean, timeTo10ClustMean, timeTo10ClustStd , timeTo10ClustStd , nucPosStd, nucPosStd);
p10.Marker = 'v';
p10.MarkerSize = 2;
p10.MarkerEdgeColor = [0.5 0.5 0.5];
p10.LineWidth = 1.5;
p10.Color =  [0.5 0.5 0.5];
p10.LineStyle = 'none';
p10.CapSize = 0;
xlabel('x/L')
ylabel('T_{en} (s)')
xlim([0.15 0.65]);

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  

% 
D_clus = 1; %unit: mu m^2/s
a_en = 0.34*10*10^(-3);

% cluster conc calculation
valMean = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@mean,x), valBin, 'un', 0)); 
valStd = cellfun(@mean, cellfun(@(x) bootstrp(nBoot,@std,x), valBin, 'un', 0)); 

concMean = (cNucAt0 /I_nuc_at0).*(valMean);
concStd =  (cNucAt0 /I_nuc_at0).*(valStd);

timeTo10EnMean = (100./(D_clus.*concMean.*a_en));
timeTo10EnStd = timeTo10EnMean.*(concStd./concMean);

p10 = errorbar(nucPosMean, timeTo10EnMean, timeTo10EnStd , timeTo10EnStd , nucPosStd, nucPosStd);
p10.Marker = 'v';
p10.MarkerSize = 2;
p10.MarkerEdgeColor = [0.5 0.5 0.5];
p10.LineWidth = 1.5;
p10.Color =  [0.5 0.5 0.5];
p10.LineStyle = 'none';
p10.CapSize = 0;
xlabel('x/L')
ylabel('T_en (s)')
xlim([0.15 0.65]);

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  

% clustConc = cellfun(@(x) x.*concFactor, spotVal, 'un', 0);
% 
% timeTo10En = cellfun(@(x)(100./(D_slow.*a_enhancer.*x)), clustConc, 'un', 0 );
% timeTo10EnAll = vertcat(timeTo10En{:});
% 
% vTimeTo10En = figure('Color', 'w');
% set(0, "CurrentFigure", vTimeTo10En)
% ax = gca;
% ax.Position = posCoord;
% hold on;
% for i=2:length(spotFitDia)
%     [phTimeTo10En{i}, ~] = scatterFitPlot((timeTo10En{i}), nucPos{i}, posCutOff, [130 167 146], 'lineNoFit');
% end
% hold on;
% [phTimeTo10EnAll, ~, ~] = scatterFitPlot((timeTo10EnAll), nucPosAll, posCutOff, [130 167 146], 'errNoFit');
% hold on;
% % [~, lambdaTimeTo10All, ~] = scatterFitPlot((timeTo10EnAll), nucPosAll, posCutOff, [130 167 146], 'justFit');
% xlim([0.15 0.65]);
% ylim([0 inf])
% ylabel('T_{en.} (s)');
% xlabel('x/L');


%%%%%%%%%%%%%%%%%%%%%%%%%% end here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clusterConcentration = nucMolAll./((4*pi/3).*(5.5/2).^3);%clusterMolFracAll./spotPerNucAll./((4*pi/3).*(spotMeanFitDiaAll/2).^3);
absoluteConc = clusterConcentration.*10^18;
timeToReadFromCluster = 100./((D_cluster*a_enhancer).*absoluteConc);
figure('Color', 'w')
[dh{1}, lambDaMolPerSpot, ~, ~] = binFitPlot(timeToReadFromCluster, nucPosAll, posCutOff, [245, 116, 29], bins, minLen, nBoot);
ylabel('Enhancer time(s)')

xlabel('x/L');
xlim([0.15 0.6])
% ................................
hold on;
D_nucleus = 7*10^-12;
a_cluster = spotMeanFitDiaAll.*(10^(-6));
nucConcentration = nucMolAll./((4*pi/3).*(5.5/2).^3);
absoluteConc = nucConcentration.*10^18;
timeToReadFromNuc = 100./((D_nucleus*a_cluster).*absoluteConc);

[dh{2}, lambDaMolPerSpot, ~, ~] = binFitPlot(timeToReadFromNuc, nucPosAll, posCutOff, [209, 152, 8], bins, minLen, nBoot);
ylabel('Time(s)')

xlabel('x/L');
xlim([0.15 0.6])

legTex = ["Enhancer time", "Cluster time"];
leg =legend(horzcat(dh{:}), legTex);
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vSpotTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", vSpotTotVal)
% % cMapNucTotVal = colormap(summer(length(nucTotVal)));
% hold on;
% for i=2:length(nucTotVal)
%     [vhSpotTotVal{i}, ~, rSqSpotTotVal(i,:)] = scatterFitPlot(spotTotVal{i}, nucVal{i}, valCutOff, [0, 112, 187], 'lineNoFit');
% end
% hold on;
% [vhSpotTotValAll, ~, ~] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, [0, 112, 187], 'errNoFit');
% hold on;
% [~, ~, rSqSpotTotValAll] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, [0, 112, 187], 'justFit');
% xlim([50 200]);
% ylabel('\Sigma I_{clust.}');
% xlabel('I_{nuc}');
% str=sprintf('R^2 = %1.2f', rSqSpotTotValAll);
% T = text(75, 0.5*10^6, str);%max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');







spotFitDiaAll = cellfun(@(x) x.emSpotProp.spotMeanIntFilt, struct2C, 'un', 0);
for i=1:length(spotFitDiaAll{1})
    spotFitDia{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitDiaAll{1}{i})';
end
spotFitDiaAll = vertcat(spotFitDia{:});







spotFitDiaAll = cellfun(@(x) x.emSpotProp.spotRawTotVal, struct2C, 'un', 0);
for i=1:length(spotFitDiaAll{1})
    spotFitDia{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitDiaAll{1}{i})';
end
spotFitDiaAll = vertcat(spotFitDia{:});

spotFitDiaAll = cellfun(@(x) x.emSpotProp.spotFitDiaFilt, struct2C, 'un', 0);
for i=1:length(spotFitDiaAll{1})
    spotFitDia{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitDiaAll{1}{i})';
end
spotFitDiaAll = vertcat(spotFitDia{:});

spotFitValAll = cellfun(@(x) x.emSpotProp.spotValFilt, struct2C, 'un', 0);
for i=1:length(spotFitValAll{1})
    spotFitVal{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitValAll{1}{i})';
end

spotFitTotValAll = cellfun(@(x) x.emSpotProp.spotTotVal, struct2C, 'un', 0);
for i=1:length(spotFitTotValAll{1})
    spotTotVal{i} = cellfun(@(x) mean(x, 'omitnan'), spotFitTotValAll{1}{i})';
end
spotFitTotValAll = vertcat(spotTotVal{:});

spotBgAll = cellfun(@(x) x.emSpotProp.spotBgFilt, struct2C, 'un', 0);
for i=1:length(spotBgAll{1})
    spotBg{i} = cellfun(@(x) mean(x, 'omitnan'), spotBgAll{1}{i})';
end
spotBgAll = vertcat(spotBg{:});










% % % for BCD 2x
% nucMolAll = nucMolsAt20pc/mean(nucValAll(nucPosAll>0.19 & nucPosAll<0.21)).*(nucValAll);
% % %

nucValRep = cellfun(@(x) x.emSpotProp.nucValRepBinCell, struct2C, 'un', 0);
nucValRepAll = vertcat(nucValRep{1}{:});
nucValRepMean = cellfun(@mean, nucValRep{1});
nucValRepStd = cellfun(@std, nucValRep{1});

nucValRep = cellfun(@(x) x.emSpotProp.nucValRepBinCell, struct2C, 'un', 0);
nucValRepAll = vertcat(nucValRep{1}{:});
nucValRepMean = cellfun(@mean, nucValRep{1});
nucValRepStd = cellfun(@std, nucValRep{1});

spotPerNuc = cellfun(@(x) x.emSpotProp.spotCountCell, struct2C, 'un', 0);
spotPerNucAll = vertcat(spotPerNuc{1}{:});
spotPerNucMean = cellfun(@mean, spotPerNuc{1});
spotPerNucStd = cellfun(@std, spotPerNuc{1});

clusterFrac = cellfun(@(x) x.emSpotProp.clusterFracNucCell, struct2C, 'un', 0);
clusterFracAll = vertcat(clusterFrac{1}{:});
clusterMolFracAll = nucMolAll.*clusterFracAll;

clusterFracMean = cellfun(@mean, clusterFrac{1});
clusterFracStd = cellfun(@std, clusterFrac{1});

clusterSumVal = cellfun(@(x) x.emSpotProp.spotRawTotValSumNucCell, struct2C, 'un', 0);
clusterSumValAll = vertcat(clusterSumVal{1}{:});

diffuseVal = cellfun(@(x) x.emSpotProp.spotBgMeanNucCell, struct2C, 'un', 0);
diffuseValAll = vertcat(diffuseVal{1}{:});
diffuseValMean = cellfun(@mean, diffuseVal{1});
diffuseValStd = cellfun(@std, diffuseVal{1});

spotMeanFitDia = cellfun(@(x) x.emSpotProp.spotFitDiaMeanNucCell, struct2C, 'un', 0);
spotMeanFitDiaAll = vertcat(spotMeanFitDia{1}{:});
spotMeanFitDiaMean = cellfun(@mean, spotMeanFitDia{1});
spotMeanFitDiaStd = cellfun(@std, spotMeanFitDia{1});
spotStdFitDia = cellfun(@(x) x.emSpotProp.spotFitDiaStdNucCell, struct2C, 'un', 0);
spotStdFitDiaAll = vertcat(spotStdFitDia{1}{:});
spotStdFitDiaMean = cellfun(@mean, spotStdFitDia{1});
spotStdFitDiaStd = cellfun(@std, spotStdFitDia{1});
spotFitDiaRep = cellfun(@(x) x.emSpotProp.spotFitDiaRepCell, struct2C, 'un', 0);
spotFitDiaRepAll = vertcat(spotFitDiaRep{1}{:});
spotFitDiaRepMean = cellfun(@mean, spotFitDiaRep{1});
spotFitDiaRepStd = cellfun(@std, spotFitDiaRep{1});

devFitDiaBin = cellfun(@(x, y) sqrt(sum((x-y).^2)./length(x)), spotMeanFitDia{1}, num2cell(spotFitDiaRepMean));
devFitDiaBinNorm = devFitDiaBin./spotFitDiaRepMean;

spotStdFitDia = cellfun(@(x) x.emSpotProp.spotFitDiaStdNucCell, struct2C, 'un', 0);
spotStdFitDiaAll = vertcat(spotStdFitDia{1}{:});
spotStdFitDiaMean = cellfun(@mean, spotStdFitDia{1});


spotMeanRawDia = cellfun(@(x) x.emSpotProp.spotRawDiaNucMeanCell, struct2C, 'un', 0);
spotMeanRawDiaAll = vertcat(spotMeanRawDia{1}{:});
spotMeanRawDiaMean = cellfun(@mean, spotMeanRawDia{1});
spotMeanRawDiaStd = cellfun(@std, spotMeanRawDia{1});
spotStdRawDia = cellfun(@(x) x.emSpotProp.spotRawDiaNucStdCell, struct2C, 'un', 0);
spotStdRawDiaAll = vertcat(spotStdRawDia{1}{:});
spotStdRawDiaMean = cellfun(@mean, spotStdRawDia{1});
spotStdRawDiaStd = cellfun(@std, spotStdRawDia{1});
spotRawDiaRep = cellfun(@(x) x.emSpotProp.spotRawDiaRepCell, struct2C, 'un', 0);
spotRawDiaRepAll = vertcat(spotRawDiaRep{1}{:});
spotRawDiaRepMean = cellfun(@mean, spotRawDiaRep{1});
spotRawDiaRepStd = cellfun(@std, spotRawDiaRep{1});

devRawDiaBin = cellfun(@(x, y) sqrt(sum((x-y).^2)./length(x)), spotMeanRawDia{1}, num2cell(spotRawDiaRepMean));
devRawDiaBinNorm = devRawDiaBin./spotRawDiaRepMean;

spotMeanfitVal = cellfun(@(x) x.emSpotProp.spotFitValMeanNucCell, struct2C, 'un', 0);
spotMeanfitValAll = vertcat(spotMeanfitVal{1}{:});

spotMeanFitTotVal = cellfun(@(x) x.emSpotProp.spotTotValMeanNucCell, struct2C, 'un', 0);%spotRawTotValMeanNucCell, struct2C, 'un', 0); %
% spotMeanFitTotVal = cellfun(@(x) x.emSpotProp.spotBgMeanNucCell, struct2C, 'un', 0);
spotMeanFitTotValAll = vertcat(spotMeanFitTotVal{1}{:});

corrSpotPerNuc = corr(spotPerNucAll(nucValAll<valCutOff(2)), nucValAll(nucValAll<valCutOff(2)));
corrSpotDia = corr(spotMeanFitDiaAll(nucValAll<valCutOff(2)), nucValAll(nucValAll<valCutOff(2)));
% corrSpotBg = corr(spotBgAll(nucValAll<nucValCutOff), nucValAll(nucValAll<nucValCutOff));
% corrSpotVal = corr(spotValAll(nucValAll<nucValCutOff), nucValAll(nucValAll<nucValCutOff));
% corrSpotTotVal = corr(spotTotValAll(nucValAll<nucValCutOff), nucValAll(nucValAll<nucValCutOff));

% % ~~~~~~~~~~~~~~~~~~~~~nuc value~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% %................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % [dh{1}, lambDaNuc] = scatterFitPlot(nucValAll, nucPosAll, posCutOff, [20 20 20]);
% [dh{1}, lambDaNuc, ~, rSqLambDaNuc] = binFitPlot(log(nucValAll), nucPosAll, posCutOff, [20 20 20], bins, minLen, nBoot);
% ylabel("Log of nuclear intensity");
% xlim([0.15 0.6])
% xlabel('Fraction of embryo length (x/L)')
% legTex{1} = {
%         append("{\lambda}_{I_{nuc}} = ", num2str(round(lambDaNuc(1),2)), char(177), ...
%         num2str(round(lambDaNuc(2),2)),  "R^2 = ", num2str(round(rSqLambDaNuc(1),2)))...
%         };
% 
% % leg =legend(horzcat(dh{1}), legTex);
% % set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% xlim([0.17 0.7])
% % ................................
% hold on;
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(nucValAll, nucPosAll, posCutOff, [20 20 20], bins, minLen, nBoot, f1, f2, 'pos', lambDaNuc);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% ylim([0 0.4])
% ylabel('Nuclear intensity (\sigma / \mu)')
% xlabel('Fraction of embryo length (x/L)');
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% ylabel('\sigma(x)/L')
% xlabel('Fraction of embryo length (x/L)');
% 
% % ~~~~~~~~~~~~~~~~~~~~~ cluster density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(0, "CurrentFigure", fv)%figure('Color', 'w');
% % [sc1, rSqSpotPerNuc] = scatterCorrPlot(spotPerNucAll, nucValAll, valCutOff(2), [200 20 20]);
% [sc1, ~, ~, rSqSpotPerNuc] = binFitPlot(spotPerNucAll, nucValAll, valCutOff, [200 20 20], bins, minLen, nBoot);
% xlim([0 200])
% ylabel("#_{cluster}");
% xlabel('I_{nuc}')
% 
% str=sprintf('R^{2} = %1.2f',rSqSpotPerNuc);
% T = text(50, max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% % ................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % [dh{1}, lambDaNuc] = scatterFitPlot(nucValAll, nucPosAll, posCutOff, [20 20 20]);
% [dh, lambDaSpotPerNuc, ~, rSqLambdaSpotPerNuc] = binFitPlot((spotPerNucAll), nucPosAll, posCutOff, [200 20 20], bins, minLen, nBoot);
% xlim([0.15 0.6])
% ylim([3 6])
% ylabel("#_{cluster}");
% xlabel('x/L')
% legTex = {
%         append("{\lambda} = ", num2str(round(lambDaSpotPerNuc(1),2)), char(177), ...
%         num2str(round(lambDaSpotPerNuc(2),2)),  newline, "R^2 = ", num2str(round(rSqLambdaSpotPerNuc(1),2)))...
%         };
% 
% leg =legend(horzcat(dh), legTex);
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% % % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(spotPerNucAll, nucPosAll, posCutOff, [200 20 20], bins, minLen, nBoot, f1, f2, 'pos', lambDaSpotPerNuc);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% xlim([0.1 0.25])
% ylabel('\sigma / \mu')
% xlabel('x/L');
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% 
% ylabel('\sigma_{x/L}')
% xlabel('x/L');

% % ~~~~~~~~~~~~~~~~~~~~~~~ clustered fraction ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% figure('Color', 'w');
% [sc1, ~, ~, rSqSpotPerNuc] = binFitPlot(clusterFracAll, nucValAll, valCutOff, [100 200 20], bins, minLen, nBoot);
% xlim([0 200])
% ylabel("\Sigma(mol_{clust.})/\Sigma(mol_{nuc})");
% xlabel('I_{nuc}')
% ylim([0.02 0.05])
% str=sprintf('R^{2} = %1.2f',rSqSpotPerNuc);
% T = text(50, max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% % % ................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % figure ('color', 'w')
% [dh{1}, lambDaclusterFrac, ~, rSqLambdaclusterFrac] = binFitPlot((clusterFracAll), nucPosAll, posCutOff, [100 200 20], bins, minLen, nBoot);
% ylabel("\Sigma(mol_{clust.})/\Sigma(mol_{nuc})");;%('(I_{tot})_{all clusters}/(I_{tot})_{nuc}') % ylabel("Clustered fraction");
% xlabel('x/L')
% legTex = {
%         append("{\lambda} = ", num2str(round(lambDaclusterFrac(1),2)), char(177), ...
%         num2str(round(lambDaclusterFrac(2),2)), "R^2 = ", num2str(round(rSqLambdaclusterFrac(1),2)))...
%         };
% 
% leg =legend(horzcat(dh{1}), legTex);
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% xlim([0.15 0.6])
% ylim([0.02 0.05])
% yticks([0.02 0.03 0.04 0.05])
% % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(clusterFracAll, nucPosAll, posCutOff, [100 200 20], bins, minLen, nBoot, f1, f2, 'pos', lambDaclusterFrac);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% ylabel('Clustered fraction (\sigma / \mu)')
% yticks([0.02 0.03 0.04 0.05])
% xlabel('x/L');
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% ylabel('\sigma_{x/L}')
% xlabel('x/L');
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~ clustered molecules ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(0, "CurrentFigure", fv)%figure('Color', 'w');
% % [sc1, rSqSpotPerNuc] = scatterCorrPlot(spotPerNucAll, nucValAll, valCutOff(2), [200 20 20]);
% [sc1, ~, ~, rSqMolFrac] = binFitPlot(clusterMolFracAll, nucValAll, valCutOff,  [9, 121, 105], bins, minLen, nBoot);
% xlim([0 200])
% ylabel('\Sigma(mol_{clust.})');
% xlabel('I_{nuc} (a. u.)')
% 
% str=sprintf('R^{2} = %1.2f',rSqMolFrac);
% T = text(50, max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% % ................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % [dh{1}, lambDaSpotPerNuc] = scatterFitPlot(spotPerNucAll, nucPosAll, posCutOff, [200 20 20]);
% [dh, lambDaMolFrac, ~, rSqLambdaMolFrac] = binFitPlot(log(clusterMolFracAll), nucPosAll, posCutOff, [9, 121, 105], bins, minLen, nBoot);
% ylabel('log(mol_{clust.})') %ylabel("Log of clustered molecules");
% xlabel('x/L')
% legTex = {
%         append("{\lambda} = ", num2str(round(lambDaMolFrac(1),2)), char(177), ...
%         num2str(round(lambDaMolFrac(2),2)), newline, "R^2 = ", num2str(round(rSqLambdaMolFrac(1),2)))...
%         };
% 
% % leg =legend(horzcat(dh{1}), legTex);
% % set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% xlim([0.15 0.6])
% % % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(clusterMolFracAll, nucPosAll, posCutOff,  [9, 121, 105], bins, minLen, nBoot, f1, f2, 'pos', lambDaMolFrac);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% ylabel('\sigma / \mu') %ylabel('Clustered molecules (\sigma / \mu)')
% xlabel('x/L');
% ylim([0 0.4])
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% ylim([0 0.1])
% ylabel('\sigma_{x/L}')
% xlabel('x/L');
% %



% % ~~~~~~~~~~~~~~~~~~~~~~~ diffused intensity~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(0, "CurrentFigure", fv)%figure('Color', 'w');
% % [sc1, rSqSpotPerNuc] = scatterCorrPlot(spotPerNucAll, nucValAll, valCutOff(2), [200 20 20]);
% [sc1, ~, ~, rSqDiff] = binFitPlot(diffuseValAll, nucValAll, valCutOff,  [251, 206, 177], bins, minLen, nBoot);
% xlim([0 200])
% ylabel("I_{non-clustered}");
% xlabel('I_{nuc} (a. u.)');
% 
% str=sprintf('R^{2} = %1.2f',rSqDiff);
% T = text(50, max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% % ................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % [dh{1}, lambDaSpotPerNuc] = scatterFitPlot(spotPerNucAll, nucPosAll, posCutOff, [200 20 20]);
% [dh{2}, lambDaDiff, ~, rSqLambdaDiff] = binFitPlot(log(diffuseValAll), nucPosAll, posCutOff, [251, 206, 177], bins, minLen, nBoot);
% ylabel('log(I_{non-cluster})') %ylabel("Log of clustered molecules");
% xlabel('x/L')
% legTex{2} = {
%         append("{\lambda}_{I_{non cluster}} = ", num2str(round(lambDaDiff(1),2)), char(177), ...
%         num2str(round(lambDaDiff(2),2)), "R^2 = ", num2str(round(rSqLambdaDiff(1),2)))...
%         };
% 
% % leg =legend(horzcat(dh), legTex);
% % set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% xlim([0.15 0.6])
% % % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(diffuseValAll, nucPosAll, posCutOff,  [251, 206, 177], bins, minLen, nBoot, f1, f2, 'pos', lambDaDiff);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% ylabel('\sigma / \mu') %ylabel('Clustered molecules (\sigma / \mu)')
% xlabel('x/L');
% ylim([0 inf])
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% ylim([0 0.1])
% ylabel('\sigma_{x/L}')
% xlabel('x/L');


% % ~~~~~~~~~~~~~~~~~~~~~~~ cluster total intensity sum ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(0, "CurrentFigure", fv)%figure('Color', 'w');
% % [sc1, rSqSpotPerNuc] = scatterCorrPlot(spotPerNucAll, nucValAll, valCutOff(2), [200 20 20]);
% [sc1, ~, ~, rSqClusterSum] = binFitPlot(clusterSumValAll, nucValAll, valCutOff,  [255, 87, 51], bins, minLen, nBoot);
% xlim([0 200])
% ylabel("I_{clustered}");
% xlabel('I_{nuc} (a. u.)');
% 
% str=sprintf('R^{2} = %1.2f',rSqClusterSum);
% T = text(50, max(get(gca, 'ylim')), str);
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% % ................................
% set(0, "CurrentFigure", fb)%figure('Color', 'w');
% % [dh{1}, lambDaSpotPerNuc] = scatterFitPlot(spotPerNucAll, nucPosAll, posCutOff, [200 20 20]);
% [dh{3}, lambDaClusterSum, ~, rSqClusterSum] = binFitPlot(log(clusterSumValAll), nucPosAll, posCutOff, [255, 87, 51], bins, minLen, nBoot);
% ylabel('log(I_{non-clustered})') %ylabel("Log of clustered molecules");
% xlabel('x/L')
% % legTex{3} = {
% %         append("{\lambda} = ", num2str(round(lambDaClusterSum(1),2)), char(177), ...
% %         num2str(round(lambDaClusterSum(2),2)), newline, "R^2 = ", num2str(round(rSqClusterSum(1),2)))...
% %         };
% legTex{3} = {
%         append("{\lambda}_{I_{cluster}} = ", num2str(round(lambDaClusterSum(1),2)), char(177), ...
%         num2str(round(lambDaClusterSum(2),2)), "R^2 = ", num2str(round(rSqClusterSum(1),2)))...
%         };
% 
% leg =legend(horzcat(dh{:}), horzcat(legTex{:}));
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% xlim([0.15 0.6])
% ylim([0 inf])
% % % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('Color', 'w');
% binErrPlot(clusterSumValAll, nucPosAll, posCutOff,  [255, 87, 51], bins, minLen, nBoot, f1, f2, 'pos', lambDaDiff);
% set(0, "CurrentFigure", f1)
% xlim([0.15 0.6])
% ylabel('\sigma / \mu') %ylabel('Clustered molecules (\sigma / \mu)')
% xlabel('x/L');
% ylim([0 0.4])
% set(0, "CurrentFigure", f2)
% xlim([0.15 0.6])
% ylim([0 0.1])
% ylabel('\sigma_{x/L}')
% xlabel('x/L');

% % % ~~~~~~~~~~~~~~~~~~~~~~~diameter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(0, "CurrentFigure", fv)%figure('Color', 'w');
% % [sc2, rSqDia] = scatterCorrPlot(spotMeanFitDiaAll, nucValAll, valCutOff(2), [20 20 200]);
% [sc2{1}, ~, diaSlope, rSqDia] = binFitPlot(spotMeanFitDiaAll, nucValAll, valCutOff, [20 20 200], bins, minLen, nBoot);
% ylim([0.3 0.5])
% xlim([0 200])
% ylabel("Cluster diameter ({\mu}m)");
% xlabel('Mean nuclear intensity (a. u.)')
% 
% % str=sprintf('R^{2} = %1.2f',rSqDia);
% % T = text(50, max(get(gca, 'ylim')), str); 
% % set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% 
% legTexV{1} = {append('R^{2}_{dia.} = ', num2str(round(rSqDia,2)))};
% legV = legend(horzcat(sc2{:}), horzcat(legTexV{:}));
% set(legV,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% %.................................
% 
% binErrPlot(spotMeanFitDiaAll, nucValAll, valCutOff, [20 20 200], bins, minLen, nBoot, f3, f4, 'conc', diaSlope);
% set(0, "CurrentFigure", f3)
% xlim([0 200])
% ylim([0 0.1])
% ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
% xlabel('I_{nuc}');
% set(0, "CurrentFigure", f4)
% xlim([0 200])
% ylim([0 0.4])
% ylabel('\sigma_{I_{nuc}}')
% xlabel('I_{nuc}');

% ................................
set(0, "CurrentFigure", fb)%figure('Color', 'w');
% [dh{1}, lambDaSpotMeanFitDia] = scatterFitPlot(spotMeanFitDiaAll, nucPosAll, posCutOff, [20 20 200]);
[dh{1}, lambDaSpotMeanFitDia, ~, rSqlambDaSpotMeanFitDia] = binFitPlot((spotMeanFitDiaAll), nucPosAll, posCutOff, [20 20 200], bins, minLen, nBoot);
ylim([-1.1 -0.7])
ylabel("Log of cluster diameter ({\mu}m)");
xlabel('x/L');
legTex{1} = {
        append("{\lambda}_{dia.} = ", num2str(round(lambDaSpotMeanFitDia(1),2)), char(177), ...
        num2str(round(lambDaSpotMeanFitDia(2),2)))...
        };

% leg =legend(horzcat(dh{1}), legTex);
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
xlim([0.15 0.6])
% ................................
% f1 = figure('Color', 'w');
% f2 = figure('Color', 'w');
binErrPlot(spotMeanFitDiaAll, nucPosAll, posCutOff, [20 20 200], bins, minLen, nBoot, f1, f2, 'pos', lambDaSpotMeanFitDia);
set(0, "CurrentFigure", f1)
xlim([0.15 0.6])
ylim([ 0 0.1])
ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
xlabel('x/L');
set(0, "CurrentFigure", f2)
xlim([0.15 0.6])
ylim([0 0.4])
ylabel('\sigma_{x/L}')
xlabel('x/L');

% figure('Color', 'w')
% handle = binErrPlot2(spotMeanFitDiaAll, nucValAll, valCutOff, [20 20 200], bins, minLen, nBoot);
% xlabel('Mean nuclear intensity (a. u.)')
% ylabel('Cluster diameter (\sigma / \mu)')
% xlim([0 200])

kernelPlotPosBin(nucPosAll, spotMeanFitDiaAll, "dia ({\mu}m)", '2x');
%%~~~~~~~~~~~~~~~~~~~ intensity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


set(0, "CurrentFigure", fv)%figure('Color', 'w');
% [sc2, rSqDia] = scatterCorrPlot(spotMeanFitDiaAll, nucValAll, valCutOff(2), [20 20 200]);
[sc2{2}, ~, valSlope, rSqMeanfitVal] = binFitPlot(spotMeanfitValAll, nucValAll, valCutOff, [20 200 200], bins, minLen, nBoot);
% ylim([0.3 0.5])
xlim([0 200])
ylabel("Cluster intensity");
xlabel('Mean nuclear intensity (a. u.)')

% str=sprintf('R^{2} = %1.2f',rSqMeanfitVal);
% T = text(50, max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

legTexV{2} = {append('R^{2}_{int.} = ',num2str(rSqMeanfitVal, 2))};
legV = legend(horzcat(sc2{:}), horzcat(legTexV{:}));
set(legV,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');

%.................................

binErrPlot(spotMeanfitValAll, nucValAll, valCutOff, [20 200 200], bins, minLen, nBoot, f3, f4, 'conc', valSlope);
set(0, "CurrentFigure", f3)
xlim([0 200])
ylim([0 0.1])
ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
xlabel('I_{nuc}');
set(0, "CurrentFigure", f4)
xlim([0 200])
ylim([0 0.4])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');
% ................................
set(0, "CurrentFigure", fb)%figure('Color', 'w');
% [dh{1}, lambDaSpotMeanFitDia] = scatterFitPlot(spotMeanFitDiaAll, nucPosAll, posCutOff, [20 20 200]);
[dh{2}, lambDaSpotMeanFitVal, ~, ~] = binFitPlot(log(spotMeanfitValAll), nucPosAll, posCutOff, [20 200 200], bins, minLen, nBoot);
ylabel("Log of cluster intensity");
xlabel('x/L');
legTex{2} = {
        append("{\lambda}_{I} = ", num2str(round(lambDaSpotMeanFitVal(1),2)), char(177), ...
        num2str(round(lambDaSpotMeanFitVal(2),2)))};

% leg =legend(horzcat(dh{1}), legTex);
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
xlim([0.15 0.6])
% ................................
% f1 = figure('Color', 'w');
% f2 = figure('color', 'w');
binErrPlot(spotMeanfitValAll, nucPosAll, posCutOff, [20 200 200], bins, minLen, nBoot, f1, f2, 'pos', lambDaSpotMeanFitVal);
set(0, "CurrentFigure", f1)
ylim([0 0.4])
xlim([0.15 0.6])
ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
xlabel('x/L');
set(0, "CurrentFigure", f2)
xlim([0.15 0.6])
ylim([0 0.4])
ylabel('\sigma_{x/L}')
xlabel('x/L');
%
%%~~~~~~~~~~~~~~~~~~~~~total intensity~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
set(0, "CurrentFigure", fv)%figure('Color', 'w');
% [sc2, rSqDia] = scatterCorrPlot(spotMeanFitDiaAll, nucValAll, valCutOff(2), [20 20 200]);
[sc2{3}, ~, totValSlope, rSqMeanfitTotVal] = binFitPlot(spotMeanFitTotValAll, nucValAll, valCutOff, [200 20 200], bins, minLen, nBoot);
ylim([0 6])
xlim([0 200])
ylabel("y");
xlabel('I_{nuc}')

% str=sprintf('R^{2} = %1.2f',rSqMeanfitTotVal);
% T = text(50, max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

legTexV{3} = {append('R^{2}_{tot. int.} = ',num2str(round(rSqMeanfitTotVal, 2)))};
legV = legend(horzcat(sc2{:}), horzcat(legTexV{:}));
set(legV,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
%.................................

binErrPlot(spotMeanFitTotValAll, nucValAll, valCutOff, [200 20 200], bins, minLen, nBoot, f3, f4, 'conc', totValSlope);
set(0, "CurrentFigure", f3)
xlim([0 200])
ylim([0 1])
ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
xlabel('I_{nuc}');
set(0, "CurrentFigure", f4)
xlim([0 200])
ylim([0 1])
ylabel('\sigma_{I_{nuc}}')
xlabel('I_{nuc}');



% ................................
set(0, "CurrentFigure", fb)%figure('Color', 'w');
% [dh{1}, lambDaSpotMeanFitDia] = scatterFitPlot(spotMeanFitDiaAll, nucPosAll, posCutOff, [20 20 200]);
[dh{3}, lambDaSpotMeanFitTotVal, ~, ~] = binFitPlot(log(spotMeanFitTotValAll), nucPosAll, posCutOff, [200 20 200], bins, minLen, nBoot);
ylabel('log(y)/log(y_{x_{0}})')
% ylabel("Log of cluster total intensity");
xlabel('x/L');
legTex{3} = {
        append("{\lambda}_{I_{tot.}} = ", num2str(round(lambDaSpotMeanFitTotVal(1),2)), char(177), ...
        num2str(round(lambDaSpotMeanFitTotVal(2),2)))};

leg = legend(horzcat(dh{:}), horzcat(legTex{:}));
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
% ................................
% f1 = figure('Color', 'w');
% f2 = figure('color', 'w');
binErrPlot(spotMeanFitTotValAll, nucPosAll, posCutOff, [200 20 200], bins, minLen, nBoot, f1, f2, 'pos', lambDaSpotMeanFitTotVal);
set(0, "CurrentFigure", f1)
xlim([0.15 0.6])
ylim([0 0.4])
ylabel('\sigma / \mu');%ylabel('Cluster diameter (\sigma / \mu)')
xlabel('x/L');
set(0, "CurrentFigure", f2)
xlim([0.15 0.6])
ylim([0 0.4])
ylabel('\sigma_{x/L}')
xlabel('x/L');

%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% %~~~~~~~~~~~~~~~~~~~~~mols per cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% figure('Color', 'w');
% % [sc2, rSqDia] = scatterCorrPlot(spotMeanFitDiaAll, nucValAll, valCutOff(2), [20 20 200]);
% [sc2, ~, ~, rSqMolPerSpot] = binFitPlot(clusterMolFracAll./spotPerNucAll, nucValAll, valCutOff, [220 20 100], bins, minLen, nBoot);
% % ylim([0.3 0.5])
% xlim([0 200])
% ylabel("N_{mol}/clust.");
% xlabel('I_{nuc}')
% 
% str=sprintf('R^{2} = %1.2f',rSqMolPerSpot);
% T = text(50, max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
% % ................................
% 
% figure('Color', 'w')
% [dh{1}, lambDaMolPerSpot, ~, ~] = binFitPlot(clusterMolFracAll./spotPerNucAll, nucPosAll, posCutOff, [220 20 100], bins, minLen, nBoot); %binFitPlot(clusterMolFracAll./spotPerNucAll, nucPosAll, posCutOff, [220 20 100], bins, minLen, nBoot);%
% ylabel('Molecules per cluster');
% ylabel('N_{mol}/clust.')
% % ylabel(append("Cluster concentration", newline, "(molecules / {\mu}m^3)"));
% xlim([0.15 0.6])
% xlabel('x/L');
% % ................................
% 
% 
% figure('Color', 'w')
% clusterConcentration = clusterMolFracAll./spotPerNucAll./((4*pi/3).*(spotMeanFitDiaAll/2).^3);
% molsPerCluster = clusterMolFracAll./spotPerNucAll;
% [dh{1}, lambDaMolPerSpot, ~, ~] = binFitPlot(clusterConcentration, nucPosAll, posCutOff, [220 20 100], bins, minLen, nBoot); %binFitPlot(clusterMolFracAll./spotPerNucAll, nucPosAll, posCutOff, [220 20 100], bins, minLen, nBoot);%
% ylabel('Molecules per cluster');
% ylabel('N_{mol}/clust.')
% ylabel('C_{mol_{clust.}} (molecules / {\mu}m^3)')
% % ylabel(append("Cluster concentration", newline, "(molecules / {\mu}m^3)"));
% xlim([0.15 0.6])
% xlabel('x/L');
% % 
% % % ................................
% % f1 = figure('Color', 'w');
% % f2 = figure('color', 'w');
% % 
% % binErrPlot(clusterMolFracAll./spotPerNucAll, nucPosAll, posCutOff, [220 20 100], bins, minLen, nBoot, f1, f2, 'pos', lambDaMolPerSpot);
% % set(0, "CurrentFigure", f1)
% % xlim([0.15 0.6])
% % ylim([0 0.4])
% % ylabel('Molecules per cluster (\sigma / \mu)')
% % xlabel('Fraction of embryo length (x/L)');
% % set(0, "CurrentFigure", f2)
% % xlim([0.15 0.6])
% % ylim([0 0.4])
% % ylabel('\sigma(x)/L')
% % xlabel('Fraction of embryo length (x/L)');
% % %
% 
% %~~~~~~~~~~~~~~~~~~~~~integration time to read precisely~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % ................................
% D_cluster = 7*10^-12;%0.3*10^-12;
% a_enhancer = 0.34*10^(-9)*10;
% clusterConcentration = nucMolAll./((4*pi/3).*(5.5/2).^3);%clusterMolFracAll./spotPerNucAll./((4*pi/3).*(spotMeanFitDiaAll/2).^3);
% absoluteConc = clusterConcentration.*10^18;
% timeToReadFromCluster = 100./((D_cluster*a_enhancer).*absoluteConc);
% figure('Color', 'w')
% [dh{1}, lambDaMolPerSpot, ~, ~] = binFitPlot(timeToReadFromCluster, nucPosAll, posCutOff, [245, 116, 29], bins, minLen, nBoot);
% ylabel('Enhancer time(s)')
% 
% xlabel('x/L');
% xlim([0.15 0.6])
% % ................................
% hold on;
% D_nucleus = 7*10^-12;
% a_cluster = spotMeanFitDiaAll.*(10^(-6));
% nucConcentration = nucMolAll./((4*pi/3).*(5.5/2).^3);
% absoluteConc = nucConcentration.*10^18;
% timeToReadFromNuc = 100./((D_nucleus*a_cluster).*absoluteConc);
% 
% [dh{2}, lambDaMolPerSpot, ~, ~] = binFitPlot(timeToReadFromNuc, nucPosAll, posCutOff, [209, 152, 8], bins, minLen, nBoot);
% ylabel('Time(s)')
% 
% xlabel('x/L');
% xlim([0.15 0.6])
% 
% legTex = ["Enhancer time", "Cluster time"];
% leg =legend(horzcat(dh{:}), legTex);
% set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 10, 'Box', 'off');
  
%
end

function [xBin, yBin] = binMaker(xData, yData, nBins, minLen)
[groupLen, valBin, idx] = histcounts(xData, nBins);
valBin = (valBin(1:end-1) + diff(valBin) / 2)';
[valBinSort, sortIdx] = sort(valBin);

tempIdx = idx;
for i = 1:length(groupLen)
    idx(tempIdx==sortIdx(i)) = (i);    
end
groupLenSort = groupLen';
groupLenSort(sortIdx) = groupLenSort;
valBinSort(groupLenSort<minLen) = [];


for i = 1:length(unique(idx))
    yTemp = yData(idx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xData(idx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));
end

function [ph, lambda, slope, rSq] = binFitPlot(yVal, xVal, cutOff, color, bins, minLen, nBoot)
color = color./255;

yVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
% yVal = log(yVal);
xVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
% dataTable = table(xVal, yVal);

[nucValBinMeanSort, allIdx, groupLengthSort] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));
%%%%%%% use for 6x %%%%%%%%
% xBin = xBin(1:7);
% yBin = yBin(1:7);

xMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xBin, 'un', 0); % mean of bins
xMeanMean = cellfun(@mean, xMean);
xStd = cellfun(@(x) bootstrp(nBoot,@std,x), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@mean, xStd);

yMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yBin, 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yMean);
yStd = cellfun(@(x) bootstrp(nBoot,@std,x), yBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, yStd);


kernelPlotPosBin(arrayfun(@(x) num2str(x, 2), round(xMeanMean, 2), 'un', 0), yBin, "dia ({\mu}m)", '2x');


%%%% normalization
yMeanMeanNorm = yMeanMean./yMeanMean(1); %yMeanMean; % 
yStdMeanNorm = yStdMean./yMeanMean(1); % yStdMean; %

pe = errorbar(xMeanMean, yMeanMeanNorm, yStdMeanNorm, yStdMeanNorm, xStdMean, xStdMean);
pe.Marker = 'o';
pe.MarkerFaceColor = color;
pe.LineWidth = 1;
pe.Color = color;
pe.LineStyle = 'none';
pe.CapSize = 0;
hold on;

dataTable = table(xMeanMean', yMeanMeanNorm');
mdl = fitlm(dataTable);
rSq = mdl.Rsquared.Adjusted;
hold on;
fitVar = mdl.Coefficients.Variables;
xArr = linspace(0, max(xVal));
yArr = fitVar(2, 1).*xArr + fitVar(1, 1);
pl = plot(xArr, yArr);
pl.LineStyle = '-';
pl.LineWidth = 1.2;
pl.Color = color;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

ph = pe;%pl;%
slope = [fitVar(2,1), fitVar(2,2)];
lambda = [-(1/(fitVar(2,1).*yMeanMean(1))), fitVar(2,2)/(fitVar(2,1).*yMeanMean(1))]; % if normalized
% lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
rSq = corr(xVal, yVal)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pl = binErrPlot(yVal, xVal, cutOff, color, bins, minLen, nBoot, f1, f2, type, factor)
yVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
xVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
color = color./255;

[nucValBinMeanSort, allIdx, groupLengthSort] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

%   Bin data
xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));

yMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yBin, 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yMean);
yMeanStd = cellfun(@std, yMean);
yStd = cellfun(@(x) bootstrp(nBoot,@std,x), yBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, yStd);
yStdStd = cellfun(@std, yStd);

xMeanMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xBin, 'un', 0); % mean of bins
xMeanMean = cellfun(@mean, xMeanMean);
xStd = cellfun(@(x) bootstrp(nBoot,@std,x), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@mean, xStd);
xStdStd = cellfun(@std, xStd);





% [xMeanMean, xStd, yMean, yStd] = getBinMeans(xVal, yVal, bins, minLen);

%%%%%%%% normalization %%%%%%%%%%%%
normErrMean = yStdMean./yMeanMean;
normErrErr = normErrMean.*sqrt((sqrt(nBoot).*(yStdStd./yStdMean).^2) + (sqrt(nBoot).*(yMeanStd./yMeanMean).^2));

% xMeanMean = cellfun(@mean, xBin);

set(0, "CurrentFigure", f1)
pl = simpleErrPlotter(xMeanMean', normErrMean, normErrErr, color);

set(0, "CurrentFigure", f2)
if strcmp(type, 'conc')
    errConcEstMean = yStdMean./factor(1);
    errConcEstStd = (yStdMean./factor(1)).*sqrt((yStdStd./yStdMean).^2+(factor(2)./factor(1)).^2);% yStdStd./factor(1);
    simpleErrPlotter(xMeanMean', errConcEstMean./xMeanMean, errConcEstStd./xMeanMean, color);
elseif strcmp(type, 'pos')
%     yMean = cellfun(@(x) bootstrp(nBoot,@mean,log(x)), yBin, 'un', 0); % mean of bins
%     yMeanMean = cellfun(@mean, yMean);
%     yMeanStd = cellfun(@std, yMean);
%     yStd = cellfun(@(x) bootstrp(nBoot,@std,log(x)), yBin, 'un', 0); % std of bins
%     yStdMean = cellfun(@mean, yStd);
%     yStdStd = cellfun(@std, yStd);

    errPosEstMean = factor(1).*yStdMean./yMeanMean;
    errPosEstStd = factor(1).*(yStdMean./yMeanMean).*sqrt((yStdStd./yStdMean).^2+(yMeanStd./yMeanMean).^2 + (factor(2)./factor(1)).^2);
    simpleErrPlotter(xMeanMean', errPosEstMean, errPosEstStd, color);
    hold on;
%     simpleErrPlotter(xMeanMean', xStdMean, xStdStd, [0.6 0.6 0.6]); %
% %     this draws the gret binning erro line

end
end

function pl = simpleErrPlotter(xMean, yMean, yErr, color)
%%%%%% use for 6x %%%%%%
% xMean = xMean(1:7);
% yMean = yMean(1:7);
% yErr = yErr(1:7);
%%%%%%%%%%%%%%%%%%
pl = errorbar(xMean', yMean, yErr);
pl.Marker = 'o';
pl.MarkerFaceColor = color;
pl.MarkerEdgeColor = color;
pl.MarkerSize = 3;
pl.LineWidth = 1;
pl.Color = color;
pl.LineStyle = '-';
pl.CapSize = 0;
hold on;
pt = patch([0 max(xMean) max(xMean) 0], [(mean(yMean) + std(yMean)) (mean(yMean) + std(yMean)) (mean(yMean) - std(yMean)) (mean(yMean) - std(yMean))], color);
pt.FaceAlpha = 0.3;
pt.EdgeAlpha = 0;

% plot([0 max(xMean)], [mean(yMean), mean(yMean)], 'Color', color, 'LineStyle','-');
% hold on;
% plot([0 max(xMean)], [mean(yMean) + std(yMean), mean(yMean) + std(yMean)], 'Color', color, 'LineStyle','--');
% hold on;
% plot([0 max(xMean)], [mean(yMean) - std(yMean), mean(yMean) - std(yMean)], 'Color', color, 'LineStyle','--');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 150;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  
end

function kernelPlotPosBin(binNames, data, xAxisLabel, plotTitle, color)
color = color./255;
figure('Color', [1, 1, 1]);
% Compute axes positions with contigunous edges
n = length(binNames); 
margins = [.13 .13 .12 .15]; %left, right, bottom, top
height = (1-sum(margins(3:4)))/n; % height of each subplot
width = 1-sum(margins(1:2)); %width of each sp
vPos = linspace(margins(3),1-margins(4)-height,n); %vert pos of each sp
subHand = gobjects(1,n);
histHand = gobjects(2,n);

for i = 1:n
    subHand(i) = axes('position',[margins(1),vPos(i),width,height]); 
    if length(data{i})>5
        histHand(:,i) = histfit(data{i}, 50, 'kernel');
    end
end

% Link the subplot x-axes
linkaxes(subHand,'x')

% Extend density curves to edges of xlim and fill.
% This is easier, more readable (and maybe faster) to do in a loop. 
xl = xlim(subHand(end));
colors = repmat(color, n, 1);% jet(n); % Use any colormap you want
for i = 1:n
     if length(data{i})>5
        x = [xl(1),histHand(2,i).XData,xl([2,1])]; 
        y = [0,histHand(2,i).YData,0,0]; 
        fillHand = fill(subHand(i), x, y, colors(i, :), ...
            'FaceAlpha', 0.9, 'EdgeColor' ,'k', 'LineWidth', 1, 'LineStyle','none');
        % Add vertical ref lines at xtick of bottom axis
%         arrayfun(@(t)xline(subHand(i),t),subHand(1).XTick); %req. >=r2018b
        % Add y axis labels
        ylh = ylabel(subHand(i),binNames{i}); 
        set(ylh,'Rotation',0,'HorizontalAlignment',...
            'right','VerticalAlignment','middle', 'FontSize', 10);
        subHand(i).LineWidth = 1;
     end
end
% Cosmetics
% Delete histogram bars & original density curves 
delete(histHand)
% remove axes (all but bottom) and 
% add vertical ref lines at x ticks of bottom axis
arrayfun(@(i)set(subHand(i),'Box','off'), 1:n);
arrayfun(@(i)set(subHand(i).XAxis,'Visible','off'),2:n)
set(subHand,'YTick',[])
set(subHand,'XLim',xl)
% ylabel('%EL');
subHand(1).XLabel.String = xAxisLabel;
subHand(1).XLabel.FontSize = 10;
subHand(1).FontSize = 10;

x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
grid off;
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w'); 
end

function kernelPlotPosBin2(binNames, data, xAxisLabel, plotTitle, color)
color = color./255;
figure('Color', [1, 1, 1]);
% Compute axes positions with contigunous edges
n = length(binNames); 
margins = [.13 .13 .12 .15]; %left, right, bottom, top
height = (1-sum(margins(3:4)))/n; % height of each subplot
width = 1-sum(margins(1:2)); %width of each sp
vPos = linspace(margins(3),1-margins(4)-height,n); %vert pos of each sp
subHand = gobjects(1,n);
histHand = gobjects(2,n);
edges = min(vertcat(data{:})):(max(vertcat(data{:}))-min(vertcat(data{:})))/50:max(vertcat(data{:}));
for i = 1:n
    subHand(i) = axes('position',[margins(1),vPos(i),width,height]); 
    if length(data{i})>5
        % histHand(:,i) = histogram(data{i}, 50, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [0 0 0]); % use 50 for all values
        % hold on;
        % histHand(:,i) = histogram(data{i}, 50, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', color, 'EdgeColor', 'none');

        histHand(:,i) = histogram(data{i}, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [0 0 0]); % use 50 for all values
        hold on;
        histHand(:,i) = histogram(data{i}, edges, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', color, 'EdgeColor', 'none');

    end
end

% Link the subplot x-axes
linkaxes(subHand,'x')

% Extend density curves to edges of xlim and fill.
% This is easier, more readable (and maybe faster) to do in a loop. 
xl = xlim(subHand(end));
colors = repmat(color, n, 1);% jet(n); % Use any colormap you want
for i = 1:n
     if length(data{i})>5
        % x = [xl(1),histHand(2,i).XData,xl([2,1])]; 
        % y = [0,histHand(2,i).YData,0,0]; 
        % fillHand = fill(subHand(i), x, y, colors(i, :), ...
        %     'FaceAlpha', 0.9, 'EdgeColor' ,'k', 'LineWidth', 1, 'LineStyle','none');
        % % Add vertical ref lines at xtick of bottom axis
        % arrayfun(@(t)xline(subHand(i),t),subHand(1).XTick); %req. >=r2018b
        % % Add y axis labels
        ylh = ylabel(subHand(i),binNames{i}); 
        set(ylh,'Rotation',0,'HorizontalAlignment',...
            'right','VerticalAlignment','middle', 'FontSize', 10);
        subHand(i).LineWidth = 1;
     end
end
% Cosmetics
% Delete histogram bars & original density curves 
% delete(histHand)
% remove axes (all but bottom) and 
% add vertical ref lines at x ticks of bottom axis
arrayfun(@(i)set(subHand(i),'Box','on'), 1:n);
arrayfun(@(i)set(subHand(i).XAxis,'Visible','on'),2:n)
arrayfun(@(i)set(subHand(i),'XTick',[]),2:n)
set(subHand,'YTick',[])
set(subHand,'XLim',xl)
% ylabel('%EL');
subHand(1).XLabel.String = xAxisLabel;
subHand(1).XLabel.FontSize = 10;
subHand(1).FontSize = 10;

x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
grid off;
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w'); 
end

function [slope, incpt] = mySlope(X, Y)
    P = polyfit(X,Y,1);    
    slope = P(1);
    incpt = P(2);
end

function [ph, lambda, rSq, slope] = scatterFitPlot(yVal, xVal, cutOff, color, fitFlag)
xVal = xVal(~isinf(yVal));
yVal = yVal(~isinf(yVal));
yVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
xVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];

color = color./255;

% [bootSlopes, bootIncpt] = bootstrp(100,@mySlope, xVal, yVal);
% bootSlopeLowCI = prctile(bootSlopes, 95);
% bootSlopeHiCI = prctile(bootSlopes, 5);
% slopeMean = mean(bootSlopes);
% slopeStd = (bootSlopeHiCI - bootSlopeLowCI)./2;
% slope = [slopeMean, slopeStd];
% lambda = [-(1/slopeMean), slopeStd/slopeMean^2];
% 
% bootIncptLowCI = prctile(bootIncpt, 95);
% bootIncptHiCI = prctile(bootIncpt, 5);
% 
% fitVar = zeros(2);
% fitVar(2,1) = slopeMean;
% fitVar(2,2) = slopeStd;

% dataTable = table(xVal, yVal);
% mdl = fitlm(dataTable);
% rSq = mdl.Rsquared.Adjusted;
% hold on;
% fitVar = mdl.Coefficients.Variables;
% lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
% slope = [fitVar(2, 1), fitVar(2,2)];


%%%%%% use for line fit r^2
if strcmp(fitFlag, 'scatNoFit')
    s1 = scatter(xVal, log(yVal), 'filled');
    s1.MarkerEdgeAlpha = 0;
    s1.MarkerFaceColor = color;
    s1.MarkerFaceAlpha = 0.2;
    s1.SizeData = 8;
    ph = s1;
end

bins = 17;
minLen = 3;
alpha = 50;

[xMean, xStd, yMean, yStd] = getBinMeans(xVal, yVal, bins, minLen);

% [xMean, xStd, yMean, yStd] = getBinMedian(xVal, yVal, bins, minLen);

dataTable = table(xMean', yMean');

%%%% use for lambda for error
mdl = fitlm(dataTable);
rSq = mdl.Rsquared.Adjusted;
hold on;
fitVar = mdl.Coefficients.Variables;
lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
slope = [fitVar(2, 1), fitVar(2,2)];


if strcmp(fitFlag, 'lineNoFit') 
    pe = plot(xMean, yMean);    
    pe.LineWidth = 1;
    pe.Color = [color, alpha/255];
    pe.LineStyle = '-';
    hold on;
    ph = pe;
end

if strcmp(fitFlag, 'errNoFit')
    pe = errorbar(xMean, yMean, yStd, yStd, xStd, xStd);    
    pe.Marker = 'none';
    pe.MarkerFaceColor = color;
    pe.LineWidth = 1;
    pe.Color = color;
    pe.LineStyle = 'none';
    pe.CapSize = 0;
    set([pe.Bar, pe.Line], 'ColorType', 'truecoloralpha', 'ColorData', [pe.Line.ColorData(1:3); 255])
    hold on;
    ph = pe;
end

if strcmp(fitFlag, 'justFit')
    xArr = linspace(0, max(xVal));
    yArr = fitVar(2, 1).*xArr + fitVar(1, 1);
    pl = plot(xArr, yArr);
    pl.LineStyle = '-';
    pl.LineWidth = 1;
    pl.Color = color;
    hold on;
    ph = pl;
end

% str=sprintf('R^{2} = %1.2f',rSq);
% % T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
% T = text(50, max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  
end

function [xMeanMean, xStdMean, yMeanMean, yStdMean] = getBinMeans(xVal, yVal, bins, minLen)
nBoot = 100;
[~, allIdx, ~] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));

xMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xBin, 'un', 0); % mean of bins
xMeanMean = cellfun(@mean, xMean);
xStd = cellfun(@(x) bootstrp(nBoot,@std,x), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@mean, xStd);

yMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yBin, 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yMean);
yStd = cellfun(@(x) bootstrp(nBoot,@std,x), yBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, yStd);
end

function [xMedianMean, xStdMean, yMedianMean, yStdMean] = getBinMedian(xVal, yVal, bins, minLen)
nBoot = 100;
[~, allIdx, ~] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));

xMedian = cellfun(@(x) bootstrp(nBoot,@median,x), xBin, 'un', 0); % mean of bins
xMedianMean = cellfun(@mean, xMedian);

xStd = cellfun(@(x) bootci(nBoot,{@median,x}, 'Alpha',0.34), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@(x) (x(2)-x(1)/2), xStd);

yMedian = cellfun(@(x) bootstrp(nBoot,@median,x), yBin, 'un', 0); % mean of bins
yMedianMean = cellfun(@mean, yMedian);
xStd = cellfun(@(x) bootci(nBoot,{@median,x}, 'Alpha',0.34), xBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, xStd);
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


function sh = histScatPlot(xVal, yVal, color)
color = color./255;
% yVal2 = 1./yVal.^3;
% figure('color', 'w')
sh = scatterhist(xVal,yVal, 'Direction','out','Kernel','on','LineWidth',1.5,'Marker','none','Color', color);
hold on;
s1 = scatter(xVal, yVal, 'filled');
s1.MarkerEdgeAlpha = 0;
s1.MarkerFaceColor = color;
s1.MarkerFaceAlpha = 0.2;
s1.SizeData = 8;
dataTable = table(xVal, yVal);
mdl = fitlm(dataTable);
fitVar = mdl.Coefficients.Variables;
slope = [fitVar(2, 1), fitVar(2,2)];
xArr = linspace(0, max(xVal));
yArr = fitVar(2, 1).*xArr + fitVar(1, 1);
pl = plot(xArr, yArr);
pl.LineStyle = '--';
pl.LineWidth = 2;
pl.Color = [0 0 0 0.8];
hold on;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  

end



function benPlot(X, Y)
sx = 100;
sy = 100;
% N = 40;
% f = ksdensity([X,Y],[X,Y],'Bandwidth',[sx,sy]/N,'Function','pdf');

%normalization
% here I actually normalize by the max for vizualization purpose.
%x1,x2 range for X
% [~,I] = histc(X,linspace(x1,x,N+1));
N = 40;
f = ksdensity([X,Y],[X,Y],'Bandwidth',[sx,sy]/N,'Function','pdf');

%normalization
% here I actually normalize by the max for vizualization purpose.
%x1,x2 range for X
[~,I] = histc(X,linspace(min(X),max(X),N+1));
for i=1:N
    Ib = I==i;
    f(Ib) = (f(Ib)-min(f(Ib)))/(max(f(Ib))-min(f(Ib)));
    % f(Ib) = f(Ib)/sum(f(Ib)); %this should give P(Y|X)
end
f(isnan(f)) = 1; %%%% delete if not needed
%plot scatter
Nc = 64;
cmap = cool(Nc);

[f,Ib] = sort(f);
scatter((X(Ib)),Y(Ib),5,cmap(round(1+f*(Nc/2-1)),:),'filled')
hold on;

end