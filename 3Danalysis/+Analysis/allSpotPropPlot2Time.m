function allSpotPropPlot2Time(txtFilePath, geneName)  
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
        c1c2SpotStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1c2SpotPropDS.mat'));
        nucStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1NucPropDS.mat'));
        c1SpotPropStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1SpotPropDS.mat'));
        c2SpotPropStruct{i} = load(append(dataFiles.(fileID), filesep, 'c2SpotPropDS.mat'));
        timeFrames{i} = length(c1c2SpotStruct{i}.c1c2SpotProp);          
        c1SpotVal{i} = cell(1, timeFrames{i});
        
        for t = 1:timeFrames{i}
            if ~isempty(nucStruct{i}.c1NucProp{t})
                nucVol{i}{t} = vertcat(nucStruct{i}.c1NucProp{t}.vol);   
                cellC1C2SpotDistUM{i}{t} = c1c2SpotStruct{i}.c1c2SpotProp{t}.C1C2SpotDistUM;
                c1SpotsPerNuc{i}{t} = cellfun(@length, cellC1C2SpotDistUM{i}{t})';
                temp1 = (nucVol{i}{t}./c1SpotsPerNuc{i}{t}).^(1/3);
                temp1(isinf(temp1)) = 0;
                c1SpotsDist{i}{t} = temp1;            
                [sortC1C2spotDistUM{i}{t}, sortDistIdx{i}{t}] = cellfun(@sort, cellC1C2SpotDistUM{i}{t}, 'Uni', 0);       
                cellC1C2SpotVolUM{i}{t} = c1c2SpotStruct{i}.c1c2SpotProp{t}.C1C2SpotVolUM;        
                distSortC1C2SpotVolUM{i}{t} = cellfun(customFun1, sortDistIdx{i}{t}, cellC1C2SpotVolUM{i}{t}, 'Uni', 0);
                
                for j = 1:length(c1SpotPropStruct{i}.c1SpotProp{t})
                    if ~isempty(c1SpotPropStruct{i}.c1SpotProp{t}(j).bb)
                        c1SpotVal{i}{t}{j} = cellfun(@mean, c1SpotPropStruct{i}.c1SpotProp{t}(j).voxVal);
                    else
                        c1SpotVal{i}{t}{j} = NaN;
                    end
                end        
                
                distSortC1SpotVal{i}{t} = cellfun(customFun1, sortDistIdx{i}{t}, c1SpotVal{i}{t}, 'Uni', 0);
                
                for k=1:length(sortC1C2spotDistUM{i}{t}) % # nuclei  
                    if ~isempty(sortC1C2spotDistUM{i}{t}{k})                        
                        nucVal{i}{t}(k, 1) = vertcat(mean(nucStruct{i}.c1NucProp{t}(k).voxVal{:}));                             
                        closeDistIdx = find(distSortC1C2SpotVolUM{i}{t}{k}>minVol, 2);
                        c1c2Close1SpotDistUM{i}{t}(k,1) = sortC1C2spotDistUM{i}{t}{k}(closeDistIdx(1));
                        c1c2Close2SpotDistUM{i}{t}(k,1) = sortC1C2spotDistUM{i}{t}{k}(closeDistIdx(2));
                        c1c2Close1SpotVolUM{i}{t}(k,1) = distSortC1C2SpotVolUM{i}{t}{k}(closeDistIdx(1));
                        c1c2Close2SpotVolUM{i}{t}(k,1) = distSortC1C2SpotVolUM{i}{t}{k}(closeDistIdx(2));
                        c1c2Close1SpotVal{i}{t}(k,1) = distSortC1SpotVal{i}{t}{k}(closeDistIdx(1));
                        c1c2Close2SpotVal{i}{t}(k,1) = distSortC1SpotVal{i}{t}{k}(closeDistIdx(2));
                        
                        c2SpotVol{i}{t}(k,1) = c2SpotPropStruct{i}.c2SpotProp{t}(k).volUM;
                        c2SpotVal{i}{t}(k,1) = mean(c2SpotPropStruct{i}.c2SpotProp{t}(k).voxVal{:});
                        c2SpotMol{i}{t}(k,1) = c2SpotVal{i}{t}(k,1).*c2SpotVol{i}{t}(k,1);                    
                    else
                        nucVal{i}{t}(k, 1) = NaN;
                        c1c2Close1SpotDistUM{i}{t}(k,1) = NaN;
                        c1c2Close2SpotDistUM{i}{t}(k,1) = NaN;
                        c1c2Close1SpotVolUM{i}{t}(k,1) = NaN;
                        c1c2Close2SpotVolUM{i}{t}(k,1) = NaN;
                        c1c2Close1SpotVal{i}{t}(k,1) = NaN;
                        c1c2Close2SpotVal{i}{t}(k,1) = NaN;
                        c2SpotVol{i}{t}(k,1) = NaN;
                        c2SpotMol{i}{t}(k,1) = NaN;

                        
                    end
%                     if length(sortC1C2spotDistUM{i}{t}{k})>1
%                         c1c2Close2SpotDistUM{i}{t}(k,1) = sortC1C2spotDistUM{i}{t}{k}(find(distSortC1C2SpotVolUM{i}{t}{k}>minVol, 1 ));  
%                     else
%                         c1c2Close2SpotDistUM{i}{t}(k,1) = NaN;
%                     end
                end     
                allNucVal = vertcat(allNucVal, vertcat(nucVal{i}{t}));
                allC1C2SpotDistUM = vertcat(allC1C2SpotDistUM, vertcat(c1c2SpotStruct{i}.c1c2SpotProp{t}.C1C2SpotDistUM{:}));
                allC1C2SpotVolUM = vertcat(allC1C2SpotVolUM, vertcat(cellC1C2SpotVolUM{i}{t}{:}));
                allCIC2meanInterC1Dist = vertcat(allCIC2meanInterC1Dist, c1SpotsDist{i}{t});
                allC1C2Close1SpotDistUM = vertcat(allC1C2Close1SpotDistUM, c1c2Close1SpotDistUM{i}{t});   
                allC1C2Close2SpotDistUM = vertcat(allC1C2Close2SpotDistUM, c1c2Close2SpotDistUM{i}{t});
                allC1C2Close1SpotVolUM = vertcat(allC1C2Close1SpotVolUM,  c1c2Close1SpotVolUM{i}{t});
                allC1C2Close2SpotVolUM = vertcat(allC1C2Close2SpotVolUM,  c1c2Close2SpotVolUM{i}{t});
                allC1C2Close1SpotVal = vertcat(allC1C2Close1SpotVal,  c1c2Close1SpotVal{i}{t});
                allC1C2Close2SpotVal = vertcat(allC1C2Close2SpotVal,  c1c2Close2SpotVal{i}{t});                
                try
                allC2SpotVal = vertcat(allC2SpotVal, c2SpotVal{i}{t});
                catch
                    aaa = 1;
                end
                allC2SpotVol = vertcat(allC2SpotVol, c2SpotVol{i}{t});
                
                if t == 1
                    close1DistNucTime{i} = (c1c2Close1SpotDistUM{i}{t});
                    close1VolNucTime{i} = (c1c2Close1SpotVolUM{i}{t});
                    close2VolNucTime{i} = (c1c2Close2SpotVolUM{i}{t});
                    close1ValNucTime{i} = (c1c2Close1SpotVal{i}{t});
                    close2ValNucTime{i} = (c1c2Close2SpotVal{i}{t});
                    c2SpotVolNucTime{i} = c2SpotVol{i}{t};
                    c2SpotMolNucTime{i} = c2SpotMol{i}{t};
                else
                    try
                    close1DistNucTime{i} = horzcat(close1DistNucTime{i}, c1c2Close1SpotDistUM{i}{t});
                    close1VolNucTime{i} = horzcat(close1VolNucTime{i}, c1c2Close1SpotVolUM{i}{t});
                    close2VolNucTime{i} = horzcat(close2VolNucTime{i}, c1c2Close2SpotVolUM{i}{t});
                    close1ValNucTime{i} = horzcat(close1ValNucTime{i}, c1c2Close1SpotVal{i}{t});
                    close2ValNucTime{i} = horzcat(close2ValNucTime{i}, c1c2Close2SpotVal{i}{t});
                    c2SpotVolNucTime{i} = horzcat(c2SpotVolNucTime{i}, c2SpotVol{i}{t});
                    c2SpotMolNucTime{i} = horzcat(c2SpotMolNucTime{i}, c2SpotMol{i}{t});
                    catch
                        aaa = 1
                    end
                end
            end   
            for p = 1:size(c2SpotVolNucTime{i}, 1)
                if p == 1                                  
                    close1DistEm{i} = close1DistNucTime{i}(p,:)';
                    close1VolEm{i} = close1VolNucTime{i}(p,:)';
                    close2VolEm{i} = close2VolNucTime{i}(p,:)';
                    close1ValEm{i} = close1ValNucTime{i}(p,:)';                    
                    close2ValEm{i} =  close2ValNucTime{i}(p,:)';
                    c2SpotVolEm{i} = c2SpotVolNucTime{i}(p,:)';
                    c2SpotMolEm{i} = c2SpotMolNucTime{i}(p,:)';      
                else                    
                    close1DistEm{i} = vertcat(close1DistEm{i}, close1DistNucTime{i}(p,:)');
                    close1VolEm{i} = vertcat(close1VolEm{i}, close1VolNucTime{i}(p,:)');
                    close2VolEm{i} = vertcat(close2VolEm{i}, close2VolNucTime{i}(p,:)');
                    close1ValEm{i} = vertcat(close1ValEm{i}, close1ValNucTime{i}(p,:)');                    
                    close2ValEm{i} = vertcat(close2ValEm{i}, close2ValNucTime{i}(p,:)');
                    c2SpotVolEm{i} = vertcat(c2SpotVolEm{i}, c2SpotVolNucTime{i}(p,:)');
                    c2SpotMolEm{i} = vertcat(c2SpotMolEm{i}, c2SpotMolNucTime{i}(p,:)');
                end
            end
        end
    end    
        
end

% xArr = 93*(1:length(c2SpotMolNucTime{1}(5,:)));
% plotDouble(xArr, c2SpotMolNucTime{1}(5,:), close1DistNucTime{1}(5,:), xArr.*0)
% plotDouble(xArr, c2SpotMolNucTime{1}(1,:), close1DiaNucTime{1}(1,:), xArr.*0)


allNucVal = nonzeros(allNucVal(~isnan(allNucVal)));
allC1C2SpotDistUM = nonzeros(allC1C2SpotDistUM(~isnan(allC1C2SpotDistUM)));
allC1C2SpotVolUM = nonzeros(allC1C2SpotVolUM(~isnan(allC1C2SpotVolUM)));
allCIC2meanInterC1Dist = nonzeros(allCIC2meanInterC1Dist(~isnan(allCIC2meanInterC1Dist)));
allC1C2Close1SpotDistUM = nonzeros(allC1C2Close1SpotDistUM(~isnan(allC1C2Close1SpotDistUM)));
allC1C2Close2SpotDistUM = nonzeros(allC1C2Close2SpotDistUM(~isnan(allC1C2Close2SpotDistUM)));
allC1C2Close1SpotVolUM = nonzeros(allC1C2Close1SpotVolUM(~isnan(allC1C2Close1SpotVolUM)));
allC1C2Close2SpotVolUM = nonzeros(allC1C2Close2SpotVolUM(~isnan(allC1C2Close2SpotVolUM)));
allC1C2SpotDiaUM = power((2*allC1C2SpotVolUM), (1/3));
allC1C2Close1SpotDiaUM = power(2*(allC1C2Close1SpotVolUM), (1/3));
allC1C2Close2SpotDiaUM = power(2*(allC1C2Close2SpotVolUM), (1/3));
allC1C2Close1SpotVal = nonzeros(allC1C2Close1SpotVal(~isnan(allC1C2Close1SpotVal)));
allC1C2Close2SpotVal = nonzeros(allC1C2Close2SpotVal(~isnan(allC1C2Close2SpotVal)));
allC2SpotVal = nonzeros(allC2SpotVal(~isnan(allC2SpotVal)));
allC2SpotVol = nonzeros(allC2SpotVol(~isnan(allC2SpotVol)));

TFSpotProp.geneName = geneName; %" example: p2 all weak";
fprintf("\ngene name = %s\n\n", TFSpotProp.geneName);
TFSpotProp.nucVal = allNucVal;
TFSpotProp.close1DistUM = allC1C2Close1SpotDistUM;
TFSpotProp.close2DistUM = allC1C2Close2SpotDistUM;
TFSpotProp.close1DiaUM = allC1C2Close1SpotDiaUM;
TFSpotProp.close2DiaUM = allC1C2Close2SpotDiaUM;
TFSpotProp.close1Val = allC1C2Close1SpotVal;
TFSpotProp.close2Val = allC1C2Close2SpotVal;
TFSpotProp.c2Val = allC2SpotVal;
TFSpotProp.c2Vol = allC2SpotVol;

resultFolder = 'G:\Dropbox (Princeton)\bcd_ss_paper\data\aa_val\geneCombineSpot2';
fileName = append('\', TFSpotProp.geneName, '_DS.mat');
save(append(resultFolder, fileName),'TFSpotProp');

plotColor = [0 1 1];
plotHist(allC1C2Close1SpotDistUM, geneName, plotColor);

% Analysis.barPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1);
% Analysis.barPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2);

% Analysis.boxPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1, 'light');
% Analysis.boxPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2, 'light');

% Analysis.scatterBoxPlotter(allC1C2SpotDiaUM, allC1C2Close1SpotDiaUM, 1, 'light');
% Analysis.scatterBoxPlotter(allC1C2Close2SpotDistUM, allC1C2Close1SpotDistUM, 2, 'light');
end

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