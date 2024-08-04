function combine2CSpotPlotter4Mod(folderPath)
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x)load(append(folderPath, filesep, x)), fileNames, 'un', 0);

% distLim =  [0.45 0.45 0.45 0.45 0.45 0.45 0.45];
distLim = [0.39 0.37];
nucValUB = 400;
nucValLB = 100;

minVolFactor = 3*3*1;%4*4*3;
minVol = minVolFactor*0.043*0.043*0.2;

nucValCutOff = [200 400]; % Lower, Upper

totalKMeans = 10;
kMeansDistLimit = 1000;
pruneBins = [1, 2, 10];
colorStruct = cell(1, 9);

% colorStruct{1} =  [132, 36, 48; 0, 0, 0];
% colorStruct{2} =  [32, 12, 158; 0, 0, 0];

colorStruct{1} =  [237, 130, 177; 0, 0, 0]; % bottleneck
colorStruct{2} =  [32, 214, 208; 0, 0, 0]; % hunchback
% colorStruct{1} =  [32, 214, 208; 0, 0, 0]; % hunchback
% 
colorStruct{1} =  [32, 214, 208; 0, 0, 0];
colorStruct{2} =  [219, 122, 103; 0, 0, 0];
colorStruct{3} =  [242, 189, 148; 0, 0, 0];
colorStruct{4} =  [133, 128, 209; 0, 0, 0];
colorStruct{5} =  [150, 150, 150; 0, 0, 0];

colorStruct{1} =  [83, 185, 234; 0, 0, 0]; % strong p2
colorStruct{2} =  [192, 95, 211; 0, 0, 0]; % weak p2

%------- Manual Reordering of data Sequence --------%

groupOrder = [1 2];
% groupOrder = [3 2 5 4 1];
groupOrder = 1:length(struct2C); % default
%------------------------------------------------------------------%

geneName = strings(1, length(struct2C));
% geneNameStruct = {'Strong', 'Weak'};
groupColor = zeros(length(struct2C), 3);

c1ValAll = cell(1, length(struct2C));
c1ValMean = cell(1, length(struct2C));
c1VolAll = cell(1, length(struct2C));
c1VolMean = cell(1, length(struct2C));
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
    c1Val{i} = vertcat(struct2C{j}.TFSpotProp.c1ValSort{:}); % split to embryos  
    c1TotVal{i} = vertcat(struct2C{j}.TFSpotProp.c1ValSort{:}); % split to embryos  
    c1ValMean{i} = cell(1, length(c1Val{i}));
    c1VolMean{i} = cell(1, length(c1Vol{i}));
    c1TotValMean{i} = cell(1, length(c1Vol{i}));
    c1Dist{i} = vertcat(struct2C{j}.TFSpotProp.distSortUM{:}); % split to embryos
    for m=1:length(c1Val{i}) % total nuclei in each embryo
        c1volUnfiltAll{i} = vertcat(c1volUnfiltAll{i}, vertcat(c1Vol{i}{m}{:}));
        c1Vol{i}{m} = cellfun(@(x,y) x(y), c1Vol{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0);   
        c1VolMean{i}{m} = cellfun(@(x) mean(x, 'omitnan'), c1Vol{i}{m}, 'un', 0);
        c1VolAll{i} = vertcat(c1VolAll{i}, vertcat(c1Vol{i}{m}{:}));
        c1VolClose1{i} = vertcat(c1VolClose1{i}, cellfun(@(x) x(1), c1Vol{i}{m}(cellfun(@(x) ~isempty(x), c1Vol{i}{m})))');
        
        c1Val{i}{m} = cellfun(@(x,y) x(y), c1Val{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0);        
        c1ValMean{i}{m} = cellfun(@(x) mean(x, 'omitnan'), c1Val{i}{m}, 'un', 0);        
        c1ValAll{i} = vertcat(c1ValAll{i}, vertcat(c1Val{i}{m}{:}));
        c1ValClose1{i} = vertcat(c1ValClose1{i}, cellfun(@(x) x(1), c1Val{i}{m}(cellfun(@(x) ~isempty(x), c1Val{i}{m})))');    

        c1TotVal{i}{m} = cellfun(@(x,y)  x.*y, c1Val{i}{m},c1Vol{i}{m}, 'un', 0);
        c1TotValMean{i}{m} = cellfun(@(x) mean(x), c1TotVal{i}{m}, 'un', 0); 
        
        c1Dist{i}{m} = cellfun(@(x,y) x(y), c1Dist{i}{m}, cellfun(@(x) x >= minVol, c1Vol{i}{m}, 'un', false), 'un', 0); 
        distAll{i} = vertcat(distAll{i}, vertcat(c1Dist{i}{m}{:}));
        c1DistClose1{i} = vertcat(c1DistClose1{i}, cellfun(@(x) x(1), c1Dist{i}{m}(cellfun(@(x) ~isempty(x), c1Dist{i}{m})))');
        
        numSpotsAll{i} = vertcat(numSpotsAll{i}, nonzeros(cellfun(@length, c1Val{i}{m}')));
        if ~isempty(vertcat(nucVal{i}{m}{:}))
            try
                nucRepTemp = arrayfun(@(x,y) repelem(x, y), vertcat(nucVal{i}{m}{:}), nonzeros(cellfun(@length, c1Val{i}{m})'), 'un', 0);
            catch
                aaa = 1
            end
            nucRepTemp = cellfun(@transpose, nucRepTemp, 'un', 0);
            nucRepTemp = vertcat(nucRepTemp{:});
            nucValRepAll{i} = vertcat(nucValRepAll{i}, nucRepTemp);
        else 
            nucValRepAll{i} = NaN;
        end
        
    end

    nucVal{i} = horzcat(nucVal{i}{:});
    nucVal{i} = vertcat(nucVal{i}{:});
    c2Val{i} = horzcat(c2Val{i}{:});
    c2Val{i} = vertcat(c2Val{i}{:});
    c1ValMean{i} = horzcat(c1ValMean{i}{:});
    c1ValMean{i} = vertcat(c1ValMean{i}{:});
    c1ValMean{i}(isnan(c1ValMean{i})) = [];
    c1ValMeanTrim{i} = c1ValMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));

    c1VolMean{i} = horzcat(c1VolMean{i}{:});
    c1VolMean{i} = vertcat(c1VolMean{i}{:});
    c1VolMean{i}(isnan(c1VolMean{i})) = [];
    c1VolMeanTrim{i} = c1VolMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));

    c1TotValMean{i} = horzcat(c1TotValMean{i}{:});
    c1TotValMean{i} = vertcat(c1TotValMean{i}{:});
    c1TotValMean{i}(isnan(c1TotValMean{i})) = [];
    c1TotValMeanTrim{i} = c1TotValMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));

    nucValTrim{i} = nucVal{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c2ValTrim{i} = c2Val{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    c1ValMeanTrim{i} = c1ValMean{i}(nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    try
    coupledDist{i} = c1DistClose1{i}(c1DistClose1{i}<=distLim(i) & nucVal{i}>nucValCutOff(1) & nucVal{i}<=nucValCutOff(2));
    catch
        aaa = 1
    end
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

    nucValRepAllTrim{i} = nucValRepAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));
    c1ValAllTrim{i} = c1ValAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));
    c1VolAllTrim{i} = c1VolAll{i}(nucValRepAll{i}>nucValCutOff(1) & nucValRepAll{i}<=nucValCutOff(2));

end


% figure('Color', 'w');
yyaxis right;

data = c1DistClose1{i};
% data = cellfun(@(x,y) x./y, c2Val, nucVal, 'un', 0);
totalBins =25;
colorTemp = groupColor(1,:)./255;
set(gca, 'YColor',colorTemp)
ff = fitgmdist(data,2);
h1 = histogram(data, 0:0.1:2,'Normalization','pdf');
hold on;
h2 = histogram(data, 0:0.1:2, 'Normalization','pdf', 'DisplayStyle', 'stairs', EdgeColor=colorTemp);
hold on;
h1.FaceColor = colorTemp;
h1.FaceAlpha = 0.25;
h2.EdgeColor = colorTemp;
h1.EdgeAlpha = 0.0;
[xMin, xMax] = bounds(data);
binEdges = linspace(xMin, xMax, totalBins)';
% hold on; plot(binEdges,pdf(ff,binEdges),'r-'); hold off

n1 = makedist('normal',ff.mu(1),sqrt(ff.Sigma(1)));
n2 = makedist('normal',ff.mu(2),sqrt(ff.Sigma(2)));
pp = ff.ComponentProportion;
yy = pp(1)*pdf(n1,binEdges) + pp(2)*pdf(n2,binEdges);
hold on; 
p1 = plot(binEdges,pp(1)*pdf(n1,binEdges),"Color", colorTemp, "LineStyle",'-'); 
hold on; 
p2 = plot(binEdges,pp(2)*pdf(n2,binEdges),"Color", colorTemp, "LineStyle",'-'); 
hold on;
plot(binEdges,yy,"Color", colorTemp, "LineStyle",':', "Marker", "none", "LineWidth",1.5); 
hold on;
plot([0.45 0.45], [0 2.5], 'k--');
hold on;
plot([0.45-0.02 0.45-0.02], [0 2.5], 'k:');
hold on;
plot([0.45+0.02 0.45+0.02], [0 2.5], 'k:');
xlim([0 1])
xlabel('Distance from mRNA center ({\mu}m)')
ylabel('Distance of nearest cluster ({\mu}m)')
hold off

x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);



%-------------------------------------------------------
figure('Color', 'w');
dataCombine = vertcat(c1DistClose1{:});
groupLength = cellfun(@numel, c1DistClose1);
axLabel = "Distance of nearest cluster ({\mu}m)";
plotHist2(dataCombine, geneName, groupLength, groupColor, axLabel, distLim(1));
%-------------------------------------------------------
%-------------------------------------------------------
figure('Color', 'w');
data = cellfun(@(x,y) x./y, c2Val, nucVal, 'un', 0);

boxPlotter1(data{1}, data{2}, groupColor, "light");

dataCombine = vertcat(data{:});
groupLength = cellfun(@numel, data);
axLabel = "Transcription hotspot intensity";
plotHist2(dataCombine, geneName, groupLength, groupColor, axLabel, distLim(1));

% dataCombine = vertcat(data{3}, data{4});
% groupLength = [numel(data{3}), numel(data{4})];
% axLabel = "Transcription hotspot intensity";
% plotHist2(dataCombine, geneName(3:4), groupLength, groupColor(3:4, :), axLabel, distLim(1));
%-------------------------------------------------------

%=======================================================
% figure('Color', 'w')
% tempColor = [37, 189, 219; 16, 86, 99];
% h1 = simpleHist(distAll{1}, tempColor(1,:));
% hold on;
% h2 = simpleHist(coupledDist{1}, tempColor(2,:));
% 
% legTex = ["All clusters"; "Coupled clusters"];
% legend([h1;h2], legTex, 'Box', 'off');
% 
% ylabel('Probability');
% xlabel('Distance from mRNA center ({\mu}m)');
% ylim([0 0.2]);
% xlim([0 4]);
%=======================================================


%=======================================================
% figure('Color', 'w')
% tempColor = [37, 189, 219; 16, 86, 99];
% h1 = simpleHist(c1ValAll{1}, tempColor(1,:));
% hold on;
% h2 = simpleHist(coupledVal{1}, tempColor(2,:));
% 
% legTex = ["All clusters"; "Coupled clusters"];
% legend([h1;h2], legTex, 'Box', 'off');
% 
% ylabel('Probability');
% xlabel('Cluster intensity (a.u.)');
% ylim([0 0.2]);
%=======================================================

%=======================================================
% figure('Color', 'w')
% tempColor = [37, 189, 219; 16, 86, 99];
% h1 = simpleHist(((6/pi)*c1VolAll{i}).^(1/3), tempColor(1,:));
% hold on;
% h2 = simpleHist(((6/pi)*coupledVol{i}).^(1/3), tempColor(2,:));
% 
% legTex = ["All clusters"; "Coupled clusters"];
% legend([h1;h2], legTex, 'Box', 'off');
% 
% ylabel('Probability');
% xlabel('Cluster diameter ({\mu}m)');
% ylim([0 0.2]);
%=======================================================



%=======================================================
coupledFraction = cellfun(@(x,y) 100.*numel(x)/(numel(x)+numel(y)), coupledDist, uncoupledDist);
figure('color', 'w')
X = categorical({'Coupled','Uncoupled'});
X = reordercats(X,{'Coupled','Uncoupled'});
X = categorical(geneName);
X = reordercats(X, geneName);

b = bar(X, coupledFraction);
b.BarWidth = 0.1;
b.FaceColor = 'flat';
b.FaceAlpha = 0.8;
b.BarWidth = 0.6;
b.LineStyle = 'none';
for i = 1:length(struct2C)
    b.CData(i,:) = groupColor(i,:)./255;
end
title('% coupled clusters', 'FontWeight', 'normal');
ylabel('%');
ylim([0 100]);
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
% ylabel('Relative TF cluster intensity (a.u.)');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 150;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
%=======================================================


% figure('Color', 'w');
% axLabel = 'Coupled NN cluster intensity';    
% dataCombine = vertcat(coupledVal{:});
% groupLength = cellfun(@numel, coupledVal);
% plotHist(dataCombine, geneName, groupLength, groupColor, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Uncoupled NN cluster intensity';    
% dataCombine = vertcat(uncoupledVal{:});
% groupLength = cellfun(@numel, uncoupledVal);
% plotHist(dataCombine, geneName, groupLength, groupColor, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Hb NN cluster intensity';    
% dataCombine = vertcat(coupledVal{1}, uncoupledVal{1});
% groupLength = [numel(coupledVal{1}), numel(uncoupledVal{1})];
% groupName = ["Coupled", "Uncoupled"];
% color = [132, 36, 48;32, 12, 158];
% plotHist(dataCombine, groupName, groupLength, color, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Bnk NN cluster intensity';    
% dataCombine = vertcat(coupledVal{2}, uncoupledVal{2});
% groupLength = [numel(coupledVal{2}), numel(uncoupledVal{2})];
% groupName = ["Coupled", "Uncoupled"];
% color = [132, 36, 48; 32, 12, 158];
% plotHist(dataCombine, groupName, groupLength, color, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Coupled NN cluster volume';    
% dataCombine = vertcat(coupledVol{:});
% groupLength = cellfun(@numel, coupledVal);
% plotHist(dataCombine, geneName, groupLength, groupColor, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Uncoupled NN cluster volume';    
% dataCombine = vertcat(uncoupledVol{:});
% groupLength = cellfun(@numel, uncoupledVol);
% plotHist(dataCombine, geneName, groupLength, groupColor, axLabel);
% hold on;
% 
% figure('Color', 'w');
% axLabel = 'Coupled NN cluster mol';    
% dataCombine = vertcat(coupledVol{:}).*vertcat(coupledVal{:});
% groupLength = cellfun(@numel, coupledVal);
% plotHist(dataCombine, geneName, groupLength, groupColor, axLabel);
% hold on;

% figure('Color', 'w');
% % fitScatter2(vertcat(nucValRepAllTrim{:}), vertcat(c1ValAllTrim{:}), [0.8 0.8 0.8]);
% % hold on;
% xData = vertcat(coupledNucVal{:});
% yData = vertcat(coupledVal{:});
% groupLength = cellfun(@numel, coupledVal);
% fitScatter(xData, yData, groupLength, groupColor)
% ylabel('Coupled NN cluster intensity (a.u.)')
% xlabel('Nuclear intensity (a.u.)')
% xlim([200 400 ])
% ylim([100 600 ])
% 
% figure('Color', 'w');
% fitMean(xData, yData, groupLength, groupColor, geneName);

% % xData2 = vertcat(coupledNucVal{:});
% % yData2 = vertcat(coupledVal{:});
% % 
% % groupLength2 = horzcat(length(coupledVal{1}), length(uncoupledDist{1}));
% % groupColor2 = vertcat(groupColor(1,:), groupColor(1,:)./2, [13, 110, 36]);
% % geneName2 = ["hb coupled", "hb uncoupled", "all"];
% % figure('Color', 'w');
% % fitMean(xData2, yData2, groupLength2, groupColor2, geneName2);
% % ylabel('NN distance \mu m')
% % xlabel('Nuclear intensity (a.u.)')
% 
% 
% xData2 = vertcat(vertcat(coupledNucVal{1}, uncoupledNucVal{1}));
% yData2 = vertcat(vertcat(coupledDist{1}, uncoupledDist{1}));
% % groupLength2 = horzcat(cellfun(@numel, coupledVal), cellfun(@numel, uncoupledVal), length(vertcat(nucValTrim)));
% groupLength2 = horzcat(length(coupledDist{1}), length(uncoupledDist{1}));
% groupColor2 = vertcat(groupColor(1,:), groupColor(1,:)./2, [13, 110, 36]);
% geneName2 = ["hb coupled", "hb uncoupled", "all"];
% figure('Color', 'w');
% fitMean(xData2, yData2, groupLength2, groupColor2, geneName2);
% ylabel('NN distance \mu m')
% xlabel('Nuclear intensity (a.u.)')
% 
% xData2 = vertcat(nucVal{:});
% yData2 = vertcat(c1DistClose1{:});
% groupLength = cellfun(@numel, c1DistClose1);
% geneName2 = ["hb NN", "bnk NN"];
% figure('Color', 'w');
% fitMean(xData2, yData2, groupLength, groupColor, geneName2);
% ylabel('NN distance \mu m')
% xlabel('Nuclear intensity (a.u.)')
% 
% xData2 = vertcat(nucVal{:});
% yData2 = vertcat(c1DistClose1{:});
% yData2(yData2>distLim(1)) = 0;
% yData2(yData2~=0) = 1;
% groupLength = cellfun(@numel, c1DistClose1);
% geneName2 = ["hb frac on", "bnk frac on"];
% figure('Color', 'w');
% fitMean(xData2, yData2, groupLength, groupColor, geneName2);
% ylabel('Fraction of clustering')
% xlabel('Nuclear intensity (a.u.)')
% 
% xData2 = vertcat(vertcat(coupledNucVal{1}, uncoupledNucVal{1}), coupledNucVal{1});
% yData2 = vertcat(vertcat(coupledVal{1}, uncoupledVal{1}), vertcat(c1ValMeanTrimCoupled{1}));
% % groupLength2 = horzcat(cellfun(@numel, coupledVal), cellfun(@numel, uncoupledVal), length(vertcat(nucValTrim)));
% groupLength2 = horzcat(length(coupledVal{1}), length(uncoupledVal{1}), length(vertcat(coupledNucVal{1})));
% groupColor2 = vertcat(groupColor(1,:), groupColor(1,:)./2, [13, 110, 36]);
% geneName2 = ["hb coupled", "hb uncoupled", "all"];
% figure('Color', 'w');
% fitMean(xData2, yData2, groupLength2, groupColor2, geneName2);
% -----------------------------------------------------------------
% figure('Color', 'w');
% fitScatter2(vertcat(nucValRepAllTrim{:}), vertcat(c1ValAllTrim{:}), [0.8 0.8 0.8]);
% hold on;
% xData = vertcat(uncoupledNucVal{:});
% yData = vertcat(uncoupledVal{:});
% groupLength = cellfun(@numel, uncoupledVal);
% fitScatter(xData, yData, groupLength, groupColor)
% ylabel('Uncoupled NN cluster intensity (a.u.)')
% xlabel('Nuclear intensity (a.u.)')
% xlim([200 400 ])
% ylim([100 600 ])
% hold off;

for i = 1:length(struct2C)
    nucTemp =   vertcat(coupledNucVal{i});

    coupValTemp = vertcat(coupledVal{i})./nucTemp;
    coupValTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    coupleValCrop{i} = coupValTemp;

    coupVolTemp = vertcat(coupledVol{i});
    coupVolTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    coupleVolCrop{i} = coupVolTemp;
    coupleDiaCrop{i} = ((6/pi)*coupVolTemp).^(1/3); 
    coupleTotValCrop{i} = coupleValCrop{i}.*coupleVolCrop{i};

    c1AllCoupValTemp = vertcat(c1ValMeanTrimCoupled{i})./nucTemp;
    c1AllCoupValTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    c1AllCoupValCrop{i} = c1AllCoupValTemp;

    c1AllCoupVolTemp = vertcat(c1VolMeanTrimCoupled{i});
    c1AllCoupVolTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    c1AllCoupVolCrop{i} = c1AllCoupVolTemp;
    c1AllCoupDiaCrop{i}= ((6/pi)*c1AllCoupVolTemp).^(1/3);

    c1AllCoupTotValTemp = vertcat(c1TotValMeanTrimCoupled{i})./nucTemp;
    c1AllCoupTotValTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    c1AllCoupTotValCrop{i} = c1AllCoupTotValTemp;

    nucTemp = vertcat(uncoupledNucVal{i});
    uncoupValTemp = vertcat(uncoupledVal{i})./nucTemp;
    uncoupValTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    uncoupleValCrop{i} = uncoupValTemp;  

    uncoupVolTemp = vertcat(uncoupledVol{i});
    uncoupVolTemp(nucTemp>nucValUB | nucTemp<nucValLB) = [];
    uncoupleVolCrop{i} = uncoupVolTemp;  
    uncoupleDiaCrop{i} = ((6/pi)*uncoupVolTemp).^(1/3);
    uncoupleTotValCrop{i} = uncoupleValCrop{i}.*uncoupleVolCrop{i};
end

figure('Color', 'w');
groupBarPlot(coupleValCrop, uncoupleValCrop, geneName, groupColor, '');
title('NN cluster intensity (a.u.)', 'FontWeight', 'normal')
ylim([0 inf])
figure('Color', 'w');
groupBarPlot(coupleDiaCrop, uncoupleDiaCrop, geneName, groupColor, '');
title('NN cluster diameter ({\mu} m)', 'FontWeight', 'normal');
ylim([0.35 0.5])
figure('Color', 'w');
groupBarPlot(coupleTotValCrop, uncoupleTotValCrop, geneName, groupColor, '');
title('NN cluster total intensity (a.u.)', 'FontWeight', 'normal');
ylim([-inf inf])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w');
groupBarPlot2(coupleValCrop, c1AllCoupValCrop, groupColor, '');
title('Cluster intensity (a.u.)', 'FontWeight', 'normal')
ylim([0.8 1.2])
figure('Color', 'w');
groupBarPlot2(coupleDiaCrop, c1AllCoupDiaCrop, groupColor, '');
title('Cluster diameter ({\mu} m)', 'FontWeight', 'normal');
ylim([0.4 0.45])
figure('Color', 'w');
groupBarPlot2(coupleTotValCrop, c1AllCoupTotValCrop, groupColor, '');
title('Cluster total intensity (a.u.)', 'FontWeight', 'normal');
ylim([0.03 0.07])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% figure('Color', 'w');
% fitScatter2(vertcat(nucValRepAllTrim{:}), vertcat(c1VolAllTrim{:}), [0.8 0.8 0.8]);
% hold on;
% xData = vertcat(coupledNucVal{:});
% yData = vertcat(coupledVol{:});
% groupLength = cellfun(@numel, coupledVal);
% fitScatter(xData, yData, groupLength, groupColor)
% ylabel('Coupled cluster volume ({\mu} m \^3)')
% hy = ylabel('Coupled cluster volume (${\mu}m^3$)  ','interpreter','latex');
% xlabel('Nuclear intensity (a.u.)')
% xlim([200 400 ])
% ylim([0 inf ])
% 
% 
figure('Color', 'w');

dataX = nucVal{1};%nucValTrim{1};
dataY = c1ValMean{1};%c1ValMeanTrim{1};
color = [180 180 180];
cutOff = [200 400];
% benPlot(dataX, dataY, color);
hold on;
[vhSpotFitDiaAll, ~, ~] = scatterFitPlot(dataY, dataX, cutOff, [0 0 0], 'scatNoFit');
hold on;
[vhSpotFitDiaAll, ~, ~] = scatterFitPlot(dataY, dataX, cutOff, [0 0 0], 'errNoFit');
hold on;
[~, ~, rSqSpotFitDiaAll, slopeSpotFitDiaAll] = scatterFitPlot(dataY, dataX, cutOff, [0 0 0], 'justFit');
xlim(cutOff);
hold on;

dataX = nucVal{1};%coupledNucVal{1};
dataY = c1ValClose1{1};%coupledVal{1};
color = [237, 130, 177];
cutOff = [200 400];
% benPlot(dataX, dataY, color);
hold on;
[vhSpotFitDiaAll, ~, ~] = scatterFitPlot(dataY, dataX, cutOff, [224, 121, 41], 'scatNoFit');
hold on;
[vhSpotFitDiaAll, ~, ~] = scatterFitPlot(dataY, dataX, cutOff, [224, 121, 41], 'errNoFit');
hold on;
[~, ~, rSqSpotFitDiaAll, slopeSpotFitDiaAll] = scatterFitPlot(dataY, dataX, cutOff, [224, 121, 41], 'justFit');
xlim(cutOff);
hold on;


% benPlot(vertcat(nucValRepAllTrim{:}), vertcat(c1ValAllTrim{:}).*vertcat(c1VolAllTrim{:}),[150 150 150])
% fitScatter2(vertcat(nucValRepAllTrim{:}), vertcat(c1ValAllTrim{:}), [0.8 0.8 0.8]);

hold on;
for i = 1:length(struct2C)
    xData = vertcat(coupledNucVal{i});
    yData = vertcat(coupledVal{i});%.*vertcat(coupledVol{i});
    fitScatter2(xData, yData, groupColor(i,:))
end
% xData = vertcat(coupledNucVal{:});
% yData = vertcat(coupledVal{:}).*vertcat(coupledVol{:});
% groupLength = cellfun(@numel, coupledVal);
% fitScatter(xData, yData, groupLength, groupColor)
ylabel('Coupled cluster mol (a.u.)')
xlabel('Nuclear intensity (a.u.)')
xlim([200 400 ])
ylim([0 inf ])
set(gca, 'YScale', 'log')
% 
% % % figure; histogram(coupledVol{1}.*coupledVal{1}, 'Normalization','probability'); hold on; histogram(coupledVol{2}.*coupledVal{2}, 'Normalization', 'probability');
% % figure; histogram(c1DistClose1{1}, 10, 'Normalization','probability'); hold on; histogram(c1DistClose1{2}, 10, 'Normalization', 'probability');
% % figure;
% % scatter(nucValRepAll{1}, c1ValAll{1},'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeAlpha', 0,'Marker', 'o')
% % hold on;
% % scatter(nucValRepAll{2}, c1ValAll{2},'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeAlpha', 0,'Marker', 'o')
% % hold on;
% % scatter(coupledNucVal{1}, coupledVal{1},'MarkerFaceColor', [0.0 0.0 0.7], 'MarkerEdgeAlpha', 0,'Marker', 'o')
% % hold on;
% % scatter(coupledNucVal{2}, coupledVal{2},'MarkerFaceColor', [0.7 0.7 0.0], 'MarkerEdgeAlpha', 0,'Marker', 'o')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fitScatter(xData, yData, groupLength, color)

color = color./255;
groupCombine = repelem(1:length(groupLength), groupLength);

hold on;
for i = 1:length(groupLength)
    p1 = scatter(xData(groupCombine==i), yData(groupCombine==i), 10, 'filled');
    p1.MarkerFaceColor = color(i,:);
    p1.MarkerFaceAlpha = 0.4;
    hold on;
    fitParam{i} = regressionFit(xData(groupCombine==i), yData(groupCombine==i), color(i,:));
end
hold on ;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
% ylabel('Relative TF cluster intensity (a.u.)');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend('off')
title('')
hold on;
end

function fitScatter2(xData, yData, color)
color = color./255;
totalKMeans = 5;
[groupLength, nucMeanBin, allIdx] = histcounts(xData, totalKMeans);
nucMeanBin = (nucMeanBin(1:end-1) + diff(nucMeanBin) / 2)';
[nucMeanBinSort, sortIdx] = sort(nucMeanBin);
tempIdx = allIdx;
for i = 1:totalKMeans
    allIdx(tempIdx==sortIdx(i)) = (i);
end
groupIdx = mat2cell(allIdx', 1, length(yData));
groupIdx = groupIdx';
    
for p = 1:length(nucMeanBin)
    yDataBin{p} = yData(allIdx==p);
    yDataBin{p} = yDataBin{p}(~isnan(yDataBin{p}));
    yDataMean{p} = mean(yDataBin{p});
    yDataSem{p} = std(yDataBin{p})./sqrt(length(yDataBin{p}));
end
yDataMean = vertcat(yDataMean{:});
yDataSem = vertcat(yDataSem{:});

p2 = errorbar(nucMeanBin, yDataMean, yDataSem);
p2.Marker = 'o';
p2.Color = color;
% p2.MarkerEdgeColor = 'k';
% p2.MarkerFaceColor = 'k';
p2.MarkerSize = 4;
p2.LineStyle = 'none';
p2.LineWidth = 1;
p2.CapSize = 1;
hold on;
% p1 = scatter(xData, yData, 10, 'filled');
% p1.MarkerFaceColor = color;
% p1.MarkerFaceAlpha = 0.2;
% hold on;
[~, pl] = regressionFit(vertcat(xData), vertcat(yData), color);
hold on;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
% ylabel('Relative TF cluster intensity (a.u.)');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend('off')
title('')
hold on;
% legend(pl, convertStringsToChars(groupNames), "Box","off")
end

function plotHist(dataCombine, groupNames, groupLength, colorPalette, labelText) 
% Initialization parameters
binType = 0; % (0 = bin numbers | 1 = fix bin position | 2 = auto)
totalBins = 25; % for option 1
binSize = 0.1; % for option 1
binMax = 2;
fitType = 'gamma'; % spline or poisson or gaussian or halfnormal or gamma
%..................................................................................................................
cdf = cell(1, length(groupLength));

colorPalette = colorPalette./255;
groupCombine = repelem(1:length(groupLength), groupLength);
namedGroups = categorical(groupCombine, 1:length(groupLength), groupNames);

binCenters = cell(1, length(groupLength));
binValues = cell(1, length(groupLength));
legendText = cell(1, length(groupLength));
binEdges = cell(1, length(groupLength));

for i = 1:length(groupLength)
    data = dataCombine(groupCombine==i);
    if binType == 0 % fix bin numbers
        [xMin, xMax] = bounds(data);
        binEdges{i} = linspace(xMin, xMax, totalBins);
        h(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
    elseif binType == 1 % fix bin position
        binEdges{i} = 0:binSize:binMax; 
        h(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
    else % auto
        h(i) = histogram(data, 'Normalization', 'probability');
        binEdges{i} = h(i).BinEdges;
    end
    
    h(i).FaceColor = colorPalette(i,:);
    h(i).FaceAlpha = 0.3;
    h(i).EdgeAlpha = 0.0;
    
    binCenters{i} = h(i).BinEdges(2:end)' - (h(i).BinWidth/2);
    binValues{i} = h(i).Values';

    hold on;
    
    %% For spline fit
    if strcmp(fitType, 'spline')
%         fitCurve = csaps(binCenters{i}, binValues{i}, 1);
%         fnplt(fitCurve);

        fitCurve = fit(binCenters{i}, binValues{i}, 'smoothingspline', 'SmoothingParam',0.95);
        ff(i) = plot(fitCurve);
        ff(i).LineWidth = 2;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 
        
%         dt = diff(binEdges{i});
%         Fvals = cumsum([0;binValues{i}.*dt]);
%         fitFun = spline(binEdges{i}, [0; Fvals; 0]);
%         DF = fnder(fitFun);  % computes its first derivative
%         fnplt(DF, colorPalette(i,:), 2, [0, 2]);
        hold on;

    %% For gaussian mixture fit
    elseif strcmp(fitType, 'gaussian')
        trials = 1; % change max trails here
        pd = cell(1,trials);
        AIC = zeros(1,trials);
        options = statset('MaxIter',500);
        hold on;
        for k = 1:trials
            pd{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
            AIC(k)= pd{k}.AIC;
        end
        [minAIC,numComponents] = min(AIC);
        bestModel = pd{numComponents};
        pd = bestModel;
        % pdIn = fitgmdist(data,3);
        pdf = pdf(pd,data);
        pdf = pdf*sum(h(i).Values * h(i).BinWidth); %pdf times area under the histogram curve
        [data, idx] = sort(data); %sort raw data, det indices
        pdf = pdf(idx); %sort y per those indices
        fitMean = sort(pd.mu);
        ff(i) = plot(data,pdf,'-', 'linewidth', 1);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 

        hold on;

        %% For poisson fit
    elseif strcmp(fitType, 'poisson')
        pd = cell(1,2);
        pd = fitdist(data,'Poisson');
        if binType == 0
            p = ceil(xMin):ceil(xMax);
            pdf = poisspdf(p, pd.lambda);
            pdf = pdf*sum(h(i).Values * h(i).BinWidth);
        elseif binType == 1
            p = binCenters{i};
            increment = (1/(binEdges{i}(2) - binEdges{i}(1)));
            pdf = poisspdf((p.*increment), increment*pd.lambda);
        end   
    %     pIn = scatter(p,pdfIn,20, plotColor, 'filled', 'o', 'MarkerEdgeColor',[0 0 0]);
        ff(i) = plot(binEdges{i}, pdf, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor',colorPalette(i,:));
        hold on;
    elseif strcmp(fitType, 'halfnormal')
        pd = fitdist(data,'hn');
        pdf = pdf(pd, binEdges{i});
        pdf = pdf*sum(h(i).Values * h(i).BinWidth);
        ff(i) = plot(binEdges{i},pdf);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 
    elseif strcmp(fitType, 'gamma')
        pd = gamfit(data);
        pdf = gampdf(binEdges{i}, pd(1), pd(2));
        cdf{i} = gamcdf(binEdges{i}, pd(1), pd(2));
        pdf = pdf*sum(h(i).Values * h(i).BinWidth);
%         cdf = cdf*sum(h(i).Values * h(i).BinWidth);
        ff(i) = plot(binEdges{i},pdf);
        hold on;
        fc(i) = plot(binEdges{i},cdf{i});
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = [colorPalette(i,:), 0.3]; 
        fc(i).LineWidth = 1;
        fc(i).LineStyle = '-';
        fc(i).Color = colorPalette(i,:); 
    end
    legendText{i} = groupNames(i);
end

% hold on; plot([radLim radLim], [0 1], '--', 'Color', [0.5 0.5 0.5 0.75], 'LineWidth', 1);
% hold off;


% global plot properties
axes1 = gca;
set(axes1,'FontSize',10,'LineWidth',1,'XColor',...
[0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
% xLabel = '{Distance from mRNA hotspot {\mu}m}';
yLabel = 'Cdf';
ylabel(yLabel);
xlabel(labelText);
% title (labelText);
gca;
grid off;
box off;
if binType==0
    xlim([xMin, xMax]);
elseif binType==1
    xlim([0, binMax]);
else
    xlim([0, inf]);
end
ylim([0 1])

x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 150;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend(fc, convertStringsToChars(groupNames), "Box","off")
end

function plotHist2(dataCombine, groupNames, groupLength, colorPalette, labelText, radLim) 
% Initialization parameters
binType = 0; % (0 = bin numbers | 1 = fix bin position | 2 = auto)
totalBins = 25; % for option 0

binSize = 25; % for option 1
binMax = 500; % for option 1

binSize = 0.1; % for option 1
binMax = 2; % for option 1
%..................................................................................................................
colorPalette = colorPalette./255;
groupCombine = repelem(1:length(groupLength), groupLength);
binEdges = cell(1, length(groupLength));

fitMean = cell(1, length(groupLength));

for i = 1:length(groupLength)
    data = dataCombine(groupCombine==i);
    if binType == 0 % fix bin numbers
        [xMin, xMax] = bounds(data);
        binEdges{i} = linspace(xMin, xMax, totalBins);
        h1(i) = histogram(data, binEdges{i}, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
        hold on;
        h2(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
        hold on;
        % h3(i) = histogram(data, binEdges{i}, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
    elseif binType == 1 % fix bin position
        binEdges{i} = 0:binSize:binMax; 
        h1(i) = histogram(data, binEdges{i}, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
        hold on;
        h2(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
        hold on;
        % h3(i) = histogram(data, binEdges{i}, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
    else % auto
        h1(i) = histogram(data, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
        binEdges{i} = h1(i).BinEdges;
        hold on;
        h2(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
        hold on;
        % h3(i) = histogram(data, binEdges{i}, 'Normalization', 'cdf', 'DisplayStyle', 'stairs');
    end
    
%     h1(i).FaceColor = colorPalette(i,:);
%     h1(i).FaceAlpha = 0.3;
    h1(i).EdgeColor = colorPalette(i,:);
    h1(i).EdgeAlpha = 1.0;
    
    h2(i).FaceColor = colorPalette(i,:);
    h2(i).FaceAlpha = 0.3;
    h2(i).EdgeColor = colorPalette(i,:);
    h2(i).EdgeAlpha = 0.0;
    h3(i).EdgeColor = colorPalette(i,:);
    h3(i).FaceAlpha = 0;
    hold on;
    
    trials = 2; % change max trails here
    pd = cell(1,trials);
    AIC = zeros(1,trials);
    options = statset('MaxIter',500);
    hold on;
    for k = 1:trials
        pd{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
        AIC(k)= pd{k}.AIC;
    end
    [minAIC,numComponents] = min(AIC);
    bestModel = pd{numComponents};
    pd = bestModel;
    % pdIn = fitgmdist(data,3);
    pdf1 = pdf(pd,data);
    cdf1 = ecdf(data);
    pdf1 = pdf1*sum(h2(i).Values * h2(i).BinWidth); %pdf times area under the histogram curve
    % cdf1 = cdf1*sum(h2(i).Values * h2(i).BinWidth); %pdf times area under the histogram curve
    [data, idx] = sort(data); %sort raw data, det indices
    pdf1 = pdf1(idx); %sort y per those indices
    % cdf1 = cdf1(idx);
    fitMean{i} = sort(pd.mu);

    % ff(i) = plot(data,pdf1,'-', 'linewidth', 1);
    % ff(i).LineWidth = 1;
    % ff(i).LineStyle = '-';
    % ff(i).Color = colorPalette(i,:); 
    % hold on;

    fc(i) = plot(sort(data), cdf1(2:end),'-', 'linewidth', 1);
    fc(i).LineWidth = 1;
    fc(i).LineStyle = '-';
    fc(i).Color = colorPalette(i,:); 
    hold on;
end 

%% ----------------------------------------------------------------------
% hold on; plot([radLim radLim], [0 max(h2(i).Values)], '--', 'Color', [0.1 0.1 0.1 0.9], 'LineWidth', 1);
hold on; plot([0 2], [0.5 0.5], '--', 'Color', [0.1 0.1 0.1 0.9], 'LineWidth', 1);
hold off;


%% global plot properties
axes1 = gca;
set(axes1,'FontSize',10,'LineWidth',1,'XColor',...
[0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
% xLabel = '{Distance from mRNA hotspot {\mu}m}';
yLabel = 'Frequency';
ylabel(yLabel);
xlabel(labelText);
% title (labelText);
gca;
grid off;
box off;
if binType==0
    xlim([xMin, xMax]);
elseif binType==1
    xlim([0, binMax]);
else
    xlim([0, inf]);
end
% ylim([0 1])

x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend(h2, convertStringsToChars(groupNames), "Box","off")
end

function [x, plm] = regressionFit(xData, yData, color)
tbl = table(xData, yData);
lm = fitlm(tbl,'linear');
%------------------------------------------------------
% p1 = plot(xData,yData,'ko');
% p1.Color = [0.3 0.3 0.3];
% p1.Marker = '.';
% p1.MarkerSize = 5;
% p1.Color = [150 150 150]./255;
% p1.LineStyle = 'none';
% hold on;
plm = plot(lm);
% plm(1).Color = [0.8 0.8 0.8];
% plm(1).Marker = '.';
% plm(1).MarkerSize = 5;
plm(1).Visible = 'off';
plm(2).Color = color;
plm(2).LineStyle = '-';
plm(2).LineWidth = 1;
plm(3).Color = color;
plm(3).LineStyle = '--';
plm(3).LineWidth = 1;
plm(4).Color = color;
plm(4).LineStyle = '--';
plm(4).LineWidth = 1;
plm(3).Visible = 'off';
plm(4).Visible = 'off';

% ax = gca;
% ax.FontSize = 8;
% ax.LineWidth = 1;
% box(ax,'off');
% % ylabel('Relative TF cluster intensity (a.u.)');
% ylabel('TF cluster diameter ( \mu m)');
% xlabel('Mean nuclear intensity (a. u.)');
% grid off;
% x0 = 100;
% y0= 100;
% plotWidth = 200;
% plotHeight = 200;
% set(gcf,'position',[x0,y0,plotWidth,plotHeight])
% legend('off')
% title('')
x = lm.Coefficients.Variables;
x = x(:,1:2);
end

function groupBarPlot(cVal, uVal, name, color, label)
color = color./255;

X = categorical({'Coupled','Uncoupled'});
X = reordercats(X,{'Coupled','Uncoupled'});

avg = [cellfun(@mean, cVal); cellfun(@mean, uVal)];
b = bar(X, avg, 'grouped');
hold on
[ngroups,nbars] = size(avg);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    b(i).BarWidth = 0.1;
    b(i).FaceColor = color(i,:);
    b(i).FaceAlpha = 0.8;
    b(i).BarWidth = 0.6;
    b(i).LineStyle = 'none';
%     for j = 1:length(name)
%         b(i).CData(j,:) = color(i,:); % Color for first data coloumn
%     end
    x(i,:) = b(i).XEndPoints;
end
err = [cellfun(@(x) std(x)./sqrt(numel(x)), cVal); cellfun(@(x) std(x)./sqrt(numel(x)), uVal)];
er = errorbar(x',avg, err,'k','linestyle','none', 'CapSize',0, 'Color',[0.2 0.2 0.2]);
ylabel(label)

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
% ylabel('Relative TF cluster intensity (a.u.)');
grid off;
x0 = 100;
y0= 100;
plotWidth = 350;
plotHeight = 100;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend('off')
title('')
hold on;
end

function fitMean(xData, yData, groupLength, groupColor, names)
totalKMeans = 10;
cutoff = 5;
% [allLength, xBin, allIdx] = histcounts(xData, totalKMeans);
% xBin = (xBin(1:end-1) + diff(xBin) / 2)';
% [xBinSort, sortIdx] = sort(xBin);
% tempIdx = allIdx;
% for i = 1:totalKMeans
%     allIdx(tempIdx==sortIdx(i)) = (i);
% end
% binIdx = mat2cell(allIdx', 1, groupLength);
% 
% groupIdx = ones(length(xData), 1);


[binCount, binEdge, binIdx] = cellfun(@(x) histcounts(x, totalKMeans), mat2cell(xData', 1, groupLength), 'un', 0);


figure('color', 'w')
hold on;
for i = 1:length(groupLength)
    yMean{i} = groupsummary(yData(sum(groupLength(1:i-1))+1:sum(groupLength(1:i))), binIdx{i}', "mean");
    yStd{i} = groupsummary(yData(sum(groupLength(1:i-1))+1:sum(groupLength(1:i))), binIdx{i}', "std");
    yLen{i} = groupsummary(yData(sum(groupLength(1:i-1))+1:sum(groupLength(1:i))), binIdx{i}', "nnz");
    ySem{i} = yStd{i}./sqrt(yLen{i});
    xMean{i} = groupsummary(xData(sum(groupLength(1:i-1))+1:sum(groupLength(1:i))), binIdx{i}', "mean");
    ph(i) = plotErr(yMean{i}(yLen{i}>cutoff), ySem{i}(yLen{i}>cutoff), xMean{i}(yLen{i}>cutoff), groupColor(i,:));
    hold on;
end
legend(ph, names, "Box", 'off');
end

function [plotHandle] = plotErr(val, err, xArr, color)
color = color./255;
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
%------------------------------------
% patchTop = val+err;
% patchTop = reshape(patchTop,1,[]);
% patchBot = val-err;
% patchBot = reshape(patchBot, 1, []);
% yPatch=[patchBot,fliplr(patchTop)];
% xPatch = [xArr',fliplr(xArr')];
% pt = patch(xPatch, yPatch, 1);
% pt.FaceColor = color(1,:);
% pt.EdgeColor = 'none';
% pt.FaceAlpha = 0.3;
% hold on;
% pl = plot(xArr, val);
% pl.LineWidth = 1;
% pl.Color = color(1,:);
% pl.LineStyle = '-';

%------------------------------------
pl = errorbar(xArr, val, err);
pl.Marker = 'o';
pl.Color = color(1,:);
pl.LineStyle = '--';
pl.CapSize = 0;
hold on;
% polyOrder = 4;
% f = polyfit(xArr, val, polyOrder);
% pf = plot(xArr, polyval(f, xArr));
%------------------------------------

% hLeg = legend(pl, geneName);
ylabel('Cluster intensity (a.u.)');
xlabel('Nuclear intensity (a.u.)');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1.0;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
% xlim([0 1.5]);
% ylim([0 0.5]);
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
plotHandle = pl(1);
end

function handle = simpleHist(data, color)
color = color./255;
h1 = histogram(data);
hold on;
h1s = histogram(data, 'Normalization', 'probability', 'DisplayStyle','stairs', 'EdgeColor',[0 0 0]);
h1.Normalization = 'probability';
h1.DisplayStyle = 'bar';
h1.FaceColor = color;
h1.FaceAlpha = 0.6;
h1.EdgeColor = 'none';
handle = h1;

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');

grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 150;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end

function groupBarPlot2(cVal, cAllVal, color, label)
allVal = vertcat(cAllVal{:});
color = color./255;
color(end+1, :) = [0.4 0.4 0.4];
avg = cellfun(@mean, cVal)';
avg = [avg; mean(allVal)];
err = cellfun(@std, cVal)./cellfun(@(x) sqrt(numel(x)), cVal);
err = err';
err = [err; std(allVal)./sqrt(numel(allVal))];
xArr = 1:length(avg); % don't change order
b = bar(xArr, diag(avg), 'stacked');
hold on
for i = 1:length(avg)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = color(i,:);
    b(i).FaceAlpha = 0.8;
    b(i).BarWidth = 0.6;
    b(i).LineStyle = '-';
end

err = errorbar(xArr,avg, err,'k','linestyle','none', 'CapSize',0, 'Color',[0.2 0.2 0.2]);

xticklabels({'Hunchback','Bottleneck','All'});

ylabel(label)

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
% ylabel('Relative TF cluster intensity (a.u.)');
grid off;
x0 = 100;
y0= 100;
plotWidth = 350;
plotHeight = 100;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend('off')
title('')
hold on;
end


function boxPlotter1(dataA, dataB, color, theme)
colorPalette = color./255;
% Boxplots of distribution fit of Spots/UnitArea/Time with wilcoxon ranksum test values
if strcmp(theme, 'light')
    figure('Color', [1 1 1]);
elseif strcmp(theme, 'dark')
        figure('Color', [0 0 0]);
else
    figure('Color', [1 1 1]);
end

dataCombine = [dataA; dataB];
g1 = repmat({'Inside'},length(dataA),1);
g2 = repmat({'Outside'},length(dataB),1);
g = [g1; g2];

x1=ones(length(dataA), 1).*(1+(rand(length(dataA), 1)-0.5)/5);
x2=ones(length(dataB), 1).*(1+(rand(length(dataB), 1)-0.5)/10);
f1=scatter(x1, dataA, 'filled');
f1.MarkerFaceColor = colorPalette(2,:);
f1.MarkerFaceAlpha = 0.1;
hold on ;
f2=scatter(x2.*2,dataB,'filled');
f2.MarkerFaceColor = colorPalette(1,:);
f2.MarkerFaceAlpha = 0.2;
hold on;

boxplot(dataCombine,g, 'Notch','off', 'Whisker',1, ...
    'ColorGroup',g, 'LabelOrientation', 'horizontal', 'Symbol', '');
set(findobj(gca,'type','line'),'linew',2);
set(findobj(gcf, 'type', 'line', 'Tag', 'Median'),'Color', [0.4, 0.4, 0.4]);

set(findobj('-regexp','Tag','(Lower|Upper) (Whisker|Adjacent Value)'),'Color',[0.4, 0.4, 0.4]);

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colorPalette(j,:),'FaceAlpha',.2);
    set(h(j),'LineWidth',1);
    set(h(j),'MarkerSize',10);
    x = get(h(j),'XData');
    y = get(h(j),'YData');
    c = get(h(j),'Color');
    l = get(h(j),'LineWidth');
    ht = y(2)-y(1);
    wd = x(3)-x(1);
    rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',colorPalette(j,:),'LineWidth',l)
end
delete(h);
hold on;

dataMean = [mean(dataA), mean(dataB)];
dataStd = [std(dataA), std(dataB)];
dataSem = [std(dataA)/sqrt(length(dataA)), std(dataB)/sqrt(length(dataB))];
eb = errorbar(dataMean, dataStd, '-','MarkerSize',20,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none', 'CapSize',10);
hold on;

ylabel('mRNA Intensity');
textString = ['p = ', num2str(ranksum(dataA, dataB))];%, '%2.3fe%05d')];
text(1.5,double(max(dataCombine)),textString,'HorizontalAlignment','center', 'FontSize',12);

ax = gca;
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end 

function benPlot(X, Y, color)
color = color./255;
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
scatter((X(Ib)),Y(Ib),10,cmap(round(1+f*(Nc/2-1)),:),'filled')
hold on;

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
    s1 = scatter(xVal, (yVal), 'filled');
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
