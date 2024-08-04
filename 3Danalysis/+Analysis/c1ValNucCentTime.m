function c1ValNucCentTime(txtFilePath, geneName)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   uses two colour data to calculate distance between the 
%   transciption spot and the bicoid spots.
%   also calculates the diameter of the bicoid spots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
nucStruct = cell(1, totalSubDirs);
centStruct = cell(1, totalSubDirs);
timePoints = cell(1, totalSubDirs);
nNuc = cell(1, totalSubDirs);
centTime = cell(1, totalSubDirs);
nucMeanVal = cell(1, totalSubDirs);
centTimeNoBg = cell(1, totalSubDirs);
c2Val = cell(1, totalSubDirs);
c2Vol = cell(1, totalSubDirs);
c2Mol = cell(1, totalSubDirs);
centNucTime = cell(1, totalSubDirs);
centDiffTime = cell(1, totalSubDirs);
centNormTime = cell(1, totalSubDirs);
centNorm = cell(1, totalSubDirs);
nucMeanTime = cell(1, totalSubDirs);
c2spotMolNuc = cell(1, totalSubDirs);
c1ValC2timeNucNoBg = cell(1, totalSubDirs);
c1ValC2timeMean = cell(1, totalSubDirs);
c1NucCentValRadMean = cell(1, totalSubDirs);
c1NucCentValRadNorm = cell(1, totalSubDirs);
c1ValC2timeMeanNoBg = cell(1, totalSubDirs);
c1ValC2timeSem = cell(1, totalSubDirs);
c1NucCentValRadSem = cell(1, totalSubDirs);
c1ValC2timeSemNoBg = cell(1, totalSubDirs);
c1ValMean = cell(1, totalSubDirs);
c1ValSem = cell(1, totalSubDirs);

radBin = cell(1, totalSubDirs);
c2SpotStruct = cell(1, totalSubDirs);
init = 1;

centNormAll = [];
for i=init:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        nucStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1NucPropDS.mat'));
        centStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1NucCentIntensityPropDS.mat'));
        nNuc{i} = length(nucStruct{i}.c1NucProp{1});
        centTime{i} = cell(1, nNuc{i});
        nucMeanVal{i} = cell(1, nNuc{i});
        centTimeNoBg{i} = cell(1, nNuc{i});
        centNucTime{i} = cell(1, nNuc{i});
        centDiffTime{i} = cell(1, nNuc{i});
        centNormTime{i} = cell(1, nNuc{i});
        centNorm{i} = cell(1, nNuc{i});
        nucMeanTime{i} = cell(1, nNuc{i});
        c1NucCentValRadMean{i} = cell(1, nNuc{i});
        c1NucCentValRadNorm{i} = cell(1, nNuc{i});
        c1NucCentValRadSem{i} = cell(1, nNuc{i});
        c1ValC2timeSemNoBg{i} = cell(1, nNuc{i});
        timePoints{i} = length(nucStruct{i}.c1NucProp);
        
        for n = 1:nNuc{i}       
            centTime{i}{n} = cell(1, timePoints{i});
            centTimeNoBg{i}{n} = cell(1, timePoints{i});
            for t = 1:timePoints{i}
                 if ~isempty(centStruct{i}.c1nucCentIntensityProp{n}.valDistProp{t}) && ~isempty(nucStruct{i}.c1NucProp{t})
                    centTime{i}{n}{t} = (centStruct{i}.c1nucCentIntensityProp{n}.valDistProp{t}.distValAv);
                    centTimeNoBg{i}{n}{t} = rescale(centStruct{i}.c1nucCentIntensityProp{n}.valBgSubDistProp{t}.distValAv); 
                    nucMeanVal{i}{n}{t} = nucStruct{i}.c1NucProp{t}(n).meanVal;
                 end
            end
            centNucTime{i}{n} = horzcat(centTime{i}{n}{:});
            nucMeanTime{i}{n} = horzcat(nucMeanVal{i}{n}{:});
            centDiffTime{i}{n} = centNucTime{i}{n} - nucMeanTime{i}{n};         
            centNormTime{i}{n} = bsxfun(@rdivide, centDiffTime{i}{n}, nucMeanTime{i}{n});
            centNorm{i}{n} =  mean(centNormTime{i}{n}, 2); % average over all radius bins (columns are times)            
            centNormAll = horzcat(centNormAll, centNorm{i}{n});
        end
    end
end

centNormMean = mean(centNormAll, 2, 'omitnan');
centNormSem = std(centNormAll, 0, 2, 'omitnan');
centNormSem = centNormSem./sqrt(sum(~isnan(centNormAll),2));

radAdd = (0.1*(1:size(centNormMean, 1)) - 0.1)'; % hard coded
handle = plotErr(centNormMean, centNormSem, radAdd,'center',[20, 20, 20; 0, 0, 0]);

leg = legend(vertcat(handle), append('center   n=', num2str(size(centNormAll, 2))));
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 8);
ylim([0, 0.5]);
hold off;
end

function plotHandle = plotErr(val, err, xArr, geneName, color)
color = color./255;
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
% figure('color', 'w');
% patchTop = val+err;
% patchTop = reshape(patchTop,1,[]);
% patchBot = val-err;
% patchBot = reshape(patchBot, 1, []);
% yPatch=[patchBot,fliplr(patchTop)];
% xPatch = [xArr',fliplr(xArr')];
% pt = patch(xPatch, yPatch, 1);
% pt.FaceColor = color(1,:)./255;
% pt.EdgeColor = 'none';
% pt.FaceAlpha = 0.6;
% hold on;
% pl = plot(xArr, val);
% pl.LineWidth = 1;
% pl.Color = color(1,:)./255;
% pl.LineStyle = '-';
%---------------------------------------------
pl = errorbar(xArr, val, err);
pl.Marker = 'none';
pl.Color = color(1,:);
pl.LineStyle = '-';
pl.CapSize = 0;
hold on;


% hold on;
% pd = fitdist(val,'HalfNormal');
% yFit = pdf(pd, xArr);
% plot(xArr, yFit);
%---------------------------------------------
% hLeg = legend(geneName);
ylabel('{(F - F_{nuc})}/{F_{nuc}}');%('^{F-{\F_nuc}}/_{\F_nuc}');
xlabel('r ({\mu}m)');
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
hold off;
ylim([0.05, 0.5]); xlim([0,1]);

set(gcf, 'color', 'w');  
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])

plotHandle = pl;
end


function plotErr2(val, err, xArr)
figure('color', 'w');
errorbar(xArr, val, err, 'o','MarkerSize',3,...
    'MarkerEdgeColor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3], 'LineStyle', 'none', 'LineWidth', 1.5);
hold on;
% pd = fitdist(val,'hn');
% pdFun = pdf(pd, xArr);
% pFit = pdFun/sum(0.1*pdFun);
f = fit(xArr, val, 'poly7');%vertcat(xCellNan{:})<0.2);
pf = plot(f);
pf(1).Color = [0.3, 0.3, 0.3];
pf.LineStyle = '--';
pf.LineWidth = 2;
% pl = plot(xArr, val);
% pl.LineWidth = 1.5;
% pl.Color = [0.4 0.4 0.4];
% pl.LineStyle = '--';
%---------------------------------------------

% hold on;
% pd = fitdist(val,'HalfNormal');
% yFit = pdf(pd, xArr);
% plot(xArr, yFit);
%---------------------------------------------
ylabel('Relative TF intensity (a. u.)');
xlabel('Distance from the mRNA center ({\mu}m) ');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=350;
plotHeight=350;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function plotDouble(xArr, y1, y2Mean, y2Dev)

fig = figure('Color', 'w');
leftColor = [0.7 0.3 0.4];
rightColor = [0.0 0.0 0.4];

pp = plot(xArr, y1);
pp.LineWidth = 1.5;
pp.Color = leftColor;
pp.LineStyle = '--';

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
pl.LineStyle = '--';

ylabel('Bicoid Intensity');
set(gca,'YTick',[]);
set(gca,'ycolor',rightColor) ;
hold on;
xlabel('Time (s)');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
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


