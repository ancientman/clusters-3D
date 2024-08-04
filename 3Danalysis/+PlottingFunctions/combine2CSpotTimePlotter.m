function combine2CSpotTimePlotter(folderPath)
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x)load(append(folderPath, filesep, x)), fileNames, 'un', 0);

distLim =  [0.45 0.45];
% distLim = [0.67 0.67 0.67];
% distLim = [1 1];
nucValUB = 400;
nucValLB = 100;

minVolFactor = 2*1*1;
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

geneName = struct2C{1}.TFSpotProp.geneName;

c2Val = horzcat(struct2C{1}.TFSpotProp.c2Val{:}); % each cell is a nucleus (with all time info inside)
c1Val = horzcat(struct2C{1}.TFSpotProp.c1ValSort{:});
c1Vol = horzcat(struct2C{1}.TFSpotProp.c1VolSort{:});
c1Dist = horzcat(struct2C{1}.TFSpotProp.distSortUM{:});
nucVal = horzcat(struct2C{1}.TFSpotProp.nucVal{:});

c1ValClose1 = cell(1, length(c2Val));
c1VolClose1 = cell(1, length(c2Val));
c1DistClose1 = cell(1, length(c2Val));



for j = 1:length(c2Val)
    nullPoint{j} = cellfun(@(x) isempty(x), c2Val{j});
    c2Val{j}(nullPoint{j}) = {0};   
    c1Vol{j}(nullPoint{j}) = {NaN};       
    c1VolClose1{j} = cellfun(@(x) x(1), c1Vol{j}, 'un', 0);
    c1Val{j}(nullPoint{j}) = {NaN};   
    c1ValClose1{j} = cellfun(@(x) x(1), c1Val{j}, 'un', 0);
    c1Dist{j}(nullPoint{j}) = {NaN};   
    c1DistClose1{j} = cellfun(@(x) x(1), c1Dist{j}, 'un', 0);
    nucVal{j}(nullPoint{j}) = {NaN};   
end

c2ValAll = horzcat(c2Val{:});
c2ValAll = horzcat(c2ValAll{:})';

nucValAll = horzcat(nucVal{:});
nucValAll = horzcat(nucValAll{:})';

c1DistClose1All = horzcat(c1DistClose1{:});
c1DistClose1All = horzcat(c1DistClose1All{:})';

c1ValClose1All = horzcat(c1ValClose1{:});
c1ValClose1All = horzcat(c1ValClose1All{:})';

c1VolClose1All = horzcat(c1VolClose1{:});
c1VolClose1All = horzcat(c1VolClose1All{:})';

c1TotValClose1All = c1ValClose1All.*c1VolClose1All;

id = 22;
filepath = 'G:\Dropbox (Princeton)\bcd_ss_paper\thomas data\';
name = append(filepath, 'nuc_', num2str(id), '.csv');
write_mat = horzcat(50.*(1:length(vertcat(c2Val{id}{:})))', vertcat(c2Val{id}{:}), vertcat(c1DistClose1{id}{:}), vertcat(c1ValClose1{id}{:}));
writematrix(write_mat, name)
figure('color', 'w')
plot(50.*(1:length(vertcat(c2Val{id}{:}))), vertcat(c2Val{id}{:}), 'ok', 'LineStyle','--')
ylabel('mRNA intensity')
hold on;
yyaxis right
plot(50.*(1:length(vertcat(c2Val{id}{:}))), vertcat(c1DistClose1{id}{:}), 'Marker','v', 'LineStyle','--')
ylabel('Distance of nearest cluster {\mu}m')
% plot(50.*(1:length(vertcat(c2Val{id}{:}))), vertcat(c1ValClose1{id}{:}))
% ylabel('Intensity of nearest cluster')
xlabel('Time (s)')
set(gcf,'position',[100,100,500,200])
figure('color', 'w')
plot(50.*(1:length(vertcat(c2Val{id}{:}))), vertcat(c2Val{id}{:}), 'ok', 'LineStyle','--')
ylabel('mRNA intensity')
hold on;
yyaxis right
plot(50.*(1:length(vertcat(c2Val{id}{:}))), vertcat(c1ValClose1{id}{:}).*vertcat(c1VolClose1{id}{:}), 'Marker','v', 'LineStyle','--')
ylabel('Intensity of the nearest cluster')
xlabel('Time (s)')
set(gcf,'position',[100,100,500,200])
% p1 = plot(50.*(1:length(c2ValAll)), c2ValAll);
figure('color', 'w')
scatter1(c2ValAll, nucValAll, 'Nuclear intensity', [0 0 0]);
end

function [pf, lambda] = scatter1(xArr, val, label, color)
sc = scatter(xArr, val);
sc.Marker = 'o';
sc.MarkerEdgeColor = color;
sc.MarkerEdgeAlpha = 0.5;
sc.MarkerFaceColor = color;
sc.MarkerFaceAlpha = 0.5;
sc.SizeData = 2;

hold on;
f = fit(xArr(~isnan(val)), val(~isnan(val)), 'exp1', 'Exclude', []);%vertcat(xCellNan{:})<0.2);
pf = plot(f);
pf(1).Color = [color, 0.5];

pf.LineStyle = '--';
pf.LineWidth = 1;
hold on;

coeff = coeffvalues(f);
lambda = -(1/coeff(2));

ylabel(label);
xlabel('Fractional egg length x/L');
xlim([0.17, 0.6])

hLeg = legend('fit');
set(hLeg,'visible','off');

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
hold on;
end