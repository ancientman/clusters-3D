function combine2CValPlotter2(folderPath)
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 

totalKMeans = 1;
kMeansDistLimit = 10000;

colorStruct = cell(1, 5);

colorStruct{1} =  [32, 214, 208; 0, 0, 0];
% colorStruct{2} =  [150, 150, 150; 0, 0, 0];
colorStruct{2} =  [219, 122, 103; 0, 0, 0];
colorStruct{3} =  [242, 189, 148; 0, 0, 0];
colorStruct{4} =  [133, 128, 209; 0, 0, 0];
colorStruct{5} =  [237, 130, 177; 0, 0, 0];

colorStruct{5} =  [150, 150, 150; 0, 0, 0];
% colorStruct{1} = [0 0 255; 0 0 0];

%===========Manual Reordering of data Sequence =============%
%
% groupOrder =[2, 1, 3, 4]; 
% groupOrder =[5 2 3 4 1];
% groupOrder = [1 2 3];
% groupOrder = [1 2];
% groupOrder = [2, 1];
% groupOrder = [1 2];
% groupOrder = [3 2 5 4 1]; % default 5
groupOrder = 1:length(struct2C); % default

% groupOrder = 1;
%
%=========================================================%

geneName = strings(1, length(struct2C));

% % Plot the radial averages for each data set
% % -------------------------------------------------------
% f1 = figure('Color', 'w');
% hold on;
% for i = 1:length(struct2C)    
%     j = groupOrder(i);
%     geneName(i) = struct2C{j}.combine2C.gene;
%     ms2 = struct2C{j}.combine2C.c2Mol;
%     crop3D = struct2C{j}.combine2C.valCrop3D;    
%     val3D = struct2C{j}.combine2C.val3D;    
%     sem3D = struct2C{j}.combine2C.sem3D;    
%     valXY = struct2C{j}.combine2C.valXY;
%     semXY = struct2C{j}.combine2C.semXY;    
%     valZ = struct2C{j}.combine2C.valZ;
%     semZ = struct2C{j}.combine2C.semZ;    
%     rad = struct2C{j}.combine2C.rad;      
%     hold on;
%     set(0, 'CurrentFigure', f1)
%     plotHandle1(i) = plotErr(valXY(~isnan(valXY)), semXY(~isnan(valXY)), rad(~isnan(valXY)), geneName(i), colorStruct{i});
%     hold on;
% end
% set(0, 'CurrentFigure', f1)
% leg = legend(plotHandle1, geneName);
% set(leg,'color','none', 'TextColor',colorStruct{10}(1,:)./255, 'FontSize', 12);
% hold off;
%% -------------------------------------------------------


%% Plot boxes for (radial) center values
%% -------------------------------------------------------
% groupData = [];
% groupName = [];
% groupColor = [];
% for i = 1:length(struct2C)    
%     j = groupOrder(i);
%     groupLength(i) = length(struct2C{j}.combine2C.valCropXY);
%     groupData = vertcat(groupData, struct2C{j}.combine2C.valCropXY);    
%     groupName = vertcat(groupName, repmat({geneName(i)},length(struct2C{j}.combine2C.valCropXY),1));
%     groupMean(i) = mean( struct2C{j}.combine2C.valCropXY, 'omitnan');
%     groupSem(i) = std( struct2C{j}.combine2C.valCropXY, 1, 'omitnan')/sqrt(length(struct2C{j}.combine2C.valCropXY));
%     groupColor = vertcat(groupColor, colorStruct{j}(1,:));
% end
% f2 = figure('Color', 'w');
% set(0, 'CurrentFigure', f2)
% plotBox(groupData, geneName, groupMean, groupSem, groupLength, groupColor, 'XY');
% hold off;
%%---------------------------------------------------------------



% Plot data sorted by nuc values
% -------------------------------------------------------

nucMeanAll = [];
groupColor = [];

for i = 1:length(struct2C)    
    j = groupOrder(i);
    geneName(i) = struct2C{j}.combine2C.gene;
%     geneName(1) = "hb";
    groupLength(i) = length(struct2C{j}.combine2C.valCrop3D);
    nucMeanAll = vertcat(nucMeanAll, struct2C{j}.combine2C.nucMean);
    rad = struct2C{j}.combine2C.rad;
end
%~~~~~~~~~~~~~~~~~~~~~~
[allIdx, nucMeanBin] =kmeans(nucMeanAll, totalKMeans);
D = sqrt(sum((nucMeanAll - nucMeanBin(allIdx,:)).^2,2));
k = D <=kMeansDistLimit;
allIdx(~k) = NaN;
[nucMeanBinSort, sortIdx] = sort(nucMeanBin);
nucMeanBin = nucMeanBinSort;
tempIdx = allIdx;
for i = 1:totalKMeans
    allIdx(tempIdx==sortIdx(i)) = (i);
end

% [groupLength, nucMeanBin, allIdx] = histcounts(nucMeanAll, totalKMeans);
% nucMeanBin = (nucMeanBin(1:end-1) + diff(nucMeanBin) / 2)';
% [nucMeanBinSort, sortIdx] = sort(nucMeanBin);
% tempIdx = allIdx;
% for i = 1:totalKMeans
%     allIdx(tempIdx==sortIdx(i)) = (i);
% end


groupIdx = mat2cell(allIdx', 1, groupLength);

normXYGroupMeanBin =  cell(1, length(struct2C));
normXYGroupSemBin =  cell(1, length(struct2C));
diffXYGroupMeanBin =  cell(1, length(struct2C));
diffXYGroupSemBin =  cell(1, length(struct2C));
diff3DGroupMeanBin =  cell(1, length(struct2C));
diff3DGroupSemBin =  cell(1, length(struct2C));

for p = 1:length(nucMeanBin)
%     fn(p) = figure('color', 'w');
    fd(p) = figure('color', 'w');
    fc(p) = figure('color', 'w');
    fh(p) = figure('color', 'w');
end

groupColor = [];
for i = 1:length(struct2C)    
    xBreak = [10 10 10 10 10];
    j = groupOrder(i);
    groupColor = vertcat(groupColor, colorStruct{j}(1,:));
    normXY = struct2C{j}.combine2C.valAllNormXY';
    diffXY = struct2C{j}.combine2C.valAllDiffXY';
    diff3D = struct2C{j}.combine2C.valAllDiff3D';
    crop3D = struct2C{j}.combine2C.valCrop3D;
    cropXY = struct2C{j}.combine2C.valCropXY;
    nucVal = struct2C{j}.combine2C.nucMean;
    
    groupIdx{i} = groupIdx{i}';
   
    for k = 1:size(normXY,2)
        normXYGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), normXY(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        normXYGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), normXY(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        normXYGroupSemBin{i}(:, k) = normXYGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));
                
        diffXYGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diffXY(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        diffXYGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diffXY(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        diffXYGroupSemBin{i}(:, k) = diffXYGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));       
        
        diff3DGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diff3D(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        diff3DGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diff3D(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        diff3DGroupSemBin{i}(:, k) = diff3DGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));  
    end
    
     fillIdx = setdiff(1:totalKMeans, 1:size(normXYGroupMeanBin{i}, 1));
    for q = 1:length(fillIdx)
        normXYGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        normXYGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));

        diffXYGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        diffXYGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));

        diff3DGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        diff3DGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
    end
    
    for p = 1:length(nucMeanBin)
        crop3DGroupBin{i}{p}(:,1) = crop3D(groupIdx{i}==p);
        crop3DGroupBin{i}{p} = crop3DGroupBin{i}{p}(~isnan(crop3DGroupBin{i}{p}));
        
        cropXYGroupBin{i}{p}(:,1) = cropXY(groupIdx{i}==p);
        cropXYGroupBin{i}{p} = cropXYGroupBin{i}{p}(~isnan(cropXYGroupBin{i}{p}));
        
        nucGroupBin{i}{p}(:,1) = nucVal(groupIdx{i}==p);
        
        cropXYNormGroupBin{i}{p}(:,1) = cropXY(groupIdx{i}==p)./nucGroupBin{i}{p}(:,1);
        cropXYNormGroupBin{i}{p} = cropXYNormGroupBin{i}{p}(~isnan(cropXYNormGroupBin{i}{p}));
    end
    
    
%     crop3DGroupMeanBin{i} = accumarray(groupIdx{i}(~isnan(groupIdx{i})), crop3D(~isnan(groupIdx{i}), 1),[],@nanmean);
%     crop3DGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), crop3D(~isnan(groupIdx{i}), 1),[],@nanstd);
%     crop3DGroupSemBin{i} = crop3DGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));
        
    for p = 1:length(nucMeanBin)
        if any(unique(groupIdx{i}(~isnan(groupIdx{i}))) == p)
%         set(0, 'CurrentFigure', fn(p))
%         [handleBinNorm{p}(i), sigma{i}(p), valPeak{i}(p)]= plotErr(normXYGroupMeanBin{i}(p,:)', normXYGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), groupColor(i,:), xBreak(i));
%         ylim([-inf inf])
%         hold on;   
        %***********************
            set(0, 'CurrentFigure', fd(p))
            [handleBinNorm{p}(i), sigma{i}(p, :), valPeak{i}(p)]= plotErr(normXYGroupMeanBin{i}(p,:)', normXYGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), colorStruct{i}(1,:), xBreak(i));
%         [handleBinDiff{p}(i), sigma{i}(p,:), valPeak{i}(p)] = plotErr(diffXYGroupMeanBin{i}(p,:)', diffXYGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), groupColor(i,:), xBreak(i));
%         [handleBinDiff{p}(i), sigma{i}(p), valPeak{i}(p)] = plotErr2(diffXYGroupMeanBin{i}(p,:)', diffXYGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), groupColor(i,:), xBreak(i));
%             [handleBinDiff{p}(i), sigma{i}(p,:), valPeak{i}(p)] = plotErr2(diff3DGroupMeanBin{i}(p,:)', diff3DGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), groupColor(i,:), xBreak(i));
            ylim([-inf inf])
            set(gca, 'layer', 'top')
            hold on; 
        end
    end
end

% figure;
% plot(rad, struct2C{1}.combine2C.valAllNormXY, 'Color', [groupColor(1,:), 70]./255);
% hold on; plot(rad, struct2C{2}.combine2C.valAllNormXY, 'Color', [groupColor(2,:), 70]./255);
% hold on; plot(rad, normXYGroupMeanBin{1},  'Color', [groupColor(1,:), 255]./255, 'LineWidth', 1.5);
% hold on; plot(rad, normXYGroupMeanBin{2},  'Color', [groupColor(2,:), 255]./255, 'LineWidth', 1.5);
% xlim([0 1.2]);
% ylim([-0.5 2]);
% ylabel('Normalized TF intensity')
% xlabel('Distance from mRNA locus (\mu m)')
% figure;
% environmentDia = cellfun(@(x) 2.33.*x, sigma, 'un', 0);
environmentDia = cellfun(@(x) 2.*x, sigma, 'un', 0); % when the output is half od fwhm
environmentDia = vertcat(environmentDia{:});

groupColor = horzcat(colorStruct{:});
groupColor = reshape(groupColor(1,:), 3, []);
groupColor = permute(groupColor, [2,1]);

%------------------------------------------
b1 = figure('color', 'w'); 
plotBar(geneName(:), environmentDia(:,1), environmentDia(:,2), groupColor(:,:), '');
set(0, 'CurrentFigure', b1)
hold on;

btStrpEnvDia = bootstrp(100,@mean,environmentDia(2:5,1));%bootstrp(100,@mean,environmentDia(1:length(struct2C),1));
meanEnvDia = mean(btStrpEnvDia);
stdEnvDia = std(btStrpEnvDia);
plot([0, length(groupOrder)], [meanEnvDia, meanEnvDia], 'k--');
hold on;
plot([0, length(groupOrder)], [meanEnvDia+stdEnvDia, meanEnvDia+stdEnvDia], 'k:');
hold on;
plot([0, length(groupOrder)], [meanEnvDia-stdEnvDia, meanEnvDia-stdEnvDia], 'k:');
xlim([0.5, 4.5])
ylim([0, 0.8])
set(gca, 'XTick', 1:length(geneName),'XTickLabel',(geneName));
hold off;
%------------------------------------------


legText = cell(1, totalKMeans);

for i = 1:length(struct2C)  
    totalNucNames = cellfun(@num2str, cellfun(@length, crop3DGroupBin{i}, 'un', 0), 'un', 0);
    nucValNames = arrayfun(@num2str, fix(nucMeanBin), 'un', 0);
    for p = 1:totalKMeans        
%         if(isgraphics(handleBinDiff{p}(i))==0)
        if(isgraphics(handleBinNorm{p}(i))==0)
%             delete(handleBinDiff{p}(i));
            delete(handleBinNorm{p}(i));
        end
        if any(unique(groupIdx{i}(~isnan(groupIdx{i}))) == p)
%             legText{p} = [legText{p};  convertCharsToStrings([geneName{i}, '.   nuc val=', nucValNames{p}])];
            legText{p} = [legText{p};  convertCharsToStrings([geneName{i}, '   n=', totalNucNames{p}])];
%             legText{p} = convertCharsToStrings(legText{p});
        else
            legText{p} = [legText{p};  ''];
        end
    end
end

% for i = 1:length(struct2C)
%     figure('color', 'w');
%     color = groupColor(i,:);
%     color = repmat(color, totalKMeans, 1);
%     plotBar(geneName(i), environmentDia{i}(:,1), environmentDia{i}(:,2), color, nucValNames);
%     title(geneName(i));
% end


for p = 1:totalKMeans
%     set(0, 'CurrentFigure', fn(p))
%     leg = legend(handleBinNorm{p}, geneName);
%     set(leg,'color','none', 'TextColor','k', 'FontSize', 8);
%     title(num2str(fix(nucMeanBin(p))));
%     hold off;    
    %***********************    
    set(0, 'CurrentFigure', fd(p))
    set(gca, 'layer', 'top')
    handles = [];
    legs = [];
    hold on;
    pv = plot([meanEnvDia, meanEnvDia], [0 1], 'k--');
    for i = 1:length(struct2C)
%         if isvalid(handleBinDiff{p}(i))
%             handles = vertcat(handles, handleBinDiff{p}(i));
%             legs = vertcat(legs, legText{p}(i,:));
%         end
        if isvalid(handleBinNorm{p}(i))
            handles = vertcat(handles, handleBinNorm{p}(i));
            legs = vertcat(legs, legText{p}(i,:));
        end
    end
%             leg = legend(handleBinDiff{p}(i), legText{p}(i,:));
    [leg, icons] = legend(handles, legs);
    set(leg,'color','none', 'TextColor','k', 'FontSize', 8, 'Box', 'off');
    title(num2str(fix(nucMeanBin(p))));    
    hold off;
end

% set(0, 'CurrentFigure', fd(p))
% shg
%~~~~~~~~~~~~~~~~~~~~~~

%**********************

for p = 1:length(nucMeanBin)    
    groupData = [];
    groupMean = [];
    groupSem = [];
    groupLength = [];
    for i = 1:length(struct2C)
%         groupData = vertcat(groupData, crop3DGroupBin{i}{p});
%         groupMean = vertcat(groupMean, mean(crop3DGroupBin{i}{p}, 'omitnan'));
%         groupSem = vertcat(groupSem, std(crop3DGroupBin{i}{p}, 1, 'omitnan')./sqrt(crop3DGroupBin{i}{p}(~isnan(crop3DGroupBin{i}{p}))));
%         groupLength = vertcat(groupLength, length(crop3DGroupBin{i}{p}(~isnan(crop3DGroupBin{i}{p}))));
        
%         groupData = vertcat(groupData, cropXYGroupBin{i}{p});
%         groupMean = vertcat(groupMean, mean(cropXYGroupBin{i}{p}, 'omitnan'));
%         groupSem = vertcat(groupSem, std(cropXYGroupBin{i}{p}, 1, 'omitnan')./sqrt(cropXYGroupBin{i}{p}(~isnan(cropXYGroupBin{i}{p}))));
%         groupLength = vertcat(groupLength, length(cropXYGroupBin{i}{p}(~isnan(cropXYGroupBin{i}{p}))));
        
        groupData = vertcat(groupData, cropXYNormGroupBin{i}{p});
        groupMean = vertcat(groupMean, mean(cropXYNormGroupBin{i}{p}, 'omitnan'));
        groupSem = vertcat(groupSem, std(cropXYNormGroupBin{i}{p}, 1, 'omitnan')./sqrt(cropXYNormGroupBin{i}{p}(~isnan(cropXYNormGroupBin{i}{p}))));
        groupLength = vertcat(groupLength, length(cropXYNormGroupBin{i}{p}(~isnan(cropXYNormGroupBin{i}{p}))));
    end
%         set(0, 'CurrentFigure', fc(p))
%         plotBox(groupData, geneName, groupMean, groupSem, groupLength, groupColor, num2str(fix(nucMeanBin)));
%         title(num2str(fix(nucMeanBin(p))));
%         hold off;
        
        set(0, 'CurrentFigure', fh(p))
        axLabel = 'F_{0}';    
        plotHist(groupData, geneName, groupLength, groupColor, axLabel);
end

% -------------------------------------------------------

crop3DAll = [];
nucMeanAll = [];
ms2All = [];
fullXYNormAll = [];
fullXYDiffAll = [];
groupName = [];
groupColor = [];

for i = 1:length(struct2C)    
    j = groupOrder(i);
    geneName(i) = struct2C{j}.combine2C.gene;
    groupLength(i) = length(struct2C{j}.combine2C.valCrop3D);
    nucMeanAll = vertcat(nucMeanAll, struct2C{j}.combine2C.nucMean);    
    fullXYNormAll = horzcat(fullXYNormAll, struct2C{j}.combine2C.valAllNormXY);
    fullXYDiffAll = horzcat(fullXYDiffAll, struct2C{j}.combine2C.valAllDiffXY);
    crop3DAll = vertcat(crop3DAll, struct2C{j}.combine2C.valCrop3D);  
    ms2All = vertcat(ms2All, struct2C{j}.combine2C.c2Mol);
    groupName = vertcat(groupName, repmat({geneName(i)},length(struct2C{j}.combine2C.valCrop3D),1));    
    groupColor = vertcat(groupColor, colorStruct{i}(1,:));
end
fullXYDiffAll = fullXYDiffAll';
fullXYNormAll = fullXYNormAll';

binCounts(:,1) = hist(allIdx(~isnan(allIdx)),unique(allIdx(~isnan(allIdx))));
for i = 1:size(fullXYNormAll,2)
    fullXYDiffMeanBin(:, i) = accumarray(allIdx(~isnan(allIdx)), fullXYDiffAll(~isnan(allIdx), i),[],@nanmean);
    fullXYDiffStdBin(:, i) = accumarray(allIdx(~isnan(allIdx)), fullXYDiffAll((~isnan(allIdx)), i),[],@nanstd);
    fullXYDiffSemBin(:, i) = fullXYDiffStdBin(:, i)./sqrt(binCounts);
    
    fullXYNormMeanBin(:, i) = accumarray(allIdx(~isnan(allIdx)), fullXYNormAll((~isnan(allIdx)), i),[],@nanmean);
    fullXYNormStdBin(:, i) = accumarray(allIdx(~isnan(allIdx)), fullXYNormAll((~isnan(allIdx)), i),[],@nanstd);
    fullXYNormSemBin(:, i) = fullXYNormStdBin(:, i)./sqrt(binCounts);
end

figure('color', 'w');
for i =  1:size(fullXYNormMeanBin, 1)   
    handleBinNorm{i} = plotErr(fullXYNormMeanBin(i, :)', fullXYNormSemBin(i, :)', rad, num2str(nucMeanBin(i)), colorStruct{i}, xBreak(i));
    handleBinDiff{i} = plotErr(fullXYDiffMeanBin(i, :)', fullXYDiffSemBin(i, :)', rad, num2str(nucMeanBin(i)), colorStruct{i}, xBreak(i));
    hold on;
end

leg = legend(vertcat(handleBinNorm{:}), num2str(fix(nucMeanBin)));
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 8);
ylim([0, 0.5]);
hold off;
leg = legend(vertcat(handleBinDiff{:}), num2str(fix(nucMeanBin)));
set(leg,'color','none', 'TextColor', [0.1 0.1 0.1], 'FontSize', 8);
% legend(plotHandle, num2str(nucMeanValBin));
ylim([-inf, inf]);
hold off;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% f4 = figure('Color', 'w');
% set(0, 'CurrentFigure', f4)
% plotHandle4(i) = plotScat(nucMeanAll, ms2All, groupColor, groupOrder);
% % plotHandle4(i) = plotScat(crop3DAll, ms2All, groupColor, groupOrder);
% leg = legend(plotHandle4, geneName);
% set(leg,'color','none', 'TextColor', [0.5 0.5 0.5], 'FontSize', 12);
% hold off;
% % hold off;

groupData = [];
groupName = [];
groupColor = [];
for i = 1:length(struct2C)    
     j = groupOrder(i);
    groupLength(i) = length(struct2C{j}.combine2C.valCrop3D);
    groupData = vertcat(groupData, struct2C{j}.combine2C.valCrop3D);    
    groupName = vertcat(groupName, repmat({geneName(i)},length(struct2C{j}.combine2C.valCrop3D),1));
    groupMean(i) = mean(struct2C{j}.combine2C.valCrop3D, 'omitnan');
    groupSem(i) = std(struct2C{j}.combine2C.valCrop3D, 1, 'omitnan')/sqrt(length(struct2C{j}.combine2C.valCrop3D));
    groupColor = vertcat(groupColor, colorStruct{i}(1,:));
end
% p = num2str(ranksum(struct2C{1}.combine2C.valCrop3D, struct2C{2}.combine2C.valCrop3D));
% f5 = figure('Color', 'w');
% set(0, 'CurrentFigure', f5)
% plotBox(groupData, geneName, groupMean, groupSem, groupLength, groupColor, 'XY');
% % plotBar(geneName, groupMean, groupSem, groupColor, '3D');
% anovaFun(groupName, groupData, groupOrder);
% legend(plotHandle4, geneName);
hold off;
end

function y = fitFun(coeff, x)
y =  (1/(sqrt(2)*pi*coeff(2))).*exp(-0.5.*((x-coeff(1)).^2)./(2*coeff(2)^2));
% % y =  sqrt(2/pi).*(1/coeff(2)).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + sqrt(2/pi).*(1/coeff(4)).*exp(-0.5.*((x-coeff(3))./coeff(4)).^2);%coeff(3) + coeff(4).*x;%
% % y =  coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + coeff(4) + coeff(5).*x + coeff(6).*x.^2  + coeff(7).*x.^3 ;%coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2);%coeff(4) + coeff(5).*x + coeff(6).*x.^2 
y =  coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2);%coeff(4) + coeff(5).*x + coeff(6).*x.^2 
end

function [plotHandle, width, valPeak] = plotErr(val, err, xArr, ~, color, xBreak)
% xBreak = 4;
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
% patchTop = val+err;
% patchTop = reshape(patchTop,1,[]);
% patchBot = val-err;
% patchBot = reshape(patchBot, 1, []);
% yPatch=[patchBot,fliplr(patchTop)];
% xPatch = [xArr',fliplr(xArr')];
% pt = patch(xPatch, yPatch, 1);
% pt.FaceColor = color(1,:)./255;
% pt.EdgeColor = 'none';
% pt.FaceAlpha = 0.3;
% hold on;
pl = errorbar(xArr, val, err);
pl.LineWidth = 1;
pl.CapSize = 0;
pl.Color = color(1,:)./255;
pl.LineStyle = '-';
%---------------------------------------------
% hold on;
% valFit = fit(xArr, val,'smoothingspline');
% ff = plot(valFit, xArr, val);
% ff(1).LineStyle = 'none';
% ff(1).Marker = 'x';
% ff(1).MarkerEdgeColor = 'none';
% ff(1).MarkerFaceColor = color(1,:)./255;
% ff(1).MarkerSize = 12;
% ff(2).Color = color(1,:)./255;
% ff(2).LineStyle = '-';
% ff(2).LineWidth = 1;
%---------------------------------------------
% hold on;
% pd = fitdist(val,'hn');
% % pd2 = makedist('HalfNormal','mu',pd.mu,'sigma',pd.sigma);
% hnPdf = pdf('hn', xArr, pd.mu, pd.sigma);
% ff = plot(xArr,0.1.*hnPdf);
% ff.LineWidth = 1;
% ff.LineStyle = '-';
% ff.Color = color./255; 
%---------------------------------------------
hold on;
xInterp = linspace(xArr(1), xArr(end));
valInterp = interp1(xArr, val, xInterp);
fun = @fitFun;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
lb =  [];%[0 0.2 0 0.5];
ub = [];%[0 0.2 0 0.5];
x0 = [0 0.2 val(1) val(xBreak) -2 0 0];
% x0 = [0 0.2];
% x0 = [0 0.2 val(1) val(xBreak) -2 0 ];
[coeff,resnorm,residual,exitflag,output, lambda, jacob] = lsqcurvefit(fun,x0,xInterp(1:end),valInterp(1:end), lb, ub, options);
ci = nlparci(coeff,residual,'jacobian',jacob);
coeffStd = (ci(:,2) - coeff')/2;
fff1 = arrayfun(@(x) coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2), xInterp);
fff2 = arrayfun(@(x) coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2), xInterp);%coeff(4).*x+coeff(3), xInterp);%

% sigma = [coeff(2), coeffStd(2)];
sigma = [coeff(2)/sqrt(2), coeffStd(2)/sqrt(2)];
valPeak = fff1(1);
sigmaFWHM = xInterp(find(abs(fff1-max(fff1)/2)==min(abs(fff1-max(fff1)/2))));
sigmaFWHMerr = sqrt(2*log(2))*coeffStd(2)/sqrt(2);

width = [sigmaFWHM, sigmaFWHMerr];
%--------------------------------------------
% for i = 1:2
%     if i==1
%         x0 = [0 0.2 0 0.5];
%         [coeff,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(fun,x0,xInterp(55:end),valInterp(55:end), lb, ub, options);
%         valBase = fun(coeff, xInterp);
%         valMod = valInterp-valBase;
% %         valMod = val;
%     elseif i ==2
%         x0 = [0 0.2 0 0.5];
%         [coeff,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(fun,x0,xInterp(1:25),valMod(1:25),lb, ub, options);
%         sigma = coeff(2);
%         valPeak = valMod(1);
%     end
%     fitData(:,1) = linspace(xArr(1), xArr(end));
%     fitData(:,2) = fun(coeff, fitData(:,1));
%     ff = plot(fitData(:,1),fitData(:,2));
%--------------------------------------------
    fitData(:,1) = linspace(xInterp(1), xInterp(end));
    fitData(:,2) = fun(coeff, fitData(:,1));
    ff = plot(fitData(:,1),fitData(:,2));
    ff.LineWidth = 1;
    ff.LineStyle = 'none';
    ff.Color = color(1,:)./255; 
    hold on;
% end
%---------------------------------------------
% tempFitX = fitData(:,1);
% sigmaFWHM = tempFitX(find(abs(fitData(:,2)-max(fitData(:,2))/2)==min(abs(fitData(:,2)-max(fitData(:,2))/2))));

% hLeg = legend(pl, geneName);
ylabel('{(I - I_{nuc})}/{I_{nuc}}'); % ('^{F-{\F_nuc}}/_{\F_nuc}');
xlabel('r ({\mu}m)');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1.0;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
xlim([0 1.5]);
ylim([0 0.5]);

set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
% ylim([0, 0.4]); xlim([0,1.5]); % change here
% plotHandle = ff(2);
plotHandle = ff(1);
% plotHandle = pt;
plotHandle = pl;

% set(gcf, 'color', 'none');  
set(gca, 'color', 'none');
set(gca,'ycolor', [25, 25, 25]./255)
set(gca,'xcolor', [25, 25, 25]./255)
end

function [plotHandle, sigma, valPeak] = plotErr2(val, err, xArr, geneName, color, xBreak)
% xBreak = 4;
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
er = errorbar(xArr, val, err);
er.Color = color./255;
er.LineStyle = 'none';
hold on;
% pl = plot(xArr, val);
% pl.LineWidth = 1.5;
% pl.Color = color(1,:)./255;
% pl.LineStyle = '-';
%---------------------------------------------
hold on;
valFit = fit(xArr, val,'smoothingspline');
ff = plot(valFit, xArr, val);
ff(1).LineStyle = 'none';
ff(1).Marker = 'x';
ff(1).MarkerEdgeColor = 'none';
ff(1).MarkerFaceColor = color(1,:)./255;
ff(1).MarkerSize = 12;
ff(2).Color = color(1,:)./255;
ff(2).LineStyle = ':';
ff(2).LineWidth = 1;
%---------------------------------------------
% hold on;
% pd = fitdist(val,'hn');
% % pd2 = makedist('HalfNormal','mu',pd.mu,'sigma',pd.sigma);
% hnPdf = pdf('hn', xArr, pd.mu, pd.sigma);
% ff = plot(xArr,0.1.*hnPdf);
% ff.LineWidth = 1;
% ff.LineStyle = '-';
% ff.Color = color./255; 
%---------------------------------------------
hold on;
xInterp = linspace(xArr(1), xArr(end));
valInterp = interp1(xArr, val, xInterp);

[breakX,breakIndex] = min(abs(xInterp-xArr(xBreak)));
vq = interp1(xInterp(breakIndex:end),valInterp(breakIndex:end),xInterp,'makima','extrap');
hold on; 
pp = plot(xInterp, vq);
pp.LineWidth = 1;
pp.LineStyle = '-';
pp.Color = color./255; 
hold on;
roiVal = valInterp - vq;

% roiVal = valInterp - vq(breakIndex);
% hold on; plot(xInterp, vq(breakIndex).*ones(length(xInterp), 1));

fun = @fitFun;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
lb =  [];%[0 0.2 0 0.5];
ub = [];%[0 0.2 0 0.5];
% x0 = [0 0.2 roiVal(1)];
% [coeff,resnorm,residual,exitflag,output, lambda, jacob] = lsqcurvefit(fun,x0,xInterp(1:breakIndex),roiVal(1:breakIndex), lb, ub, options);

 x0 = [0 0.2 val(1) val(xBreak) -2 0 0];
[coeff,resnorm,residual,exitflag,output, lambda, jacob] = lsqcurvefit(fun,x0,xInterp(1:end),valInterp(1:end), lb, ub, options);

ci = nlparci(coeff,residual,'jacobian',jacob);
coeffStd = (ci(:,2) - coeff')/2;
fff1 = arrayfun(@(x) coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2), xInterp);

% fff1 = arrayfun(@(x) coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2), xInterp);
% fff2 = arrayfun(@(x) coeff(4) + coeff(5).*x + coeff(6).*x.^2  + coeff(7).*x.^3, xInterp);%coeff(4).*x+coeff(3), xInterp);%

sigma = [coeff(2), coeffStd(2)];
valPeak = fff1(1);
% for i = 1:2
%     if i==1
%         x0 = [0 0.2 0 0.5];
%         [coeff,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(fun,x0,xInterp(55:end),valInterp(55:end), lb, ub, options);
%         valBase = fun(coeff, xInterp);
%         valMod = valInterp-valBase;
% %         valMod = val;
%     elseif i ==2
%         x0 = [0 0.2 0 0.5];
%         [coeff,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(fun,x0,xInterp(1:25),valMod(1:25),lb, ub, options);
%         sigma = coeff(2);
%         valPeak = valMod(1);
%     end
%     fitData(:,1) = linspace(xArr(1), xArr(end));
%     fitData(:,2) = fun(x, fitData(:,1));
%     ff = plot(fitData(:,1),fitData(:,2));

hold on;
fitData(:,1) = linspace(xInterp(1), xInterp(end));
fitData(:,2) = fun(coeff, fitData(:,1));
ff = plot(fitData(:,1),fitData(:,2));
ff.LineWidth = 1;
ff.LineStyle = '-';
ff.Color = color./255; 
hold on;
% end
%---------------------------------------------


% hLeg = legend(pl, geneName);
ylabel('Relative TF enrichment');
xlabel('Distance from mRNA center ({\mu}m)');
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1.0;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 150;

set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
% ylim([0, 0.4]); xlim([0,1.5]); % change here
% plotHandle = ff(2);
plotHandle = ff(1);

% set(gcf, 'color', 'none');  
set(gca, 'color', 'none');
set(gca,'ycolor', [25, 25, 25]./255)
set(gca,'xcolor', [25, 25, 25]./255)
end

function plotBox(dataCombine, groupNames, dataMean, dataSem, groupLength, colorPalette, titleText)
colorPalette = colorPalette./255;
groupCombine = repelem(1:length(groupLength), groupLength);
t = tiledlayout(1, length(groupLength));

for i = 1:length(groupNames)
    dataMean(i) = mean(dataCombine(groupCombine==i), 'omitnan');
    dataSem(i) = std(dataCombine(groupCombine==i), 1, 'omitnan')/sqrt(groupLength(i));    
    ax(i) = nexttile;    
    b(i) = boxchart(ax(i),dataCombine(groupCombine==i),'GroupByColor',groupCombine(groupCombine==i));
    b(i).BoxWidth = 0.3;
    b(i).Notch = 'off';
    b(i).BoxFaceColor = colorPalette(i,:);
    b(i).BoxFaceAlpha = 0.4;
    b(i).WhiskerLineColor = [0.7 0.7 0.7];
    b(i).LineWidth = 1;
    b(i).MarkerColor = [1 1 1];
    b(i).MarkerStyle = '.';
    b(i).JitterOutliers = 'off';
    hold on;
    
    x = ones(groupLength(i), 1).*(1+(rand(groupLength(i), 1)-0.5)/(5));
    f(i) = scatter(x,dataCombine(groupCombine == i), 3, 'filled');
    f(i).MarkerFaceColor = colorPalette(i,:);
    f(i).MarkerFaceAlpha = 0.2;
    hold on;
    
    e(i) = errorbar(dataMean(i), dataSem(i), '.','MarkerSize',10,...
        'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
    hold on;
    
    if i ==1
        ylabel({'Fold enrighment of TF'});
        ax(i).FontSize = 8;
        ax(i).LineWidth = 1.0;
    end
    
    if i>1
        ax(i).YTick = [];
        ax(i).YAxis.Visible = 0;   
    end
    
    ax(i).YLim = [min(dataCombine) max(dataCombine)];    
    ax(i).XAxis.Visible = 0;       
    ax(i).XAxis.TickLabels = groupNames(i);       
    tx = text(1, min(dataCombine), {groupNames(i)}, 'HorizontalAlignment', 'center');
    tx.Rotation = 45;
    tx.FontSize = 8;
    xline(0.05,'b')
    set(gca,'Color','none')
    
    box off
end

% title(t, titleText, 'FontWeight', 'bold');
set(gca,'ycolor',[0.1 0.1 0.1])
x0 = 100;
y0= 100;
plotWidth = 150;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

[p, t, stats] = anova1(dataCombine, groupCombine, 'off');
p
% 
% boxplot(dataCombine, groupCombine, 'Notch','off', 'Whisker',1, ...
%     'LabelOrientation', 'horizontal', 'Symbol', '');
% set(findobj(gca,'type','line'),'linew',2);
% set(findobj(gcf, 'type', 'line', 'Tag', 'Median'),'Color', [0.5 0.5 0.5]);
% 
% set(findobj('-regexp','Tag','(Lower|Upper) (Whisker|Adjacent Value)'),'Color',[0.5, 0.5, 0.5]);
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),colorPalette(j,:),'FaceAlpha',.2);
%     set(h(j),'LineWidth',2);
%     set(h(j),'MarkerSize',10);
%     x = get(h(j),'XData');
%     y = get(h(j),'YData');
%     c = get(h(j),'Color');
%     l = get(h(j),'LineWidth');
%     ht = y(2)-y(1);
%     wd = x(3)-x(1);
%     rectangle('position',[x(1),y(1),wd,ht],'EdgeColor',colorPalette(j,:),'LineWidth',l)
%     
%     hold on;
%     x = ones(groupLength(j), 1).*(j+(rand(groupLength(j), 1)-0.5)/(5));
%     f(j) = scatter(x,dataCombine(groupCombine == j), 'filled');
%     f(j).MarkerFaceColor = colorPalette(j,:);
%     f(j).MarkerFaceAlpha = 0.2;
%     hold on;
% end
% 
% delete(h);
% hold on;
% 
% 
%     
% errorbar(dataMean, dataSem, '.','MarkerSize',10,...
%     'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
% hold on;


% 
% title(titleText, 'FontWeight', 'bold');
% ax = gca;
% ax.FontSize = 12;
% ax.LineWidth = 1.5;
% box(ax,'on');
% grid off;
% x0 = 100;
% y0= 100;
% plotWidth=350;
% plotHeight=400;

set(gcf,'position',[x0,y0,plotWidth,plotHeight])
% set(gcf, 'color', 'none');  


% ylim([-0.5,1]);

ylabel({'Fold enrichment of TF',''});
% title({'TF enrichment at transcription site ','{', titleText, '}'}, 'FontWeight', 'bold');

end 

function plotBar(dataNames, dataMean, dataSem, color, titleText)
color = color./255;
x = categorical(cellstr(dataNames));
x = reordercats(x,string(x));

% b = bar(x, dataMean);
b = bar(1:length(dataMean), dataMean);
for i = 1:length(b)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = 'flat';
    b(i).FaceAlpha = 0.4;
    b(i).BarWidth = 0.6;
    b(i).LineStyle = 'none';
    for j = 1:length(dataNames)
        b(i).CData(j,:) = color(j,:); % Color for first data coloumn
    end
end
hold on;
% er = errorbar(x,dataMean,dataSem); 
er = errorbar(1:length(dataMean), dataMean,dataSem); 
er.Color = [0.1 0.1 0.1];                            
er.LineStyle = 'none';
er.LineWidth = 1;
er.CapSize = 0;

ylabel({'Enrichment FWHM ({\mu}m)'});

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
xticklabels(titleText);

set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w');  
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
ylim([0 1]);
hold off;
end

function plotHandle = plotScat(valX, valY, color, dataOrder)
scatter(valX, valY, 3, [0.8 0.8 0.8], 'filled');
hold on;
%------------------------------------------------------
tbl = table(valX, valY);
lm = fitlm(tbl,'linear');

plm = plot(lm);
plm(1).Color = [0.3 0.3 0.3];
plm(1).Marker = '.';
plm(1).MarkerSize = 5;
plm(2).Color = 'r';
plm(2).LineStyle = '-';
plm(2).LineWidth = 2;

plm(3).Color = [0.6 0.6 0.6];
plm(3).LineStyle = '--';
plm(3).LineWidth = 1.5;
plm(4).Color = [0.6 0.6 0.6];
plm(4).LineStyle = '--';
plm(4).LineWidth = 1.5;
%------------------------------------------------------


%------------------------------------------------------
%%% Hill Fitting
%%% MAPPING: Emax = b(1),  EC50 = b(2)
% hill_fit = @(b,x)  b(1).*x./(b(2)+x);
%  
% b0 = [0.2; 0.1];    % Initial Parameter Estimates
% B = lsqcurvefit(hill_fit, b0, valX, valY);
% AgVct = linspace(min(valX), max(valX));   % Plot Finer Resolution
% 
% % plot(valX, valY, 'bp');
% hold on
% plm = plot(AgVct, hill_fit(B,AgVct));
% plm.LineStyle = '-';
% plm.LineWidth = 2;
% plm.Color = [1 0 0];
%------------------------------------------------------

ylabel('mRNA signal (a. u.)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel('Mean nuclear intensity (a.u.)');

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
box(ax,'off');

set(gcf, 'color', 'none');   
set(gca, 'color', 'none');
set(gca,'ycolor',[0.5 0.5 0.5])
set(gca,'xcolor',[0.5 0.5 0.5])

% xlim([150 450]);

grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

plotHandle = plm(1); 
end

function anovaFun(groupName, groupData, groupOrder)

groupIdx = findgroups(vertcat(groupName{:}));

[p, t, stats] = anova1(groupData, vertcat(groupName{:}), 'off');
[c,m,h] = multcompare(stats);

nSamples = length(unique(vertcat(groupName{:})));
elPerms = num2cell(nchoosek(1:nSamples, 2), 2);
pWilcox = cell(length(elPerms), 1);
[x, y] = meshgrid(1:nSamples, 1:nSamples);
imTemp = zeros(nSamples);
for i = 1:length(elPerms)
pWilcox{i} = ranksum(groupData(groupIdx == elPerms{i}(1)), groupData(groupIdx == elPerms{i}(2)));
imTemp(elPerms{i}(1), elPerms{i}(2)) = pWilcox{i};
end
im = rescale(imTemp);
imR = repelem(im, 100, 100);
figure;
imagesc(imR);
% im(imTemp==0) = NaN;
% im(imTemp<0.05) = 1;
% im(imTemp>0.05) = 0;
% pWil = cellfun(@(x, y) ranksum(x{y(1)}, x{y(2)}), yCell, elPerms, 'un', 0);
end

function plotHist(dataCombine, groupNames, groupLength, colorPalette, labelText) 
% Initialization parameters
binType = 2; % (0 = bin numbers | 1 = fix bin position | 2 = auto)
totalBins = 10; % for option 1
binSize = 0.1; % for option 1
binMax = 2;
fitType = 'spline'; % spline or poisson or gaussian or halfnormal or gamma
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
    h(i).FaceAlpha = 0.5;
    h(i).EdgeAlpha = 1;
    h(i).EdgeColor = colorPalette(i,:);
    
    binCenters{i} = h(i).BinEdges(2:end)' - (h(i).BinWidth/2);
    binValues{i} = h(i).Values';
    
    %% Plot 1
%     axes1 = gca;
%     set(axes1,'FontSize',12,'LineWidth',1.5,'XColor',...
%     [0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
%     set(h,'Parent',axes1, 'LineWidth',1.5,...
%         'EdgeColor',colorPalette(i,:), 'FaceColor',colorPalette(i,:), 'FaceAlpha', 0.4);

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

radLim = 0.5;

hold on; plot([radLim radLim], [0 1], '--', 'Color', [0.5 0.5 0.5 0.75], 'LineWidth', 2);
hold off;


%% global plot properties
axes1 = gca;
set(axes1,'FontSize',8,'LineWidth',1,'XColor',...
[0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
% xLabel = '{Distance from mRNA hotspot {\mu}m}';
yLabel = 'Normalized probability';
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
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend(ff, convertStringsToChars(groupNames))
end

