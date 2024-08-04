function [fim, totalSpots] = adapThres3(im, nucLabel)
maxIter = 4; % nuclear threshold
technique = 1; %%%%%%%%%%%%%%%
byPass = 0; %%%%%%%%%%%%%%%
imF= im;
imF(nucLabel == 0) = 0;
imF= medfilt3(rescale(imF));
nNuc = max(nucLabel, [], 'all');

thres = cell(1, nNuc);
imThresNuc = cell(1, nNuc);
% imThresLoc = cell(1, max(nucLabel, [], 'all'));
for i = 1:nNuc 
    thres{i} = zeros(maxIter, 1);
    imThresNuc{i} = zeros([size(im), maxIter]);
%     imThresLoc{i} = zeros([size(im), maxIter]);
end

if technique == 1
%___________________________________________________
%   Calculate  multi-thresholds. technique#1: auto-otsu
imThresNuc = nucThres(nucLabel, im, maxIter, thres, imThresNuc);
if byPass == 0
    [imThresLoc, totalSpots] = locThres(nucLabel, imThresNuc); 
else
    imThresLoc = imThresNuc{1}{maxIter};
    CC = bwconncomp(imThresLoc, 8);
    totalSpots = CC.NumObjects;
end
    
%___________________________________________________
elseif technique == 2
%___________________________________________________
%   Calculate  multi-thresholds. technique#2: iterative local std        
%         imAvTempSq =  conv2(imTemp.^2, ones(3)/3^2, 'same');
%         imVar = imAvTempSq - imAvTemp.^2;
%         imVar(imVar<0) = 0;
%         imVar = rescale(imVar);
%         imDiff = imVar - imStdTemp.^2;
%         imDiff(imDiff<0) = 0;
%         imDiff = rescale(imDiff);        
end
%___________________________________________________
fim = imThresLoc;
end

%___________________________________________________
function imThres = nucThres(nucLabel, im, maxIter, thres, imThres)
imThresTemp = cell(1, maxIter);
nNuc = max(nucLabel, [], 'all');
for i= 1:nNuc
    imTemp = im;
    imTemp(nucLabel~=i) = 0;
    imTemp = rescale(imTemp);
    if maxIter>1 && maxIter<=21
        maxMultiLevel = maxIter;
        thres{i}(2:maxIter) = ...
        multithresh(imTemp, maxMultiLevel-1);
    elseif maxIter>21
        maxMultiLevel = 21;
        thresTemp(2:maxMultiLevel) = ...
            multithresh(imTemp, maxMultiLevel-1);
        thres{i}(2:maxIter) = interp1(linspace(1, 2, maxMultiLevel-1), ...
            thresTemp(2:maxMultiLevel), linspace(1, 2, maxIter-1));
    end
%     figure('Color', 'w');
%     [row, col] = find(nucLabel==i);
%     nucCrop = imcrop(imTemp, [min(col), min(row), ...
%                 (max(col) - min(col)), (max(row) - min(row))]);
    
    %   Plot histogram
%     plotImgHist(nucCrop, thres{i});
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     figure('color', 'w');
    for j=1:maxIter
        imTemp = imTemp - thres{i}(j);
        imTemp(imTemp<0) = 0;
        imTemp = rescale(imTemp);
        imThresTemp{j} = imTemp;        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %   Visualization
%         if j<=9
%             thresCrop = imcrop(imTemp, [min(col), min(row), ...
%                     (max(col) - min(col)), (max(row) - min(row))]);
%             subplot(3, 3, j)
%             imshow(thresCrop,[]);
%         end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end    
    imThres{i} = imThresTemp;
%     ax = gca;
%     ax = gca;
%     ax.FontSize = 12;
%     ax.LineWidth = 1.5;
%     box(ax,'on');
%     grid off;
%     x0 = 100;
%     y0= 100;
%     plotWidth=500;
%     plotHeight=500;
%     set(gcf,'position',[x0,y0,plotWidth,plotHeight])
%     hold off;
end
end

function    [imThresLoc, totalSpots] = locThres(nucLabel, im)

nNuc = max(nucLabel, [], 'all');
nThres = size(im{1}, 2);
imThresLoc = zeros(size(nucLabel));
filtWinSize = 25; % 15
if size(nucLabel, 3)>1
    filtWinZ = 11;
else
    filtWinZ = 1;
end

maxIter = 2; % use 2, 3 or 4 
maxIter = maxIter^2; % local threshold

thresInd = 10;    %%%%%%%%% change here %%%%%%%%%
jMin = fix(thresInd/maxIter)+1;
if jMin == nThres
    jMax = nThres;
else
    jMax = jMin+1;
end

jMax = 4;
jMin = 4;

imTempStack = cell(1,nNuc);
cropSpotValStack = cell(1,nNuc);
cropSpotBinStack = cell(1,nNuc);
spotBin = cell(1, nNuc);
bgMean = cell(1, nNuc);
totalSpots = zeros(1, nNuc);

ssimVal = cell(1,nNuc);
corrVal = cell(1,nNuc);
ssimBin = cell(1,nNuc);
corrBin = cell(1,nNuc);
numSpots = cell(1,nNuc);
spotAreaAv = cell(1,nNuc);
spotAreaStd = cell(1,nNuc);

for i= 1:nNuc  % Nuclei    
    imTempStack{i} = cell(1, nThres);    
    spotBin{i} = cell(1, nThres);
    cropSpotValStack{i} = cell(1, nThres);
    cropSpotBinStack{i} = cell(1, nThres);
    [row, col, zz] = ind2sub(size(nucLabel), find(nucLabel==i));
    kk = 1;
%     cropSpotStack = zeros([max(row)-min(row)+1, max(col)-min(col)+1, maxIter]);
    for j = jMax:jMax %    Nuclear Thresholds
        imTemp = im{i}{j};
        imTempStack{i}{j} = cell(1, maxIter);
        cropSpotValStack{i}{j} = cell(1, maxIter);
        cropSpotBinStack{i}{j} = cell(1, maxIter);
        spotBin{i}{j} = cell(1, maxIter);
        imTemp(nucLabel~=i) = 0;
        imTemp = rescale(imTemp);
%         figure('Color', 'w');
        for p= maxIter:maxIter  % Local threshold iterations
            imAvTemp = convn(imTemp, ones(filtWinSize,filtWinSize, filtWinZ)./((filtWinSize^2)*filtWinZ), 'same');
            imStdTemp = stdfilt(imTemp, ones(filtWinSize,filtWinSize, filtWinZ));
            %~~~~~~~~~~~~~~~~~~~~~~~
%             subplot(sqrt(maxIter), sqrt(maxIter), p)
%             noisePlotter(imAvTemp, imStdTemp);
            %~~~~~~~~~~~~~~~~~~~~~~~
            thres = imAvTemp + imStdTemp;
            imTemp = imTemp - thres;
            imTemp(imTemp<0) = 0;
            imTemp = rescale(imTemp);
            imTempStack{i}{j}{p} = imTemp;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Cut one nucleus out
            nucCrop = imcrop3(imTempStack{i}{j}{p}, [min(col), min(row), min(zz)...
                (max(col) - min(col)), (max(row) - min(row)), (max(zz) - min(zz))]);
            cropSpotValStack{i}{j}{p} = nucCrop;
%             imTempStack{i}{j}{p} = imbinarize(nucCrop);
            cropSpotBinStack{i}{j}{p} = imbinarize(nucCrop);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Calculate background pixels
            if p==1 && j==1
                spotBin{i}{j}{p} = imbinarize(imTemp);
                imTemp2 = im{i}{j};
                bgMean{i} = mean(imTemp2(nucLabel==i & spotBin{i}{j}{p}==0), 'all');
            end
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if p==1 && j==1
                cropSpotVal1 = cropSpotValStack{i}{j}{p};
                cropSpotBin1 = cropSpotBinStack{i}{j}{p};
            end
                     
%             [ssimVal{i}(kk),~] = ssim(cropSpotVal1,  cropSpotValStack{i}{j}{p});
%             corrVal{i}(kk) =  corr2(cropSpotVal1,  cropSpotBinStack{i}{j}{p});%ifftn(fftn(cropSpotVal1).*conj(fftn(cropSpotBinStack{i}{j}{p})),'symmetric');
%             [ssimBin{i}(kk),~] = ssim(single(cropSpotBin1),  single(cropSpotBinStack{i}{j}{p}));
%             corrBin{i}(kk) =  corr2(cropSpotBin1,  cropSpotBinStack{i}{j}{p});%ifftn(fftn(single(cropSpotBin1)).*conj(fftn(single(cropSpotBinStack{i}{j}{p}))),'symmetric');
            
            CC = bwconncomp(cropSpotBinStack{i}{j}{p}, 8);
            numSpots{i}(kk) = CC.NumObjects;
            A = regionprops(CC, 'area');
            spotAreaAll = cat(1, A.Area);
            spotAreaAv{i}(kk) = mean(spotAreaAll);
            spotAreaStd{i}(kk) = std(spotAreaAll);
            
            kk = kk+1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   Visualization   
%             subplot(sqrt(maxIter), sqrt(maxIter), p)
%             imshow(cropSpotValStack{i}{j}{p}(:,:,5),[]);
%             histogram(cropSpotValStack{i}{j}{p},0.05:0.05:1);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
    end
    
%------------------------------------------------------------
%   Visualization   
%     figure; plot (numSpots);
%     ylabel ('# Spots per nucleus');
%     yyaxis right;
%     errorbar(spotAreaAv, spotAreaStd);
%     ylabel('Mean spot area (pixels)');            
%     xlabel('Threshold iteration');
%     
%     ax = gca;
%     ax = gca;
%     ax.FontSize = 12;
%     ax.LineWidth = 1.5;
%     box(ax,'on');
%     grid off;
%     x0 = 100;
%     y0= 100;
%     plotWidth=350;
%     plotHeight=350;
%     set(gcf,'position',[x0,y0,plotWidth,plotHeight])
%     hold off;    

%------------------------------------------------------------
%   Visualization
%     figure; plot(ssimVal{i}); title('ssim val');
%     figure; plot(ssimBin); title('ssim bin');
%     figure; plot(corrVal); title('corr val');
%     figure; plot(corrBin); title('corr bin');
%     figure; plot(numSpots{i}); title('total spots');

%------------------------------------------------------------
%   Visualization
%     figure('color', 'w'); 
%     plot(numSpots{i}, 'k'); findpeaks(numSpots{i}); ylabel('#spots per nucleus');
% %     plot(ssimBin{i}, 'k'); findpeaks(ssimBin{i}); ylabel('ssim spot pix');
%     xlabel('Threshold iterations');
%     hold on; 
%     yyaxis right;
%         plot(ssimBin{i}, 'k'); findpeaks(ssimBin{i}); ylabel('ssim spot pix');
% %     plot(corrBin); findpeaks(corrBin); ylabel('Corr coeff. spot pix');%     
%     ylim([0 1]);    
%     ax = gca;
%     ax = gca;
%     ax.FontSize = 12;
%     ax.LineWidth = 1.5;
%     box(ax,'on');
%     grid off;
%     x0 = 100;
%     y0= 100;
%     plotWidth=350;
%     plotHeight=350;
%     set(gcf,'position',[x0,y0,plotWidth,plotHeight])
%     hold off;
%------------------------------------------------------------
    
%     [~, thresInd] = max(corrVal(2:end));    
%     thresInd = 24;    %%%%%%%%%%%%% change here %%%%%%%%%
%     
%     jj = fix(thresInd/maxIter)+1;
%     if jj==0
%         jj=1;
%     end

%     if rem(thresInd,  maxIter) >0
%         pp = rem(thresInd,  maxIter);
%     else 
%         pp = 2;
%     end
    
%     if rem(thresInd,  maxIter) >0
%         pp = rem(thresInd,  maxIter);
%     else 
%         pp = maxIter;
%     end

    pp = maxIter;
    imThresLoc = imThresLoc + imTempStack{i}{jMax}{pp};
%     figure; imshow(imcrop(imThresLoc, [min(col), min(row), (max(col) - min(col)), (max(row) - min(row))]),[]);
%     totalSpots(1, i) = numSpots{i}(pp);
totalSpots(1, i) = numSpots{i}(1);
end
end

%_____________________________________________________________
%   Plots pixel distribution and marks the bins of threshold
function plotImgHist(im, thres)
im = rescale(im);
[p1, e1] = histcounts(im, 'Normalization','probability');
binC = e1 + (e1(2) - e1(1))/2;
plot(binC(1:end-1), p1, 'r-');
hold on;
H1=area(binC(1:end-1), p1);
set(H1(1),'FaceColor',[0.9 0.9 0.7]);
hold on;
% t1 = graythresh(im);
% idx=binC<t1;
% Ht1=area(binC(idx),p1(idx));
% set(Ht1(1),'FaceColor',[0.7 0.7 0.5]);
for i = 1:length(thres)-1
    idx= find(binC>thres(i) & binC<thres(i+1));
    Ht1=area(binC(idx),p1(idx));
    
    if ~isempty(Ht1)
        set(Ht1(1),'FaceColor',[1/i 0.2 0.4]);
    end
    hold on;
end
xlabel('Normalized counts')
ylabel('Probability density')
title('Multilevel Otsu thresholding')
axis('square')
xlim([0 1]);
% legend('nucleus', 'background');
hold on;
end

function noisePlotter(imAvTemp, imStdTemp)
meanArr = reshape(imAvTemp, [numel(imAvTemp), 1]);
varArr = reshape(imStdTemp.^2, [numel(imStdTemp), 1]);
scatter(meanArr, varArr, '.r');
[sortMean,indSort] = sort(meanArr);
sortVar = varArr(indSort);
[lB, uB]= bounds(meanArr(2:end));
binMeanInd = discretize(sortMean,linspace(lB, uB, 10));
end