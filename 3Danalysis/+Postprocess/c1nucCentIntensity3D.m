function c1nucCentIntensityProp = c1nucCentIntensity3D(c1Stack, c1NucLabelStack, c1NucProp, metaDataDS)
c1nucCentIntensityProp = struct([]);
imXMax = metaDataDS.imagingInfo.X;
imYMax = metaDataDS.imagingInfo.Y;
imZMax = metaDataDS.imagingInfo.Zslices;
% imXMax = metaDataDS.imagingInfo.stackSizeX;
% imYMax = metaDataDS.imagingInfo.stackSizeY;
% imZMax = metaDataDS.analysisInfo.nZslices;
maxcorrPlane = 1;
distXYMax = 25; %%% hard code
distZZMax = 3; %%% hard code
im = c1Stack;
timePoints = size(im, 4);
nNuc = length(nonzeros(HelperFunctions.count_unique(c1NucLabelStack(:,:,:,:))));

bgMean = double(zeros(timePoints, nNuc));
bgStd = double(zeros(timePoints, nNuc));

timePoints = size(c1Stack, 4);
nNuc = max(c1NucLabelStack, [], 'all');
if nNuc == 0
    nNuc = 1;
end
nucCent = cell(1, nNuc);
peakLoc =  cell(1, nNuc);
corrCoeff = cell(1, nNuc);
mockIm = zeros(size(c1Stack(:,:,:,1)));
for i = 1:nNuc
    nucCent{i} = zeros(timePoints, 3);
end

for t=1:timePoints
    for i = 1:size(c1NucProp{t}, 2)
        if ~isempty(c1NucProp{t}(i).center)
            nucCent{i}(t,:) = round(c1NucProp{t}(i).center);
        end
    end
end

radValAv= cell(1, nNuc);
radValSem= cell(1, nNuc);

for i = 1:nNuc
    radValAv{i} = cell(1, timePoints);
    radValSem{i} = cell(1, timePoints);
    for t = 1:timePoints
        
        if nucCent{i}(t,3) ~=0       
            %--------------------------------------------------------------
            nucCentTemp = nucCent{i}(t,:);         
            
            xEdgeLeft = (nucCentTemp(1) - distXYMax);
            if xEdgeLeft<1
                xEdgeLeft = 1;
                distXMax = nucCentTemp(1) - 1;
            else
                xEdgeLeft = (nucCentTemp(1) - distXYMax);
                distXMax = distXYMax;
            end
            xEdgeRight = (nucCentTemp(1) + distXYMax);
            if xEdgeRight>imXMax
                xEdgeRight = imXMax;
                distXMax = imXMax - nucCentTemp(1);
            else
                xEdgeRight = (nucCentTemp(1) + distXYMax);
                distXMax = distXYMax;
            end
            yEdgeLeft = (nucCentTemp(2) - distXYMax);
            if yEdgeLeft<1
                yEdgeLeft = 1;
                distYMax = nucCentTemp(2) - 1;
            else
                yEdgeLeft = (nucCentTemp(2) - distXYMax);
                distYMax = distXYMax;
            end
            yEdgeRight = (nucCentTemp(2) + distXYMax);
            if yEdgeRight>imYMax
                yEdgeRight = imYMax;
                distYMax = imYMax - nucCentTemp(2);
            else
                yEdgeRight = (nucCentTemp(2) + distXYMax);
                distYMax = distXYMax;
            end
            zEdgeLeft = (nucCentTemp(3) - distZZMax);
            if zEdgeLeft<1
                zEdgeLeft = 1;
                distZMax = nucCentTemp(3) - 1;
            else
                zEdgeLeft = (nucCentTemp(3) - distZZMax);
                distZMax = distZZMax;
            end
            zEdgeRight = (nucCentTemp(3) + distZZMax);
            if zEdgeRight>imZMax
                zEdgeRight = imZMax;
                distZMax = imZMax - nucCentTemp(3);
            else
                zEdgeRight = (nucCentTemp(3) + distZMax);
                distZMax = distZZMax;
            end      
            
            mockIm = zeros([(2*distXYMax+1), (2*distXYMax+1), (2*distZZMax+1)]);           
            [pixDistMat] = pixDistCalculator(mockIm, metaDataDS);
            %--------------------------------------------------------------
            try
                c1CropTemp = imcrop3(c1Stack(:,:,:,t), [xEdgeLeft, yEdgeLeft, zEdgeLeft, 2*distXMax, 2*distXMax, 2*distZMax]);
            catch
                warning('failed to make crop for time = %d, nuc = %d\n', t, i);
            end
            
%             if bgSub == 0
                c1Crop = (c1CropTemp);
                valDistProp = avIntensityDist(pixDistMat, c1Crop);
%             elseif bgSub == 1
                c1Temp2 = c1CropTemp;
                bgMean(t, i) = mean(nonzeros(c1Temp2),'all', 'omitnan');
                bgStd(t, i) = std(nonzeros(c1Temp2),0,'all', 'omitnan');
                imTempNew = c1Temp2 - (bgMean(t, i) + 2*bgStd(t, i));
                imTempNew(imTempNew<0) = 0;
                imTempNew(isnan(imTempNew)) = 0;     
                c1Crop = (imTempNew);
                valBgSubDistProp = avIntensityDist(pixDistMat, c1Crop);
%             end
            %--------------------------------------------------------------
            c1nucCentIntensityProp{i}.valDistProp{t} = valDistProp;
            c1nucCentIntensityProp{i}.valBgSubDistProp{t} = valBgSubDistProp;
            %--------------------------------------------------------------
        else
            c1nucCentIntensityProp{i}.valDistProp{t} = [];
            c1nucCentIntensityProp{i}.valBgSubDistProp{t} = [];
        end
    end
end
%     Visualizations%%%%%%%%%%%%%%
%     linePlot1(radValBin{i});
%     linePlot2(0.042*peakLoc{i});
%     linePlot2(corrCoeff{i});
end


function [pixDistMat] = pixDistCalculator(imCrop, metaDataDS)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
imRows = size(imCrop, 2);
imCols = size(imCrop, 1);
imPages = size(imCrop, 3);

if mod(imRows, 2) == 1
    rowCent = idivide(int16(imRows), 2);
else 
    rowCent = idivide(int16(imRows), 2);
    warning('crop center is shifted by one pixel \n');
end

if mod(imCols, 2) == 1
    colCent = idivide(int16(imCols), 2);
else 
    colCent = idivide(int16(imCols), 2);
    warning('crop center is shifted by one pixel \n');
end

if mod(imRows, 2) == 1
    pageCent = idivide(int16(imPages), 2);
else 
    pageCent = idivide(int16(imPages), 2);
    warning('crop center is shifted by one pixel \n');
end

pixDistMat = zeros(size(imCrop));
for i = 1:numel(imCrop)
    [row, col, page] = ind2sub(size(imCrop), i);
    distSqr =(double(row-rowCent)*xPixUM)^2 + (double(col-colCent)*yPixUM)^2 + (double(page-pageCent)*zPixUM)^2;
    pixDistMat(i) = sqrt(distSqr);
end

end

function valDistlProp = avIntensityDist(pixDistMat, imCrop)
distMax = max(pixDistMat, [], 'all'); %distance in microns
distBinSize = 0.1; %distance in microns
totalBins = max(1, floor(distMax/distBinSize)); 
distBin = 0:distBinSize:(distBinSize*totalBins);
distVal = cell(1, totalBins);
for k = 1:numel(pixDistMat)
    for i = 1:totalBins
        if (pixDistMat(k) > distBin(i)) && (pixDistMat(k) <= distBin(i+1))
            try
                distVal{i} = vertcat(distVal{i}, imCrop(k));
            catch
                aaa = 0;
            end
        end
    end
end

valDistlProp.distBin = distBin(2:end)';
valDistlProp.distValAv = (cellfun(@mean, distVal))';
valDistlProp.distValSem = (cellfun(@(x) std(x)/sqrt(numel(x)), distVal)');
end





function [valRadProp] = avIntensityRad(rMax, rZMax, imCrop, metaDataDS)
valRadProp = struct([]);
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;

distBin = 0.1; %distance in microns
colPoints = (-rMax:rMax)*xPixUM;
rowPoints = (-rMax:rMax)*yPixUM;
zPoints = (-rZMax:rZMax)*zPixUM;
[X, Y, Z] = ndgrid(colPoints, rowPoints, zPoints);
pixRadDist = sqrt(X.^2 + Y.^2 + Z.^2);
distMax = min(pixRadDist(rMax+1, rMax+1, 1), pixRadDist(1, rMax+1, rZMax+1));
radDistBin = 0:distBin:distMax;
shellValAll = cell(1, length(radDistBin)-1);
for i = 1:(2*rMax+1)
    for j = 1:(2*rMax+1)
        for k = 1:(2*rZMax+1)
            for p = 1:length(radDistBin)-1
                if pixRadDist(i,j,k) > radDistBin(p) && pixRadDist(i,j,k) < radDistBin(p+1)
                    shellValAll{p} = vertcat(shellValAll{p}, imCrop(i, j, k));
                end
            end
        end
    end
end

for p = 1: length(shellValAll)
    shellValAll{p}(1) = [];
end
shellValAvg = cellfun(@mean, shellValAll);
shellValStd = cellfun(@std, shellValAll);
shellValSem = cellfun(@(x) std(x)/sqrt(numel(x)), shellValAll);

valRadProp{1}.radDistBin = radDistBin';
valRadProp{1}.shellValAvg = shellValAvg';
valRadProp{1}.shellValStd = shellValStd';
valRadProp{1}.shellValSem = shellValSem';
end


function linePlot1(radValBin)
figure('color', 'w');
for i=1:size(radValBin, 2)
    p1= plot(0.042*(1:size(radValBin, 1)), radValBin(:,i));
    p1.LineWidth = 1.5;
    p1.Color = [0.7 0.7 0.7];
    p1.LineStyle = '-';    
    hold on;
end

nonZeroRows = any(radValBin, 1);
radValBinAvg = mean(radValBin(:,nonZeroRows), 2);
radValBinErr = std(radValBin(:,nonZeroRows), 0, 2);

p2= plot(0.042*(1:size(radValBin, 1)), radValBinAvg);
p2.LineWidth = 2;
p2.Color = [0.8 0.3 0.3];
p2.LineStyle = '-';
ylabel('Relative bcd intensity (a. u.)');
xlabel('Apparent radius ({\mu}m)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;

plotErr(radValBinAvg, radValBinErr);
end

function plotErr(val, err)
figure('color', 'w');
pixScale = 0.042; 
patchTop = val+err;
patchTop = reshape(patchTop,1,[]);
patchBot = val-err;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch=[pixScale*(1:length(val)),fliplr(pixScale*(1:length(val)))];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = [0.3, 0.7, 0.6];
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.6;
hold on;
pl = plot(pixScale*(1:length(val)), val);
pl.LineWidth = 1.5;
pl.Color = [0.4 0.4 0.4];
pl.LineStyle = '--';
ylabel('Relative bcd intensity (a. u.)');
xlabel('Apparent radius (r)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function linePlot2(radValBin)
figure('color', 'w');

p1= plot(1.1*(1:size(radValBin, 1)), radValBin);
p1.LineWidth = 1.5;
p1.Color = [0.3 0.3 0.3];
p1.LineStyle = '-';   
ylabel('Correlation coefficient');
% ylabel('Brightest Bicoid from hb ({\mu}m');
xlabel('Time (s)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function corrMat = corr2D(cropIm)
cropZeroPad = 100;
centCropPad = zeros(size(cropIm, 1) + 2*cropZeroPad,  size(cropIm, 2) + 2*cropZeroPad);
centCropPad(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = cropIm;
mask = zeros(size(centCropPad));
mask(cropZeroPad+1: size(centCropPad, 1)-cropZeroPad, cropZeroPad+1: size(centCropPad, 2)-cropZeroPad) = 1;
nParticles = sum(cropIm, 'all');  % number of particles within mask
nPixMask = sum(mask, 'all');      % area of mask
NP = real(fftshift(ifft2(abs(fft2(mask)).^2))); % Normalization for correct boundary conditions
G1 = nPixMask^2/nParticles^2*real(fftshift(ifft2(abs(fft2(centCropPad)).^2)))./NP; % 2D G(r) with proper normalization
corrMat = imcrop(G1, [floor(size(G1, 2)/2+1)-cropZeroPad, floor(size(G1, 1)/2+1)-cropZeroPad, 2*cropZeroPad, 2*cropZeroPad]);  %only return valid part of G
[valAvg, valSem] = cart2PolCorr(cropZeroPad, corrMat);
corrCentRadAvg = valAvg;
corrCentRadSem = valSem;
plotErr(corrCentRadAvg, corrCentRadSem);
end

function [valRadBinAvg, valRadBinSem] = cart2PolCorr(rMax, corrCut)
xVals = ones(1, 2*rMax+1)'*(-rMax:rMax);
yVals = (-rMax:rMax)'*ones(1, 2*rMax+1);
[theta, rho, valTemp] = cart2pol(xVals,yVals,  corrCut);  % convert x, y to polar coordinates
aR = reshape(rho,1, (2*rMax+1)^2);
aValTemp = reshape(valTemp,1, (2*rMax+1)^2);
[rrTemp,ind] = sort(aR);
corrRadVal = aValTemp(ind);
rad= 0:floor(max(rrTemp));
[n, bin] = histcounts(rrTemp, rad);
valRadBinAvg = zeros(length(bin)-1, 1);
valRadBinSem = zeros(length(bin)-1, 1);
for j = 1:length(bin)-1
    valRadBinAvg(j) = sum(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sum(rrTemp>bin(j) & rrTemp<=bin(j+1));
    valRadBinSem(j) = std(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sqrt(sum(rrTemp>bin(j) & rrTemp<=bin(j+1)));                        
end
end