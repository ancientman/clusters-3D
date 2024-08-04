function c1c2IntensityProp = c1IntensityC2(c1Stack, c2SpotLabelStack, c2SpotProp, c1NucLabelStack, c1SpotLabelStack)
bgSub = 1;
maxcorrPlane = 1;
distMax = 25;
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
spotCentNuc = cell(1, nNuc);
peakLoc =  cell(1, nNuc);
corrCoeff = cell(1, nNuc);
mockIm = zeros(size(c1Stack(:,:,1,1)));
for i = 1:nNuc
    spotCentNuc{i} = zeros(timePoints, 3);
end

for t=1:timePoints
    for i = 1:size(c2SpotProp{t}, 2)
        if ~isempty(c2SpotProp{t}(i).center)
            spotCentNuc{i}(t,:) = round(c2SpotProp{t}(i).center);
        end
    end
end

radValBin= cell(1, nNuc);
for i = 1:nNuc
    radValBin{i} = zeros(distMax, timePoints);
    peakLoc{i} = zeros(timePoints, 1);
    corrCoeff{i} = zeros(timePoints, 1);
    iter = 0;
    flag = 0;
    for t = 1:timePoints
        if spotCentNuc{i}(t,3) ~=0
            zTemp = spotCentNuc{i}(t,3);
            labelMatTemp = c1NucLabelStack(:,:,zTemp,t);
            labelMatTemp(labelMatTemp~=i) = 0;
            c1Temp = c1Stack(:,:,zTemp,t).*labelMatTemp;           
            % For subtraction
            %--------------------------------------------------------------
            if bgSub == 1
                c1Temp2 = c1Temp;          
                c1Temp2(c1SpotLabelStack(:,:,zTemp,t)==i) = 0;       
                bgMean(t, i) = mean(nonzeros(c1Temp2),'all', 'omitnan');
                bgStd(t, i) = std(nonzeros(c1Temp2),0,'all', 'omitnan');
                imTempNew = c1Temp - (bgMean(t, i) + 2*bgStd(t, i));
                imTempNew(imTempNew<0) = 0;
                imTempNew(isnan(imTempNew)) = 0;     
                c1Temp = rescale(imTempNew);
            end
            %--------------------------------------------------------------
            spotCentTemp = spotCentNuc{i}(t,1:2);
            
            xEdge = (spotCentTemp(1) - distMax) : (spotCentTemp(1) + distMax);
            yEdge = (spotCentTemp(2) - distMax) : (spotCentTemp(2) + distMax);
            [xNuc,yNuc] = meshgrid(xEdge, yEdge);
            if  ~isempty(setdiff(xNuc,1:size(mockIm, 1))) || ~isempty(setdiff(yNuc,1:size(mockIm, 2)))
                continue
            else            
                mockIm(sub2ind(size(mockIm), yNuc, xNuc)) = 1;
                c1Temp(~mockIm) = 0;
                c1Crop = imcrop(c1Temp, [(spotCentTemp(1) - distMax), (spotCentTemp(2) - distMax), 2*distMax, 2*distMax]);
                c1Crop = medfilt2(c1Crop);
                if bgSub == 1
                    c1Crop = rescale(c1Crop); % For subtraction             
                else
                    c1Crop = c1Crop/c1Crop(distMax+1, distMax+1); % For no subtraction
                end               
%                 corrMat = corr2D(c1Crop);
                
                [valRadBinAvg, valRadBinSem] = cart2PolAvg(distMax, c1Crop);
%                 radValBin{i} =  horzcat(radValBin{i},valRadBinAvg);
                radValBin{i}(:,t) =  valRadBinAvg;
                iter = iter + 1;
            end
        end
        if any(radValBin{i}(:,t))
            if flag == 0
                c1CropPrime = c1Crop;
            end
            corrCoeff{i}(t) = corr2(c1CropPrime, c1Crop);
%             [pks,loc] = findpeaks(radValBin{i}(:,t));
            [maxVal,loc] = max(radValBin{i}(:,t));
            peak1Loc = min(loc);
            if peak1Loc == 0
                peak1Loc = 1;
            end
            peakLoc{i}(t) = peak1Loc;
            flag = flag+1;
        end
    end

     %  Visualizations%%%%%%%%%%%%%%
%     linePlot1(radValBin{i});
%     linePlot2(0.042*peakLoc{i});
%     linePlot2(corrCoeff{i});
end

c1c2IntensityProp.radValBin = radValBin;
c1c2IntensityProp.corrCoeff = corrCoeff;
c1c2IntensityProp.peakLoc = peakLoc;

end


function [valRadBinAvg, valRadBinSem] = cart2PolAvg(rMax, imCrop)
xVals = ones(1, 2*rMax+1)'*(-rMax:rMax);
yVals = (-rMax:rMax)'*ones(1, 2*rMax+1);
[theta, rho, valTemp] = cart2pol(xVals,yVals,  imCrop);  % convert x, y to polar coordinates
aR = reshape(rho,1, (2*rMax+1)^2);
aValTemp = reshape(valTemp,1, (2*rMax+1)^2);
[rrTemp,ind] = sort(aR);
corrRadVal = aValTemp(ind);
% rad= 0:floor(max(rrTemp));
rad= 0:rMax;
[n, bin] = histcounts(rrTemp, rad);
valRadBinAvg = zeros(length(bin)-1, 1);
valRadBinSem = zeros(length(bin)-1, 1);
for j = 1:length(bin)-1
    valRadBinAvg(j) = sum(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sum(rrTemp>bin(j) & rrTemp<=bin(j+1));
    valRadBinSem(j) = std(corrRadVal(rrTemp>bin(j) & rrTemp<=bin(j+1)))./sqrt(sum(rrTemp>bin(j) & rrTemp<=bin(j+1)));                        
end
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