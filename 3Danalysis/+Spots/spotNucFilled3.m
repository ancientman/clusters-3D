function [bwSpot, nucLabelSpot] = spotNucFilled3(im, nucLabel, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Gets called if "deviationofGradient" is used
%   Is actually Deviation from mean   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nNuc = max(nucLabel,[],'all');
elementSize = metaDataDS.analysisInfo.elementSize;
imSpot = zeros(size(im));
imNucLabelSpot = zeros(size(im));
thresNuc = zeros(nNuc, 1);

% spotFilter = metaDataDS.analysisInfo.spotFilter;
% meanNuc = zeros(nNuc,1);
% sigmaNuc = zeros(nNuc,1);
% imBlur1 = medfilt3(imScale,elementSize*ones(1, 3));
% imBlur2 = imgaussfilt3(imBlur1);
% thresholdFraction = 0.99;

thresholdFraction = 0.9;
el = elementSize;
se = strel('disk', el, 8);
imBlur = imtophat(im,se); 

%% Technique based on deviation from max
for i = 1:nNuc
    thresNuc(i) = thresholdFraction*prctile(imBlur(nucLabel==i), (99),'all'); 
    imSpotTemp = imBlur>thresNuc(i);
    imSpotTemp(nucLabel ~= i) = 0;
    imSpotTemp = imfill(imSpotTemp, 'holes');
    imSpot = imSpot + imSpotTemp;
    try
    imNucLabelSpot = imNucLabelSpot + double(i).*imSpotTemp;
    catch
        aaa = 1
    end
end

%% Technique based on deviation from the mean
% for i = 1:nNuc
%     meanNuc(i) = mean2(im(nucLabel==i));
%     sigmaNuc(i) = std2(im(nucLabel==i));
%     thresNuc(i) = prctile(im(nucLabel==i & im>(sigmaNuc(i)+meanNuc(i))),...
%         (thresholdFraction*100),'all');    
%     imSpotTemp = im>thresNuc(i);    
%     imSpotTemp(nucLabel ~= i) = 0;
%     imSpotTemp = imfill(imSpotTemp, 'holes');
%     imSpot = imSpot + imSpotTemp;
%     imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
% end

%% Watershed the spots
% imSpotW = watershedSpot(im, imSpot, metaDataDS);

bwSpot = imSpot;
nucLabelSpot = imNucLabelSpot;
end

function [fim] = watershedSpot(im, mask, metaDataDS)
el = metaDataDS.analysisInfo.elementSize;
% im = rescale(im);
im(mask==0) = 0;

imBlur = imgaussfilt(im);        
se = strel('disk', el, 8);
imBlur = imtophat(im,se); 

bw = imBlur>0;
D = bwdist(~bw);
D = -D;
D(~bw) = -Inf;
L = watershed(D);
bw2 = mask;
bw2(L==0) = 0;
bw3 = bw2;
bw3 = imclearborder(bw3);
fim = bw3;
end