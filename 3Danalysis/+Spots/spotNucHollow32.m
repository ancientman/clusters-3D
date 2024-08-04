function [bwSpot, nucLabelSpot] = spotNucHollow32(im, nucLabel, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imUseType = metaDataDS.analysisInfo.imUseType;
howManySigmas = 4; % for hunchback use 7
thresholdFraction = 0.98; %for hunchback use 0.98
labels = sort(nonzeros(HelperFunctions.count_unique(nucLabel)));
nNuc = length(labels);
meanNuc = zeros(nNuc,1);
sigmaNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc, 1);
minMRNAspotVol = metaDataDS.analysisInfo.minMRNAspotVol;
imSpot = zeros(size(im));
% imFilt = imgaussfilt3(medfilt3(im));
% im = imFilt;
imNucLabelSpot = zeros(size(im));

% for hunchback use
% H1 = fspecial3('gaussian',[21, 21, size(im, 3) ], 7);
% H2 = fspecial3('gaussian',[21, 21, size(im, 3) ], 5);

% for bottleneck use
H1 = fspecial3('gaussian',[15, 15, size(im, 3) ], 7);
H2 = fspecial3('gaussian',[15, 15, size(im, 3) ], 5);

dog = H2 - H1;
imDogFilt = convn(im,dog,'same');
im = rescale(imDogFilt);
for i = 1:nNuc
    meanNuc(i) = mean(im(nucLabel==labels(i)),'all');
    sigmaNuc(i) = std(im(nucLabel==labels(i)), 0, 'all');
    thresNuc(i) = prctile(im(nucLabel==labels(i) & im>(howManySigmas*sigmaNuc(i)+meanNuc(i))),...
        (thresholdFraction),'all');
end
    bwSpot = im>mean(thresNuc, 'omitnan');
    bwSpot = imfill(bwSpot, 'holes');
    bwSpot = bwareaopen(bwSpot, round(minMRNAspotVol/3));
    
%     imSpotTemp(nucLabel~=labels(i)) = 0; % use carefully, only for convex hull
%     imSpot = imSpot + imSpotTemp;
% %     imNucLabelSpot = imNucLabelSpot + i.*imSpotTemp;
% end

% bw1 = imbinarize(imSpot, 0.9);
% bw2 = bwareaopen(bw1, minMRNAspotVol/3);
% bw2 = bwareaopen(bw1, minMRNAspotVol);


% Assuming transcription spots are twice as large and appears over twice as many frames
% if imUseType==3 % use "sharp"
%     bw2 = bwareaopen(bw1, minMRNAspotVol/3);
% else
%     bw2 = bwareaopen(bw1, minMRNAspotVol);
% end
% bwSpot = bw2;

nucLabelSpot = [];
% nucLabelSpot = bw2.*imNucLabelSpot;

end
