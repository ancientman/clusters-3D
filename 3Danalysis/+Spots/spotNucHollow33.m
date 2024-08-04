function [bwSpotTemp, nucLabelSpot] = spotNucHollow33(im, nucLabel, metaDataDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
howManySigmas = 3; % for hunchback use 7
thresholdFraction = 0.98; %for hunchback use 0.98
labels = sort(nonzeros(HelperFunctions.count_unique(nucLabel)));
nNuc = length(labels);
meanNuc = zeros(nNuc,1);
sigmaNuc = zeros(nNuc,1);
thresNuc = zeros(nNuc, 1);
minMRNAspotVol = metaDataDS.analysisInfo.minMRNAspotVol;
bwSpot = zeros(size(im));

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
    bwSpotTemp = im;
    bwSpotTemp(nucLabel~=labels(i)) = 0;
    bwSpotTemp = bwSpotTemp>thresNuc(i);
    bwSpotTemp = imfill(bwSpotTemp, 'holes');
    bwSpotTemp = bwareaopen(bwSpotTemp, minMRNAspotVol/3);
    bwSpotTemp = imclearborder(bwSpotTemp);
    bwSpot = or(bwSpot,  bwSpotTemp); 
end
nucLabelSpot = [];
end
