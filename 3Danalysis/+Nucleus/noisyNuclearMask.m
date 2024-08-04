function [fim, t1] = noisyNuclearMask(im, mask, blurControl, nucSize, clearBorder)
mask2 = activecontour(imcomplement(im), mask, 10, 'edge', 'SmoothFactor', 0.5, 'ContractionBias', 1);
imScale = rescale(im);
imFilt = imgaussfilt(imadjust(imScale));
t1 = adaptthresh(imFilt, blurControl, 'ForegroundPolarity', 'bright', 'Statistic', 'gaussian'); 
bw1 = imcomplement(imbinarize(imFilt, t1));
bw1 = imfill(bw1, 'holes');
bw1 = bwareaopen(bw1, ceil(nucSize/3));
bw2 = mask & mask2;
bw3 = imfill(bw2, 'holes');
bw3 = bwareaopen(bw3, ceil(nucSize/3));
D = -bwdist(~bw3);
D(bw3)= -inf;
DL = watershed(D);
bgm = DL == 0;
% gmag = imimposemin(imgradient(imScale), bw3 | bgm);
% DL2 = watershed(gmag);
% S = regionprops(DL2,'Area');
% Area = vertcat(S.Area);
% [~,bgValue] = max(Area);
% bw4 = any((DL2~=bgValue & DL2~=0),3);
bw4 = any((bw1==1 & bgm==0),3);
bw5 = imfill(bw4, 'holes');
bw5 = bwareaopen(bw5, nucSize);
if clearBorder
    bw5 = imclearborder(bw5);
end
fim = bw5;
end