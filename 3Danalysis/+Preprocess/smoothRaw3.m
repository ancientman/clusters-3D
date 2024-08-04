function [fim] = smoothRaw3(im, metaDataDS)
filterParam = metaDataDS.analysisInfo.smoothingParam;

if rem(filterParam,2)==1
    sx = filterParam;
    sy = filterParam;
else
    sx = filterParam+1;
    sy = filterParam+1;
end
sz = (2*filterParam) + 1;
imSmooth = imgaussfilt3(im, [0.5,0.5,0.7], 'FilterSize', [sx, sy, sz]);   
fim = imSmooth;
end