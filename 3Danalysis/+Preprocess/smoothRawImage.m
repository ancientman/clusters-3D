function [fim] = smoothRawImage(im, metaDataDS)
filterType = metaDataDS.analysisInfo.nucSmoothFilter;
filterParam = metaDataDS.analysisInfo.smoothingParam;
imTemp = im;
sx = filterParam;
sy = filterParam;
switch filterType
    case 1
    imSmooth = HelperFunctions.gaussianFilter2D(imTemp,sx,sy);    
    case 2
        imSmooth = HelperFunctions.wienerFilter2D(imTemp, sx, sy);
    case 3
        imSmooth = HelperFunctions.meanFilter2D(imTemp, sx, sy);
    case 4
        imSmooth = HelperFunctions.medianFilter2D(imTemp, sx, sy);
    case 5
        imSmooth = HelperFunctions.bilateralFilter2D(imTemp, sx, sy);
    case 6
        imSmooth = HelperFunctions.nonlocalFilter2D(imTemp, sx, sy);
    case 7
        %using abandoned filter type: DoG
        Sxy = [19, 25];
        imSmooth = HelperFunctions.dogFilter2D(imTemp, Sxy);
    otherwise
        warning ('no nucleus smooth type mentioned, using non local')
        imSmooth = HelperFunctions.wienerFilter2D(imTemp, sx, sy);
end
% imDiff = im - imSmoothClean;
% figure; imshowpair(im, imSmoothClean, 'montage'); title (strcat('smooth as ', {' '}, 'ass'));
fim = imSmooth;
end