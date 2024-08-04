function [deconvIm] = sharpRaw3(im,metaDataDS)
xPixUM = metaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.analysisInfo.yPixUM;
xyPsfUM = metaDataDS.analysisInfo.xyPsfUM;
zPixUM = metaDataDS.analysisInfo.zPixUM;
zPsfUM = metaDataDS.analysisInfo.zPsfUM;
sigma(1) = (xyPsfUM/xPixUM)/(2*sqrt(2*log(2)));
sigma(2) = (xyPsfUM/yPixUM)/(2*sqrt(2*log(2)));
sigma(3) = (zPsfUM/zPixUM)/(2*sqrt(2*log(2)));
hSize = sigma;
hSize(:) = 2*ceil(max(sigma));
psf = fspecial3('gaussian',hSize,sigma)/sum(ones(hSize),'all');
psf2 = fspecial('gaussian',hSize(1:2),sigma(1))/sum(ones(hSize(1:2)),'all');
I = zeros(size(im));
for i=1:size(im, 3)
   I(:,:,i) = edgetaper(im(:,:,i),psf2);
end
J = deconvlucy(I, psf, 20);
% JJ = deconvblind({J}, {psf});
deconvIm = (J);
end