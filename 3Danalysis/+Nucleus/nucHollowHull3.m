function [fim] = nucHollowHull3(mask, channelFilled, metaDataDS)
stackLen = size(mask, 3);
maskHull = zeros(size(mask));
% if channelFilled==0
    for i=1:stackLen
        mask(:,:,i) = imfill(mask(:,:,i), 'holes');
        maskHull(:,:,i) = bwconvhull(mask(:,:,i), 'objects');        
    end
% end
fim= maskHull;
end