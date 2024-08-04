function [imNucRemoved, mask] = removeNuc3(im, iter, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allows to interactively mask unwanted nuclei one by one
% Roi masking happens only in the first frame
% in all other frames this function is called 
% but the same roi, drawn in
% the first frame is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iter==1
%     polyVert = {};
    prompt = "another nucleus?";
    i = 1;
    imZproj = max(im,[],3);
    mask = zeros(size(imZproj));
    fprintf("\ndraw polygons around nuclei to reject, one by one\n");
    while i>=1
        imshow(imZproj,[]);
        p = drawpolygon('LineWidth',5,'Color','magenta');
        mask = mask + createMask(p, size(mask,1), size(mask, 2));
%         polyVert{i} = p.Position;        
        another = input(prompt);% 1 if another nucleus is needed to be masked, 0 otherwise
        if another == 1
            i = i+1;
        else
            i = 0;
        end
    end
%     J = im;
%     for j=1:length(polyVert)
%         J = regionfill(J,polyVert{j}(:,1), polyVert{j}(:,2));
%     end
end
J = double(im).*(~repmat(mask, [1, 1, size(im, 3)]));
imNucRemoved = J;
end
        
  