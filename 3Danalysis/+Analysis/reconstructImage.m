function reconstructImage(folderPath)
dirInfo = dir(folderPath);
if size(dirInfo, 1)==0
    error('data folder is empty')
end
c1NucStruct = load(append(folderPath, filesep, 'c1NucPropDS.mat'));
c1SpotStruct = load(append(folderPath, filesep, 'c1SpotPropDS.mat'));
metaData = load(append(folderPath, filesep, 'metaDataDS.mat'));
imSizeX = metaData.metaDataDS.imagingInfo.X;
imSizeY = metaData.metaDataDS.imagingInfo.Y;
imSizeZ = metaData.metaDataDS.imagingInfo.Zslices;
binSpotIm = zeros(imSizeY, imSizeX, imSizeZ);
binNucIm = zeros(imSizeY, imSizeX, imSizeZ);
binSpotIm(cell2mat(vertcat(c1SpotStruct.c1SpotProp{1}(:).voxIdx))) = 1;
binNucIm(cell2mat(vertcat(c1NucStruct.c1NucProp{1}(:).voxIdx))) = 1;
HelperFunctions.imshow3D(or(binSpotIm, bwperim(binNucIm)));
end