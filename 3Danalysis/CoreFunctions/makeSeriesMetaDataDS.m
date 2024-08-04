function [serieMetaDataDS] = makeSeriesMetaDataDS(analysisFolder)
% makeExpInfo Stores all relevant information/data that will be used
% in the following analyses
% correct folder
expInfoFileName = append(analysisFolder, filesep,'expInfoCZI.txt');
expInfoFileChar = convertStringsToChars(expInfoFileName);
if exist(expInfoFileName, 'file')
    info = HelperFunctions.readtext(expInfoFileChar,'\t');
    userInfo = cell2struct(info(:,2),info(:,1),1);
else
    error('User must supply expInfo.txt in the analysis directory')
end
if isfield(userInfo,'rawDataFilePath')
    serieMetaDataDS.expInfo.rawDataFilePath = userInfo.rawDataFilePath;
    [rawFilePath,rawFileName,rawFileExt] = fileparts(userInfo.rawDataFilePath);
    serieMetaDataDS.expInfo.analysisFolder = analysisFolder;
    serieMetaDataDS.expInfo.rawImagesFolderName = strcat(analysisFolder, filesep, 'out');
    serieMetaDataDS.expInfo.procImagesFolderName = strcat(analysisFolder, filesep, 'DS_', rawFileName);
else
    error('expInfo.txt should include rawDataFilePath')
end
% psf in the z direction
if isfield(userInfo,'Zpsf')
    Zpsf = userInfo.Zpsf;
else
    Zpsf = '500';
    warning('no feature size reference found, using 500 nanometers')
end
%psf along the horizontal plane
if isfield(userInfo,'XYpsf')
    XYpsf = userInfo.XYpsf;
else
    XYpsf = '200';
    warning('no feature size reference found, using 200 nanometers')
end
%what is the starting z for analysis
if isfield(userInfo,'startZ')
    startZ = userInfo.startZ;
else
    startZ = 0;
    warning('using starting z slice = 1');
end
%what is the ending z for analysis
if isfield(userInfo,'endZ')
    endZ = userInfo.endZ;
else
    endZ = 0;
    warning('using ending z slice = %d',Zslices);
end
%which is the protein channel
if isfield(userInfo,'bicoidChannel')
    bicoidChannel = userInfo.bicoidChannel;
else
    bicoidChannel = 1;
end
% Flag to determine which image to use for all analyses
if isfield(userInfo,'imageToUse')
    imageToUse = userInfo.imageToUse;
else
    warning('no use image type specified, using raw images for all analyses')
end
switch imageToUse 
    case 'raw'
        imUseType = 1;
    case 'smooth'
        imUseType = 2;
    case 'sharp'
        imUseType = 3;
    otherwise
        imUseType = 1;
end
% if the nuclei touching the borders need to be cleared
if isfield(userInfo,'imClearBorder')
    imClearBorder = userInfo.imClearBorder;
else
    warning('no valid input smoothing filter type found, not clearing border')
end
switch imClearBorder 
    case 'yes'
        clearBorder = 1;
    case 'no'
        clearBorder = 0;
    otherwise
        clearBorder = 0;
end
% what kind of smoothing filter should be used on the iamges
if isfield(userInfo,'smoothFilterType')
    smoothFilterType = userInfo.smoothFilterType;
else
    warning('no valid input smoothing filter type found, using wiener')
end
switch smoothFilterType 
    case 'gaussian'
        nucSmoothFilter = 1;
    case 'wiener'
        nucSmoothFilter = 2;    
    case 'mean'
        nucSmoothFilter = 3;
    case 'median'
        nucSmoothFilter = 4;
    case 'bilateral'
        nucSmoothFilter = 5;
    case 'nonlocal'
        nucSmoothFilter = 6;
    case 'dog'
        nucSmoothFilter = 7;    %Don't use
    otherwise
        nucSmoothFilter = 5;
end
% what kind of sharpening filter should be used on the iamges
if isfield(userInfo,'deconvolutionFilterType')
    deconvolutionFilterType = userInfo.deconvolutionFilterType;
else
    warning('no valid input deconvolution filter found, using wiener')
end
switch deconvolutionFilterType 
    case 'lucy'
        deconvolutionType = 1;
    case 'wiener'
        deconvolutionType = 2;    
    case 'regular'
        deconvolutionType = 3;    
    otherwise
        deconvolutionType = 2;
end
%how to detect spots (most common would be deviation from mean). devation
%from max is better paired with "imClearBorder" set to yes.
if isfield(userInfo,'spotDetectType')
    spotDetectType = userInfo.spotDetectType;
else
    warning('no valid input spot filter type found, deviationFromMax')
end
switch spotDetectType 
    case 'deviationFromMean'
        spotDetect = 1;
    case 'deviationFromMax'
        spotDetect = 2;    
    case 'deviationofGradient'
        spotDetect = 3;
    otherwise
        spotDetect = 2;
end
%what kind of filtering to use in the spot dete analysis. (most common
%would be top-hat
if isfield(userInfo,'spotFilterType')
    spotFilterType = userInfo.spotFilterType;
else
    warning('no valid input spot filter type found, using top hat')
end
switch spotFilterType 
    case 'gaussian'
        spotFilter = 1;
    case 'bilateral'
       spotFilter = 2;    
    case 'tophat'
        spotFilter = 3;
    case 'median'
        spotFilter = 4;
    otherwise
        spotFilter = 4;
end
%options for detecting nucleus (most common would be watershed)
if isfield(userInfo,'nucleusDetectMethod')
    nucleusDetectMethod = userInfo.nucleusDetectMethod;
else
    nucleusDetectMethod = 'morphclose';
    warning('no nucleus detection type specified, using morphclose')
end
switch nucleusDetectMethod
    case 'reconstruct'
        nucleusDetectMode = 1;
    case 'morphclose'        
        nucleusDetectMode = 2;
    case 'morpherode'
        nucleusDetectMode = 3;
    case 'watershed'
        nucleusDetectMode = 4;
    otherwise
        nucleusDetectMode = 1;
end
% parameter for image thresholding for nucleus detection
if isfield(userInfo,'nucleusIntensityThreshold')
    nucThres = userInfo.nucleusIntensityThreshold;
else
    nucThres = '0.85';
    warning('no threshold value found, using 0.85')
end
% parameter for image smoothing
if isfield(userInfo,'smoothingParam')
    smoothingParam = userInfo.smoothingParam;
else
    smoothingParam = '5';
    warning('no smoothing parameter found, using 5')
end
% cutoff for detected spots in nm (0 = Auto)
if isfield(userInfo,'minSpotSize')
    minSpotSize = userInfo.minSpotSize;
else
    minSpotSize = '0';
    warning('no sharp parameter found, using auto mode = 0')
end
% desired size of the structural elements for image morphological
% transformation
if isfield(userInfo,'elementSize')
    elementSize = userInfo.elementSize;
else
    elementSize = '3';
    warning('no structural element size found, using 3')
end
% size of boundary pad on the edges to be blackened out
if isfield(userInfo,'padSize')
    padSize = userInfo.padSize;    
else
    padSize = '0';    
end
% how many time points to analyze
if isfield(userInfo,'endTimePoint')
    endTimePoint = userInfo.endTimePoint;
else
    endTimePoint = '0';
    warning('using all time points from the file')
end
%starting time frame of analysis
if isfield(userInfo,'startTimePoint')
    startTimePoint = userInfo.startTimePoint;
else
    startTimePoint = '1';
    warning('using starting time point = 1')
end
%rough cutoff for nucleus diameter in [um]
if isfield(userInfo,'nucleusFeatureSize')
    nucleusFeatureSize = userInfo.nucleusFeatureSize;
else
    nucleusFeatureSize = '3';
    warning('no feature size reference found, using 3 microns')
end

%analyze 3d from the start?
if isfield(userInfo,'analyze3D')
    analyze3D = userInfo.analyze3D;
else
    analyze3D = 'No';
    warning('no 3d analysis instruction found, forcing 2d')
end
switch analyze3D
    case 'Yes'
        analyze3D = 1;
    case 'No'        
        analyze3D = 0;    
    otherwise
        analyze3D = 0;
end

serieMetaDataDS.imagingInfo.XYpsf = XYpsf;
serieMetaDataDS.imagingInfo.bicoidChannel = bicoidChannel;
serieMetaDataDS.imagingInfo.Zpsf = Zpsf;
serieMetaDataDS.analysisInfo.imUseType = imUseType;
serieMetaDataDS.analysisInfo.nucSmoothFilter = nucSmoothFilter;
serieMetaDataDS.analysisInfo.deconvolutionType = deconvolutionType;
serieMetaDataDS.analysisInfo.nucleusDetectMode = nucleusDetectMode;
serieMetaDataDS.analysisInfo.nucThres = nucThres;
serieMetaDataDS.analysisInfo.smoothingParam = smoothingParam;
serieMetaDataDS.analysisInfo.minSpotSize = minSpotSize;
serieMetaDataDS.analysisInfo.clearBorder = clearBorder;
serieMetaDataDS.analysisInfo.elementSize = elementSize;
serieMetaDataDS.analysisInfo.padSize = padSize;
serieMetaDataDS.analysisInfo.nucleusFeatureSize = nucleusFeatureSize;
serieMetaDataDS.analysisInfo.endTimePoint = endTimePoint;
serieMetaDataDS.analysisInfo.startTimePoint = startTimePoint;

totalTimePoints = endTimePoint - startTimePoint + 1; 
serieMetaDataDS.analysisInfo.totalTimePoints = totalTimePoints; 

serieMetaDataDS.analysisInfo.spotDetect = spotDetect;
serieMetaDataDS.analysisInfo.spotFilter = spotFilter;
serieMetaDataDS.imagingInfo.timeResolution = DeltaT;
serieMetaDataDS.analysisInfo.analyze3D = analyze3D;
% by how much the centroid shift of each nucleus between consecutive frames
% be allowed (in pixels)
if isfield(userInfo,'centroidShiftConstraint')
    centroidShiftConstraint = userInfo.centroidShiftConstraint;
else
    centroidShiftConstraint = ceil(nucleusFeatureSize/(2*userInfo.XpixelSize));
    warning('no centroid shift constraint found, using calculated = %d',centroidShiftConstraint)
end
serieMetaDataDS.analysisInfo.centroidShiftConstraint = centroidShiftConstraint;
% total time for which the hotspots are to be analyzed (in seconds)
if isfield(userInfo,'analysisTime')
    analysisTime = userInfo.analysisTime;
else
    analysisTime = floor(endTimePoint/DeltaT);    
    warning('no hotspot analysis time span found, using %d seconds', analysisTime);
end
% total time for which the hotspots are to be analyzed (in seconds)
if isfield(userInfo,'timeWindow')
    timeWindow = userInfo.timeWindow;
else
    timeWindow = floor(analysisTime/3);    
    warning('no moving time window for hotspot analysis found, using %d seconds', timeWindow);
end
% slide side of moving time window for hotspot analysis (in seconds)
if isfield(userInfo,'timeWindow')
    timeSlide = userInfo.timeSlide;
else
    timeSlide = '1';    
    warning('no time slide for hotspot analysis found, using %d seconds', timeSlide);
end
serieMetaDataDS.analysisInfo.analysisTime = analysisTime;
serieMetaDataDS.analysisInfo.timeWindow = timeWindow;
serieMetaDataDS.analysisInfo.timeSlide = timeSlide;

if startZ==0 && endZ==0
    startZ = 1;
    endZ = Zslices;
    nZslices = Zslices;
elseif startZ==0 && endZ>0
    startZ = 1;
    if endZ<Zslices
        nZslices = endZ;
    else
    endZ = Zslices;
    nZslices = Zslices;
    end
elseif startZ>0 && endZ>0
    if startZ<Zslices && startZ<endZ && endZ<=Zslices
        nZslices = endZ-startZ+1;
    else
        error('check starting and or ending z slice');
    end
end
serieMetaDataDS.imagingInfo.startZ = startZ;
serieMetaDataDS.imagingInfo.endZ = endZ;
serieMetaDataDS.analysisInfo.nZslices = nZslices; 

nucleusFeatureSize = serieMetaDataDS.analysisInfo.nucleusFeatureSize;% already in microns
minNucAreaPix = ceil(3*(nucleusFeatureSize^2)/(xPixUM*yPixUM)); %all in microns
serieMetaDataDS.analysisInfo.minNucSize = minNucAreaPix;%%%%%%%% Assignment in pixels here

if (2*nucleusFeatureSize)<=(nZslices*zPixUM)
    minNucVolPix = ceil(4*nucleusFeatureSize^3/(xPixUM*yPixUM*zPixUM)); %all in microns
elseif (2*nucleusFeatureSize)>(nZslices*zPixUM)
    minNucVolPix = ceil(minNucAreaPix*(2*nucleusFeatureSize/zPixUM));
end
serieMetaDataDS.analysisInfo.minNucVol = ceil(minNucVolPix);

minSpotSizeUM = (minSpotSize/1000);
serieMetaDataDS.analysisInfo.minSpotSizeUM = minSpotSizeUM;

mrnaSpotSizeFactor = 1.5;%%%%%%%%%%%%%%%%%%%%%%
if minSpotSizeUM==0
    minSpotAreaPix = 0;
    minMRNAspotAreaPix = mrnaSpotSizeFactor*ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
elseif minSpotSizeUM>0 && minSpotSizeUM<xyPsfUM% Convert to pixels here
    minSpotAreaPix = ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
else
    minSpotAreaPix = ceil((pi*(minSpotSizeUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
end
serieMetaDataDS.analysisInfo.minSpotArea = ceil(minSpotAreaPix); %%%%%%%% assignment in pixels
serieMetaDataDS.analysisInfo.minMRNAspotArea = ceil(minMRNAspotAreaPix); %%%%%%%% assignment in pixels

if zPixUM<zPsfUM
    minSpotVolPix = ceil(minSpotAreaPix)*ceil(zPsfUM/zPixUM);%all in microns
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
else
    minSpotVolPix = 2*ceil(minSpotAreaPix);
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
end
serieMetaDataDS.analysisInfo.minSpotVol = ceil(minSpotVolPix); %%%%%%%% assignment in pixels
serieMetaDataDS.analysisInfo.minMRNAspotVol = ceil(minMRNAspotVolPix); %%%%%%%% assignment in pixels
end