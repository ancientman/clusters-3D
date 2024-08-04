function [metaDataDS] = makeMetaDataDSNew(analysisFolder)

% Where is the input text file
expInfoFileName = append(analysisFolder, filesep,'expInfo.txt');
expInfoFileChar = convertStringsToChars(expInfoFileName);
% Read the input text file
if exist(expInfoFileName, 'file')
    info = HelperFunctions.readtext(expInfoFileChar,'\t');
    userInfo = cell2struct(info(:,2),info(:,1),1);
else
    error('User must supply expInfo.txt in the analysis directory')
end
% Where is the raw file stored (read from the input)
if isfield(userInfo,'rawDataFilePath')
    metaDataDS.expInfo.rawDataFilePath = userInfo.rawDataFilePath;
    [rawFilePath,rawFileName,rawFileExt] = fileparts(userInfo.rawDataFilePath);
else
    error('expInfo.txt should include rawDataFilePath')
end
% Assign all the folder names
if isfield(userInfo,'DSFilePath')
    metaDataDS.expInfo.analysisFolder = userInfo.DSFilePath;
    metaDataDS.expInfo.rawImagesFolderName = strcat(metaDataDS.expInfo.analysisFolder, filesep, 'out');
    metaDataDS.expInfo.procImagesFolderName = strcat(metaDataDS.expInfo.analysisFolder, filesep, 'DS_', rawFileName);
    else
    error('expInfo.txt should include valid DSFilePath')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read OME metadata
reader = bfGetReader(metaDataDS.expInfo.rawDataFilePath);
omeMeta = reader.getMetadataStore();
omeGlobalMeta = reader.getGlobalMetadata();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% what is the bit depth of the acquisition
bitDepth = str2double(omeGlobalMeta.get(strcat('BitsPerPixel')));
% what is the time resolution
DeltaT = str2double(omeGlobalMeta.get(strcat('Information|Image|Channel|LaserScanInfo|FrameTime #1')));
% what is the frameTime

%how many x sections in the image
X = str2double(omeGlobalMeta.get(strcat('SizeX')));
%how many y sections in the image
Y = str2double(omeGlobalMeta.get(strcat('SizeY')));
%how many z sections in the image
imageZslices = str2double(omeGlobalMeta.get(strcat('SizeZ')));
%how many t sections in the image
T = str2double(omeGlobalMeta.get(strcat('SizeT')));
%how many color channels are used
colorChannels = str2double(omeGlobalMeta.get(strcat('SizeC')));

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

%which is the filled channel
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
    endTimePoint = 0;
    warning('using all time points from the file')
end
%starting time frame of analysis
if isfield(userInfo,'startTimePoint')
    startTimePoint = userInfo.startTimePoint;
else
    startTimePoint = 1;
    warning('using starting time point = 1')
end
% how many time points to analyze
if isfield(userInfo,'overrideTimePoints')
    overrideTimePoints = userInfo.overrideTimePoints;
else
    overrideTimePoints = 0;
    warning('using all time points from the file')
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

% Data from image metadata
metaDataDS.imagingInfo.bitDepth = bitDepth;
metaDataDS.imagingInfo = struct('X',X,'Y',Y,'T',T,'Zslices',imageZslices,'DeltaT',DeltaT);
metaDataDS.imagingInfo.timeResolution = DeltaT;
xPixUM = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER)); % in µm
metaDataDS.analysisInfo.xPixUM = xPixUM;
yPixUM = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER)); % in µm
metaDataDS.analysisInfo.yPixUM = yPixUM;
if imageZslices>1 && ~isempty(omeMeta.getPixelsPhysicalSizeZ(0))
    zPixUM =double( omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER)); % in µm
else
    zPixUM = 1e6*str2double(omeGlobalMeta.get(strcat('Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1')));
end
metaDataDS.analysisInfo.zPixUM = zPixUM;

metaDataDS.imagingInfo.colorChannels = colorChannels;
metaDataDS.analysisInfo.startTimePoint = startTimePoint;
if endTimePoint == 0 || endTimePoint>T
    endTimePoint = T;
end
metaDataDS.analysisInfo.endTimePoint = endTimePoint;
totalTimePoints = endTimePoint - startTimePoint + 1; 
if overrideTimePoints~=0 && overrideTimePoints<totalTimePoints
    totalTimePoints = overrideTimePoints;
end
metaDataDS.analysisInfo.totalTimePoints = totalTimePoints; 

%data from input text file
metaDataDS.imagingInfo.XYpsf = XYpsf;
metaDataDS.imagingInfo.Zpsf = Zpsf;
xyPsfUM = metaDataDS.imagingInfo.XYpsf/1000;% changed to microns
zPsfUM = metaDataDS.imagingInfo.Zpsf/1000;% changed to microns
metaDataDS.analysisInfo.xyPsfUM = xyPsfUM;
metaDataDS.analysisInfo.zPsfUM = zPsfUM;
if bicoidChannel<=colorChannels
    metaDataDS.imagingInfo.bicoidChannel = bicoidChannel;
else
    error('check bicoid channel assignment')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis parameters
metaDataDS.analysisInfo.imUseType = imUseType;
metaDataDS.analysisInfo.nucSmoothFilter = nucSmoothFilter;
metaDataDS.analysisInfo.deconvolutionType = deconvolutionType;
metaDataDS.analysisInfo.nucleusDetectMode = nucleusDetectMode;
metaDataDS.analysisInfo.nucThres = nucThres;
metaDataDS.analysisInfo.smoothingParam = smoothingParam;
metaDataDS.analysisInfo.minSpotSize = minSpotSize;
metaDataDS.analysisInfo.clearBorder = clearBorder;
metaDataDS.analysisInfo.elementSize = elementSize;
metaDataDS.analysisInfo.padSize = padSize;
metaDataDS.analysisInfo.nucleusFeatureSize = nucleusFeatureSize;
metaDataDS.analysisInfo.spotDetect = spotDetect;
metaDataDS.analysisInfo.spotFilter = spotFilter;
metaDataDS.analysisInfo.analyze3D = analyze3D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by how much the centroid shift of each nucleus between consecutive frames
% be allowed (in pixels)
if isfield(userInfo,'centroidShiftConstraint')
    centroidShiftConstraint = userInfo.centroidShiftConstraint;
else
    centroidShiftConstraint = ceil(nucleusFeatureSize/(2*userInfo.XpixelSize));
    warning('no centroid shift constraint found, using calculated = %d',centroidShiftConstraint)
end
metaDataDS.analysisInfo.centroidShiftConstraint = centroidShiftConstraint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
metaDataDS.analysisInfo.analysisTime = analysisTime;
metaDataDS.analysisInfo.timeWindow = timeWindow;
metaDataDS.analysisInfo.timeSlide = timeSlide;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    warning('using ending z slice = %d',imageZslices);
end

if startZ==0 && endZ==0
    startZ = 1;
    endZ = imageZslices;
    useZslices = imageZslices;
elseif startZ==0 && endZ>0
    startZ = 1;
    if endZ<imageZslices
        useZslices = endZ;
    else
    endZ = imageZslices;
    useZslices = imageZslices;
    end
elseif startZ>0 && endZ==0
    if startZ<=imageZslices 
        endZ = imageZslices;
        useZslices = endZ-startZ+1;
    end
elseif startZ>0 && endZ>0
    if startZ<imageZslices && startZ<endZ && endZ<=imageZslices
        useZslices = endZ-startZ+1;
    elseif startZ==endZ && startZ==imageZslices
        useZslices = imageZslices;
    else        
        error('check starting and or ending z slice');
    end
end

metaDataDS.imagingInfo.startZ = startZ;
metaDataDS.imagingInfo.endZ = endZ;
metaDataDS.analysisInfo.nZslices = useZslices; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum nuclear volume to be considered
nucleusFeatureSize = metaDataDS.analysisInfo.nucleusFeatureSize;% already in microns
minNucAreaPix = ceil(3*(nucleusFeatureSize^2)/(xPixUM*yPixUM)); %all in microns
metaDataDS.analysisInfo.minNucSize = minNucAreaPix;%%%%%%%% Assignment in pixels here

if (2*nucleusFeatureSize)<=(useZslices*zPixUM)
    minNucVolPix = ceil(4*nucleusFeatureSize^3/(xPixUM*yPixUM*zPixUM)); %all in microns
elseif (2*nucleusFeatureSize)>(useZslices*zPixUM)
%     minNucVolPix = ceil(minNucAreaPix*(2*nucleusFeatureSize/zPixUM));
    minNucVolPix = ceil(minNucAreaPix*(useZslices));
end
metaDataDS.analysisInfo.minNucVol = ceil(minNucVolPix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum spot size to be considered
minSpotSizeUM = (minSpotSize/1000);
metaDataDS.analysisInfo.minSpotSizeUM = minSpotSizeUM;

if isfield(userInfo,'mrnaSpotSizeFactor')
    mrnaSpotSizeFactor = userInfo.mrnaSpotSizeFactor;
else
    mrnaSpotSizeFactor = 3;
    warning('using mrnaSpotSizeFactor = 3');
end

if minSpotSizeUM==0
    minSpotAreaPix = 0;
    minMRNAspotAreaPix = mrnaSpotSizeFactor*ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
elseif minSpotSizeUM>0 && minSpotSizeUM<0%xyPsfUM% Convert to pixels here
    minSpotAreaPix = ceil((pi*(xyPsfUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
else
    minSpotAreaPix = ceil((pi*(minSpotSizeUM/2)^2)/(xPixUM*yPixUM));
    minMRNAspotAreaPix = mrnaSpotSizeFactor*minSpotAreaPix;
end
metaDataDS.analysisInfo.minSpotArea = ceil(minSpotAreaPix); %%%%%%%% assignment in pixels
metaDataDS.analysisInfo.minMRNAspotArea = ceil(minMRNAspotAreaPix); %%%%%%%% assignment in pixels

if zPixUM<zPsfUM
    zSections = min([floor(zPsfUM/zPixUM), useZslices]);
    minSpotVolPix = ceil(minSpotAreaPix)*zSections;%all in microns
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
else
    minSpotVolPix = 2*ceil(minSpotAreaPix);
    minMRNAspotVolPix = mrnaSpotSizeFactor*minSpotVolPix;
end
metaDataDS.analysisInfo.minSpotVol = ceil(minSpotVolPix); %%%%%%%% assignment in pixels
metaDataDS.analysisInfo.minMRNAspotVol = ceil(minMRNAspotVolPix); %%%%%%%% assignment in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if useZslices>1 && analyze3D==0
    error('turn on analyze 3D');
elseif useZslices==1 && analyze3D==1
    error('turn off analyze 3D');    
end
end