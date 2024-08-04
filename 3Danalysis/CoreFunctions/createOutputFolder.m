function [metaDataDS] = createOutputFolder(metaDataDS)
rawImagesFolderName = metaDataDS.expInfo.rawImagesFolderName;
procfname = metaDataDS.expInfo.procImagesFolderName;
rawDataFilePath = append(metaDataDS.expInfo.rawDataFilePath);
endTimePoint = metaDataDS.analysisInfo.endTimePoint;
startTimePoint = metaDataDS.analysisInfo.startTimePoint;
colorChannels = metaDataDS.imagingInfo.colorChannels;
zSlices = metaDataDS.imagingInfo.Zslices;
startZ = metaDataDS.imagingInfo.startZ;
endZ = metaDataDS.imagingInfo.endZ;
analyze3D = metaDataDS.analysisInfo.analyze3D;

c1ImagesFolderName = append(rawImagesFolderName,filesep, 'c1');
metaDataDS.expInfo.c1ImagesFolderName = c1ImagesFolderName;%%% assignment
if colorChannels>1
    c2ImagesFolderName = append(rawImagesFolderName,filesep, 'c2');    
    metaDataDS.expInfo.c2ImagesFolderName = c2ImagesFolderName;%%% assignment
elseif colorChannels>2
    c3ImagesFolderName = append(rawImagesFolderName,filesep, 'c3');    
    metaDataDS.expInfo.c3ImagesFolderName = c3ImagesFolderName;%%% assignment
elseif colorChannels>3
    c4ImagesFolderName = append(rawImagesFolderName,filesep, 'c4');    
    metaDataDS.expInfo.c4ImagesFolderName = c4ImagesFolderName;%%% assignment
end

if ~exist(procfname,'file')
    mkdir(procfname);
else
    dirInfo = dir(procfname);
    dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
    dirInfo([dirInfo.isdir]) = [];% skip subdirectories
    fileNames = fullfile(procfname, {dirInfo.name});
    if size(dirInfo, 1)>0
        delete (fileNames{:})
    end
end
if ~exist(rawImagesFolderName,'file')
    mkdir(rawImagesFolderName);
    mkdir(c1ImagesFolderName);
    if colorChannels>1
        mkdir(c2ImagesFolderName);
    elseif colorChannels>2
        mkdir(c3ImagesFolderName);
    elseif colorChannels>3
        mkdir(c4ImagesFolderName);
    end
else
    if ~exist(c1ImagesFolderName,'file')
        mkdir(c1ImagesFolderName);
    else
        dirInfo = dir(c1ImagesFolderName);
        dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
        fileNames = fullfile(c1ImagesFolderName, {dirInfo.name});
        subdirList = dirInfo([dirInfo.isdir]);
        subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
        for iDir = 1:numel(subdirList)
            delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
            rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
        end
        if size(dirInfo, 1)>0
            cellfun(@delete,  fileNames(~ismember({dirInfo.name}, {'.', '..'})));
%             delete (fileNames{:});            
        end        
    end    
    
    if colorChannels>1 && ~exist(c2ImagesFolderName,'file')
        mkdir(c2ImagesFolderName);
    elseif colorChannels>1
        dirInfo = dir(c2ImagesFolderName);
        dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));        
        fileNames = fullfile(c2ImagesFolderName, {dirInfo.name});
        subdirList = dirInfo([dirInfo.isdir]);
        subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
        for iDir = 1:numel(subdirList)
            delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
            rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
        end
        if size(dirInfo, 1)>0
           cellfun(@delete,  fileNames(~ismember({dirInfo.name}, {'.', '..'})));
        end
    end
    
    if colorChannels>2 && ~exist(c3ImagesFolderName,'file')
        mkdir(c3ImagesFolderName);
    elseif colorChannels>2
        dirInfo = dir(c3ImagesFolderName);
        dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
        fileNames = fullfile(c3ImagesFolderName, {dirInfo.name});
        subdirList = dirInfo([dirInfo.isdir]);
        subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
        for iDir = 1:numel(subdirList)
            delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
            rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
        end
        if size(dirInfo, 1)>0
           cellfun(@delete,  fileNames(~ismember({dirInfo.name}, {'.', '..'})));
        end
    end
        
    if colorChannels>3 && ~exist(c4ImagesFolderName,'file')
        mkdir(c4ImagesFolderName);
    elseif colorChannels>3
        dirInfo = dir(c4ImagesFolderName);
        dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
        fileNames = fullfile(c4ImagesFolderName, {dirInfo.name});
        subdirList = dirInfo([dirInfo.isdir]);
        subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
        for iDir = 1:numel(subdirList)
            delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
            rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
        end
        if size(dirInfo, 1)>0
           cellfun(@delete,  fileNames(~ismember({dirInfo.name}, {'.', '..'})));
        end
    end    
end

if analyze3D==0 % Single plane
    options.append = 0;
    options.message = 0;
    for t = startTimePoint:endTimePoint
        timePoint = startTimePoint-1+t;
        for z = startZ:endZ
            zPoint = startZ-1+z;
            framePoint = colorChannels*((timePoint-1)*zSlices)+ colorChannels*(zPoint-1)+1;
            c1Frame = imread(rawDataFilePath, framePoint);
            fname = append(c1ImagesFolderName, filesep,'c1_t',...
                num2str((t-startTimePoint+1),'%04.f'), '.tif');            
            HelperFunctions.saveTiff((c1Frame),fname, options);
            if colorChannels==2 || colorChannels==3 || colorChannels==4
                framePoint = framePoint+1;
                c2Frame = imread(rawDataFilePath, framePoint);            
                fname = append(c2ImagesFolderName, filesep,'c2_t',...
                    num2str((t-startTimePoint+1),'%04.f'), '.tif');
                HelperFunctions.saveTiff((c2Frame),fname, options);
            end
            if colorChannels==3 || colorChannels==4
                framePoint = framePoint+2;
                c3Frame = imread(rawDataFilePath, framePoint);
                fname = append(c3ImagesFolderName, filesep,'c3_t',...
                    num2str((t-startTimePoint+1),'%04.f'), '.tif');
                HelperFunctions.saveTiff((c3Frame),fname, options);
            end 
            if colorChannels==4
                framePoint = framePoint+3;
                c4Frame = imread(rawDataFilePath, framePoint);
                fname = append(c4ImagesFolderName, filesep,'c4_t',...
                    num2str((t-startTimePoint+1),'%04.f'), '.tif');
                HelperFunctions.saveTiff((c4Frame),fname, options);
            end 
        end
    end
else % 3D image (works for single plane too)
    stackInfo = imfinfo(rawDataFilePath);
    c1Frame = zeros(stackInfo(1).Height, stackInfo(1).Width, (endZ - startZ+1));
    if colorChannels>1
        c2Frame = zeros(stackInfo(1).Height, stackInfo(1).Width, (endZ - startZ+1));
    elseif colorChannels>2
        c3Frame = zeros(stackInfo(1).Height, stackInfo(1).Width, (endZ - startZ+1));
    elseif colorChannels>3
        c4Frame = zeros(stackInfo(1).Height, stackInfo(1).Width, (endZ - startZ+1));
    end
    options.append = 1;
    options.message = 0;
    for t = startTimePoint:endTimePoint
        timePoint = startTimePoint-1+t;
        for z = startZ:endZ
            zPoint = startZ-1+z;
            framePoint = colorChannels*((timePoint-1)*zSlices) + colorChannels*(z-1)+1;
            c1Frame(:,:,(z-startZ+1)) = imread(rawDataFilePath, framePoint);
            if colorChannels>1
                framePoint = framePoint+1;
                c2Frame(:,:,(z-startZ+1)) = imread(rawDataFilePath, framePoint);            
            end
            if colorChannels>2
                framePoint = framePoint+1;
                c3Frame(:,:,(z-startZ+1)) = imread(rawDataFilePath, framePoint);
            end 
            if colorChannels>3
                framePoint = framePoint+1;
                c4Frame(:,:,(z-startZ+1)) = imread(rawDataFilePath, framePoint);
            end 
        end
        fname = append(c1ImagesFolderName, filesep,'c1_t',...
            num2str((t-startTimePoint+1),'%04.f'),  '.tif');
        HelperFunctions.saveTiff(uint16((c1Frame)),fname, options);

        if colorChannels>1
            fname = append(c2ImagesFolderName, filesep,'c2_t',...
                num2str((t-startTimePoint+1),'%04.f'), '.tif');
            HelperFunctions.saveTiff(uint16((c2Frame)),fname, options);
        end

        if colorChannels>2
            fname = append(c3ImagesFolderName, filesep,'c3_t',...
                num2str((t-startTimePoint+1),'%04.f'), '.tif');
            HelperFunctions.saveTiff(uint16(c3Frame),fname, options);
        end    
        
        if colorChannels>3
            fname = append(c3ImagesFolderName, filesep,'c4_t',...
                num2str((t-startTimePoint+1),'%04.f'), '.tif');
            HelperFunctions.saveTiff(uint16(c4Frame),fname, options);
        end    
    end
end
end