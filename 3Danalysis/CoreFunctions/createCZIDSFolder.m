function createCZIDSFolder(seriesMetaDataDS)
% rawImagesFolderName = seriesMetaDataDS.expInfo.rawImagesFolderName;
% c1ImagesFolderName = append(rawImagesFolderName,filesep, 'c1');
% c2ImagesFolderName = append(rawImagesFolderName,filesep, 'c2');
% c3ImagesFolderName = append(rawImagesFolderName,filesep, 'c3');
% seriesMetaDataDS.expInfo.c1ImagesFolderName = c1ImagesFolderName; %%% assignment
% seriesMetaDataDS.expInfo.c2ImagesFolderName = c2ImagesFolderName; %%% assignment
% seriesMetaDataDS.expInfo.c3ImagesFolderName = c3ImagesFolderName; %%% assignment
procfname = seriesMetaDataDS.expInfo.procImagesFolderName;

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
% if ~exist(rawImagesFolderName,'file')
%     mkdir(rawImagesFolderName);
%     mkdir(c1ImagesFolderName);
%     mkdir(c2ImagesFolderName);
%     mkdir(c3ImagesFolderName);
% else
%     if ~exist(c1ImagesFolderName,'file')
%         mkdir(c1ImagesFolderName);
%     else
%         dirInfo = dir(c1ImagesFolderName);
%         dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
%         fileNames = fullfile(c1ImagesFolderName, {dirInfo.name});
%         subdirList = dirInfo([dirInfo.isdir]);
%         subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
%         for iDir = 1:numel(subdirList)
%             delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
%             rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
%         end
%         if size(dirInfo, 1)>0
%             delete (fileNames{:});
%         end        
%     end
%     if ~exist(c2ImagesFolderName,'file')
%         mkdir(c2ImagesFolderName);
%     else
%         dirInfo = dir(c2ImagesFolderName);
%         dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));        
%         fileNames = fullfile(c2ImagesFolderName, {dirInfo.name});
%         subdirList = dirInfo([dirInfo.isdir]);
%         subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
%         for iDir = 1:numel(subdirList)
%             delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
%             rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
%         end
%         if size(dirInfo, 1)>0
%             delete (fileNames{:})
%         end
%     end
%     if ~exist(c3ImagesFolderName,'file')
%         mkdir(c3ImagesFolderName);
%     else
%         dirInfo = dir(c3ImagesFolderName);
%         dirInfo = dirInfo(~cellfun('isempty', {dirInfo.name}));
%         fileNames = fullfile(c2ImagesFolderName, {dirInfo.name});
%         subdirList = dirInfo([dirInfo.isdir]);
%         subdirList = subdirList(~ismember({subdirList.name}, {'.', '..'}));
%         for iDir = 1:numel(subdirList)
%             delete(fullfile(subdirList(iDir).folder, subdirList(iDir).name, '\*'));
%             rmdir(fullfile(subdirList(iDir).folder, subdirList(iDir).name));
%         end
%         if size(dirInfo, 1)>0
%             delete (fileNames{:})
%         end
%     end
%         
% end

end