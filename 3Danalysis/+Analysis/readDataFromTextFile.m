function readDataFromTextFile(textFileFolder)
dirInfo = dir(textFileFolder);
dirInfo([dirInfo.isdir]) = [];
totalFiles = size(dirInfo, 1);
fileNames = cell(1, totalFiles);
dataCell = cell(1, totalFiles);
for i = 1:totalFiles
    fileNames{i} = fullfile(dirInfo(i).folder, dirInfo(i).name);
    fileID = fopen(fileNames{i},'rt');
    data = textscan(fileID,'%f%f','Delimiter',{'\t'}, 'CollectOutput',1);
    data = data{1};
    data1 = data(:,1);
    data1 = data1(~isnan(data1));
    data2 = data(:,2);
    data2 = data2(~isnan(data2));
    data = horzcat(data1, data2);
    [dataMin, minIdx] = min(data2([1: ceil(length(data2)/2)]));
    data1 = data1 - data1(minIdx);
    fclose(fileID);
    dataCell{i} = horzcat(data1, data2);
end

figure;
hold on;
cellfun(@(x) plot(x(:,1), x(:,2)), dataCell, 'un', 0);

txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalFiles = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
end

% fid = fopen('D:\Tyrone_analysis\FRAP\20210228_bcd2xa_nc13_frap_1.txt','rt');
% D = textscan(fid,'%f32%f32%f32%f32','Delimiter',{'\t'}, 'CollectOutput',1);
% fclose(fid);
% Data = cell2mat(D);
% Data1 = Data(:,[1 2]);                          % Columns 1 & 2
% Data2 = Data(:,[3 4]);                          % Columns 3 & 4