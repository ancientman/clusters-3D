function c1ValC2SpotTime(txtFilePath, geneName)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   uses two colour data to calculate distance between the 
%   transciption spot and the bicoid spots.
%   also calculates the diameter of the bicoid spots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalKMeans = 1;
kMeansDistLimit = 20;

cropLen = 3; % use 3 closest points hardcoded

colorStruct = cell(1, 9);
colorStruct{1} =  [120, 120, 120; 0, 0, 0];
colorStruct{2} = [31,120,180; 0, 0, 0];
colorStruct{3} = [49,54,149; 0, 0, 0];
colorStruct{4} = [51,160,44;0, 0, 0];
colorStruct{5} = [251,154,153; 0, 0, 0];
colorStruct{6} = [253,191,111; 0, 0, 0];
colorStruct{7} = [166,206,227; 0, 0, 0];
colorStruct{8} = [202,178,214; 0, 0, 0];
colorStruct{9} = [178,223,138;0, 0, 0];

txtFilePathChar = convertStringsToChars(txtFilePath);
if exist(txtFilePath, 'file')
    info = HelperFunctions.readtext(txtFilePathChar,'\t');
    totalSubDirs = size(info, 1);
    dataFiles = cell2struct(info(:,2), info(:,1),1);
else
    error('data folder is empty')
end
nucStruct = cell(1, totalSubDirs);
c1c2ValStruct = cell(1, totalSubDirs);
timePoints = cell(1, totalSubDirs);
nNuc = cell(1, totalSubDirs);
c1ValC2time = cell(1, totalSubDirs);
c1ValC2timeXY = cell(1, totalSubDirs);
c1ValC2timeZ = cell(1, totalSubDirs);
c1ValC2timeNoBg = cell(1, totalSubDirs);
nucMeanVal = cell(1, totalSubDirs);
c2SpotVal = cell(1, totalSubDirs);
c2SpotVol = cell(1, totalSubDirs);
c2SpotMol = cell(1, totalSubDirs);
c2MolTimeNuc = cell(1, totalSubDirs);
c1ValC2timeNuc = cell(1, totalSubDirs);
c1ValC2timeNucXY = cell(1, totalSubDirs);
c1ValC2timeNucZ = cell(1, totalSubDirs);
c1ValC2timeNoBgNuc = cell(1, totalSubDirs);
nucMeanTime = cell(1, totalSubDirs);
c2spotMolNuc = cell(1, totalSubDirs);
c1ValC2timeNucNoBg = cell(1, totalSubDirs);
c1ValC2timeMean = cell(1, totalSubDirs);
c1ValC2RadMean = cell(1, totalSubDirs);
c1ValC2NoBgRadNorm = cell(1, totalSubDirs);

c1ValC2RadDiff = cell(1, totalSubDirs);
c1ValC2RadDiffXY = cell(1, totalSubDirs);
c1ValC2RadDiffZ = cell(1, totalSubDirs);

c1ValC2RadNorm = cell(1, totalSubDirs);
c1ValC2RadNormXY = cell(1, totalSubDirs);
c1ValC2RadNormZ = cell(1, totalSubDirs);
c1ValC2timeMeanNoBg = cell(1, totalSubDirs);
c1ValC2timeSem = cell(1, totalSubDirs);
c1ValC2RadSem = cell(1, totalSubDirs);
c1ValC2timeSemNoBg = cell(1, totalSubDirs);
c2Strength = cell(1, totalSubDirs);
c1ValMean = cell(1, totalSubDirs);
c1ValSem = cell(1, totalSubDirs);
radBin = cell(1, totalSubDirs);
c2SpotStruct = cell(1, totalSubDirs);
c1ValC2Av = cell(1, totalSubDirs);
c1ValC2Sem = cell(1, totalSubDirs);

init = 1;

for i=init:totalSubDirs
    fileID = append('file_', num2str(i,'%03d'));
    fileName = dataFiles.(fileID);
    dirInfo = dir(fileName);
    dirInfo([dirInfo.isdir]) = [];
    if ~isempty({dirInfo.name})
        c1c2ValStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1c2IntensityPropDS.mat'));
        c2SpotStruct{i} = load(append(dataFiles.(fileID), filesep, 'c2SpotPropDS.mat'));
        nucStruct{i} = load(append(dataFiles.(fileID), filesep, 'c1NucPropDS.mat'));
        try
        nNuc{i} = length(nucStruct{i}.c1NucProp{1}); 
        catch
            aaa = 1
        end
        c1ValC2time{i} = cell(1, nNuc{i});
        c1ValC2timeXY{i} = cell(1, nNuc{i});
        c1ValC2timeZ{i} = cell(1, nNuc{i});
        c1ValC2timeNoBg{i} = cell(1, nNuc{i});
        nucMeanVal{i} = cell(1, nNuc{i});
        c2SpotVal{i} = cell(1, nNuc{i});
        c2SpotVol{i} = cell(1, nNuc{i});
        c2SpotMol{i} = cell(1, nNuc{i});
        c1ValC2timeNoBgNuc{i} = cell(1, nNuc{i});
        nucMeanTime{i} = cell(1, nNuc{i});
        c2MolTimeNuc{i} = cell(1, nNuc{i});
        c1ValC2timeNuc{i} = cell(1, nNuc{i});
        c1ValC2timeNucXY{i} = cell(1, nNuc{i});
        c1ValC2timeNucZ{i} = cell(1, nNuc{i});
        c2spotMolNuc{i} = cell(1, nNuc{i});
        c1ValC2timeNucNoBg{i} = cell(1, nNuc{i});
        c1ValC2timeMean{i} = cell(1, nNuc{i});
        c1ValC2RadMean{i} = cell(1, nNuc{i});
        c1ValC2NoBgRadNorm{i} = cell(1, nNuc{i});        
        c1ValC2RadDiff{i} = cell(1, nNuc{i});
        c1ValC2RadDiffXY{i} = cell(1, nNuc{i});
        c1ValC2RadDiffZ{i} = cell(1, nNuc{i});        
        c1ValC2RadNorm{i} = cell(1, nNuc{i});
        c1ValC2RadNormXY{i} = cell(1, nNuc{i});
        c1ValC2RadNormZ{i} = cell(1, nNuc{i});
        c1ValC2timeMeanNoBg{i} = cell(1, nNuc{i});
        c1ValC2timeSem{i} = cell(1, nNuc{i});
        c1ValC2RadSem{i} = cell(1, nNuc{i});
        c1ValC2timeSemNoBg{i} = cell(1, nNuc{i});
        timePoints{i} = length(c2SpotStruct{i}.c2SpotProp);
        
        for n = 1:nNuc{i}       
%             if ~isempty(c1c2ValStruct{i}.c1c2IntensityProp{n})
%                 timePoints{i} = length(c1c2ValStruct{i}.c1c2IntensityProp{n}.valNoBgDistProp);                    
%             end            
            c1ValC2time{i}{n} = cell(1, timePoints{i});
            c1ValC2timeXY{i}{n} = cell(1, timePoints{i});
            c1ValC2timeZ{i}{n} = cell(1, timePoints{i});
            c1ValC2timeNoBg{i}{n} = cell(1, timePoints{i});
            nucMeanVal{i}{n} = cell(1, timePoints{i});
            c2SpotVal{i}{n} = cell(1, timePoints{i});
            c2SpotVol{i}{n} = cell(1, timePoints{i});
            c2SpotMol{i}{n} = cell(1, timePoints{i});
            for t = 1:timePoints{i}
                if ~isempty(c1c2ValStruct{i}.c1c2IntensityProp{n})
                    if ~isempty(c1c2ValStruct{i}.c1c2IntensityProp{n}.valDistProp3D{t}) && ~isempty(nucStruct{i}.c1NucProp{t})
                        nucMeanVal{i}{n}{t} = mean(nucStruct{i}.c1NucProp{t}(n).voxVal{:});
                        c2SpotVal{i}{n}{t} = mean(c2SpotStruct{i}.c2SpotProp{t}(n).voxVal{:});
                        c2SpotVol{i}{n}{t} = c2SpotStruct{i}.c2SpotProp{t}(n).vol;
                        c2SpotMol{i}{n}{t} = c2SpotVal{i}{n}{t}.*c2SpotVol{i}{n}{t};
                        c1ValC2time{i}{n}{t} = (c1c2ValStruct{i}.c1c2IntensityProp{n}.valDistProp3D{t}.distValAv);
                        c1ValC2timeXY{i}{n}{t} = (c1c2ValStruct{i}.c1c2IntensityProp{n}.valDistPropXY{t}.distValAv);
                        c1ValC2timeZ{i}{n}{t} = (c1c2ValStruct{i}.c1c2IntensityProp{n}.valDistPropZ{t}.distValAv);
%                         c1ValC2time{i}{n}{t} = (c1c2ValStruct{i}.c1c2IntensityProp{n}.valDistProp3D{t}.distValAv);
                        c1ValC2timeNoBg{i}{n}{t} = rescale(c1c2ValStruct{i}.c1c2IntensityProp{n}.valBgSubDistProp3D{t}.distValAv);
                    else
                        c2SpotMol{i}{n}{t} = 0;
                        c1ValC2time{i}{n}{t} = zeros(20, 1); % hardcoded %%%%%%%%%%%%%%%%%%%%%%%%%%
                        c1ValC2timeNoBg{i}{n}{t} =  zeros(20, 1); % hardcoded %%%%%%%%%%%%%%%%%%%%%%%%%%;
                    end
                end
            end
%             c2spotMolNuc{i}{n} = horzcat(c2SpotMol{i}{n}{:}{:}); %% incomplete
            nullFlag = vertcat(c2SpotMol{i}{n}{:});
            kk = 0;
            for ff = 1:length(nullFlag)
                if nullFlag(ff) ~=0
                    kk = kk+1;
                    if kk ==1
                        c1ValC2timeNoBgNuc{i}{n} = c1ValC2timeNoBg{i}{n}{ff};
                        
%                         c2MolTimeNuc{i}{n} = c2SpotVal{i}{n}{ff};
                        c2MolTimeNuc{i}{n} = c2SpotMol{i}{n}{ff};
                        c1ValC2timeNuc{i}{n} = c1ValC2time{i}{n}{ff};
                        c1ValC2timeNucXY{i}{n} = c1ValC2timeXY{i}{n}{ff};
                        c1ValC2timeNucZ{i}{n} = c1ValC2timeZ{i}{n}{ff};
                        nucMeanTime{i}{n} = nucMeanVal{i}{n}{ff};                        
                    else
                        c1ValC2timeNoBgNuc{i}{n} = horzcat(c1ValC2timeNoBgNuc{i}{n}, c1ValC2timeNoBg{i}{n}{ff});                        
                        c2MolTimeNuc{i}{n} = horzcat(c2MolTimeNuc{i}{n}, c2SpotMol{i}{n}{ff});
                        c1ValC2timeNuc{i}{n} = horzcat(c1ValC2timeNuc{i}{n}, c1ValC2time{i}{n}{ff});
                        c1ValC2timeNucXY{i}{n} = horzcat(c1ValC2timeNucXY{i}{n}, c1ValC2timeXY{i}{n}{ff});
                        c1ValC2timeNucZ{i}{n} = horzcat(c1ValC2timeNucZ{i}{n}, c1ValC2timeZ{i}{n}{ff});
                        nucMeanTime{i}{n} = horzcat(nucMeanTime{i}{n}, nucMeanVal{i}{n}{ff});
                    end
%                 c1ValC2timeNuc{i}{n} = horzcat(c1ValC2time{i}{n}{:});
%                 nucMeanTime{i}{n} = horzcat(nucMeanVal{i}{n}{:});
                end
            end
            if ~isempty(nucMeanTime{i}{n})
                c1ValC2NoBgRadNorm{i}{n} = c1ValC2timeNoBgNuc{i}{n};
%                 c1ValC2NoBgRadNorm{i}{n} = bsxfun(@rdivide, c1ValC2timeNoBgNuc{i}{n}, nucMeanTime{i}{n});
%                 c1ValC2RadNorm{i}{n} = c1ValC2timeNuc{i}{n};
                
                c1ValC2RadDiff{i}{n} = c1ValC2timeNuc{i}{n} - nucMeanTime{i}{n};
                c1ValC2RadNorm{i}{n} = bsxfun(@rdivide, c1ValC2RadDiff{i}{n}, nucMeanTime{i}{n});
                c1ValC2RadDiffXY{i}{n} = c1ValC2timeNucXY{i}{n} - nucMeanTime{i}{n};
                c1ValC2RadNormXY{i}{n} = bsxfun(@rdivide, c1ValC2RadDiffXY{i}{n}, nucMeanTime{i}{n});
                c1ValC2RadDiffZ{i}{n} = c1ValC2timeNucZ{i}{n} - nucMeanTime{i}{n};
                c1ValC2RadNormZ{i}{n} = bsxfun(@rdivide, c1ValC2RadDiffZ{i}{n}, nucMeanTime{i}{n});
                
%                 c1ValC2RadNorm{i}{n} = bsxfun(@rdivide, c1ValC2timeNuc{i}{n}, max(c1ValC2timeNuc{i}{n},[],1));             
%                 c1ValC2RadNorm{i}{n} = bsxfun(@rdivide, c1ValC2timeNuc{i}{n}, c1ValC2timeNuc{i}{n}(1,:)); % average over all radius bins (columns are times)
            end

            c1ValC2RadMean{i}{n} = mean(c1ValC2timeNuc{i}{n}, 1, 'omitnan'); % average over all radius bins (columns are times)
            c1ValC2RadSem{i}{n} = std(c1ValC2timeNuc{i}{n}, 0, 1, 'omitnan')./sqrt(size(c1ValC2timeNuc{i}{n}, 1));    
            c1ValC2timeMean{i}{n} = mean(nonzeros(c1ValC2timeNuc{i}{n}), 2, 'omitnan'); % Averaged over all time 
            c1ValC2timeSem{i}{n} = std(nonzeros(c1ValC2timeNuc{i}{n}), 0, 2, 'omitnan')./sqrt(size(nonzeros(c1ValC2timeNuc{i}{n}), 2));    
            c1ValC2timeNucNoBg{i}{n} = horzcat(c1ValC2timeNoBg{i}{n}{:});
            c1ValC2timeMeanNoBg{i}{n} = mean(c1ValC2timeNucNoBg{i}{n}, 2, 'omitnan');
            c1ValC2timeSemNoBg{i}{n} = std(c1ValC2timeNucNoBg{i}{n}, 0, 2, 'omitnan')./sqrt(size(c1ValC2timeNucNoBg{i}{n}, 2));    
            
            if i==1&&n==1
%                 c1ValC2RadNormAll = c1ValC2RadNorm{i}{n};
                nucMeanValCell{1} = nucMeanTime{i}{n};
                c2MolTimeNucAllCell{1} = c2MolTimeNuc{i}{n};
                c1ValC2NoBgRadNormAllCell{1} =  c1ValC2NoBgRadNorm{i}{n};
                
                c1ValC2RadDiffAllCell{1} =  c1ValC2RadDiff{i}{n};
                c1ValC2RadDiffAllCellXY{1} =  c1ValC2RadDiffXY{i}{n};
                c1ValC2RadDiffAllCellZ{1} =  c1ValC2RadDiffZ{i}{n};
                
                c1ValC2RadNormAllCell{1} =  c1ValC2RadNorm{i}{n};
                c1ValC2RadNormAllCellXY{1} =  c1ValC2RadNormXY{i}{n};
                c1ValC2RadNormAllCellZ{1} =  c1ValC2RadNormZ{i}{n};
            else
%                 c1ValC2RadNormAll = horzcat(c1ValC2RadNormAll, c1ValC2RadNorm{i}{n});
                nucMeanValCell{end+1} = nucMeanTime{i}{n};
                c2MolTimeNucAllCell{end+1} = c2MolTimeNuc{i}{n};
                c1ValC2NoBgRadNormAllCell{end+1} =  c1ValC2NoBgRadNorm{i}{n};
                
                c1ValC2RadDiffAllCell{end+1} =  c1ValC2RadDiff{i}{n};
                c1ValC2RadDiffAllCellXY{end+1} =  c1ValC2RadDiffXY{i}{n};
                c1ValC2RadDiffAllCellZ{end+1} =  c1ValC2RadDiffZ{i}{n};
                
                c1ValC2RadNormAllCell{end+1} =  c1ValC2RadNorm{i}{n};
                c1ValC2RadNormAllCellXY{end+1} =  c1ValC2RadNormXY{i}{n};
                c1ValC2RadNormAllCellZ{end+1} =  c1ValC2RadNormZ{i}{n};
            end
        end
    end
    
    %______________________________________________________________
%     figure('Color', 'w');
%     xArr = 93*(1:length(c2MolTimeNuc{i}{5})); % hard coded 62 seconds is assumed to be the frame time
%     plotDouble(xArr, c2MolTimeNuc{i}{5},  mean(c1ValC2RadNormXY{i}{5}(1:2, :), 'omitnan'), std(c1ValC2RadNormXY{i}{5}(1:2, :), 1, 'omitnan'));
%     letsPlot1(c1ValC2RadNormXY{i}{5});
%     letsPlot2(c2MolTimeNuc{i}{5})
    %______________________________________________________________
    if ~isempty(c2spotMolNuc{i})
        if ~isempty(vertcat(c2spotMolNuc{i}{:}))
            c2Strength{i} = vertcat(c2spotMolNuc{i}{:}); % columns: time; rows: nuc
            c1ValMean{i} = vertcat((c1ValC2RadMean{i}{:})); % columns: time; rows: nuc
            c1ValSem{i} = vertcat((c1ValC2RadSem{i}{:})); % columns: time; rows: nuc
        end
    end
end

maxLen = min(nonzeros(cellfun(@length, c1ValC2RadNormAllCell)));

j = 0;
for i = 1:length(c1ValC2RadNormAllCell)
    if ~isempty(c1ValC2RadNormAllCell{i})
        j = j+1;
        if j == 1 
            nucMeanValAll = nucMeanValCell{i}';
            c2MolTimeNucAll = c2MolTimeNucAllCell{i}';   
            
            c1ValC2RadDiffAll = c1ValC2RadDiffAllCell{i}(1:maxLen, :);
            c1ValC2RadDiffCrop3D = mean(c1ValC2RadDiffAllCell{i}(1:cropLen, :))'; % initialize the crop
            c1ValC2RadDiffAllXY = c1ValC2RadDiffAllCellXY{i}(1:maxLen, :);
            c1ValC2RadDiffCropXY = mean(c1ValC2RadDiffAllCellXY{i}(1:cropLen, :))'; % initialize the crop
            c1ValC2RadDiffAllZ = c1ValC2RadDiffAllCellZ{i}(1:maxLen, :);
            
            c1ValC2RadNormAll = c1ValC2RadNormAllCell{i}(1:maxLen, :);
            c1ValC2RadNormCrop3D = mean(c1ValC2RadNormAllCell{i}(1:cropLen, :))'; % initialize the crop
            c1ValC2RadNormAllXY = c1ValC2RadNormAllCellXY{i}(1:maxLen, :);
            c1ValC2RadNormCropXY = mean(c1ValC2RadNormAllCellXY{i}(1:cropLen, :))'; % initialize the crop
            c1ValC2RadNormAllZ = c1ValC2RadNormAllCellZ{i}(1:maxLen, :);
            c1ValC2NoBgRadNormAll = c1ValC2NoBgRadNormAllCell{i}(1:maxLen, :);
        else
            nucMeanValAll = vertcat(nucMeanValAll, nucMeanValCell{i}');
            c2MolTimeNucAll = vertcat(c2MolTimeNucAll, c2MolTimeNucAllCell{i}');     
            
            c1ValC2RadDiffAll = horzcat(c1ValC2RadDiffAll, c1ValC2RadDiffAllCell{i}(1:maxLen, :));
            c1ValC2RadDiffCrop3D = vertcat(c1ValC2RadDiffCrop3D, mean(c1ValC2RadDiffAllCell{i}(1:cropLen, :), 1, 'omitnan')');
            c1ValC2RadDiffAllXY = horzcat(c1ValC2RadDiffAllXY, c1ValC2RadDiffAllCellXY{i}(1:maxLen, :));
            c1ValC2RadDiffCropXY = vertcat(c1ValC2RadDiffCropXY, mean(c1ValC2RadDiffAllCellXY{i}(1:cropLen, :), 1, 'omitnan')');
            c1ValC2RadDiffAllZ = horzcat(c1ValC2RadDiffAllZ, c1ValC2RadDiffAllCellZ{i}(1:maxLen, :));
            
            c1ValC2RadNormAll = horzcat(c1ValC2RadNormAll, c1ValC2RadNormAllCell{i}(1:maxLen, :));
            c1ValC2RadNormCrop3D = vertcat(c1ValC2RadNormCrop3D, mean(c1ValC2RadNormAllCell{i}(1:cropLen, :), 1, 'omitnan')');
            c1ValC2RadNormAllXY = horzcat(c1ValC2RadNormAllXY, c1ValC2RadNormAllCellXY{i}(1:maxLen, :));
            c1ValC2RadNormCropXY = vertcat(c1ValC2RadNormCropXY, mean(c1ValC2RadNormAllCellXY{i}(1:cropLen, :), 1, 'omitnan')');
            c1ValC2RadNormAllZ = horzcat(c1ValC2RadNormAllZ, c1ValC2RadNormAllCellZ{i}(1:maxLen, :));
            
            c1ValC2NoBgRadNormAll = horzcat(c1ValC2NoBgRadNormAll, c1ValC2NoBgRadNormAllCell{i}(1:maxLen, :));
        end
    end
end

c1ValC2RadDiffMean = mean(c1ValC2RadDiffAll, 2, 'omitnan');
c1ValC2RadDiffSem = std(c1ValC2RadDiffAll, 0, 2, 'omitnan');
c1ValC2RadDiffSem = c1ValC2RadDiffSem./sqrt(sum(~isnan(c1ValC2RadDiffAll),2));

c1ValC2RadDiffMeanXY = mean(c1ValC2RadDiffAllXY, 2, 'omitnan');
c1ValC2RadDiffSemXY = std(c1ValC2RadDiffAllXY, 0, 2, 'omitnan');
c1ValC2RadDiffSemXY = c1ValC2RadDiffSemXY./sqrt(sum(~isnan(c1ValC2RadDiffAllXY),2));

c1ValC2RadDiffMeanZ = mean(c1ValC2RadDiffAllZ, 2, 'omitnan');
c1ValC2RadDiffSemZ = std(c1ValC2RadDiffAllXY, 0, 2, 'omitnan');
c1ValC2RadDiffSemZ = c1ValC2RadDiffSemZ./sqrt(sum(~isnan(c1ValC2RadDiffAllZ),2));


c1ValC2RadNormMean = mean(c1ValC2RadNormAll, 2, 'omitnan');
c1ValC2RadNormSem = std(c1ValC2RadNormAll, 0, 2, 'omitnan');
c1ValC2RadNormSem = c1ValC2RadNormSem./sqrt(sum(~isnan(c1ValC2RadNormAll),2));

c1ValC2RadNormMeanXY = mean(c1ValC2RadNormAllXY, 2, 'omitnan');
c1ValC2RadNormSemXY = std(c1ValC2RadNormAllXY, 0, 2, 'omitnan');
c1ValC2RadNormSemXY = c1ValC2RadNormSemXY./sqrt(sum(~isnan(c1ValC2RadNormAllXY),2));

c1ValC2RadNormMeanZ = mean(c1ValC2RadNormAllZ, 2, 'omitnan');
c1ValC2RadNormSemZ = std(c1ValC2RadNormAllXY, 0, 2, 'omitnan');
c1ValC2RadNormSemZ = c1ValC2RadNormSemZ./sqrt(sum(~isnan(c1ValC2RadNormAllZ),2));

c1ValC2NoBgRadNormMean = mean(c1ValC2NoBgRadNormAll, 2, 'omitnan');
c1ValC2NoBgRadNormSem = std(c1ValC2NoBgRadNormAll, 0, 2, 'omitnan');
c1ValC2NoBgRadNormSem = c1ValC2NoBgRadNormSem./sqrt(sum(~isnan(c1ValC2NoBgRadNormAll),2));

radAdd = (0.1*(1:size(c1ValC2RadNormMean, 1)) - 0.1)'; % hard coded

combine2C.gene = geneName;
combine2C.nucMean = nucMeanValAll;
combine2C.rad = radAdd;
combine2C.c2Mol = c2MolTimeNucAll;
combine2C.valCrop3D = c1ValC2RadDiffCrop3D;
combine2C.valCropXY = c1ValC2RadDiffCropXY;
combine2C.valAllNorm3D = c1ValC2RadNormAll;
combine2C.valAllNormXY = c1ValC2RadNormAllXY;
combine2C.valAllDiff3D = c1ValC2RadDiffAll;
combine2C.valAllDiffXY = c1ValC2RadDiffAllXY;
combine2C.val3D = c1ValC2RadNormMean;
combine2C.sem3D = c1ValC2RadNormSem;
combine2C.valXY = c1ValC2RadNormMeanXY;
combine2C.semXY = c1ValC2RadNormSemXY;
combine2C.valZ = c1ValC2RadNormMeanZ;
combine2C.semZ = c1ValC2RadNormSemZ;

letsPlot3D(c1ValC2RadNormXY{1}{7});

fileName = append('combine','_',geneName,'_','DS');

combineFolder = 'D:\Tyrone_analysis\test\val_combine';
save([combineFolder, filesep, fileName, '.mat'], 'combine2C');
figure('color', 'w');
% plotErr(c1ValC2RadNormMean, c1ValC2RadNormSem, radAdd, geneName);
plotErr(c1ValC2RadNormMeanXY, c1ValC2RadNormSemXY, radAdd, geneName, colorStruct{1});
hold off;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% nucValBinEdge =[150, 250, 350, 450]; % hard coded
% nucValBin = discretize(nucMeanValAll, nucValBinEdge);
% dataIdx = nucValBin;
% % dataIdx(isnan(nucValBin)) = [];
% nucMeanVals = nucMeanValAll;
% % nucMeanVals(isnan(nucValBin)) = [];
% % binCounts(:,1) = hist(dataIdx,unique(dataIdx));
% [binCounts, ~] = hist(dataIdx(~isnan(dataIdx)),unique(dataIdx(~isnan(dataIdx))));
% binCounts = binCounts';
% nucMeanValBin = accumarray(dataIdx(~isnan(dataIdx)), nucMeanVals(~isnan(dataIdx)),[],@nanmean);
% nucStdValBin = accumarray(dataIdx(~isnan(dataIdx)), nucMeanVals(~isnan(dataIdx)),[],@nanstd);
% nucSemValBin = nucStdValBin./sqrt(binCounts);
% fullXYDiff = c1ValC2RadDiffAllXY';
% % fullXYDiff(isnan(nucValBin), :) = [];
% fullXYNorm = c1ValC2RadNormAllXY';
% % fullXYNorm(isnan(nucValBin), :) = [];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[dataIdx, nucMeanValBin] =kmeans(nucMeanValAll, totalKMeans);
D = sqrt(sum((nucMeanValAll - nucMeanValBin(dataIdx,:)).^2,2));
k = D <=kMeansDistLimit;
dataIdx(~k) = NaN;

binCounts(:,1) = hist(dataIdx(~isnan(dataIdx)),unique(dataIdx(~isnan(dataIdx))));
fullXYDiff = c1ValC2RadDiffAllXY';
fullXYNorm = c1ValC2RadNormAllXY';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i = 1:size(fullXYDiff,2)
    fullXYDiffMeanBin(:, i) = accumarray(dataIdx(~isnan(dataIdx)), fullXYDiff(~isnan(dataIdx), i),[],@nanmean);
    fullXYDiffStdBin(:, i) = accumarray(dataIdx(~isnan(dataIdx)), fullXYDiff((~isnan(dataIdx)), i),[],@nanstd);
    fullXYDiffSemBin(:, i) = fullXYDiffStdBin(:, i)./sqrt(binCounts);
    
    fullXYNormMeanBin(:, i) = accumarray(dataIdx(~isnan(dataIdx)), fullXYNorm((~isnan(dataIdx)), i),[],@nanmean);
    fullXYNormStdBin(:, i) = accumarray(dataIdx(~isnan(dataIdx)), fullXYNorm((~isnan(dataIdx)), i),[],@nanstd);
    fullXYNormSemBin(:, i) = fullXYNormStdBin(:, i)./sqrt(binCounts);
end

figure('color', 'w');
for i =  1:size(fullXYDiffMeanBin, 1)   
    plotHandle(i) = plotErr(fullXYDiffMeanBin(i, :)', fullXYDiffSemBin(i, :)', radAdd, num2str(binCounts(i)), colorStruct{i});
    hold on;
end

leg = legend(plotHandle, num2str(round(nucMeanValBin)));
set(leg,'color','none', 'TextColor', [0.5 0.5 0.5], 'FontSize', 12);
% legend(plotHandle, num2str(nucMeanValBin));
ylim([0, inf]);
title(geneName);
hold off;


figure('color', 'w');
hold on;
for i =  1:size(fullXYNormMeanBin, 1)   
    plotHandle(i) = plotErr(fullXYNormMeanBin(i, :)', fullXYNormSemBin(i, :)', radAdd, num2str(binCounts(i)), colorStruct{i});
    hold on;
end

leg = legend(plotHandle, num2str(round(nucMeanValBin)));
set(leg,'color','none', 'TextColor', [0.5 0.5 0.5], 'FontSize', 10);
% legend(plotHandle, num2str(nucMeanValBin));
ylim([0, inf]);
title(geneName);
hold off;


cropXYNorm = c1ValC2RadNormCropXY;
% cropXYNorm(isnan(nucValBin)) = [];
cropXYNormMeanBin = accumarray(dataIdx(~isnan(dataIdx)), cropXYNorm(~isnan(dataIdx)),[],@nanmean);
cropXYNormStdBin = accumarray(dataIdx(~isnan(dataIdx)), cropXYNorm(~isnan(dataIdx)),[],@nanstd);
cropXYNormSemBin = cropXYNormStdBin./sqrt(binCounts);

cropXYDiff = c1ValC2RadDiffCropXY;
% cropXYDiff(isnan(nucValBin)) = [];
cropXYDiffMeanBin = accumarray(dataIdx(~isnan(dataIdx)), cropXYDiff(~isnan(dataIdx)),[],@nanmean);
cropXYDiffStdBin = accumarray(dataIdx(~isnan(dataIdx)), cropXYDiff(~isnan(dataIdx)),[],@nanstd);
cropXYDiffSemBin = cropXYDiffStdBin./sqrt(binCounts);

crop3DDiff = c1ValC2RadDiffCrop3D;
% cropXYDiff(isnan(nucValBin)) = [];
crop3DDiffMeanBin = accumarray(dataIdx(~isnan(dataIdx)), crop3DDiff(~isnan(dataIdx)),[],@nanmean);
crop3DDiffStdBin = accumarray(dataIdx(~isnan(dataIdx)), crop3DDiff(~isnan(dataIdx)),[],@nanstd);
crop3DDiffSemBin = crop3DDiffStdBin./sqrt(binCounts);


[~, ~, vv] = unique(dataIdx(~isnan(dataIdx)));
totalSamples = accumarray(vv,1);
% errorbar(nucMeanValBin, cropXYNormMeanBin, cropXYNormSemBin);

plotLine(nucMeanValBin, crop3DDiffMeanBin, crop3DDiffSemBin);
plotLine(nucMeanValBin, cropXYDiffMeanBin, cropXYDiffSemBin);

groupColor = cellfun(@(x) vertcat(x(1,:)), colorStruct, 'un', 0);
groupColor = vertcat(groupColor{:});
plotBox(nucMeanValBin, cropXYDiff, dataIdx, cropXYDiffMeanBin, cropXYDiffSemBin, binCounts, groupColor, 'XY');
plotBox(nucMeanValBin, crop3DDiff, dataIdx, crop3DDiffMeanBin, crop3DDiffSemBin, binCounts, groupColor, '3D');
% plotBox(nucMeanValBin, cropXYNorm, dataIdx, cropXYNormMeanBin, cropXYNormSemBin, binCounts, groupColor, 'XY');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(nucMeanValBin)
    valRadDiffMeanBin(:,i) = mean(c1ValC2RadDiffAll(:, dataIdx==i), 2, 'omitnan');
    valRadDiffStdBin(:,i) = std(c1ValC2RadDiffAll(:, dataIdx==i), 0, 2, 'omitnan');
    valRadDiffSemBin(:,i) = valRadDiffStdBin(:,i)./sqrt(binCounts(i));
    
    valRadNormMeanBin(:,i) = mean(c1ValC2RadNormAll(:, dataIdx==i), 2, 'omitnan');
    valRadNormStdBin(:,i) = std(c1ValC2RadNormAll(:, dataIdx==i), 0, 2, 'omitnan');
    valRadNormSemBin(:,i) = valRadNormStdBin(:,i)./sqrt(binCounts(i));
    
    
    valRadXYDiffMeanBin(:,i) = mean(c1ValC2RadDiffAllXY(:, dataIdx==i), 2, 'omitnan');
    valRadXYDiffStdBin(:,i) = std(c1ValC2RadDiffAllXY(:, dataIdx==i), 0, 2, 'omitnan');
    valRadXYDiffSemBin(:,i) = valRadXYDiffStdBin(:,i)./sqrt(binCounts(i));
    
    valRadXYNormMeanBin(:,i) = mean(c1ValC2RadNormAllXY(:, dataIdx==i), 2, 'omitnan');
    valRadXYNormStdBin(:,i) = std(c1ValC2RadNormAllXY(:, dataIdx==i), 0, 2, 'omitnan');
    valRadXYNormSemBin(:,i) = valRadXYNormStdBin(:,i)./sqrt(binCounts(i));
end

errorbar(repmat(radAdd, 1, length(nucMeanValBin)), valRadDiffMeanBin, valRadDiffSemBin);
errorbar(repmat(radAdd, 1, length(nucMeanValBin)), valRadNormMeanBin, valRadNormSemBin);
errorbar(repmat(radAdd, 1, length(nucMeanValBin)), valRadXYDiffMeanBin, valRadXYDiffSemBin);
errorbar(repmat(radAdd, 1, length(nucMeanValBin)), valRadXYNormMeanBin, valRadXYNormSemBin);
% xlim([0, 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c1ValC2NoBgAll = horzcat(c1ValC2timeMeanNoBg{:});
% c1ValC2NoBgAll = horzcat(c1ValC2NoBgAll{:});
% c1ValC2NoBgAll(isnan(c1ValC2NoBgAll)) = 0;
% c1ValC2NoBgMeanAll = mean(c1ValC2NoBgAll, 2);
% c1ValC2NoBgSemAll = std(c1ValC2NoBgAll, 0, 2)./sqrt(size(c1ValC2NoBgAll, 2));

% c1ValC2All = horzcat(c1ValC2timeMean{:});
% c1ValC2All = horzcat(c1ValC2All{:});
% c1ValC2All(isnan(c1ValC2All)) = 0;
% c1ValC2MeanAll = mean(c1ValC2All, 2);
% c1ValC2SemAll = std(c1ValC2All, 0, 2)./sqrt(size(c1ValC2All, 2));

% xArr = 93*(1:length(c2Strength{1}(5,:))); % hard coded 62 seconds is assumed to be the frame time
% plotDouble(xArr, c2Strength{1}(8,:), c1ValMean{1}(5,:), c1ValSem{1}(5,:));
% 
% radAdd = (0.1*(1:size(c1ValC2NoBgAll, 1)) - 0.1)'; % hard coded
% % plotErr(c1ValC2NoBgMeanAll, c1ValC2NoBgSemAll, radAdd);
% plotErr(c1ValC2MeanAll, c1ValC2SemAll, radAdd);

end

function plotHandle = plotErr(val, err, xArr, geneName, color)
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
% figure('color', 'w');
patchTop = val+err;
patchTop = reshape(patchTop,1,[]);
patchBot = val-err;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch = [xArr',fliplr(xArr')];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = color(1,:)./255;
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.6;
hold on;
pl = plot(xArr, val);
pl.LineWidth = 1.5;
pl.Color = color(1,:)./255;
pl.LineStyle = '-';
%---------------------------------------------

% hold on;
% pd = fitdist(val,'HalfNormal');
% yFit = pdf(pd, xArr);
% plot(xArr, yFit);
%---------------------------------------------
% hLeg = legend(geneName);
ylabel('Relative bcd intensity (a. u.)');
xlabel('Apparent radius (r)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
ylim([0, 0.4]); xlim([0,1.5]);

set(gcf, 'color', 'none');  
set(gca, 'color', 'none');
set(gca,'ycolor',[0.5 0.5 0.5])
set(gca,'xcolor',[0.5 0.5 0.5])

plotHandle = pl;
end

function plotErr2(val, err, xArr, geneName)
figure('color', 'w');
errorbar(xArr, val, err, 'o','MarkerSize',3,...
    'MarkerEdgeColor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3], 'LineStyle', 'none', 'LineWidth', 1.5);
hold on;
% pd = fitdist(val,'hn');
% pdFun = pdf(pd, xArr);
% pFit = pdFun/sum(0.1*pdFun);
nanIdx = find(isnan(val));
ff = fit(xArr, val, 'poly7', 'Exclude', nanIdx);
pf = plot(ff);
pf(1).Color = [0.3, 0.3, 0.3];
pf.LineStyle = '--';
pf.LineWidth = 2;
% pl = plot(xArr, val);
% pl.LineWidth = 1.5;
% pl.Color = [0.4 0.4 0.4];
% pl.LineStyle = '--';
%---------------------------------------------

% hold on;
% pd = fitdist(val,'HalfNormal');
% yFit = pdf(pd, xArr);
% plot(xArr, yFit);
%---------------------------------------------
hLeg = legend(geneName);
% set(hLeg,'visible','off')
ylabel('Relative TF intensity (a. u.)');
xlabel('Distance from the mRNA center ({\mu}m) ');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=350;
plotHeight=350;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
hold off;
end

function plotDouble(xArr, y1, y2Mean, y2Dev)

leftColor = [0.9 0.3 0.3];
rightColor = [0.3 0.7 0.5];

yyaxis left;
pp = plot(xArr, y1);
pp.LineWidth = 1.5;
pp.Color = leftColor;
pp.LineStyle = '-';

ylabel('MS2 Intensity');
set(gca,'YTick',[]);
set(gca,'ycolor',leftColor) ;
hold on;
yyaxis right;
patchTop = y2Mean+y2Dev;
patchTop = reshape(patchTop,1,[]);
patchBot = y2Mean-y2Dev;
patchBot = reshape(patchBot, 1, []);
yPatch=[patchBot,fliplr(patchTop)];
xPatch = [xArr,fliplr(xArr)];
pt = patch(xPatch, yPatch, 1);
pt.FaceColor = rightColor;
pt.EdgeColor = 'none';
pt.FaceAlpha = 0.4;
hold on;
pl = plot(xArr, y2Mean);
pl.LineWidth = 1.5;
pl.Color = rightColor;
pl.LineStyle = '-';

ylabel('Bicoid Intensity');
set(gca,'YTick',[]);
set(gca,'ycolor',rightColor) ;
hold on;
xlabel('Time (s)');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=500;
plotHeight=200;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf,'defaultAxesColorOrder',[leftColor; rightColor]);
hold off;
end

function plotBox(xNames, dataCombine, g, dataMean, dataSem, groupLength, colorPalette, titleText)
colorPalette = colorPalette./255;
xNameInt = uint16(xNames);
xNameChar = num2str(xNameInt);
xNameChar = cellstr(xNameChar);
xNameStr = convertCharsToStrings(xNameChar);
% xNames = ["aa", 'bb', 'cc', 'dd', 'ee', 'ff'];%num2str(xNames);
xNamesCat = categorical(g, 1:length(xNames), xNameStr);
% b = boxchart(dataCombine, 'GroupByColor', g);
b = boxchart(xNamesCat, dataCombine);
boxIndex = b.SeriesIndex;
for j = 1:length(b)
    b(j).BoxWidth = 0.3;
    b(j).Notch = 'off';
    b(j).BoxFaceColor = colorPalette(j,:);
    b(j).WhiskerLineColor = colorPalette(j,:);
    b(j).LineWidth = 2;
    b(j).MarkerStyle = 'none';
    b(j).MarkerColor = colorPalette(j,:);    
end

hold on;

errorbar(dataMean, dataSem, '.','MarkerSize',10,...
    'MarkerEdgeColor','k', 'Color', 'k', 'LineStyle', 'none');
hold on;

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
box(ax,'on');
grid off;
x0 = 100;
y0= 100;
plotWidth=300;
plotHeight=300;
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
end

function letsPlot1(data)
timePoints = size(data, 2);
% dataXTemp1 = 0.1.*(1:size(data, 1)) - 0.1;
% dataXTemp2 = -1.*fliplr(0.1.*(1:size(data, 1)));
% dataX = [dataXTemp2, dataXTemp1]';

dataXTemp1 = 0.1.*(1:size(data, 1)) - 0.1;
dataX = dataXTemp1;

yLims = zeros(timePoints, 2);

figure;
hold on;
for i = 1:timePoints
%     dataYTemp1 = data(:,timePoints-i+1);
%     dataYTemp2 = flipud(dataYTemp1);
%     dataY = [dataYTemp2;  dataYTemp1];

    dataY = data(:,end-i+1);
    subplot(timePoints,1, i);    
    p(i) = plot(dataX, dataY, '-g', 'LineWidth',2);
    ax(i) = gca;
    yLims(i,:) = ax(i).YLim;
    ax(i).FontSize = 8;
    ax(i).LineWidth = 2;
    box(ax(i),'off');    
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    set(gcf, 'color', 'none');   
    set(gca, 'color', 'none');
    set(gca,'ycolor',[0.1 0.1 0.1])
    set(gca,'xcolor',[0.1 0.1 0.1])
%     ylim([0 1]);
    if i==timePoints
        xticks('auto');
        xlabel('Distance ({\mu}m)')
    end
    ylabel([num2str((timePoints-i+1)*93), 's']);
end
yLimAll = [min(yLims(:,1)), max(yLims(:,2))];
for i = 1:timePoints
    ax(i).YLim = yLimAll;
end

grid off;
x0 = 100;
y0= 100;
plotWidth = 200;
plotHeight = 500;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])

end

function letsPlot2(data)
timePoints = size(data, 2);
figure;
hold on;
for i = 1:timePoints
    dataY = data(timePoints-i+1);
    subplot(timePoints,1, i);    
    bar(dataY, 'r');
    
    ax = gca;
    set(gcf, 'color', 'none');   
    set(gca, 'color', 'none');
    set(gca,'ycolor',[0.5 0.5 0.5])
    set(gca,'xcolor',[0.5 0.5 0.5])
    
    ax.FontSize = 12;
    ax.LineWidth = 2;
    box(ax,'off');
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    
    ylim([0 max(data)]);
    if i==timePoints
        xlabel('MS2 hotspot signal')
    end
    ylabel([num2str((timePoints-i+1)*90), 's']);
end

grid off;
x0 = 100;
y0= 100;
plotWidth = 200;
plotHeight = 500;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end

function letsPlot3D(data1, data2)
timePoints = size(data1, 2);
figure;
hold on;

x = (0:0.1:1.9).';
xMat = repmat(x, 1, timePoints); %// For plot3
y = 0:93:(timePoints-1)*93;
yMat = repmat(y, numel(x), 1); %//For plot3

% x = 0:93:(timePoints-1)*93;
% xMat = repmat(x, numel(x), 1); %//For plot3
% y = (0:0.1:1.9).';
% yMat = repmat(y, 1, timePoints); %// For plot3


% plot3(xMat, yMat, data, 'b'); %// Make all traces blue
plot3(yMat, xMat, data1, 'b'); %// Make all traces blue
grid;
xlabel('x'); ylabel('y'); zlabel('z');
view(30, 60); %// Adjust viewing angle so you can clearly see data
ZL = zlim(gca);
DZ = 0.1*(ZL(2)-ZL(1));
for k=1:size(xMat,2)
%     hPatch(k) = patch( ...
%         [xMat(:,k);    flipud(xMat(:,k))   ], ...
%         [yMat(:,k);    flipud(yMat(:,k))   ], ...
%         [data(:,k);    flipud(data(:,k))-DZ], ...
%         'w');
        hPatch(k) = patch( ...
        [yMat(:,k);    flipud(yMat(:,k))   ], ...
        [xMat(:,k);    flipud(xMat(:,k))   ], ...        
        [data1(:,k);    flipud(data1(:,k))-DZ], ...
        'w');
    
    set(hPatch(k), 'EdgeColor','none', 'FaceColor','w', 'FaceAlpha',0.9 );
end
hold on;


axis tight

zlim([min(data1, [], 'all') max(data1, [], 'all')]);
ylabel('Distance from MS2 hotspot ({\mu}m)');
xlabel('Time (s)');
zlabel('TF intensity (a. u.)');


ylh = get(gca,'ylabel');
gyl = get(ylh);   % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',30, 'Position',ylp, 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
% grid off;
grid on;
% x0 = 100;
% y0= 100;
% plotWidth = 200;
% plotHeight = 500;
% set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end

function plotLine(xVal, yVal, yErr)

  arr = [161.9795   26.2661    4.3253
  323.7421   40.3925   11.0859
  255.1990   28.4189    7.7454
  413.3192   53.6688   20.4751
   161.9795   26.2661    4.3253
  323.7421   40.3925   11.0859
  255.1990   28.4189    7.7454
  413.3192   53.6688   20.4751];

errorbar(arr(:,1), arr(:,2), arr(:,3), '-g', 'LineWidth',2);

errorbar(xVal, yVal, yErr, '-g', 'LineWidth',2);
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
box(ax,'off');

set(gca,'ytick',[])
set(gca,'xtick',[])
set(gcf, 'color', 'none');   
set(gca, 'color', 'none');
set(gca,'ycolor',[0.5 0.5 0.5])
set(gca,'xcolor',[0.5 0.5 0.5])
end

function btstrpMean = btstrpFun(x)
btstrpMean = bootstrp(500, @nanmean, x);
end

