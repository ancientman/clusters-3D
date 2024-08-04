function C1SpotPropCZIplot1(DSfolder)
%______________________________________________________________________________________
%	Works for single color multiposition data generated from CZI file format.
%   Calculates various properties based on position bins and saves the
%   datastructures within a specified folder.
%	Plots total spots per nucleus etc. versus EL position
%   Input#1:  DS folder
%   Input#2: (string) data identifier e.g. 'bcd 2xa'
%______________________________________________________________________________________


%   Hardcoded parameters.
%______________________________________________________________________________________
positionBin = 0.00:0.05:1; % hard coded position bins. 
posDSName = 'spotPosPropDS'; % save params in this DS
fitDSName = 'spotFitPropDS'; % save spot fit params in this DS
fitFlag = 'on'; %  options 'on' and 'off' . Off turns off all fitting.

spotXYpixelLim = 5; % The number of pixels to consider around the intensity centroid
spotZpixelLim = 2; % The number of pixels to consider around the intensity centroid
spotWindow = 5; % 2D window size for intensity analysis around the spot centroid
%______________________________________________________________________________________

if exist(append(DSfolder, filesep, 'c1SpotPropDS.mat'), 'file')
    spotDS = load(append(DSfolder, filesep, 'c1SpotPropDS.mat')); 
else
      error('no relevant spot prop DS found')
end

if exist(append(DSfolder, filesep, 'c1NucPropDS.mat'), 'file')
    nucDS = load(append(DSfolder, filesep, 'c1NucPropDS.mat')); 
else
      error('no relevant nuc prop DS found')
end

if exist(append(DSfolder, filesep, 'positionListDS.mat'), 'file')
     posDS = load(append(DSfolder, filesep, 'positionListDS.mat')); 
else
      error('no relevant position prop DS found')
end

if exist(append(DSfolder, filesep, 'embryoDimDS.mat'), 'file')
     emDimDS = load(append(DSfolder, filesep, 'embryoDimDS.mat')); 
else
      error('no relevant embryo DS found')
end

if exist(append(DSfolder, filesep, 'seriesMetaDataDS.mat'), 'file')
     metaDataDS = load(append(DSfolder, filesep, 'seriesMetaDataDS.mat')); 
else
      error('no relevant metaData DS found')
end

if metaDataDS.seriesMetaDataDS.analysisInfo.xPixUM == ...
        metaDataDS.seriesMetaDataDS.analysisInfo.yPixUM
    xyPixSize = metaDataDS.seriesMetaDataDS.analysisInfo.xPixUM;
    zPixSize = metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
    voxSize = (metaDataDS.seriesMetaDataDS.analysisInfo.xPixUM)^2 ...
    *metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
else
    error('x pixel size is not equal to y pixel size');
end

sceneCenterXpix = round(metaDataDS.seriesMetaDataDS.imagingInfo.stackSizeX/2);
sceneCenterYpix = round(metaDataDS.seriesMetaDataDS.imagingInfo.stackSizeY/2);

p = 1; % counter for all nuclei deteted in the czi file
spotPropPerNuc = {}; % combines all nuclei from every positions
fitDS = {}; % collects all the spot fit paramters

if length(nucDS.c1NucProp) == length(posDS.positionList)
%     for i = 1:6 % short run over n scenes only
    for i = 1:length(nucDS.c1NucProp) % full run, over all scenes
        sceneCenterXum = posDS.positionList{i}.xPosition;
        sceneCenterYum = posDS.positionList{i}.yPosition;
        k = 1; % counter for total nuclei in one CZI scene
        while(k <= length(nucDS.c1NucProp{i}{1}))
            %These two would accumulate values from the fit function
            temp2DfitTotVal = 0;
            temp2DfitDiaUM = 0;           
            temp3DfitTotVal = 0;
            temp3DfitDiaUM = 0;
            temp2DfitDeconvDiaUM = 0;
            temp3DfitDeconvDiaUM = 0;
            %______________________________________________________________________________________________
            %   Find absolute EL position of each nuclei
            nucDS.c1NucProp{i}{1}(k).nucCentUM = [sceneCenterXum, sceneCenterYum] + ...
                 (nucDS.c1NucProp{i}{1}(k).center(1:2) - ...
                 [sceneCenterXpix, sceneCenterYpix])*xyPixSize;
             
             nucDS.c1NucProp{i}{1}(k).nucAntDist = pdist([(emDimDS.emDim.ant); ...
                 (nucDS.c1NucProp{i}{1}(k).nucCentUM)]);
             
            xDiffUM = (nucDS.c1NucProp{i}{1}(k).nucCentUM(1))-(emDimDS.emDim.ant(1));
            yDiffUM = (nucDS.c1NucProp{i}{1}(k).nucCentUM(2))-(emDimDS.emDim.ant(2));
            nucDS.c1NucProp{i}{1}(k).nucAxisAngle = atan(yDiffUM/xDiffUM)...
                - emDimDS.emDim.angle;

            nucDS.c1NucProp{i}{1}(k).nucDistAP = nucDS.c1NucProp{i}{1}(k).nucAntDist*...
                 cos(nucDS.c1NucProp{i}{1}(k).nucAxisAngle);        

             nucDS.c1NucProp{i}{1}(k).nucFracEL =  nucDS.c1NucProp{i}{1}(k).nucDistAP/...
                 emDimDS.emDim.len;
   
%______________________________________________________________________________________________
            spotMeanVal = cellfun(@mean, spotDS.c1SpotProp{i}{1}(k).voxVal);
 
%______________________________________________________________________________________________
%
%                           Properties of each spot within the nuclei      
%______________________________________________________________________________________________
            j = 1; % j = total spots in one nucleus
            while(j <= size(spotDS.c1SpotProp{i}{1}(k).paLen, 1))
                
                majALen = spotDS.c1SpotProp{i}{1}(k).paLen(j, 1); % Major axis length
                midALen = spotDS.c1SpotProp{i}{1}(k).paLen(j, 2); % Middle axis length
                minALen = spotDS.c1SpotProp{i}{1}(k).paLen(j, 3); % Minor axis length
                
                %   Calculate spot properties directly from regionprops
                %   parameters
                %__________________________________________________________________________________
                %   Calculate effective spot dia, from PAs
                %-------------------------------------------------

%                 spotDS.c1SpotProp{i}{1}(k).spotDiaPA(j) = xyPixSize*sqrt(2).* ...
%                     (det(diag([(midALen^2 + minALen^2), ...
%                     (majALen^2 + minALen^2), ...
%                     (majALen^2 + midALen^2)]))^0.5);

                %   Dia (based on PA) = sqrt(2)*(maj^2 + min^2)^0.5
                %------------------------------------------------
                spotDS.c1SpotProp{i}{1}(k).spotDiaPA(j) = xyPixSize*sqrt(2).* ...
                    ((majALen^2 + minALen^2)^0.5);

                %   Dia (based on spot volume) = (6/pi)*(volume)^1/3
                %------------------------------------------------
                volDia = ((6/pi)*spotDS.c1SpotProp{i}{1}(k).volUM(j))^(1/3);
                spotDS.c1SpotProp{i}{1}(k).spotDiaVol(j) = volDia;

                %   Spot intensity mean
                %------------------------------------------------
                spotDS.c1SpotProp{i}{1}(k).spotMeanVal(j) = spotMeanVal(j);               
                
                %   Condition to select only spots with independent critera
                %   on the 2D area and the z span of the bounding box 
                %------------------------------------------------
                if spotDS.c1SpotProp{i}{1}(k).bb(j, 4)*spotDS.c1SpotProp{i}{1}(k).bb(j, 5) < spotXYpixelLim || ...
                        spotDS.c1SpotProp{i}{1}(k).bb(j, 6)< spotZpixelLim % only spots above the limit    
                    
                    spotDS.c1SpotProp{i}{1}(k).spotDiaVol(j) = 0;
                    spotDS.c1SpotProp{i}{1}(k).spotDiaPA(j) = 0;
                    spotDS.c1SpotProp{i}{1}(k).spotMeanVal(j) = 0;                          
                end        
                %____________________________________________________________________________________

                
                % Calculate sigma and amp of 2D gauss fit of spots
                %____________________________________________________________________________________                
                if spotDS.c1SpotProp{i}{1}(k).bb(j, 4)*spotDS.c1SpotProp{i}{1}(k).bb(j, 5) < spotXYpixelLim || ...
                        spotDS.c1SpotProp{i}{1}(k).bb(j, 6)< spotZpixelLim || ...
                        strcmp(fitFlag, 'off')% only spots above the limit    
                    
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j} = [];
                else                    
                %____________________________________________________________________________________ 
                
                    %   Send for 3D deconvolution all raw pixels within a
                    %   window around the max val center
                    %____________________________________________________________________________________
%                     spotDS.c1SpotProp{i}{1}(k).deconvVal{j} = spotBBdeconv(...
%                         spotDS.c1SpotProp{i}{1}(k).voxValCenter(j, :), ...
%                         nucDS.c1NucProp{i}{1}(k).voxList{1}, nucDS.c1NucProp{i}{1}(k).voxVal{1}, ...
%                         spotWindow, metaDataDS);
                    %____________________________________________________________________________________
                    
                    %	Starting from the value weighted cenrtroid, it considers xy pixels 
                    %   within the specified window and 2D gaussian fits the intensity profile
                    %   Then calculates fit properties including x and y sigma
                    %----------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j} = spotFitXY(...
                        spotDS.c1SpotProp{i}{1}(k).voxValCenter(j,:), ...
                        nucDS.c1NucProp{i}{1}(k).voxList, ...
                        nucDS.c1NucProp{i}{1}(k).voxVal, spotWindow);
                    %----------------------------------------------------------------
                    
                    %	Starting from the value weighted cenrtroid, it considers xz pixels 
                    %   within the specified window and 2D gaussian fits the intensity profile
                    %   Then calculates fit properties including x and z
                    %   sigma. This is mainly for the sigma in the z
                    %   direction.
                    %----------------------------------------------------------------                    
%                     spotDS.c1SpotProp{i}{1}(k).spotFitProp{j} = spotFitXZ(...
%                         spotDS.c1SpotProp{i}{1}(k).voxValCenter(j,:), ...
%                         nucDS.c1NucProp{i}{1}(k).voxList, ...
%                         nucDS.c1NucProp{i}{1}(k).voxVal, spotWindow);
                    %----------------------------------------------------------------
                    
                    
                    %	Does 1D gaussian fitting along z to calculate z sigma
                    %----------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.zSigma = ...
                        spotFitZ(spotDS.c1SpotProp{i}{1}(k).voxValCenter(j,:), ...
                        nucDS.c1NucProp{i}{1}(k).voxList, ...
                        nucDS.c1NucProp{i}{1}(k).voxVal, spotWindow, metaDataDS);
                    %----------------------------------------------------------------
                    

                    %______________________________________________________________________________________
                    %   Calculate radii from fits (contrary to the
                    %   regionprop values.
                    %______________________________________________________________________________________
                    
                    %   2D equivalent radius = x_sigma*y_sigma/(x_sigma+y_sigma)  
                    %----------------------------------------------------------------
                    temp2Drad = (spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.xSigma*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.ySigma)/...
                        (spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.xSigma + ...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.ySigma);                    
                    
                    %   Z equivalent radius = z_sigma 
                    %-----------------------------------------------------------------
                    tempZrad = spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.zSigma;
                    
                    %   2D dia accumulate um (equivalent dia = 2*pix*eq_rad) 
                    %-----------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DfitDiaUM = ...
                        2*(xyPixSize)*(   temp2Drad   );
                    
                    %	3D equivalent radius: 
                    %   NOTE: the pix size along x and Z are different.
                    %   Hence this way of dinfing equivalent 3D
                    %-----------------------------------------------------------------
                    temp3DradUM = sqrt(((xyPixSize*temp2Drad)^2 + (zPixSize*tempZrad)^2)/2);
                    
                    %   3D dia accumulate um (equivalent dia = 2*pix*sqrt((2Drad^2 + z_sigma^2)/2))           
                    %-----------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DfitDiaUM = ...
                        2*(	temp3DradUM  );
                    %________________________________________________________________________________________
                    
                    
                    %________________________________________________________________________________________
                    %   Send for 3D deconvolution all raw pixels within a
                    %   window around the max val center
                    %________________________________________________________________________________________
                     spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.deconvPostFitProp =...
                         spotPostFitDeconv(spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}, ...
                        spotWindow, metaDataDS);                    
                    %_________________________________________________________________________________________
                    %   All deconvolved stuff. here
                    
                    
                    %   2D equivalent deconv radius = x_de_sigma*y_de_sigma/(x_de_sigma+y_de_sigma)  
                    %-----------------------------------------------------------------
                    tempXYDeconvRad = ...
                        (spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.deconvPostFitProp.xSigmaDeconv*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.deconvPostFitProp.ySigmaDeconv)/...
                        (spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.deconvPostFitProp.xSigmaDeconv + ...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.deconvPostFitProp.ySigmaDeconv);
                    
                    
                    %   Z equivalent deconv radius = z_sigma - zPsfSigma  
                    %-----------------------------------------------------------------
                    zPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
                    zPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPsfUM;
                    zPsfSigma = (zPsfUM/zPixUM)/(2*sqrt(2*log(2)));
                    tempZDeconvRad = spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.zSigma - zPsfSigma;
                    
                    %   2D Deconv dia accumulate um (equivalent dia = 2*pix*eq_rad) 
                    %------------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DfitDeconvDiaUM = ...
                        2*(xyPixSize)*(   tempXYDeconvRad   );
                    
                    %	3D deConv equivalent radius = ((sigma_2d^2 +
                    %	sigma_z^2)/2)^1/2
                    %------------------------------------------------------------------
                    temp3DdeconvRadUM = sqrt(((xyPixSize*tempXYDeconvRad)^2 +...
                        (zPixSize*tempZDeconvRad)^2)/2);                    
                    
                    %   3D dia accumulate um (equivalent dia = 2*pix*sqrt((2Drad^2 + z_sigma^2)/2))    
                    %------------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DfitDeconvDiaUM = ...
                        2*(	 temp3DdeconvRadUM  );
                    %______________________________________________________________________________________________
                    
                   
                    %______________________________________________________________________________________________
                    %   Intensities based on fits
                    %______________________________________________________________________________________________
                    %   2D total intensity accumulate =
                    %   2*pi*x_sigma*y_sigma*I_0
                    %-------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DTotVal = ...
                        2*pi*spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spotAmp*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.xSigma*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.ySigma;
                    
                    %   3D total intensity accumulate =
                    %   intensity within a sphere of radius 3D equivalent
                    %   2*pi*I_0_x_sigma*y_sigma*z_sigma
                    %--------------------------------------------------------------
                    spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DTotVal = ...
                        (2*pi)*spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spotAmp*...
                         spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.xSigma*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.ySigma*...
                        spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.zSigma;
                    
                    %   Intensity based on volume of spehere = 4pi/3*r_3D^3
                    %---------------------------------------------------------------
%                     spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DTotVal = ...
%                         (4*pi/3)*spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spotAmp*...
%                         temp3DradUM^3;
                    %_____________________________________________________________________________________________

                    
                    if (~isempty( spotDS.c1SpotProp{i}{1}(k).spotFitProp{j})) % redundant but there
                        temp2DfitDiaUM = vertcat(temp2DfitDiaUM, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DfitDiaUM);
                        temp2DfitTotVal = vertcat(temp2DfitTotVal, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DTotVal);            
                        
                        temp3DfitDiaUM = vertcat(temp3DfitDiaUM, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DfitDiaUM);
                        temp3DfitTotVal = vertcat(temp3DfitTotVal, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DTotVal);  
                        
                        temp2DfitDeconvDiaUM = vertcat(temp2DfitDeconvDiaUM, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot2DfitDeconvDiaUM);
                        temp3DfitDeconvDiaUM = vertcat(temp3DfitDeconvDiaUM, spotDS.c1SpotProp{i}{1}(k).spotFitProp{j}.spot3DfitDeconvDiaUM);  
                    end
                end               
                
                %____________________________________________________________________________________
                fitDS{i}.nuc{k}.spotFitProp{j} = spotDS.c1SpotProp{i}{1}(k).spotFitProp{j};
                %____________________________________________________________________________________        

                j = j+1;
            end



            %_______________________________________________________________________________________
            % collect values for spots for all scenes, each cell is a nucleus
            %_______________________________________________________________________________________
            if(~isempty(spotDS.c1SpotProp{i}{1}(k))) 
                spotPropPerNuc{p}.nucPosEL = nucDS.c1NucProp{i}{1}(k).nucFracEL;
                spotPropPerNuc{p}.nucSpotMeanVals = nonzeros(spotDS.c1SpotProp{i}{1}(k).spotMeanVal);
                spotPropPerNuc{p}.nucSpotMeanVolDias = nonzeros(spotDS.c1SpotProp{i}{1}(k).spotDiaVol);  
                spotPropPerNuc{p}.nucSpotMeanPADias = nonzeros(spotDS.c1SpotProp{i}{1}(k).spotDiaPA);  
                spotPropPerNuc{p}.nucSpot2DFitDiasUM = nonzeros(temp2DfitDiaUM);
                spotPropPerNuc{p}.nucSpotFit2DTotVals = nonzeros(temp2DfitTotVal);         
                
                spotPropPerNuc{p}.nucSpot3DFitDiasUM = nonzeros(temp3DfitDiaUM);
                spotPropPerNuc{p}.nucSpotFit3DTotVals = nonzeros(temp3DfitTotVal); 
                
                spotPropPerNuc{p}.nucSpot2DFitDeconvDiasUM = nonzeros(temp2DfitDeconvDiaUM);
                spotPropPerNuc{p}.nucSpot3DFitDeconvDiasUM = nonzeros(temp3DfitDeconvDiaUM); 

                p = p+1;
            end

             fprintf('scene\t#%d\t|\tnuc\t#%d\tdone.\tTotal Spots = %d\n',i,k,(length(spotPropPerNuc{p-1}.nucSpotMeanVals)));  % counter
             k = k+1;
        end
    end
end

%________________________________________________________________________
% Compute spot properties by EL position of the nucleus
%________________________________________________________________________
spotPropPosBinDS = cell(1, length(positionBin));
spotMeanVals = cell(1, length(positionBin));
spotMeanVolDias = cell(1, length(positionBin));
spotMeanPADias = cell(1, length(positionBin));
spot2DFitDias = cell(1, length(positionBin));
spot2DFitVals = cell(1, length(positionBin));
spot3DFitDias = cell(1, length(positionBin));
spot3DFitVals = cell(1, length(positionBin));
spot2DFitDeconvDias = cell(1, length(positionBin));
spot3DFitDeconvDias = cell(1, length(positionBin));
totalSpots = cell(1, length(positionBin));
totalNucs = cell(1, length(positionBin));

for i=1:length(positionBin)
    spotPropPosBinDS{i}.spotMeanVals = 0;
    spotPropPosBinDS{i}.spotMeanVolDias = 0;
    spotPropPosBinDS{i}.spotMeanPADias = 0;
    spotPropPosBinDS{i}.spot2DFitDias = 0;
    spotPropPosBinDS{i}.spot2DFitVals = 0;    
    spotPropPosBinDS{i}.spot3DFitDias = 0;
    spotPropPosBinDS{i}.spot3DFitVals = 0;    
    spotPropPosBinDS{i}.spot2DFitDeconvDias = 0;
    spotPropPosBinDS{i}.spot3DFitDeconvDias = 0;
    spotPropPosBinDS{i}.totalSpots = 0;
end

for p = 1:length(spotPropPerNuc) % total number of nuclei
    i = 1;
    while i<length(positionBin) % total position bins
        % find which position bin a nucleus falls within
        if spotPropPerNuc{p}.nucPosEL > positionBin(i) && spotPropPerNuc{p}.nucPosEL <= positionBin(i+1) 
            
            if totalSpots{i}==0
                spotPropPosBinDS{i}.spotMeanVals = spotPropPerNuc{p}.nucSpotMeanVals;
                spotPropPosBinDS{i}.spotMeanVolDias =spotPropPerNuc{p}.nucSpotMeanVolDias;
                spotPropPosBinDS{i}.spotMeanPADias =spotPropPerNuc{p}.nucSpotMeanPADias;
                spotPropPosBinDS{i}.spot2DFitDias = spotPropPerNuc{p}.nucSpot2DFitDiasUM;
                spotPropPosBinDS{i}.spot2DFitVals = spotPropPerNuc{p}.nucSpotFit2DTotVals;    
                spotPropPosBinDS{i}.spot3DFitDias = spotPropPerNuc{p}.nucSpot3DFitDiasUM;
                spotPropPosBinDS{i}.spot3DFitVals = spotPropPerNuc{p}.nucSpotFit3DTotVals;    
                spotPropPosBinDS{i}.spot2DFitDeconvDias = spotPropPerNuc{p}.nucSpot2DFitDeconvDiasUM;
                spotPropPosBinDS{i}.spot3DFitDeconvDias = spotPropPerNuc{p}.nucSpot3DFitDeconvDiasUM;
                                
                spotPropPosBinDS{i}.totalSpots = length(spotPropPerNuc{p}.nucSpotMeanVals);  
                
                spotMeanVals{i} = spotPropPerNuc{p}.nucSpotMeanVals;
                spotMeanVolDias{i} =spotPropPerNuc{p}.nucSpotMeanVolDias;
                spotMeanPADias{i} =spotPropPerNuc{p}.nucSpotMeanPADias;
                spot2DFitDias{i} = spotPropPerNuc{p}.nucSpot2DFitDiasUM;
                spot2DFitVals{i} = spotPropPerNuc{p}.nucSpotFit2DTotVals;    
                spot3DFitDias{i} = spotPropPerNuc{p}.nucSpot3DFitDiasUM;
                spot3DFitVals{i} = spotPropPerNuc{p}.nucSpotFit3DTotVals;    
                totalSpots{i} = length(spotPropPerNuc{p}.nucSpotMeanVals);  
                spot2DFitDeconvDias{i} = spotPropPerNuc{p}.nucSpot2DFitDeconvDiasUM;
                spot3DFitDeconvDias{i} = spotPropPerNuc{p}.nucSpot3DFitDeconvDiasUM;
            else
                spotPropPosBinDS{i}.spotMeanVals = vertcat(spotPropPosBinDS{i}.spotMeanVals, spotPropPerNuc{p}.nucSpotMeanVals);
                spotPropPosBinDS{i}.spotMeanVolDias = vertcat(spotPropPosBinDS{i}.spotMeanVolDias, spotPropPerNuc{p}.nucSpotMeanVolDias);
                spotPropPosBinDS{i}.spotMeanPADias = vertcat(spotPropPosBinDS{i}.spotMeanPADias, spotPropPerNuc{p}.nucSpotMeanPADias);
                
                spotPropPosBinDS{i}.spot2DFitDias = vertcat(spotPropPosBinDS{i}.spot2DFitDias, spotPropPerNuc{p}.nucSpot2DFitDiasUM);
                spotPropPosBinDS{i}.spot2DFitVals = vertcat(spotPropPosBinDS{i}.spot2DFitVals, spotPropPerNuc{p}.nucSpotFit2DTotVals);    
                
                spotPropPosBinDS{i}.spot3DFitDias = vertcat(spotPropPosBinDS{i}.spot3DFitDias, spotPropPerNuc{p}.nucSpot3DFitDiasUM);
                spotPropPosBinDS{i}.spot3DFitVals = vertcat(spotPropPosBinDS{i}.spot3DFitVals, spotPropPerNuc{p}.nucSpotFit3DTotVals);    
                
                spotPropPosBinDS{i}.spot2DFitDeconvDias = vertcat(spotPropPosBinDS{i}.spot2DFitDeconvDias, spotPropPerNuc{p}.nucSpot2DFitDeconvDiasUM);
                spotPropPosBinDS{i}.spot3DFitDeconvDias = vertcat(spotPropPosBinDS{i}.spot3DFitDeconvDias, spotPropPerNuc{p}.nucSpot3DFitDeconvDiasUM);

                
                spotPropPosBinDS{i}.totalSpots = vertcat(totalSpots{i}, length(spotPropPerNuc{p}.nucSpotMeanVals));    
                
                spotMeanVals{i} = vertcat(spotMeanVals{i}, spotPropPerNuc{p}.nucSpotMeanVals);
                spotMeanVolDias{i} = vertcat(spotMeanVolDias{i}, spotPropPerNuc{p}.nucSpotMeanVolDias);
                spotMeanPADias{i} = vertcat(spotMeanPADias{i}, spotPropPerNuc{p}.nucSpotMeanPADias);
                
                spot2DFitDias{i} = vertcat(spot2DFitDias{i}, spotPropPerNuc{p}.nucSpot2DFitDiasUM);
                spot2DFitVals{i} = vertcat(spot2DFitVals{i}, spotPropPerNuc{p}.nucSpotFit2DTotVals);  
                
                spot3DFitDias{i} = vertcat(spot3DFitDias{i}, spotPropPerNuc{p}.nucSpot3DFitDiasUM);
                spot3DFitVals{i} = vertcat(spot3DFitVals{i}, spotPropPerNuc{p}.nucSpotFit3DTotVals);  
                
                spot2DFitDeconvDias{i} = vertcat(spot2DFitDeconvDias{i}, spotPropPerNuc{p}.nucSpot2DFitDeconvDiasUM);
                spot3DFitDeconvDias{i} = vertcat(spot3DFitDeconvDias{i}, spotPropPerNuc{p}.nucSpot3DFitDeconvDiasUM);
                
                totalSpots{i} = vertcat(totalSpots{i}, length(spotPropPerNuc{p}.nucSpotMeanVals));   
            end
            
            totalNucs{i} = totalNucs{i} + 1;
            break            
        end
        spotPropPosBinDS{i}.spotSizeCutoff = spotXYpixelLim; %mindless repetition of the same value oever all position bins
        spotPropPosBinDS{i}.binLeftEdge = positionBin(i);
        i = i+1;
    end
end

save([DSfolder, filesep, fitDSName, '.mat'], 'fitDS');
save([DSfolder, filesep, posDSName, '.mat'], 'spotPropPosBinDS');
% histogram(spotMeanVals{:});

% Plotting spot intensities per egg length position
%______________________________________________________________________________________________
figure('Color',[1 1 1]);
title(['cutoff = ', num2str(spotXYpixelLim), ' pix'])
% Intensity of spots
y = num2cell(positionBin);
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spotMeanVals, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. intensity a. u');

% Plotting spot dias per egg length position
%______________________________________________________________________________________________
figure('Color',[1 1 1]);
title(['cutoff = ', num2str(spotXYpixelLim), ' pix'])
% Size of spots: dia assuming volume within sphere
subplot(4, 1, 1)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spotMeanVolDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. dia. (vol.) ({\mu}m');

% Size of spots: dia via principle axes of regionprops
subplot(4, 1, 2)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spotMeanPADias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. dia. (PA.) ({\mu}m');

% Size of spots: dia via 3D fitting
subplot(4, 1, 3)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot3DFitDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. dia. (3D fit) ({\mu}m');

% Size of spots: 3D deconv
subplot(4, 1, 4)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot3DFitDeconvDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. dia. (3D deconv) ({\mu}m');
xlabel('%EL');

% Total spots
%_______________________________________________________________________________________________
figure('Color', 'w');
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], cellfun(@mean, totalSpots, 'un', 0), y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('Agg. per nuc');
xlabel('%EL');

%	2D fit plots
%____________________________________________________________________________________
figure('Color',[1 1 1]);
title(['2D fit | cutoff = ', num2str(spotXYpixelLim), ' pix'])
% Intensity of spots
subplot(2, 1, 1)
y = num2cell(positionBin);
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot2DFitVals, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('2D Agg. intensity a. u');
% Size of spots
subplot(2, 1, 2)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot2DFitDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('2D Agg. dia. ({\mu}m');

%	3D fit plots
%___________________________________________________________________________________
figure('Color',[1 1 1]);
title(['3D fit | cutoff = ', num2str(spotXYpixelLim), ' pix']);
% Intensity of spots
subplot(2, 1, 1)
y = num2cell(positionBin);
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot3DFitVals, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('3D Agg. intensity a. u');
% Size of spots
subplot(2, 1, 2)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot3DFitDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('3D Agg. dia. ({\mu}m');
ylim([0, 2]);

%	deconv fit plots
%__________________________________________________________________________________
figure('Color',[1 1 1]);
title(['Deconv fit | cutoff = ', num2str(spotXYpixelLim), ' pix']);
% Intensity of spots
subplot(2, 1, 1)
y = num2cell(positionBin);
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot2DFitDeconvDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('2D deconv dia ({\mu}m');
% Size of spots
subplot(2, 1, 2)
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], spot3DFitDeconvDias, y, 'un', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))
ylabel('3D deconv dia. ({\mu}m');
ylim([0, 2]);
end
%   End of main function
%______________________________________________________________________________________________


%______________________________________________________________________________________________
%
%                               Gaussian fit of individual in the z plane
%______________________________________________________________________________________________
function [zSigma] = spotFitZ(spotVoxCent, nucVoxList, nucValList, winSize, metaDataDS)
winSize = winSize;
spotCent = round(spotVoxCent);
z = -winSize:winSize;
winMat = zeros((2*winSize + 1), 1); %define z matrix around the spot centroid
winZRange = bsxfun(@plus, -winSize:winSize, spotCent(3)); % real subs of image pixels

for i = 1:length(winZRange)
        % Find the relevant voxel index from the voxel sub array
        tempIdx = find(nucVoxList{1}(:,1) == spotCent(1) &...
            nucVoxList{1}(:,2) == spotCent(2) & ...
            nucVoxList{1}(:,3) == winZRange(i));
        
        % Construct the matrix with values from the window voxels
        if length(tempIdx) ==1
            winMat(i) = nucValList{1}(tempIdx);
        end
end

F = griddedInterpolant(z, winMat);
xq = -winSize:0.2:winSize;
vq = F(xq);
A = [winMat(winSize+1), 0, winSize/2, 0];
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f1 = fit(xq',vq',gaussEqn,'Start', A);

% plot(f1, xq, vq); hold on;

zSigma = f1.c/sqrt(2);

% zPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
% zPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPsfUM;
% psfSigma = (zPsfUM/zPixUM)/(2*sqrt(2*log(2)));
% 
% zFit = exp(-((z)./f1.c).^2);
% zPsf = exp(-((z)./psfSigma).^2);

end

%______________________________________________________________________________________________
% 2D Gaussian fit of individual in x-y plane
% Works for square windows only
%______________________________________________________________________________________________
function [spotXYFitProp] = spotFitXY(spotVoxCent, nucVoxList, nucValList, winSize)

fitFunction = 'ido';	% 'rot' for angle, 'norot'  for no angle, 'ido' for ido golding scheme

spotCent = round(spotVoxCent);
winArr = -winSize:winSize;
winMat = zeros(2*winSize + 1); %	define matrix around the spot centroid
winXRange = bsxfun(@plus, -winSize:winSize, spotCent(2)); %	real subs of image pixels
winYRange = bsxfun(@plus, -winSize:winSize, spotCent(1)); % real subs of image pixels

for i = 1:length(winXRange)
    for j = 1:length(winYRange)
        %	Find the relevant voxel indez from the voxel sub array
        tempIdx = find(nucVoxList{1}(:,1) == winYRange(j) & ...
            nucVoxList{1}(:,2) == winXRange(i) & ...
            nucVoxList{1}(:,3) == spotCent(3));
        
        %	Construct the matrix with values from the window voxels
        if length(tempIdx) ==1
            winMat(j, i) = nucValList{1}(tempIdx);
        end
    end
end

%_____________________________________________________________________________________________
% Plot smooth surface from vals and add real datapoints
%_____________________________________________________________________________________________
% figure('Color', [1, 1, 1]);
% 
% %   Interpolated surface plot
% % winMatInterp = interp2(winMat,5,'bicubic');
% % ss = surf(winMatInterp);
% 
% % real surface plot
% % ss = surf(winXRange, winYRange, winMat); % smooth val profile with real pix subs
% ss = surf(winArr, winArr, winMat); % smooth val profile with center relative pix subs
% ss.EdgeColor = 'none';
% ss.FaceAlpha =  0.5000;
% ss.FaceColor = 'interp';
% ax = gca;
% ax.BoxStyle = 'full';
% % ax.XTick = [0 40 80 120 160];
% % ax.YTick = [0 40 80 120 160];
% grid off;
% shading interp;
% hold on;
% 
% 
% % allX = repmat(winXRange, length(winYRange),1 ); % with real pix subs
% % allY = repmat(winYRange', 1, length(winXRange) ); % with real pix subs
% 
% allX = repmat(winArr, length(winArr),1 ); % relative pix subs
% allY = repmat(winArr', 1, length(winArr) );% relative pix subs
% 
% %   Scatter plot of the data points 
% sc = scatter3(allX(:), allY(:), winMat(:)); % actual val data points
% sc.MarkerEdgeColor = [0 0 0];
% sc.MarkerFaceColor = [0 .75 .75];
% hold on;
%_____________________________________________________________________________________________


%_____________________________________________________________________________________________
%   Fit data to 2D gaussian.
%_____________________________________________________________________________________________
%   2D Gaussian fit of the data
%   X-mat convention: size(n,n,2) where  %	X(:,:,1) : x-coord,  %	X(:,:,2) : y-coord.
%_____________________________________________________________________________________________
% 	2D Gaussian function ( A requires 5 coefs )
gauss2D = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );

%   2D Rotated Gaussian function ( A requires 6 coefs ).
%   Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DRot = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );%+ A(7);

%   Copy of golding 2015 ( A requires 7 coefs ).
%   Coeficients A convention:
%	A = [I_amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DIdo = @(A,X) A(1).*exp(-((( (cos(A(6)).^2)/(2.*A(3).^2) ) + ...
    ( (sin(A(6)).^2)/(2.*A(5).^2) )).*(X(:,:,1) - A(2)).^2 + 2.*(( -(sin(2*A(6))/(4.*A(3).^2)) ) + ...
    ( sin(2.*A(6))/(4.*A(5).^2) )).*(X(:,:,1) - A(2)).*(X(:,:,2) - A(4)) + ...
    (sin(A(6)).^2/(2.*A(3).^2) + cos(A(6)).^2/(2.*A(5).^2)).*(X(:,:,2) - A(4)).^2)) + A(7);
%______________________________________________________________________________________________

%   Work on the coefficients
%   coeffMat format :  [Amp,xo,sigmax,yo,sigmay,rotAngle]
%______________________________________________________________________________________________
ampVal = winMat(winSize, winSize); % Center voxel value
x0 = 0; % Relative center
y0 = 0; % relative center
xWidth = winSize/2; % Sigma in x direction
yWidth = winSize/2; % Sigma in y direction
rotAngle = 0;
bgVal = 100;
coeffMatInit = [ampVal, x0, xWidth, y0, yWidth, rotAngle, bgVal];   % Inital (guess) parameters

%	Define lower and upper bounds [Amp,xo,sigmax,yo,sigmay,rotAngle]
lBound = [0, -winSize, 0, -winSize, 0, 0, 0];
uBound = [realmax('double'), winSize, (winSize)^2, winSize, (winSize)^2, pi/4, realmax('double')];

%	Fit data and get the fit coefficients
%________________________________________________________________________________________________

%   Build numerical Grids
[x, y] = meshgrid(-winSize:winSize, -winSize:winSize); 
X = zeros(2*winSize + 1, 2*winSize + 1, 2); 
X(:, :, 1)=x; 
X(:, :, 2)=y;

% %	High Resolution Grid
% resFactor = 3; 
% [xh, yh] = meshgrid(-winSize:1/resFactor:winSize, -winSize:1/resFactor:winSize); 
% Xh = zeros(resFactor*2*winSize + 1, resFactor*2*winSize + 1, 2); 
% Xh(:,:,1) = xh; 
% Xh(:,:,2) = yh;

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'off');

switch fitFunction
    case 'norot', [coeffMat] = ...
            lsqcurvefit(gauss2D, coeffMatInit(1:5), X, winMat, lBound(1:5), uBound(1:5), options);
    case 'rot',  [coeffMat] = ...
            lsqcurvefit(gauss2DRot, coeffMatInit(1:7), X, winMat, lBound, uBound, options);
    case 'ido',  [coeffMat] = ...
            lsqcurvefit(gauss2DIdo, coeffMatInit(1:7), X, winMat, lBound, uBound, options);
    otherwise, error('invalid entry');
end
%________________________________________________________________________________________________

% disp(output); % display summary of LSQ algorithm

%________________________________________________________________________________________________
% Plot 3D Data and Fitted curve
%________________________________________________________________________________________________

% %   Plot the actual data surface
% ss = surface(x, y, winMat); 
% ss.EdgeColor = 'none';
% ss.FaceAlpha =  0.7000;
% ss.FaceColor = 'interp';
% colormap pink;
% hold on

% %   Scatter plot of the data points 
% allX = repmat(winArr, length(winArr),1 ); % relative pix subs
% allY = repmat(winArr', 1, length(winArr) );% relative pix subs
% sc = scatter3(allX(:), allY(:), winMat(:)); % actual val data points
% sc.MarkerEdgeColor = [0 0 0];
% sc.MarkerFaceColor = [0 .75 .75];
% hold on;
% 
% %   Plot of the fit curve
% switch fitFunction
%     case 'norot'
%         C = del2(gauss2D(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2D(coeffMat,Xh),C); 
%         hold on
%     case 'rot'
%         C = del2(gauss2DRot(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2DRot(coeffMat,Xh),C); 
%         hold on
%     case 'ido'
%         C = del2(gauss2DIdo(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2DIdo(coeffMat,Xh),C); 
%         hold on
% end
% 
% mm.EdgeColor = 'flat';
% mm.LineStyle = '-';
% mm.FaceColor = [1 1 1];
% mm.FaceLighting = 'none';
% mm.FaceAlpha = 0.5000;
% view(-60,20); 
% grid off;
% hold off;
% _____________________________________________________________________________________________
%   Turn off fit related warnings
% w = warning('query','last');
% id = w.identifier;
% warning('off',id);

% % _____________________________________________________________________________________________
% %   Return coefficient values
% % 	coeffMat format :  [Amp,xo,sigmax,yo,sigmay,rotAngle]
% % _____________________________________________________________________________________________
spotXYFitProp.spotAmp = coeffMat(1);
spotXYFitProp.x0 = coeffMat(2);
spotXYFitProp.xSigma = coeffMat(3);
spotXYFitProp.y0 = coeffMat(4);
spotXYFitProp.ySigma = coeffMat(5);
if ~strcmp(fitFunction, 'norot')
    spotXYFitProp.rotAngle = coeffMat(6);
    spotXYFitProp.spotBg = coeffMat(7);
end
%_____________________________________________________________________________________________
end

%______________________________________________________________________________________________
% 2D Gaussian fit of individual in x-z plane
% Works for square windows only
%______________________________________________________________________________________________
function [spotXZFitProp] = spotFitXZ(spotVoxCent, nucVoxList, nucValList, winSize)

fitFunction = 'ido';	% 'rot' for angle, 'norot'  for no angle, 'ido' for ido golding scheme

spotCent = round(spotVoxCent);
winArr = -winSize:winSize;
winMat = zeros(2*winSize + 1); %	define matrix around the spot centroid
winXRange = bsxfun(@plus, -winSize:winSize, spotCent(2)); %	real subs of image x pixels
winZRange = bsxfun(@plus, -winSize:winSize, spotCent(3)); % real subs of image z pixels

for i = 1:length(winXRange)
    for j = 1:length(winZRange)
        %	Find the relevant voxel index from the voxel sub array
        tempIdx = find(nucVoxList{1}(:,1) == spotCent(2) & ...
            nucVoxList{1}(:,2) == winXRange(i) & ...
            nucVoxList{1}(:,3) == winZRange(j));
        
        %	Construct the matrix with values from the window voxels
        if length(tempIdx) ==1
            winMat(j, i) = nucValList{1}(tempIdx);
        end
    end
end

%_____________________________________________________________________________________________
% Plot smooth surface from vals and add real datapoints
%_____________________________________________________________________________________________
% figure('Color', [1, 1, 1]);
% 
% %   Interpolated surface plot
% % winMatInterp = interp2(winMat,5,'bicubic');
% % ss = surf(winMatInterp);
% 
% % real surface plot
% % ss = surf(winXRange, winYRange, winMat); % smooth val profile with real pix subs
% ss = surf(winArr, winArr, winMat); % smooth val profile with center relative pix subs
% ss.EdgeColor = 'none';
% ss.FaceAlpha =  0.5000;
% ss.FaceColor = 'interp';
% ax = gca;
% ax.BoxStyle = 'full';
% % ax.XTick = [0 40 80 120 160];
% % ax.YTick = [0 40 80 120 160];
% grid off;
% shading interp;
% hold on;
% 
% 
% % allX = repmat(winXRange, length(winYRange),1 ); % with real pix subs
% % allY = repmat(winYRange', 1, length(winXRange) ); % with real pix subs
% 
% allX = repmat(winArr, length(winArr),1 ); % relative pix subs
% allY = repmat(winArr', 1, length(winArr) );% relative pix subs
% 
% %   Scatter plot of the data points 
% sc = scatter3(allX(:), allY(:), winMat(:)); % actual val data points
% sc.MarkerEdgeColor = [0 0 0];
% sc.MarkerFaceColor = [0 .75 .75];
% hold on;
%_____________________________________________________________________________________________


%_____________________________________________________________________________________________
%   Fit data to 2D gaussian.
%_____________________________________________________________________________________________
%   2D Gaussian fit of the data
%   X-mat convention: size(n,n,2) where  %	X(:,:,1) : x-coord,  %	X(:,:,2) : y-coord.
%_____________________________________________________________________________________________
% 	2D Gaussian function ( A requires 5 coefs )
gauss2D = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );

%   2D Rotated Gaussian function ( A requires 6 coefs ).
%   Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DRot = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );%+ A(7);

%   Copy of golding 2015 ( A requires 7 coefs ).
%   Coeficients A convention:
%	A = [I_amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DIdo = @(A,X) A(1).*exp(-((( (cos(A(6)).^2)/(2.*A(3).^2) ) + ...
    ( (sin(A(6)).^2)/(2.*A(5).^2) )).*(X(:,:,1) - A(2)).^2 + 2.*(( -(sin(2*A(6))/(4.*A(3).^2)) ) + ...
    ( sin(2.*A(6))/(4.*A(5).^2) )).*(X(:,:,1) - A(2)).*(X(:,:,2) - A(4)) + ...
    (sin(A(6)).^2/(2.*A(3).^2) + cos(A(6)).^2/(2.*A(5).^2)).*(X(:,:,2) - A(4)).^2)) + A(7);
%______________________________________________________________________________________________

%   Work on the coefficients
%   coeffMat format :  [Amp,xo,sigmax,yo,sigmay,rotAngle]
%______________________________________________________________________________________________
ampVal = winMat(winSize, winSize); % Center voxel value
x0 = 0; % Relative center
z0 = 0; % relative center
xWidth = winSize/2; % Sigma in x direction
zWidth = winSize/2; % Sigma in y direction
rotAngle = 0;
bgVal = 100;
coeffMatInit = [ampVal, x0, xWidth, z0, zWidth, rotAngle, bgVal];   % Inital (guess) parameters

%	Define lower and upper bounds [Amp,xo,sigmax,yo,sigmay,rotAngle]
lBound = [0, -winSize, 0, -winSize, 0, 0, 0];
uBound = [realmax('double'), winSize, (winSize)^2, winSize, (winSize)^2, pi/4, realmax('double')];

%	Fit data and get the fit coefficients
%________________________________________________________________________________________________

%   Build numerical Grids
[x, z] = meshgrid(-winSize:winSize, -winSize:winSize); 
X = zeros(2*winSize + 1, 2*winSize + 1, 2); 
X(:, :, 1)=x; 
X(:, :, 2)=z;

%	High Resolution Grid
resFactor = 3; 
[xh, yh] = meshgrid(-winSize:1/resFactor:winSize, -winSize:1/resFactor:winSize); 
Xh = zeros(resFactor*2*winSize + 1, resFactor*2*winSize + 1, 2); 
Xh(:,:,1) = xh; 
Xh(:,:,2) = yh;

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'off');

switch fitFunction
    case 'norot', [coeffMat] = ...
            lsqcurvefit(gauss2D, coeffMatInit(1:5), X, winMat, lBound(1:5), uBound(1:5), options);
    case 'rot',  [coeffMat] = ...
            lsqcurvefit(gauss2DRot, coeffMatInit(1:7), X, winMat, lBound, uBound, options);
    case 'ido',  [coeffMat] = ...
            lsqcurvefit(gauss2DIdo, coeffMatInit(1:7), X, winMat, lBound, uBound, options);
    otherwise, error('invalid entry');
end
%________________________________________________________________________________________________

% disp(output); % display summary of LSQ algorithm

%________________________________________________________________________________________________
% Plot 3D Data and Fitted curve
%________________________________________________________________________________________________

%   Plot the actual data surface
ss = surface(x, z, winMat); 
ss.EdgeColor = 'none';
ss.FaceAlpha =  0.7000;
ss.FaceColor = 'interp';
colormap pink;
hold on

%   Scatter plot of the data points 
allX = repmat(winArr, length(winArr),1 ); % relative pix subs
allZ = repmat(winArr', 1, length(winArr) );% relative pix subs
sc = scatter3(allX(:), allZ(:), winMat(:)); % actual val data points
sc.MarkerEdgeColor = [0 0 0];
sc.MarkerFaceColor = [0 .75 .75];
hold on;

%   Plot of the fit curve
switch fitFunction
    case 'norot'
        C = del2(gauss2D(coeffMat,Xh)); 
        mm = mesh(xh,yh,gauss2D(coeffMat,Xh),C); 
        hold on
    case 'rot'
        C = del2(gauss2DRot(coeffMat,Xh)); 
        mm = mesh(xh,yh,gauss2DRot(coeffMat,Xh),C); 
        hold on
    case 'ido'
        C = del2(gauss2DIdo(coeffMat,Xh)); 
        mm = mesh(xh,yh,gauss2DIdo(coeffMat,Xh),C); 
        hold on
end

mm.EdgeColor = 'flat';
mm.LineStyle = '-';
mm.FaceColor = [1 1 1];
mm.FaceLighting = 'none';
mm.FaceAlpha = 0.5000;
view(-60,20); 
grid off;
hold off;
%_____________________________________________________________________________________________
%   Turn off fit related warnings
w = warning('query','last');
id = w.identifier;
warning('off',id);

%_____________________________________________________________________________________________
%   Return coefficient values
%	coeffMat format :  [Amp,xo,sigmax,yo,sigmay,rotAngle]
%_____________________________________________________________________________________________
spotXZFitProp.spotAmp = coeffMat(1);
spotXZFitProp.x0 = coeffMat(2);
spotXZFitProp.xSigma = coeffMat(3);
spotXZFitProp.y0 = coeffMat(4);
spotXZFitProp.ySigma = coeffMat(5);
if ~strcmp(fitFunction, 'norot')
    spotXZFitProp.rotAngle = coeffMat(6);
    spotXZFitProp.spotBg = coeffMat(7);
end
%_____________________________________________________________________________________________
end

%_____________________________________________________________________________________________
% Deconvolve raw spots
%_____________________________________________________________________________________________
function [deconvIm] = spotBBdeconv(center, nucVoxList, nucValList, spotWindow, metaDataDS)
spotWindow = spotWindow*2;
xPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.yPixUM;
xyPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.xyPsfUM;
zPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
zPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPsfUM;

spotMat = zeros((2*spotWindow+1), (2*spotWindow+1), (2*spotWindow+1));
spotCent = round(center);

win1Range = bsxfun(@plus, -spotWindow:spotWindow, spotCent(1)); %	real subs of image pixels
win2Range = bsxfun(@plus, -spotWindow:spotWindow, spotCent(2)); % real subs of image pixels
win3Range = bsxfun(@plus, -spotWindow:spotWindow, spotCent(3)); % real subs of image pixels

for i = 1:(2*spotWindow+1) %    (1)
    for j = 1:(2*spotWindow+1) %    (2)
        for k = 1:(2*spotWindow+1) %    (3)
            %   Find the relevant voxel index from the voxel-sub array
            tempIdx = ...
                find(nucVoxList(:,1) == win1Range(i) &...
                nucVoxList(:,2) == win2Range(j) & ...
                nucVoxList(:,3) == win3Range(k));

            %	Construct the matrix with values from the window voxels
            if length(tempIdx) ==1
                spotMat(i, j, k) = nucValList(tempIdx);
            end
        end
    end
end

sigma(1) = (xyPsfUM/xPixUM)/(2*sqrt(2*log(2)));
sigma(2) = (xyPsfUM/yPixUM)/(2*sqrt(2*log(2)));
sigma(3) = (zPsfUM/zPixUM)/(2*sqrt(2*log(2)));

hSize = sigma;
hSize(:) = 2*ceil(max(sigma));

psf3 = fspecial3('gaussian',hSize,sigma)/sum(ones(hSize),'all');
psf2 = fspecial('gaussian',hSize(1:2),sigma(1))/sum(ones(hSize(1:2)),'all');

I = zeros(size(spotMat));
for i=1:size(spotMat, 3)
   I(:,:,i) = edgetaper(spotMat(:,:,i),psf2);
end
J = deconvlucy(I, psf3, 100);
% JJ = deconvblind({J}, {psf});
deconvIm = (J);
end

%_____________________________________________________________________________________________
%   Take spot's XY sigma value from gauss fit. Reconstruct a gaussian and
%   deconvolve it with a psf based gaussian. Then fit the deconvolved
%   gaussian again and extract the parameters.
%_____________________________________________________________________________________________
function [deconvPostFitProp] =  spotPostFitDeconv(spotFitProp, spotWindow, metaDataDS)
spotWindow = spotWindow*2;
xPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.xPixUM;
yPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.yPixUM;
xyPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.xyPsfUM;
zPixUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPixUM;
zPsfUM = metaDataDS.seriesMetaDataDS.analysisInfo.zPsfUM;
%%
psfSigma(1) = (xyPsfUM/xPixUM)/(2*sqrt(2*log(2)));
psfSigma(2) = (xyPsfUM/yPixUM)/(2*sqrt(2*log(2)));
psfSigma(3) = (zPsfUM/zPixUM)/(2*sqrt(2*log(2)));

hSize = psfSigma;
hSize(:) = 4*ceil(2*psfSigma)+1;
psf3 = fspecial3('gaussian',hSize,psfSigma)/sum(ones(hSize),'all');
psf2 = fspecial('gaussian',hSize(1:2),psfSigma(1))/sum(ones(hSize(1:2)),'all');
%%

fitSigmaX = spotFitProp.xSigma;
fitSigmaY = spotFitProp.ySigma;

x = -spotWindow:spotWindow;
y = -spotWindow:spotWindow;

% [X, Y] = meshgrid(x, y);

[x, y] = meshgrid(-spotWindow:spotWindow, -spotWindow:spotWindow); 
X = zeros(2*spotWindow + 1, 2*spotWindow + 1, 2); 
X(:, :, 1)=x; 
X(:, :, 2)=y;

% spotMat = rescale(exp(-X(:, :, 1).^2/(2*fitSigmaX^2)-X(:, :, 2).^2/(2*fitSigmaY^2)));

A = [1, 0, fitSigmaX, 0, fitSigmaY, spotFitProp.rotAngle, 0];

spotMat = A(1).*exp(-((( (cos(A(6)).^2)/(2.*A(3).^2) ) + ...
    ( (sin(A(6)).^2)/(2.*A(5).^2) )).*(X(:,:,1) - A(2)).^2 + 2.*(( -(sin(2*A(6))/(4.*A(3).^2)) ) + ...
    ( sin(2.*A(6))/(4.*A(5).^2) )).*(X(:,:,1) - A(2)).*(X(:,:,2) - A(4)) + ...
    (sin(A(6)).^2/(2.*A(3).^2) + cos(A(6)).^2/(2.*A(5).^2)).*(X(:,:,2) - A(4)).^2)) + A(7);

deconvSpotMat = deconvlucy(spotMat, psf2, 100);
%   Fit 2D gaussian
coeffMat = gaussianFit(deconvSpotMat, spotWindow, 'ido'); % use ido

deconvPostFitProp.xSigmaDeconv = coeffMat(3);
deconvPostFitProp.ySigmaDeconv = coeffMat(5);
end


%_____________________________________________________________________________________________
%   Fit data to 2D gaussian.
%_____________________________________________________________________________________________
function coeffMat = gaussianFit(im, winSize, whichFun)
im = rescale(im);
fitFunction = whichFun;
%   X-mat convention: size(n,n,2) where  %	X(:,:,1) : x-coord,  %	X(:,:,2) : y-coord.
%_____________________________________________________________________________________________
% 	2D Gaussian function ( A requires 5 coefs )
gauss2D = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );

%   2D Rotated Gaussian function ( A requires 6 coefs ).
%   Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DRot = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );%+ A(7);

%   Copy of golding 2015 ( A requires 7 coefs ).
%   Coeficients A convention:
%	A = [I_amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), I_background]
gauss2DIdo = @(A,X) A(1).*exp(-((( (cos(A(6)).^2)/(2.*A(3).^2) ) + ...
    ( (sin(A(6)).^2)/(2.*A(5).^2) )).*(X(:,:,1) - A(2)).^2 + 2.*(( -(sin(2*A(6))/(4.*A(3).^2)) ) + ...
    ( sin(2.*A(6))/(4.*A(5).^2) )).*(X(:,:,1) - A(2)).*(X(:,:,2) - A(4)) + ...
    (sin(A(6)).^2/(2.*A(3).^2) + cos(A(6)).^2/(2.*A(5).^2)).*(X(:,:,2) - A(4)).^2)) + A(7);
%______________________________________________________________________________________________

%   Work on the coefficients
%   coeffMat format :  [Amp,xo,sigmax,yo,sigmay,rotAngle]
%______________________________________________________________________________________________
ampVal = im(winSize+1, winSize+1); % Center voxel value
x0 = 0; % Relative center
y0 = 0; % relative center
xWidth = winSize/2; % Sigma in x direction
yWidth = winSize/2; % Sigma in y direction
rotAngle = 0;
bgVal = 0;
coeffMatInit = [ampVal, x0, xWidth, y0, yWidth, rotAngle, bgVal];   % Inital (guess) parameters

%	Define lower and upper bounds [Amp,xo,sigmax,yo,sigmay,rotAngle]
lBound = [0, -winSize, 0, -winSize, 0, 0, 0];
uBound = [realmax('double'), winSize, (winSize)^2, winSize, (winSize)^2, pi/4, realmax('double')];
%	Fit data and get the fit coefficients
%________________________________________________________________________________________________

%   Build numerical Grids
[x, y] = meshgrid(-winSize:winSize, -winSize:winSize); 
X = zeros(2*winSize + 1, 2*winSize + 1, 2); 
X(:, :, 1)=x; 
X(:, :, 2)=y;

%	High Resolution Grid
resFactor = 3; 
[xh, yh] = meshgrid(-winSize:1/resFactor:winSize, -winSize:1/resFactor:winSize); 
Xh = zeros(resFactor*2*winSize + 1, resFactor*2*winSize + 1, 2); 
Xh(:,:,1) = xh; 
Xh(:,:,2) = yh;

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'off');

switch fitFunction
    case 'norot', [coeffMat] = ...
            lsqcurvefit(gauss2D, coeffMatInit(1:5), X, im, lBound(1:5), uBound(1:5), options); % major axis along x
    case 'rot',  [coeffMat] = ...
            lsqcurvefit(gauss2DRot, coeffMatInit(1:7), X, im, lBound, uBound, options); % rotated with angle 
    case 'ido',  [coeffMat] = ...
            lsqcurvefit(gauss2DIdo, coeffMatInit(1:7), X, im, lBound, uBound, options); % rotated with angle use this
    otherwise, error('invalid entry');
end

%   Plot the actual data surface
%_____________________________________________________________________________
% ss = surf(x, y, im); 
% ss.EdgeColor = 'none';
% ss.FaceAlpha =  0.7000;
% ss.FaceColor = 'interp';
% colormap pink;
% hold on
% 
% winArr = -winSize:winSize;
% %   Scatter plot of the data points 
% allX = repmat(winArr, length(winArr),1 ); % relative pix subs
% allY = repmat(winArr', 1, length(winArr) );% relative pix subs
% sc = scatter3(allX(:), allY(:), im(:)); % actual val data points
% sc.MarkerEdgeColor = [0 0 0];
% sc.MarkerFaceColor = [0 .75 .75];
% hold on;

%   Plot of the fit curve
%______________________________________________________________________________
% switch fitFunction
%     case 'norot'
%         C = del2(gauss2D(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2D(coeffMat,Xh),C); 
%         hold on
%     case 'rot'
%         C = del2(gauss2DRot(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2DRot(coeffMat,Xh),C); 
%         hold on
%     case 'ido'
%         C = del2(gauss2DIdo(coeffMat,Xh)); 
%         mm = mesh(xh,yh,gauss2DIdo(coeffMat,Xh),C); 
%         hold on
% end
end


 