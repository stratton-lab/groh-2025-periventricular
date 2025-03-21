function [] = plotGradData(slabFol, tissueType, qMRIType, varargin)

%PURPOSE: Plot mean T2*/T1/QSM over geodesic bands (going ependyma-in) (and) 
%polar angle (counter-clockwise with respect to center line of ventricle) 

%INPUTS: 
%1) slabFol: string that gives name of slab 7 T MRI directory 

%2) tissueType: Type of tissue for tissue-specific T2*/T1/QSM gradients
%WM = white matter 
%cGM = cortical gray matter 

%3) qMRIType: Type of qMRI metric 
%T2star = Transverse relaxation time, T2* 
%T1 = Longitudinal relaxation time, T1 
%T2star&T1 = Plot mean T2* and T1 on same plot on different y-axes 

%OPTIONAL ARGUMENTS:
%1) varargin{1}: array giving numbers of geodesic bands to consider 
%2) varargin{2}: array giving numbers of angular bands to consider

%ASSUMPTIONS:
% -assuming unselected (inviable) slices in geodesic band map
% 'geodesic_band_map.nii' and CCW angular map 'angular_band_map.nii' have
% been zeroed out 

% -assuming geodesic distance map and CCW angular map have been masked
% with brain mask 'BrainMaskT2starW.nii'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 

%Load tissue-specific segmentation data 
mainDir = '/export02/data/risa/02_Ex-Vivo_MS_Slab_Scans/01_MEGRE_and_MP2RAGE/'; slabFol = char(slabFol); 
Seg = load_untouch_nii([mainDir slabFol '/Segmentation.nii']).img; 
segLbls = struct("WM",1,"CGM",2); 
segNames = struct("WM","White Matter","CGM","Cortical Gray Matter"); 
segLbl = getfield(segLbls,string(upper(tissueType))); 

%Load geodesic and angular band maps 
geoBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/geodesic_band_map.nii']).img; 
angBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/angular_band_map.nii']).img;

outFigPath = [mainDir slabFol '/figures/']; 
if ~isfolder(outFigPath), mkdir(outFigPath); end 

if nargin < 4 % 3 (official) arguments passed-in
    BandsGeo = 2:uint16(max(geoBandMap,[],'all','omitnan')); BandsAng = 2:uint16(max(angBandMap,[],'all','omitnan')); 
elseif nargin < 5 % 4 arguments passed-in (1 extra)
    BandsGeo = varargin{1}; BandsAng = 2:uint16(max(angBandMap,[],'all','omitnan'));
elseif nargin < 6 % 5 arguments passed-in (2 extras)
    BandsGeo = varargin{1}; BandsAng = varargin{2}; 
end 
noBandsGeo = length(BandsGeo); noBandsAng = length(BandsAng); 

%While loading qMRI data, NaN data from tissues other than the selected tissue type 
qMRIType = string(lower(qMRIType)); 

if qMRIType=="t2star&t1"
    names = ["/T2star_uncorr_AVG.nii" "/T1MAP_AVG.nii"]; 
    ylbls = ["Mean T_{2}^{*} (ms)" "Mean T_{1} (ms)"]; 
    tlbl = sprintf("Mean T2* & T1 over geodesic bands in %s",getfield(segNames,upper(tissueType))); 
    ylimits = {[10 60], [100 500]}; 
    qMRI = cell(1,2); 
    qMRI{1} = double(load_untouch_nii([mainDir slabFol char(names(1))]).img); 
    qMRI{1}(geoBandMap==0) = 0; qMRI{1}(Seg~=segLbl) = NaN;
    qMRI{2} = double(load_untouch_nii([mainDir slabFol char(names(2))]).img); 
    qMRI{2}(geoBandMap==0) = 0; qMRI{2}(Seg~=segLbl) = NaN;

    %if pia mask exists, zero out voxels corresponding to pia mater 
    piaMPath = [mainDir slabFol '/piaMask.nii']; sentinel = isfile(piaMPath); 
    if sentinel
        piaMask = double(load_untouch_nii(piaMPath).img); 
        qMRI{1}(piaMask~=0) = NaN; qMRI{2}(piaMask~=0) = NaN; 
    end 
    
    %A) Plot periventricular qMRI gradient over geodesic and angular bands 
    noFigures = uint32(ceil(noBandsAng/4)); 
    IdxAng = 0; 
    for nFig = 1:noFigures 
        ax = figure("WindowState","maximized"); 
        t = tiledlayout(2,2); 

        for nPlot = 1:4 
            IdxAng = IdxAng + 1; 
            %Compute mean and std of qMRI values over geodesic bands 
            %NOTE: We are assuming a NORMAL distribution of qMRI values within the
            %selected slice range, selected tissue type & selected band when computing 
            %mean and standard deviation values
            meanGeo = zeros(2,noBandsGeo); stdGeo = meanGeo;
            for q = 1:2
                for IdxGeo = 1:noBandsGeo
                    
                    if ~sentinel %if no pia mater mask present, use denoising on whole 'geodesic band' qMRI distribution 
                        tmp1 = qMRI{q}; tmp1(geoBandMap~=BandsGeo(IdxGeo)) = NaN; 
                        [~,rmidxs] = rmoutliers(tmp1(:),'median','ThresholdFactor',3); %remove outliers within 3 MAD from median
                        tmp1(rmidxs) = NaN; 

                        tmp = tmp1(angBandMap==BandsAng(IdxAng) & geoBandMap==BandsGeo(IdxGeo)); 
                        tmp = rmmissing(tmp); %remove missing entries 
                        %tmp = rmoutliers(tmp,'median','ThresholdFactor',1.5); %remove outliers                     
                    else %otherwise, use denoising on region ROI_ij only 
                        tmp = qMRI{q}(angBandMap==BandsAng(IdxAng) & geoBandMap==BandsGeo(IdxGeo)); 
                        tmp = rmmissing(tmp); %remove missing entries 
%                         tmp = rmoutliers(tmp,'median','ThresholdFactor',1.5); %remove outliers within 1.5 MAD from median 
                        tmp = rmoutliers(tmp,'median','ThresholdFactor',2.5); %TEMPORARY: for deep GM structures
                    end 
                    
                    meanGeo(q,IdxGeo) = mean(tmp,'all'); stdGeo(q,IdxGeo) = std(tmp,0,'all'); 
                    %fprintf("Angular Band %d, Geodesic Band %d, number of voxels = %d\n",IdxAng, IdxGeo, length(tmp)); 

                    %if segmented brain portion (according to angular band, geodesic band and tissue type) comprises less 
                    %than 100 voxels,then it is susceptible to noise -> qMRI average and std will be highly influenced by 
                    %noise, so omit! 
%                     if length(tmp) < 100, meanGeo(q,IdxGeo) = NaN; stdGeo(q,IdxGeo) = NaN; end 
                    if length(tmp) < 10, meanGeo(q,IdxGeo) = NaN; stdGeo(q,IdxGeo) = NaN; end %TEMPORARY: for deep GM structures, changed limit to 15 voxels
                end 
            end 
            plotTle = sprintf("Angular Band #%d",BandsAng(IdxAng));
            xlbl = "Geodesic band"; 
            nexttile; 

            yyaxis left; 
            errorbar(BandsGeo,meanGeo(1,:),stdGeo(1,:),'LineWidth',2.0);  
            ylim(ylimits{1}); xlabel(xlbl); ylabel(ylbls(1)); 

            yyaxis right; 
            errorbar(BandsGeo,meanGeo(2,:),stdGeo(2,:),'LineWidth',2.0);  
            ylim(ylimits{2}); xlabel(xlbl); ylabel(ylbls(2)); 
            xticks(BandsGeo); title(plotTle);         
            if (IdxAng == noBandsAng), break; end 
        end 
        title(t,tlbl);

        %Save figure
        figName = char(sprintf("/Fig_angular_bands_set_%d_%s.jpg",nFig,qMRIType)); 
        exportgraphics(ax, [outFigPath figName],'Resolution',300); 
    end 

    %B) Plot periventricular qMRI gradient over geodesic bands 
    ax = figure("WindowState","maximized");
    for q = 1:2 
        for IdxGeo = 1:noBandsGeo
            tmp = qMRI{q}(geoBandMap==BandsGeo(IdxGeo)); 
            tmp = rmmissing(tmp); %remove missing entries 
            tmp = rmoutliers(tmp,'median','ThresholdFactor',3); %remove outliers within 3 MAD from median 
            meanGeo(q,IdxGeo) = mean(tmp,'all'); stdGeo(q,IdxGeo) = std(tmp,0,'all'); 
%             if (length(tmp) < 200), meanGeo(q,IdxGeo) = NaN; stdGeo(q,IdxGeo) = NaN; end 
            if (length(tmp) < 20), meanGeo(q,IdxGeo) = NaN; stdGeo(q,IdxGeo) = NaN; end %TEMPORARY: for deep GM structures, changed limit to 30 voxels 
        end  
    end 
    xlbl = "Geodesic band"; 
    yyaxis left; 
    errorbar(BandsGeo,meanGeo(1,:),stdGeo(1,:),'LineWidth',1.2);  
    ylim(ylimits{1}); xlabel(xlbl); ylabel(ylbls(1)); 
    yyaxis right; 
    errorbar(BandsGeo,meanGeo(2,:),stdGeo(2,:),'LineWidth',1.2);  
    ylim(ylimits{2}); xlabel(xlbl); ylabel(ylbls(2)); 
    xticks(BandsGeo); title(tlbl); 

    %Save figure
    figName = char(sprintf("/Fig_geodesic_bands_%d_%s.jpg",noFigures+1,qMRIType));
    exportgraphics(ax, [outFigPath figName],'Resolution',300);

else 
    if qMRIType=="t2star"
        name = '/T2star_uncorr_AVG.nii'; 
        ylbl = "Mean T_{2}^{*} (ms)"; 
        tlbl = sprintf("Mean T2* over geodesic bands in %s",getfield(segNames,upper(tissueType))); 
        ylimit = [10 60]; 
    elseif qMRIType=="t1"
        name = '/T1MAP_AVG.nii'; 
        ylbl = "Mean T_{1} (ms)"; 
        tlbl = sprintf("Mean T1 over geodesic bands in %s",getfield(segNames,upper(tissueType))); 
        ylimit = [100 500]; 
    end 
    qMRI = double(load_untouch_nii([mainDir slabFol name]).img); 
    qMRI(geoBandMap==0) = 0; qMRI(Seg~=segLbl) = NaN; 
    
    %if pia mask exists, zero out voxels corresponding to pia mater 
    piaMPath = [mainDir slabFol '/piaMask.nii']; 
    if isfile(piaMPath)
        piaMask = double(load_untouch_nii(piaMPath).img); 
        qMRI{1}(piaMask~=0) = NaN; qMRI{2}(piaMask~=0) = NaN; 
    end 

    %A) Plot periventricular qMRI gradient over geodesic and angular bands 
    noFigures = uint32(ceil(noBandsAng/4)); 
    IdxAng = 0; 
    for nFig = 1:noFigures 
        ax = figure("WindowState","maximized"); 
        t = tiledlayout(2,2); 

        for nPlot = 1:4 
            IdxAng = IdxAng + 1; 
            %Compute mean and std of qMRI values over geodesic bands 
            %NOTE: We are assuming a NORMAL distribution of qMRI values within the
            %selected slice range, selected tissue type & selected band when computing 
            %mean and standard deviation values
            meanGeo = zeros(1,noBandsGeo); stdGeo = meanGeo;
            for IdxGeo = 1:noBandsGeo
                tmp = qMRI(angBandMap==BandsAng(IdxAng) & geoBandMap==BandsGeo(IdxGeo)); 
                tmp = rmmissing(tmp); %remove missing entries 
%                 tmp = rmoutliers(tmp,'median','ThresholdFactor',1.5); %remove outliers 
                tmp = rmoutliers(tmp,'median','ThresholdFactor',2.5); %TEMPORARY: for deep GM structures
                meanGeo(IdxGeo) = mean(tmp,'all'); stdGeo(IdxGeo) = std(tmp,0,'all'); 
                %fprintf("Angular Band %d, Geodesic Band %d, number of voxels = %d\n",IdxAng, IdxGeo, length(tmp)); 

                %if segmented brain portion (according to angular band,
                %geodesic band and tissue type) comprises less than 100 voxels,
                %then it is susceptible to noise -> qMRI average and std will
                %be highly influenced by noise, so omit! 
%                 if length(tmp) < 100, meanGeo(IdxGeo) = NaN; stdGeo(IdxGeo) = NaN; end 
                if length(tmp) < 10, meanGeo(IdxGeo) = NaN; stdGeo(IdxGeo) = NaN; end %TEMPORARY: for deep GM structures, changed limit to 15 voxels
            end 
            plotTle = sprintf("Angular Band #%d",BandsAng(IdxAng));
            xlbl = "Geodesic band"; 
            nexttile; 

            errorbar(BandsGeo,meanGeo,stdGeo,'LineWidth',1);  
            ylim(ylimit); xlabel(xlbl); ylabel(ylbl); 
            xticks(BandsGeo); title(plotTle); 
            if (IdxAng == noBandsAng), break; end 
        end 
        title(t,tlbl);

        %Save figure
        figName = char(sprintf("/Fig_angular_bands_set_%d_%s.jpg",nFig,qMRIType)); 
        exportgraphics(ax, [outFigPath figName],'Resolution',300); 
    end 

    %B) Plot periventricular qMRI gradient over geodesic bands 
    ax = figure("WindowState","maximized");
    for IdxGeo = 1:noBandsGeo
        tmp = qMRI(geoBandMap==BandsGeo(IdxGeo));
        tmp = rmmissing(tmp); %remove missing entries 
        tmp = rmoutliers(tmp,'median','ThresholdFactor',3); %remove outliers 
        meanGeo(IdxGeo) = mean(tmp,'all'); stdGeo(IdxGeo) = std(tmp,0,'all'); 
%         if (length(tmp) < 200), meanGeo(IdxGeo) = NaN; stdGeo(IdxGeo) = NaN; end   
        if (length(tmp) < 20), meanGeo(IdxGeo) = NaN; stdGeo(IdxGeo) = NaN; end %TEMPORARY: for deep GM structures, changed limit to 30 voxels 
    end  
    xlbl = "Geodesic band"; 
    errorbar(BandsGeo,meanGeo,stdGeo,'LineWidth',2.0);  
    ylim(ylimit); xlabel(xlbl); ylabel(ylbl); 
    xticks(BandsGeo); title(tlbl); 

    %Save figure
    figName = char(sprintf("/Fig_geodesic_bands_%d_%s.jpg",noFigures+1,qMRIType));
    exportgraphics(ax, [outFigPath figName],'Resolution',300); 
end 

end 