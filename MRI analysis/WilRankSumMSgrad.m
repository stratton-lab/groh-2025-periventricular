function [p,h,stats] = WilRankSumMSgrad(qMRIType,slabFols,geoBand,angBand,angBands,geoMax)

%INPUTS:
%1) qMRIType: Type of qMRI metric 
%T2star = Transverse relaxation time, T2* 
%T1 = Longitudinal relaxation time, T1 

%2) slabFols: Array of strings giving folder names of slabs. Folder name
%corresponding to MS slab must be first, followed by the HC slabs 

%4) geoBand: (number of) geodesic band to perform stats on

%5) angBand: (number of) angular band to perform stats on

%6) angBands: noSlabs-length cell array, where ith cell gives array of "angular bands to include" 
% for corresponding slab, for computing slab mean or median (order must be in agreement with 'slabFols') 

%7) geoMax: noSlabs-length cell array, where ith element gives furthest geodesic band for corresponding 
% slab, for computing slab mean or median (order must be in agreement with 'slabFols') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
mainDir = '/export02/data/risa/02_Ex-Vivo_MS_Slab_Scans/01_MEGRE_and_MP2RAGE/'; 

noSlabs = length(slabFols); qMRIdat = cell(1,noSlabs); qMRImeanHCs = qMRIdat; 
for i = 1:noSlabs
    slabFol = char(slabFols(i));
    Seg = load_untouch_nii([mainDir slabFol '/Segmentation.nii']).img;  

    %Load geodesic and angular band maps 
    geoBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/geodesic_band_map.nii']).img; 
    angBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/angular_band_map.nii']).img;
    qMRIType = string(lower(qMRIType)); 
    if qMRIType=="t2star"
        name = '/T2star_uncorr_AVG.nii'; 
    elseif qMRIType=="t1"
        name = '/T1MAP_AVG.nii';  
    end 
    qMRI = double(load_untouch_nii([mainDir slabFol name]).img); 
    
    %NaN all image voxels that are not WM 
    qMRI(Seg~=1) = NaN; 
    
    %if pia mask exists, null voxels corresponding to pia mater 
    piaMPath = [mainDir slabFol '/piaMask.nii']; sentinel = isfile(piaMPath); 
    if sentinel
        piaMask = double(load_untouch_nii(piaMPath).img); 
        qMRI(piaMask~=0) = NaN;
%     else %if no pia mater mask present, use denoising on whole 'geodesic band' qMRI distribution
%         tmp1 = qMRI; tmp1(geoBandMap~=geoBand) = NaN; 
%         [~,rmidxs] = rmoutliers(tmp1(:),'median','ThresholdFactor',3); %remove outliers within 3 MAD from median
%         qMRI(rmidxs) = NaN; 
    end 
    
    if i~=1 %HC brain samples 
        %if current angular band or geodesic band does not fall into selected set for HC slab, skip loop
        if ~any(angBands{i}==angBand) || (geoBand > geoMax{i}), qMRIdat{i} = []; continue; end        
    end 
    
    %extract QMRI voxelwise values in ROI defined by ith geodesic band and jth angular band from MS slab 
    qMRIdat{i} = qMRI(angBandMap==angBand & geoBandMap==geoBand);
    qMRIdat{i} = rmmissing(qMRIdat{i}); %remove missing entries 

    if sentinel %if pia mater exists, remove outliers in ROI_ij with 1.5 MAD of median instead 
        qMRIdat{i} = rmoutliers(qMRIdat{i},'median','ThresholdFactor',1.5); %remove outliers
    end         

    %if ROI_ij for slab comprises less than 100 voxels, omit it 
    if (length(qMRIdat{i}) < 100), qMRIdat{i} = []; end 

    if i == 1 %store qMRI mean of MS slab 
        qMRImeanMS = mean(qMRIdat{i},'all','omitnan');
    else %Normalize distribution of qMRI values for HC slab to have same qMRI mean as MS slab  
        tmp2 = []; 
        for a = 1:length(angBands{i})
            tmp = qMRI(angBandMap==angBands{i}(a) & geoBandMap <= geoMax{i});
            tmp2 = [tmp2 tmp']; 
        end
        qMRImeanHCs{i} = mean(tmp2,'all','omitnan');  
        qMRIdat{i} = qMRIdat{i} - qMRImeanHCs{i}; 
    end 
end 

%Using median of congregated qMRI WM distribution pertaining to HC slab group 
qMRIHC = []; 
for i = 2:noSlabs, qMRIHC = [qMRIHC qMRIdat{i}']; end 

if (length(qMRIHC) < 100) || isempty(qMRIdat{1}), p = NaN; h = 0; stats = NaN; return; end 

y = qMRIdat{1}; y = y(:); 
[p,h,stats] = ranksum(y,qMRIHC);  
end

