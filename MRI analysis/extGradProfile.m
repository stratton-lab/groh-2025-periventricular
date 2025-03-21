function [meanGeo,stdGeo] = extGradProfile(slabFol,angBandRange,geoBandRange,sliceRange,offset,qMRI,varargin)

%INPUTS: 
%1) slabFol: string that gives name of slab 7 T MRI directory 
%2) angBandRange
%3) geoBandRange
%4) sliceRange: 
%-slice1st = first permissible slice (corresponding to silver surface
%side), since Adam did tissue sectioning upwards from silver surface side
%-sliceLst = last permissible slice 
%NOTE #1: assuming index starts at 0 as in fsleyes
%NOTE #2: slice1st and sliceLst excludes 'slice offset'
%5) slice offset 
%6) qMRI: qMRI metric "T2star" or "T1"

%Since matlab starts indexing at 1
sliceRange = sliceRange + offset + 1; 

addpath(genpath("/export02/data/risa/NIfTI_20140122/"));
mainDir = '/export02/data/risa/02_Ex-Vivo_MS_Slab_Scans/01_MEGRE_and_MP2RAGE/'; 
slabFol = char(slabFol); 

%Load tissue-specific segmentation data 
Seg = load_untouch_nii([mainDir slabFol '/Segmentation.nii']).img; 
%Load geodesic and angular band maps 
geoBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/geodesic_band_map.nii']).img; 
angBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/angular_band_map.nii']).img;

if string(lower(qMRI))=="t2star"
    name = "/T2star_uncorr_AVG.nii"; 
elseif string(lower(qMRI))=="t1"
    name = "/T1MAP_AVG.nii"; 
end  
qMRI = double(load_untouch_nii([mainDir slabFol char(name)]).img); 
qMRI(geoBandMap==0) = 0; qMRI(Seg~=1) = NaN; %selecting cerebral WM 

%if pia mask exists, zero out voxels corresponding to pia mater 
piaMPath = [mainDir slabFol '/piaMask.nii']; sentinel = isfile(piaMPath); 
if sentinel
    piaMask = double(load_untouch_nii(piaMPath).img); 
    qMRI(piaMask~=0) = NaN; 
end 

noBands = length(geoBandRange); 
meanGeo = zeros(1,noBands); stdGeo = meanGeo;
for idxGeo = 1:noBands %going through geodesic bands 'a' and 'b'
    tmp = qMRI;  
    idxTmp = ~ismember(angBandMap,angBandRange); 
    tmp(idxTmp) = NaN; %select angular band range
    idxTmp = ~ismember(geoBandMap,geoBandRange(idxGeo)); 
    tmp(idxTmp) = NaN; %select current geodesic band

    %TESTING
    [~,rmidxs] = rmoutliers(tmp(:),'median','ThresholdFactor',1.5); %remove outliers within 1.5 MAD from median
    tmp(rmidxs) = NaN; 

    tmp = tmp(:,sliceRange,:); %select slice range    
    meanGeo(idxGeo) = mean(tmp,'all','omitnan'); stdGeo(idxGeo) = std(tmp,0,'all','omitnan'); 
    if length(tmp) < 10, meanGeo(idxGeo) = NaN; stdGeo(idxGeo) = NaN; end 
end 
meanGeo = meanGeo - meanGeo(1); 

end 


