function [T, T1maxdf, T2starmaxdf, T1del, T2stardel] = maxRateOfChange(slabFol,patholData,geoBandWdth)
%patholData: array of 4 columns
%1st column lists angular region, 
%2nd column lists corresponding geodesic band 'a', 
%3rd column lists corresponding geodesic band 'b', 
%and 4th column specifies whether gradient is 'ependyma-in' (1) or 'pia-in' (0)

tic 
noRows = size(patholData,1); slabFol = char(slabFol); 
T1dat = cell(1,noRows); T2stardat = T1dat;
T1df = cell(1,noRows); T2stardf = T1df; 
T1maxdf = zeros(1,noRows); T2starmaxdf = T1maxdf; 
T1del = zeros(1,noRows); T2stardel = T1del; 
try 
parpool(8); 
parfor row = 1:noRows
    angBand = patholData(row,1); %going through angular bands for slab
    bandA = patholData(row,2);
    bandB = patholData(row,3); 
    noBands = bandB-bandA+1; 
    T2stardat{row} = zeros(1,noBands); T1dat{row} = T2stardat{row};
    col = 1; 
    for band = bandA:bandB %going through geodesic bands 'a' and 'b' for slab
        tmp = getqMRIROIij("T2star",slabFol,band,angBand); 
        T2stardat{row}(col) = mean(tmp,'all','omitnan');
       
        tmp = getqMRIROIij("T1",slabFol,band,angBand); 
        T1dat{row}(col) = mean(tmp,'all','omitnan');    
        col = col+1; 
    end 
    T1df{row} = diff(T1dat{row})/geoBandWdth; %[ms/mm]
    T2stardf{row} = diff(T2stardat{row})/geoBandWdth; %[ms/mm]    
    if patholData(row,4) %ependyma-in gradient (negative T1/T2* change)
        T1maxdf(row) = min(T1df{row},[],'all','omitnan'); 
        T2starmaxdf(row) = min(T2stardf{row},[],'all','omitnan');
    else %pia-in gradient (positive T1/T2* change)
        T1maxdf(row) = max(T1df{row},[],'all','omitnan'); 
        T2starmaxdf(row) = max(T2stardf{row},[],'all','omitnan');        
    end 
    T1del(row) = T1dat{row}(end) - T1dat{row}(1); 
    T2stardel(row) = T2stardat{row}(end) - T2stardat{row}(1); 
end 
catch ME
    delete(gcp('nocreate'));
    rethrow(ME); 
end 
delete(gcp('nocreate'));
toc 

result = NaN*ones(noRows,7);
for n = 1:noRows
   angBand = patholData(n,1); geoBanda = patholData(n,2); geoBandb = patholData(n,3); 
   result(n,1) = angBand;
   result(n,2) = geoBanda; 
   result(n,3) = geoBandb; 
   result(n,4) = T1maxdf(n); 
   result(n,5) = T2starmaxdf(n); 
   result(n,6) = T1del(n); 
   result(n,7) = T2stardel(n); 
end 

tmp = char(slabFol); gname = string(tmp(4:end)); 
gname = replace(gname,"_"," "); 

varNames = cell(1,7); 
varNames{1} = 'Angular Region'; varNames{2} = 'Geodesic band a'; varNames{3} = 'Geodesic band b'; 
varNames{4} = 'Maximum dt/dm in T1'; varNames{5} = 'Maximum dt/dm in T2*'; 
varNames{6} = 'delT1 (ms)'; varNames{7} = 'delT2* (ms)'; 
T = array2table(result,'VariableNames',varNames); 

fprintf("\n\nFor brain sample %s: \n\n",gname); 
disp(T); 
end 

function [qMRIdat] = getqMRIROIij(qMRIType,slabFol,geoBand,angBand)

%INPUTS:
%1) qMRIType: Type of qMRI metric 
%T2star = Transverse relaxation time, T2* 
%T1 = Longitudinal relaxation time, T1 

%2) slabFols: String giving folder name of slab 

%4) geoBand: (number of) geodesic band 

%5) angBand: (number of) angular band 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
mainDir = '/export02/data/risa/02_MEGRE_and_MP2RAGE/'; 

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
end 

%extract QMRI voxelwise values in ROI defined by ith geodesic band and jth angular band 
qMRIdat = qMRI(angBandMap==angBand & geoBandMap==geoBand);
qMRIdat = rmmissing(qMRIdat); %remove missing entries 

if sentinel %if pia mater does not exists, remove outliers in ROI_ij with 1.5 MAD of median instead 
    qMRIdat = rmoutliers(qMRIdat,'median','ThresholdFactor',1.5); %remove outliers
end         

%if ROI_ij for slab comprises less than 100 voxels, omit it 
if (length(qMRIdat) < 50), qMRIdat = []; end 
% if (length(qMRIdat) < 10), qMRIdat = []; end %TEMPORARY: for deep GM structures, changed limit to 10 voxels 

end 