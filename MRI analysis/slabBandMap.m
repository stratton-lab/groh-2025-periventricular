function [noBands] = slabBandMap(slabFol, lowThres, upThres, step, opt)

%INPUTS:
%1) slabFol: name of slab MRI directory 
%2) upThres: upper threshold value of geodesic (mm) or angular angle (deg) map 
%3) lowThres: lower threshold value of geodesic (mm) or angular angle (deg) map
%4) step: increment value for geodesic (mm) or angular bands (deg)
%5) opt: string that is 'geodesic' if segmenting geodesic map OR 'angular' if segmenting
%polar angle map 

%OUTPUT:
%1) noBands: number of concentric geodesic bands 

%ASSUMPTIONS:
% -assuming unselected (inviable) slices in geodesic distance map
% 'distance_map_masked.nii' and CCW angular map 'angleCCWabtVent.nii' have
% been zeroed out 

% - assuming geodesic distance map and CCW angular map have been masked
% with brain mask 'BrainMaskT2starW.nii'

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
mainDir = '/export02/data/risa/02_MEGRE_and_MP2RAGE/'; 
slabFol = char(slabFol); 
outSegPath = [mainDir slabFol '/ang_geo_seg_output/']; 

if (string(opt)== "geodesic")
    inpName = 'distance_map_masked.nii'; 
    outName = 'geodesic_band_map.nii'; 
elseif (string(opt) == "angular")
    inpName = 'angleCCWabtVent.nii'; 
    outName = 'angular_band_map.nii'; 
end 
maskedMap = load_untouch_nii([outSegPath inpName]);
DimDat = maskedMap.hdr.dime.pixdim; Hist = maskedMap.hdr.hist; maskedMap = double(maskedMap.img);

noBands = uint16(ceil((upThres-lowThres)/step)); 
geoBands = zeros(size(maskedMap)); 

%starting value for first concentric geodesic band 'ependyma-in' or angular band counter-clockwise 
startVal = 0.0; 
clrVal = 1; %color of selected band 
for i = 1:noBands
    geoBands(maskedMap>startVal & maskedMap<=(startVal+step)) = clrVal; 
    startVal = startVal + step; 
    clrVal = clrVal + 1; 
end 

GB_nii = make_nii(geoBands); 
GB_nii.hdr.dime.pixdim = DimDat; GB_nii.hdr.hist = Hist; 
save_nii(GB_nii,[outSegPath outName]); 

end 