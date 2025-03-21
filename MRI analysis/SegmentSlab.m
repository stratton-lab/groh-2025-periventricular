function [] = SegmentSlab(slabDir, selSlices, geoStep, angStep, x_angleSft, y_orient, z_orient)

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
slabDir = char(slabDir); 
mainDir = '/export02/data/risa/02_MEGRE_and_MP2RAGE'; 
outSegPath = [mainDir '/' slabDir '/ang_geo_seg_output/']; 
if ~isfolder(outSegPath), mkdir(outSegPath); end 

VM_path = [outSegPath 'ventMask.nii']; 
if ~isfile(VM_path)
    com = ['mv ' mainDir '/' slabDir '/ventMask.nii ' outSegPath]; [~,~] = system(com); 
end

DM_path = [outSegPath 'distance_map.nii'];
if ~isfile(DM_path)
    com1 = 'source /export02/data/conda/etc/profile.d/conda.sh';  
    com2 = 'conda activate falconenv'; 
    com3 = ['falcon_math distancemap -img=' VM_path ' -out=' DM_path ' -val=300 '];
    com4 = 'conda deactivate';
    com = [com1 ';' com2 ';' com3 ';' com4]; system(com); 
end

BM_path = [mainDir '/' slabDir '/BrainMaskT2starW.nii'];
DM_masked_path = [outSegPath 'distance_map_masked.nii']; 
if ~isfile(DM_masked_path) %if you want to generate new distance map (include more slices), you need to delete it manually
    DM = load_untouch_nii(DM_path); VM = load_untouch_nii(VM_path).img; BM = load_untouch_nii(BM_path).img; 
    DimDat = DM.hdr.dime.pixdim; Hist = DM.hdr.hist; DM = double(DM.img); 
    tmp = zeros(size(DM)); tmp(:,(selSlices)+1,:) = DM(:,(selSlices)+1,:); DM = tmp; 
    DM(BM==0) = 0; DM(VM~=0)=0; DM = -1*DM; 

    DM_nii = make_nii(DM); 
    DM_nii.hdr.dime.pixdim = DimDat; DM_nii.hdr.hist = Hist; 
    save_nii(DM_nii,DM_masked_path); 
end 

DM_masked = double(load_untouch_nii(DM_masked_path).img); 
upThresGeo = max(DM_masked,[],'all','omitnan'); 

slabBandMap(slabDir,0,upThresGeo,geoStep,'geodesic')
angularCCWmap(slabDir, x_angleSft, y_orient, z_orient, selSlices); 
slabBandMap(slabDir,0,360,angStep,'angular'); 

end

