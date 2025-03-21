function [] = angularCCWmap(slabFol, x_angleSft, y_orient, z_orient, selSlices)
%INPUTS: 
%1) slabFol: name of slab 7 T MRI directory 
%2) x_angleSft: angular shift 'theta' of x-y axes in degrees (CCW positive) 
%3) y_orient: string giving orientation of y-axes where...
% - "M>>L" is medial>>lateral going +y axis to -y axis 
% - "L>>M" is lateral>>medial going +y axis to -y axis 
%4) z_orient: string giving orientation of z-axes where...
% - "R>>C" is rostral>>caudal going +z axis to -z axis 
% - "C>>R" is caudal>>rostal going +z axis to -z axis 
%5) selSlices: selected slice range along z axis 

%ASSUMPTIONS: 
% -assuming unselected (inviable) slices in slab has been zeroed out in 
% ventricular mask 'ventMask.nii' 

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
y_orient = upper(y_orient); z_orient = upper(z_orient); 
slabFol = char(slabFol); 

mainDir = '/export02/data/risa/02_MEGRE_and_MP2RAGE/'; 
outSegPath = [mainDir slabFol '/ang_geo_seg_output/']; 
VM = double(load_untouch_nii([outSegPath 'ventMask.nii']).img);
MEGRE = load_untouch_nii([mainDir slabFol '/MULTIECHOGRE_AVG.nii']); 
DimDat = MEGRE.hdr.dime.pixdim; Hist = MEGRE.hdr.hist; 
BM = double(load_untouch_nii([mainDir slabFol '/BrainMaskT2starW.nii']).img);

centMask = zeros(size(VM)); 
theta = zeros(size(VM)); dim = size(VM); x_len = int32(dim(1)); y_len = int32(dim(3)); 
for z = 1:dim(2)
    if (VM(:,z,:) == 0), continue; end 
    centroids = cat(1,regionprops(squeeze(VM(:,z,:))).Centroid);
    Cy = int32(centroids(:,1)); Cx = int32(centroids(:,2)); 
    centMask(Cx,z,Cy) = 1; 
    
    for x = 1:x_len
        for y = 1:y_len
            %Calculate according to 4-quadrant CCW angle grid 
            Ly = double(y-Cy); Lx = -1*double(x-Cx); A = 1; C = 0; 
            
            if (Ly >= 0 && Lx < 0) %2nd Quadrant
                A = -1; C = 180; 
            elseif (Ly < 0 && Lx <= 0) %3rd Quadrant
                A = 1; C = 180; 
            elseif (Ly < 0 && Lx > 0) %4th Quadrant
                A = -1; C = 360;  
            elseif (Ly == 0 && Lx == 0)
                continue; 
            end 
            
            %Accommodate for CCW four-quadrant angle calculation
            Ly = abs(Ly); Lx = abs(Lx); 
            theta(x,z,y) = A*atand(Ly/Lx)+C; 
        end 
    end 
end 

CM_nii = make_nii(centMask); 
CM_nii.hdr.dime.pixdim = DimDat; CM_nii.hdr.hist = Hist; 
save_nii(CM_nii,[outSegPath 'centMask.nii']); 

%Based on x_angleSft and y_orient, transform initially computed angle map
if (y_orient == "M>>L" && z_orient == "C>>R") || (y_orient == "L>>M" && z_orient == "R>>C")
    theta = -1.*(theta-360); %apply operation to entire angle map 
    x_angleSft = -1*x_angleSft; %negate x-axis angular shift 
end 

if x_angleSft > 0 %counter-clockwise
    tmp = theta - abs(x_angleSft); %subtract angular shift from entire angle map
    
    for x = 1:size(theta,1)
        for y = 1:size(theta,3)
            for z = 1:size(theta,2) 
                if (tmp(x,z,y)  >= -abs(x_angleSft) && tmp(x,z,y) <= 0)
                    theta(x,z,y) = tmp(x,z,y) + 360; 
                else 
                    theta(x,z,y) = tmp(x,z,y); 
                end  
            end 
        end 
    end 
     
elseif x_angleSft < 0 %clockwise 
    tmp = theta + abs(x_angleSft); %add angular shift to entire angle map 
    
    for x = 1:size(theta,1)
        for y = 1:size(theta,3)
            for z = 1:size(theta,2)
                if (tmp(x,z,y) >= 360 && tmp(x,z,y)  <= (360+abs(x_angleSft)))
                    theta(x,z,y) = tmp(x,z,y) - 360; 
                else 
                    theta(x,z,y) = tmp(x,z,y); 
                end 
            end 
        end 
    end 
    
end

theta(theta==0) = 360; %change 0 degrees to 360 degrees to avoid confusing with brain masking
tmp = zeros(size(theta)); tmp(:,(selSlices)+1,:) = theta(:,(selSlices)+1,:); theta = tmp; %zero out slices outside selected range

theta_nii = make_nii(theta); 
theta_nii.hdr.dime.pixdim = DimDat; theta_nii.hdr.hist = Hist; 
save_nii(theta_nii,[outSegPath 'angleCCWabtVent_unmasked.nii']); 

theta(BM==0) = 0; %zero out non-brain matter
theta_nii = make_nii(theta); 
theta_nii.hdr.dime.pixdim = DimDat; theta_nii.hdr.hist = Hist; 
save_nii(theta_nii,[outSegPath 'angleCCWabtVent.nii']); 
end 