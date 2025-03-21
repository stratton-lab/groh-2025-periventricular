function [] = plotGradDataSlicewise(slabFol,slice1st, sliceLst, offset, angBandNo, varargin)
%What we want function to be: 

%INPUTS: 
%1) slabFol: string that gives name of slab 7 T MRI directory 
%2) slice1st: first permissible slice (corresponding to silver surface
%side), since Adam did tissue sectioning upwards from silver surface side
%-assuming index starts at 0 as in fsleyes
%3) sliceLst: last permissible slice 
%4) angBandNo: array containing indices of angular bands to consider or
%numeric value which is angular band #

%OPTIONAL ARGUMENTS:
%1) varargin{1}: array giving numbers of geodesic bands to consider 

%Since matlab starts indexing at 1
slice1st = slice1st+1; sliceLst = sliceLst + 1; offset = offset + 1; 

addpath(genpath("/export02/data/risa/NIfTI_20140122/"));
mainDir = '/export02/data/risa/02_MEGRE_and_MP2RAGE/'; slabFol = char(slabFol); 

outFigPath = [mainDir slabFol '/slicewiseFigures/']; 
if ~isfolder(outFigPath), mkdir(outFigPath); end 

%Load tissue-specific segmentation data 
Seg = load_untouch_nii([mainDir slabFol '/Segmentation.nii']).img; 
%Load geodesic and angular band maps 
geoBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/geodesic_band_map.nii']).img; 
angBandMap = load_untouch_nii([mainDir slabFol '/ang_geo_seg_output/angular_band_map.nii']).img;

if nargin < 6 % 5 (official) arguments passed-in
    BandsGeo = 2:uint16(max(geoBandMap,[],'all','omitnan')); 
elseif nargin < 7 % 6 arguments passed-in (1 extra)
    BandsGeo = varargin{1}; 
end 
noBandsGeo = length(BandsGeo); 

names = ["/T2star_uncorr_AVG.nii" "/T1MAP_AVG.nii"]; 
ylbls = ["Mean T_{2}^{*} (ms)" "Mean T_{1} (ms)"]; 
ylimits = {[10 60], [100 500]}; 
qMRI = cell(1,2); 
qMRI{1} = double(load_untouch_nii([mainDir slabFol char(names(1))]).img); 
qMRI{1}(geoBandMap==0) = 0; qMRI{1}(Seg~=1) = NaN; 
qMRI{2} = double(load_untouch_nii([mainDir slabFol char(names(2))]).img); 
qMRI{2}(geoBandMap==0) = 0; qMRI{2}(Seg~=1) = NaN; 

if length(angBandNo) > 1
    qMRI{1}(~ismember(angBandMap,angBandNo)) = NaN;
    qMRI{2}(~ismember(angBandMap,angBandNo)) = NaN; 
    tlbl = sprintf("Mean WM T2* & T1 over geodesic bands for Angular Bands #%d to #%d",angBandNo(1),angBandNo(end)); 
else
    qMRI{1}(angBandMap~=angBandNo) = NaN;
    qMRI{2}(angBandMap~=angBandNo) = NaN;
    tlbl = sprintf("Mean WM T2* & T1 over geodesic bands for Angular Band #%d",angBandNo);
end 

% %TEMPORARY
% qMRI{1}(ismember(geoBandMap,[6 7 8])) = NaN; 
% qMRI{2}(ismember(geoBandMap,[6 7 8])) = NaN; 

%if pia mask exists, zero out voxels corresponding to pia mater 
piaMPath = [mainDir slabFol '/piaMask.nii']; sentinel = isfile(piaMPath); 
if sentinel
    piaMask = double(load_untouch_nii(piaMPath).img); 
    qMRI{1}(piaMask~=0) = NaN; qMRI{2}(piaMask~=0) = NaN; 
end 

if slice1st > sliceLst
    step = -1; 
    sliceNo = offset - slice1st; 
elseif slice1st < sliceLst
    step = 1; 
    sliceNo = slice1st - offset;
end 

tileNo = 0; nFig = 0; 
for nSlice = slice1st:step:sliceLst %going through 1st permissble slice to last permissible slice 
    meanGeo = zeros(2,noBandsGeo); stdGeo = meanGeo;
    for q = 1:2 %T1 or T2* metric 
        for IdxGeo = 1:noBandsGeo
            %Since we are doing single slice and thus have limited voxels for ROI_ij => OMITTING DENOISING PROCESS
            tmp = qMRI{q}; %extract T1 or T2* map data 
            tmp(geoBandMap~=BandsGeo(IdxGeo)) = NaN; %select current geodesic band
            
%             %TESTING
%             [~,rmidxs] = rmoutliers(tmp(:),'median','ThresholdFactor',1.5); %remove outliers within 1.5 MAD from median
%             tmp(rmidxs) = NaN; 
            
            tmp = tmp(:,nSlice,:); %select current slice   
            meanGeo(q,IdxGeo) = mean(tmp,'all','omitnan'); stdGeo(q,IdxGeo) = std(tmp,0,'all','omitnan'); 
            if length(tmp) < 10, meanGeo(q,IdxGeo) = NaN; stdGeo(q,IdxGeo) = NaN; end 
        end 
    end 
    plotTle = sprintf("Slice #%d",sliceNo);
    xlbl = "Geodesic band"; 
    
    if tileNo == 0
        ax = figure("WindowState","maximized");
        t = tiledlayout(3,3); 
        title(t,tlbl);        
        nexttile;
    elseif mod(tileNo,9) == 0
        %Save figure
        nFig = nFig + 1; 
        if length(angBandNo) > 1 
            figName = char(sprintf("/Fig_angBand_%dto%d_set_%d.jpg",angBandNo(1),angBandNo(end),nFig)); 
        else 
            figName = char(sprintf("/Fig_angBand_%d_set_%d.jpg",angBandNo,nFig));
        end 
        exportgraphics(ax, [outFigPath figName],'Resolution',300); 
    
        ax = figure("WindowState","maximized");
        t = tiledlayout(3,3); 
        title(t,tlbl);        
        nexttile;
    else 
        nexttile; 
    end 

    yyaxis left; 
    errorbar(BandsGeo,meanGeo(1,:),stdGeo(1,:),'LineWidth',1);  
    ylim(ylimits{1}); xlabel(xlbl); ylabel(ylbls(1)); 

    yyaxis right; 
    errorbar(BandsGeo,meanGeo(2,:),stdGeo(2,:),'LineWidth',1);  
    ylim(ylimits{2}); xlabel(xlbl); ylabel(ylbls(2)); 
    xticks(BandsGeo); title(plotTle); 

    sliceNo = sliceNo+1; tileNo = tileNo + 1; 
    
    if nSlice==sliceLst
        nFig = nFig + 1; 
        if length(angBandNo) > 1 
            figName = char(sprintf("/Fig_angBand_%dto%d_set_%d.jpg",angBandNo(1),angBandNo(end),nFig)); 
        else 
            figName = char(sprintf("/Fig_angBand_%d_set_%d.jpg",angBandNo,nFig));
        end 
        exportgraphics(ax, [outFigPath figName],'Resolution',300); 
    end 
end

end 