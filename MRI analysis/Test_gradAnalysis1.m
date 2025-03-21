%%
mainDir = '/export02/data/risa/02_Ex-Vivo_MS_Slab_Scans/01_MEGRE_and_MP2RAGE/14_1143_Caudal_Cl/';
LF = load_untouch_nii([mainDir '/ang_geo_seg_output/angular_band_map.nii']); 
DimDat = LF.hdr.dime.pixdim; Hist = LF.hdr.hist; LF = double(LF.img);
LF(:,40,:) = 0; 
LF_nii = make_nii(LF); 
LF_nii.hdr.dime.pixdim = DimDat; LF_nii.hdr.hist = Hist; 
save_nii(LF_nii,[mainDir '/ang_geo_seg_output/angular_band_map.nii']); 

%%
mainDir = '/export02/data/risa/02_Ex-Vivo_MS_Slab_Scans/01_MEGRE_and_MP2RAGE/05_1002_Frontal_WithVentricle_MS/'; 

LF = load_untouch_nii([mainDir 'laplacian_field.nii']); 
BM = load_untouch_nii([mainDir 'outCortMask.nii']).img; 
DimDat = LF.hdr.dime.pixdim; Hist = LF.hdr.hist; LF = double(LF.img); 
Inds = find(BM~=0); %find field voxel indices 

tmp = double(LF(Inds)); tmp = tmp(tmp~=1000); 

tmp = double(zeros(size(LF))); 
tmp(Inds) = -1*(LF(Inds)-(1000));
logLF = double(zeros(size(LF))); 
logLF(Inds) = -reallog(tmp(Inds)); 
logLF(isnan(logLF)) = 0; logLF(isinf(logLF)) = 0; 
% logLF = LF; 

%Normalize field map to between (0,1]
minlogLF = logLF(Inds); minlogLF = min(minlogLF(:)); 
maxlogLF = logLF(Inds); maxlogLF = max(maxlogLF(:)); 
logLF(Inds) = (logLF(Inds) - minlogLF)./(maxlogLF - minlogLF); 

LF_nii = make_nii(logLF); 
LF_nii.hdr.dime.pixdim = DimDat; LF_nii.hdr.hist = Hist; 
save_nii(LF_nii,[mainDir 'laplacian_field_log.nii']); 

%% GENERATE GEODESIC AND ANGULAR BAND MAPS
%ignoring slabs 1002 Frontal (no ventricle), 1090 Caudal Small, 1421
%Frontal Small, 1736 Frontal Small, 1421 Middle, 989 Middle 

ScansMS = ["01_1011_Caudal_MS", "02_1011_Middle_MS", "03_1011_Frontal_MS",...
    "05_1002_Frontal_WithVentricle_MS", "06_1002_Middle_MS", "07_1002_Caudal_MS", "08_1090_Frontal_MS",...
    "09_1090_Middle_MS", "10_1090_Caudal_Large_MS", ...
    "16_1421_Frontal_Large_MS", "18_1421_Caudal_MS",...
    "26_1736_Frontal_Large_MS", "27_1736_Middle_MS", "28_1736_Caudal_MS"]; 
IdxsMS = {24:34, 24:32, 26:38, 23:36, 22:32, 27:34, 21:41, 15:43, 17:40,...
    16:40, 17:41, 17:42, 28:32, 20:33}; 
x_angleSftMS = [-9 0 180 0 16 0 0 -4 0 -19 0 0 0 27]; 
y_orientMS = ["L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L"]; 
z_orientMS = ["C>>R" "C>>R" "R>>C" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R"]; 

ScansCl = ["12_1143_Frontal_Cl", "13_1143_Middle_Cl", "14_1143_Caudal_Cl", "19_1670_Frontal_Cl",...
    "20_1670_Middle_Cl", "21_1670_Caudal_Cl", "22_1722_Frontal_Cl", "23_1722_Middle_Cl",...
    "24_1722_Caudal_Cl", "32_989_Frontal_Cl", "33_989_Caudal_Cl"];

IdxsCl = {28:34, 25:35, 20:39, 23:35, 25:31, 23:35, 18:43, 20:37, 24:38, 27:38, 22:29};
x_angleSftCl = [-3.5 -5 0 15 -4.0 -23 0 -30.5 -3 5.5 -15.5]; 
y_orientCl = ["L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "M>>L" "M>>L" "M>>L" "L>>M" "L>>M"]; 
z_orientCl = ["C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R"];

parfor i = 1:length(ScansMS)
    SegmentSlab(ScansMS(i), IdxsMS{i}, 1.5, 15, x_angleSftMS(i), y_orientMS(i), z_orientMS(i));
end 
delete(gcp('nocreate'));

parfor i = 1:length(ScansCl)
    SegmentSlab(ScansCl(i), IdxsCl{i}, 1.5, 15, x_angleSftCl(i), y_orientCl(i), z_orientCl(i));
end 
delete(gcp('nocreate'));
SegmentSlab(ScansMS(1), IdxsMS{1}, 5, 15, x_angleSftMS(1), y_orientMS(1), z_orientMS(1));
%SegmentSlab("34_989_Middle_Cl", 27:30, 1.5, 15, 0, "L>>M", "C>>R"); 
SegmentSlab("17_1421_Middle_MS",17:43,1.2,15,0,"M>>L","C>>R"); 
SegmentSlab(ScansCl(5), IdxsCl{5}, 1.5, 15, x_angleSftCl(5), y_orientCl(5), z_orientCl(5));

%% LOAD VARIABLES FOR STATISTICAL ANALYSIS 
MS_Caudal = ["01_1011_Caudal_MS", "07_1002_Caudal_MS","10_1090_Caudal_Large_MS", "18_1421_Caudal_MS",...
    "28_1736_Caudal_MS"]; 
HC_Caudal = ["14_1143_Caudal_Cl","21_1670_Caudal_Cl","24_1722_Caudal_Cl","33_989_Caudal_Cl"];
HC_AB_Caudal = {[1:12,20:24], [1:14,19:24], [1:13,19:24], [1:13,20:24]}; 
HC_GBmax_Caudal = {36, 27, 28, 26}; 
PMI_MS_Caudal = [63.0 11.25 57.75 22.50 14.05]; 
PMI_HC_Caudal = [20.75 26.25 17.67 14.75]; 
Age_HC_Caudal = [64 51 59 79]; 
gnames_MS_Caudal = ["1011 Caudal", "1002 Caudal", "1090 Caudal Large", "1421 Caudal", "1736 Caudal"]; 
gnames_HC_Caudal = ["1143 Caudal", "1670 Caudal", "1722 Caudal", "989 Caudal"]; 

% MS Slab 1421 (Middle) and HC Slab 989 (Middle) have been omitted 
MS_Middle = ["02_1011_Middle_MS", "06_1002_Middle_MS", "09_1090_Middle_MS", "27_1736_Middle_MS"]; 
HC_Middle = ["13_1143_Middle_Cl", "20_1670_Middle_Cl","23_1722_Middle_Cl"]; 
HC_AB_Middle = {[1:7,10:11,21:24], [1:6,23:24], [1:8,15:24]}; 
HC_GBmax_Middle = {36, 26, 27}; 
PMI_MS_Middle = [63.0 11.25 57.75 14.05]; 
PMI_HC_Middle = [20.75 26.25 17.67]; 
Age_HC_Middle = [64 51 59]; 
gnames_MS_Middle = ["1011 Middle", "1002 Middle", "1090 Middle", "1736 Middle"]; 
gnames_HC_Middle = ["1143 Middle", "1670 Middle", "1722 Middle"]; 
    
MS_Frontal = ["03_1011_Frontal_MS", "05_1002_Frontal_WithVentricle_MS", "08_1090_Frontal_MS",...
    "16_1421_Frontal_Large_MS", "26_1736_Frontal_Large_MS"]; 
HC_Frontal = ["12_1143_Frontal_Cl",  "19_1670_Frontal_Cl", "22_1722_Frontal_Cl",...
     "32_989_Frontal_Cl"];  
HC_AB_Frontal = {[1:24], [1:14,23:24], [1:24], [1:14,23:24]}; 
HC_GBmax_Frontal = {24, 22, 23, 19}; 
PMI_MS_Frontal = [63.0 11.25 57.75 22.50 14.05]; 
PMI_HC_Frontal = [20.75 26.25 17.67 14.75]; 
Age_HC_Frontal = [64 51 59 79]; 
gnames_MS_Frontal = ["1011 Frontal", "1002 Frontal", "1090 Frontal", "1421 Frontal Large", "1736 Frontal Large"]; 
gnames_HC_Frontal = ["1143 Frontal", "1670 Frontal", "1722 Frontal", "989 Frontal"]; 

%% STATS ON NORMALIZED RELAXOMETRY DATA for angular segments seen as 'exhibiting pathology'
% patholData = [12 2 10; 21 2 7; 22 2 8; 23 2 13; 24 2 19]; 
% MS_AB = {[1:12,20:24]}; MS_GBmax = {24};
% T = getStatsTableMSslab("01_1011_Caudal_MS",HC_Caudal,...
% [MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal],patholData);

% patholData = [5 7 20; 6 2 5; 7 2 6; 8 2 5]; 
% MS_AB = {[1:24]}; MS_GBmax = {24};
% T = getStatsTableMSslab("03_1011_Frontal_MS",HC_Frontal,...
%     [MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal],patholData); 

% patholData = [2 2 9; 3 2 12; 4 4 14; 5 5 11; 11 2 7; 12 2 7; 20 2 5]; 
% MS_AB = {[1:24]}; MS_GBmax = {23}; 
% T = getStatsTableMSslab("05_1002_Frontal_WithVentricle_MS",HC_Frontal,...
%     [MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal],patholData);

% patholData = [6 2 19; 7 2 22; 8 2 19]; 
% MS_AB = {[1:9,21:24]}; MS_GBmax = {26};
% T = getStatsTableMSslab("06_1002_Middle_MS",HC_Middle,...
%     [MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle],patholData);

% patholData = [1 2 24; 2 2 24;3 2 24; 4 2 11; 5 2 7; 6 4 12; 7 4 16; 8 3 8; 9 3 6; 10 2 6; 11 2 9; 12 2 8]; 
% MS_AB = {[1:12,22:24]}; MS_GBmax = {25};
% T = getStatsTableMSslab("07_1002_Caudal_MS",HC_Caudal,...
%     [MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal],patholData);

% patholData = [3 4 13; 4 10 18; 7 2 14; 8 2 20; 9 2 18; 10 2 16;...
%     11 2 4; 11 4 10; 12 2 11; 13 2 12; 14 2 13; 15 2 9; 18 2 5;...
%     20 2 5; 21 2 5; 22 2 5]; 
% MS_AB = {[1:24]}; MS_GBmax = {20}; 
% T = getStatsTableMSslab("08_1090_Frontal_MS",HC_Frontal,...
%     [MS_AB, HC_AB_Frontal],[MS_GBmax, HC_GBmax_Frontal],patholData);

% patholData = [4 11 18; 5 11 19; 6 14 20; 7 2 5; 7 5 18; 8 3 6;...
%     8 6 20; 9 2 6; 9 6 22; 10 2 5; 10 5 23; 11 3 7; 23 2 6; 24 2 5]; 
% MS_AB = {[1:11,23:24]}; MS_GBmax = {24}; 
% T = getStatsTableMSslab("09_1090_Middle_MS",HC_Middle,...
%     [MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle],patholData);

% patholData = [2 2 20; 3 2 22; 4 2 21; 5 9 21; 6 12 19; 7 2 18; 8 2 15;...
%     9 2 17; 10 2 12; 23 2 12]; 
% MS_AB = {[1:24]}; MS_GBmax = {23};
% T = getStatsTableMSslab("10_1090_Caudal_Large_MS",HC_Caudal,...
%     [MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal],patholData);
    
% patholData = [1 2 13; 2 2 14; 3 2 14; 4 2 14; 5 2 13; 6 2 18; 7 2 18;...
%     8 2 4; 9 2 4; 10 2 4; 11 2 4; 12 2 7; 13 2 8; 24 4 10]; 
% MS_AB = {[1:24]}; MS_GBmax = {21};  
% T = getStatsTableMSslab("26_1736_Frontal_Large_MS",HC_Frontal,...
%     [MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal],patholData); 

% patholData = [8 5 23; 9 9 20]; 
% MS_AB = {[1:9,22:24]}; MS_GBmax = {25};
% T = getStatsTableMSslab("27_1736_Middle_MS",HC_Middle,...
%     [MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle],patholData);

patholData = [1 3 22; 2 3 21; 3 3 21; 5 4 21; 8 4 10; 9 3 7; 10 2 7; 11 2 7]; 
MS_AB = {[1:3,5,8:11]}; MS_GBmax = {22};
T = getStatsTableMSslab("28_1736_Caudal_MS",HC_Caudal,...
    [MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal],patholData);

%% RAW/UNNORMALIZED AND NORMALIZED RELAXOMETRY DATA PLOTTING     
MS_AB = {[1:12,20:24]}; MS_GBmax = {24};
plotDiffqMRIData("01_1011_Caudal_MS",HC_Caudal,[MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal]);
plotGradData("01_1011_Caudal_MS", "WM", "T2star&T1", 2:24, [1:12 20:24]); 

plotGradData("02_1011_Middle_MS", "WM", "T2star&T1", 2:32, [1:11 24]); 

MS_AB = {[1:24]}; MS_GBmax = {24};
plotDiffqMRIData("03_1011_Frontal_MS",HC_Frontal,[MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal]);
plotGradData("03_1011_Frontal_MS", "WM", "T2star&T1", 2:24, [1:24]); 

MS_AB = {[1:24]}; MS_GBmax = {23}; 
plotDiffqMRIData("05_1002_Frontal_WithVentricle_MS",HC_Frontal,[MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal]); 
plotGradData("05_1002_Frontal_WithVentricle_MS", "WM", "T2star&T1", 2:23, [1:24]); 

MS_AB = {[1:9,21:24]}; MS_GBmax = {26};
plotDiffqMRIData("06_1002_Middle_MS",HC_Middle,[MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle]);
plotGradData("06_1002_Middle_MS", "WM", "T2star&T1", 2:26, [1:9 21:24]); 

MS_AB = {[1:12,22:24]}; MS_GBmax = {25};
plotDiffqMRIData("07_1002_Caudal_MS",HC_Caudal,[MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal]);
plotGradData("07_1002_Caudal_MS", "WM", "T2star&T1", 2:25, [1:12 22:24]);

MS_AB = {[1:24]}; MS_GBmax = {20}; 
plotDiffqMRIData("08_1090_Frontal_MS",HC_Frontal, [MS_AB, HC_AB_Frontal],[MS_GBmax, HC_GBmax_Frontal]);
plotGradData("08_1090_Frontal_MS", "WM", "T2star&T1", 2:20, [1:24]);

MS_AB = {[1:11,23:24]}; MS_GBmax = {24}; 
plotDiffqMRIData("09_1090_Middle_MS",HC_Middle,[MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle]);  
plotGradData("09_1090_Middle_MS", "WM", "T2star&T1", 2:24, [1:11 23:24]);

MS_AB = {[1:24]}; MS_GBmax = {23};
plotDiffqMRIData("10_1090_Caudal_Large_MS",HC_Caudal,[MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal]);
plotGradData("10_1090_Caudal_Large_MS", "WM", "T2star&T1", 2:23, [1:24]);

MS_AB = {[1:12]}; MS_GBmax = {27};
plotDiffqMRIData("16_1421_Frontal_Large_MS",HC_Frontal,[MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal]);
plotGradData("16_1421_Frontal_Large_MS", "WM", "T2star&T1", 2:27, [1:12]);

MS_AB = {[1:4,10:11]}; MS_GBmax = {30};
plotDiffqMRIData("17_1421_Middle_MS",HC_Middle,[MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle]);
plotGradData("17_1421_Middle_MS", "WM", "T2star&T1", 2:30, [1:4 10:11]);

MS_AB = {[1:24]}; MS_GBmax = {33};
plotDiffqMRIData("18_1421_Caudal_MS",HC_Caudal,[MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal]);
plotGradData("18_1421_Caudal_MS", "WM", "T2star&T1", 2:33, [1:24]);

MS_AB = {[1:24]}; MS_GBmax = {21};  
plotDiffqMRIData("26_1736_Frontal_Large_MS",HC_Frontal,[MS_AB,HC_AB_Frontal],[MS_GBmax,HC_GBmax_Frontal]); 
plotGradData("26_1736_Frontal_Large_MS", "WM", "T2star&T1", 2:21, [1:24]);

MS_AB = {[1:9,22:24]}; MS_GBmax = {25};
plotDiffqMRIData("27_1736_Middle_MS",HC_Middle,[MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle]);
plotGradData("27_1736_Middle_MS", "WM", "T2star&T1", 2:25, [1:9 22:24]);

MS_AB = {[1:3,5,8:11]}; MS_GBmax = {22};
plotDiffqMRIData("28_1736_Caudal_MS",HC_Caudal,[MS_AB,HC_AB_Caudal],[MS_GBmax,HC_GBmax_Caudal]);
plotGradData("28_1736_Caudal_MS", "WM", "T2star&T1", 2:22, [1:3 5 8:11]);

plotGradData("12_1143_Frontal_Cl", "WM", "T2star&T1", 2:24, [1:24]);
plotGradData("13_1143_Middle_Cl", "WM", "T2star&T1", 2:36, [1:7 10:11 21:23]);
plotGradData("14_1143_Caudal_Cl", "WM", "T2star&T1", 2:36, [1:12 20:24]);
plotGradData("19_1670_Frontal_Cl", "WM", "T2star&T1", 2:22, [1:14 23:24]);
plotGradData("20_1670_Middle_Cl", "WM", "T2star&T1", 2:26, [1:6 23:24]);
plotGradData("21_1670_Caudal_Cl", "WM", "T2star&T1", 2:27, [1:14 19:24]);
plotGradData("22_1722_Frontal_Cl", "WM", "T2star&T1", 2:23, [1:24]);
plotGradData("23_1722_Middle_Cl", "WM", "T2star&T1", 2:27, [1:8 15:24]);
plotGradData("24_1722_Caudal_Cl", "WM", "T2star&T1", 2:28, [1:13 19:24]);
plotGradData("32_989_Frontal_Cl", "WM", "T2star&T1", 2:19, [1:14 23:24]);
plotGradData("33_989_Caudal_Cl", "WM", "T2star&T1", 2:26, [1:13 20:24]);

%%
MS_AB = {[1:11,23:24]}; MS_GBmax = {24}; 
plotPrctDiff("09_1090_Middle_MS",HC_Middle,[MS_AB,HC_AB_Middle],[MS_GBmax,HC_GBmax_Middle]); 

%%
patholData = [22 2 9; 23 2 13; 24 2 25]; 
slabFol = "33_989_Caudal_Cl"; 
tbl = getStatsTableanySlab(slabFol,patholData,1.50); 

%% GENERATE MAP WITH BOUNDARIES OF SELECTED ANGULAR REGIONS MARKED
%ignoring slabs 1002 Frontal (no ventricle), 1090 Caudal Small, 1421
%Frontal Small, 1736 Frontal Small, 1421 Middle, 989 Middle 

ScansMS = ["01_1011_Caudal_MS", "02_1011_Middle_MS", "03_1011_Frontal_MS",...
    "05_1002_Frontal_WithVentricle_MS", "06_1002_Middle_MS", "07_1002_Caudal_MS", "08_1090_Frontal_MS",...
    "09_1090_Middle_MS", "10_1090_Caudal_Large_MS", ...
    "16_1421_Frontal_Large_MS", "18_1421_Caudal_MS",...
    "26_1736_Frontal_Large_MS", "27_1736_Middle_MS", "28_1736_Caudal_MS"]; 
IdxsMS = {24:34, 25:32, 26:38, 23:36, 22:32, 27:34, 22:41, 15:43, 17:40,...
    16:40, 17:41, 17:42, 29:32, 20:33}; 
x_angleSftMS = [-9 0 180 0 16 0 0 -4 0 -19 0 0 0 27]; 
y_orientMS = ["L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L" "M>>L"]; 
z_orientMS = ["C>>R" "C>>R" "R>>C" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R"]; 

ScansCl = ["12_1143_Frontal_Cl", "13_1143_Middle_Cl", "14_1143_Caudal_Cl", "19_1670_Frontal_Cl",...
    "20_1670_Middle_Cl", "21_1670_Caudal_Cl", "22_1722_Frontal_Cl", "23_1722_Middle_Cl",...
    "24_1722_Caudal_Cl", "32_989_Frontal_Cl", "33_989_Caudal_Cl"];

IdxsCl = {28:34, 25:33, 20:39, 23:35, 26:31, 23:35, 18:43, 20:37, 24:38, 27:38, 22:29};
x_angleSftCl = [-3.5 -5 0 15 -4.0 -23 0 -30.5 -3 5.5 -15.5]; 
y_orientCl = ["L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "L>>M" "M>>L" "M>>L" "M>>L" "L>>M" "L>>M"]; 
z_orientCl = ["C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R" "C>>R"];

markAngSegBorders(ScansMS(7), 0, 360, 15, [22:24], [0.5 0.5 0.5]); 



