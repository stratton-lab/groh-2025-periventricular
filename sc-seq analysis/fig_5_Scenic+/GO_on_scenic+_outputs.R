#====================================================
# Outline - GO Analysis on Scenic+ TF Target Genes
#====================================================
# 1. Loading Libraries and Setting Up
# 2. Defining Master Regulator Genes
# 3. Combining Gene Lists & Mapping to Entrez
# 4. GO Enrichment Analysis
# 5. Visualizing Enrichment Results

#====================================================
# 1. Loading Libraries and Setting Up
#====================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("E:/multiOmics/integrated object/multiome/GO_Scenicplus")
set.seed(1)

#====================================================
# 2. Defining Master Regulator Genes
#====================================================

###################### Fig 5D ######################


# ---- STAT3 ----
STAT3_genes <- c(
  'SAMD4A','CD44','CSGALNACT2','RGS20','YBX3','COLEC12','RBM7','RASGEF1A',
  'ARHGEF3','DGKD','ZFYVE28','NTRK3','UGCG','BAALC','C11orf71'
)

# ---- BACH2 ----
BACH2_genes <- c(
  'CD44','NTRK3','ARHGEF3','TNC','CLIC4','SAMD1','KCNJ9','SAMD4A','GBP2','MAPRE2',
  'SLC1A2','LRIG1','UBTD2','SMARCC2','FRMD5','DGKD','PHLDA1','NTNG1','CSGALNACT2',
  'M1AP','LTBR','USP12','NTM','ARID5A','HS6ST1','KCNJ10','MAN1C1','ECE1','FLNC',
  'TAF4B','NFASC','PRKCQ','PNPLA2','WDR48','GNA12','DBI','TJP2','KCNH1','ZNF708',
  'PFKP','CNTN2','DLGAP1','MRPS5','SUOX','ABCA3','PDE4B','NEK6','MAFK','ANKRD23',
  'TRIM9','CSRNP1'
)

# ---- ZBTB20 ----
ZBTB20_genes <- c(
  'NCAM2','PKN2','PLCE1','SOX5','PSPC1','TJP1','ADGRB3','ASB3','CCDC88A','UPF3A',
  'MAPK1','KIDINS220','SIK3','PDE8A','DENND5B','RAP1GDS1','MAP4K5','NEK1','MBD2',
  'EDA','ATF6','ARMC2','IFT81','UBL3','ZNF135','AOPEP','ITPKB','MBOAT2','SGK3',
  'SUPT20H','CYTH1','ZNF592','MOK','ANKIB1','ZSCAN30','STK33','ANKRD26','ACYP2',
  'CSNK1G1','MAPRE2','IGF1R','TRIM66','TRDMT1','WASHC3','CCDC144A','SECISBP2L',
  'RHEB','DMTF1','TASOR2','PTRH2','CYSTM1','C1orf87','RNF217','DNM3','NIN','ITPR1',
  'SMIM8','MIA3','PARG','IFT74','TMEM165','CC2D2A','INO80D','ANKRD30BL','FGGY',
  'TXNDC16','GNAO1','ADHFE1','KRAS','LHFPL6','TMX3','FSTL3','LRRC23','NPC1',
  'CSGALNACT1','FCER1G','DPYD','WDR35','NAB1','INIP','SCFD2','ZNF624','PLEKHH1',
  'LMAN1','ARHGAP31','AMOTL1','ENDOD1','SCFD1','ZC3H13','SNX6','SCLT1','MOB1A',
  'RSPH4A','SAXO2','SLC10A7','ZNF345','DLD','SMC6','SPATA17','CYP51A1','RHOBTB3',
  'NAPEPLD','ZCCHC2','INTS3','PTPDC1','PCDH10','GK5','XIAP','FRYL','CUBN','OSBPL6',
  'KLHL2','ZDHHC20','TPM1','FGF14','LRP11','APBA2','LIPA','SPATS2L','SPAG6',
  'TMEM33','CLK1','PAQR3','HNRNPF','ALPK3','IQUB','ZNF487','GPR39','STX17','WDR27',
  'NDUFA5','IKBKB','APH1B','FBXL13','CCDC180','TMEM68','BZW1','SMAD1'
)

# ---- ZNF704 ----
ZNF704_genes <- c(
  'PLLP','RNF207','NFKB1','KIF13B','SNX6','PVALEF','OTUD7A','RMDN1','SGK3','SEMA5A',
  'AMZ1','ZNF652','CTNNA3','APLP1','SLC5A11','LDB3','CDC14B','TMTC2','DNM3',
  'SECISBP2L','NLK','TRPC4AP','MITF','CCDC88A','PLCL1','TMEM98','ZDHHC20','TUBG2',
  'KIZ','FAM171A1','C12orf76','PPP1R16B','ANAPC4','MAP4K5','CAMK2N1','SLC44A1',
  'SORT1','TMEM144','FAF1','WASHC3','PXK','AGTPBP1','HIPK2','FZD5','FDFT1',
  'FAM171B','SEPTIN7','PPM1N','MSMO1','MAP7','SGK1','KAT2B','ZNF536','RUNX2',
  'SEMA4D','PRKACB','ACSS2','PTPRD','UBR2','ECPAS','PTPDC1','ANKRD13A','ENDOD1',
  'GPR137C','R3HDM4','LIPA','INPP1','TTLL7','LPGAT1','JAKMIP3','KIRREL3','LGR5',
  'SHTN1','GARNL3','FGFR2','DOCK10','TRIM59','CLDN11','PEX5L','GPSM2','KLHL2',
  'CWC27','LSM10','ADGRB3','CERCAM','ANAPC5','ACTN2','GAB2','LRP2','ERI3','NHEJ1',
  'HS3ST5'
)

# ---- NFIB ----
NFIB_genes <- c(
  'ROPN1L','ARMH1','MROH1','CLYBL','CCDC6','PLEKHA7','NRG3','AOPEP','ZHX3','CFAP61',
  'PACRG','CARMIL1','CLIC6','KANSL1L','IQCK','NPAS3','BICC1','LPIN1','TRAF3',
  'SLC20A2','AGBL4','RNF169','DRC1','SORBS1','USP43','PRICKLE2','FGF2','LIMCH1',
  'SPAG9','NAV2','CCDC13','GLIS3','ICE1','PARD3','STXBP5L','ATP6V1D','NFIA',
  'C11orf58','DAZAP2','TVP23B','CEP126','DDX5','LRGUK','C21orf62','GPM6A','HSPB1',
  'SRGAP3','STK33','RAB10','KITLG','MYO9A','SLC4A8','LMNTD1','SLC30A1','SRP9',
  'CELF6','CEBPG','KDM3A','LRRC8A','SUSD6','GPC6','CFAP52','PSAP','ADCY2','TSC22D2',
  'DDB1','DLEC1','NMT1','LCA5L','GOLGA2','LRRN2','AZIN1','AMPD3','MINDY4','SLC6A16',
  'ALDH7A1','CYSTM1','FGGY','DST','GPR176','RNASE2','LUC7L3','RPGRIP1L','KIAA2012',
  'TTLL5','PATZ1','NCOR1','GET1','LDLRAD4','SLC1A2','TCF7L2','CDH11','EIF4G3',
  'RSPH14','CTNND2','RCN2','LTBP3','BAIAP2','RPL18','COLEC12','ZBTB47','WDR35',
  'COQ4','PPP3CA','AHR','PBX1','RAB18','BMP2K','SOX13','MYLK3','RGS3','PLCG1','C1S',
  'PASK'
)

#====================================================
# 3. Combining Gene Lists & Mapping to Entrez
#====================================================
gene_list <- unique(c(
  STAT3_genes, BACH2_genes, ZBTB20_genes, 
  ZNF704_genes, NFIB_genes
))

# Convert gene symbols to Entrez IDs
gene_list_go <- mapIds(
  org.Hs.eg.db,
  keys      = gene_list,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
gene_list_go <- na.omit(gene_list_go)

#====================================================
# 4. GO Enrichment Analysis
#====================================================
scenic_plus_TF_go_results <- enrichGO(
  gene          = gene_list_go,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

scenic_plus_TF_go_results <- pairwise_termsim(scenic_plus_TF_go_results)

#====================================================
# 5. Visualizing Enrichment Results
#====================================================
# Tree plot
p2 <- treeplot(scenic_plus_TF_go_results)

# cnetplot
p3 <- cnetplot(
  scenic_plus_TF_go_results,
  showCategory = 5,
  circular     = TRUE,
  colorEdge    = TRUE
)
