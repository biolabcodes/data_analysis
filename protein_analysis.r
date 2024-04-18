library(ggplot2)
library(openxlsx)
library(ComplexHeatmap)
library(DEP)
library(clusterProfiler)

data<- read.xlsx('data/raw/data.xlsx')
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")


LFQ_columns <- grep("LFQ.", colnames(data_unique))
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)


LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)
data_filt <- filter_missval(data_se, thr = 1)
data_norm <- normalize_vsn(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff <- test_diff(data_imp, type = "control", control = "WT",test = 'KO')

enrich.go <- enrichGO(gene = Gene.names  
    OrgDb, 
    ont = 'ALL', 
    pAdjustMethod = 'fdr')
