library(dplyr)
library(data.table)
library(harmony)
library(miQC)
library(Seurat)
library(SingleR)
library(Seurat)
library(ggrepel)
library(glue)
library(SCP)
library(openxlsx)
library(CellChat)
library(ggplot2)


# load seu
sample_path <- './data/'  # the sample path
samples<- list.files(sample_path)  # sample ID in the path

seu<- Read10X(glue::glue(sample_path,'/filtered_feature_bc_matrix'))
    seu<- CreateSeuratObject(seu)
    seu$percent.mito <- PercentageFeatureSet(seu,pattern = mito_pattern)
    seu$percent.ribo <- PercentageFeatureSet(seu,pattern = ribo_pattern)
    seu$precent.hb <- PercentageFeatureSet(seu,pattern = hb_pattern)

# cell cyle
seu<- CellCycleScoring(seu,  #Insufficient data values to produce 24 bins.
                       s.features = cell.cycle.genes$human$s.genes, g2m.features = cell.cycle.genes$human$g2m.genes, set.ident = TRUE)  ### cell.cycle.genes


# Cell filtration
Idents(seu)<- 'orig.ident'
VlnPlot(seu,features = c('nCount_RNA','nFeature_RNA','percent.mito','percent.ribo','precent.hb'),ncol = 5,pt.size = 0)
seu_filtered <- subset(seu,
                    percent.mito <= 15 &
                    percent.ribo <= 20 &
                    nCount_RNA >=500  &
                    nFeature_RNA >=250  & nFeature_RNA<7000)


features<-c("nCount_RNA", "nFeature_RNA", "percent.mito",'percent.ribo','precent.hb')
p_before <- SCpubr::do_ViolinPlot(sample = seu, ncol = 5,
    features = features)
p_after <- SCpubr::do_ViolinPlot(sample = seu_filtered, ncol = 5,
    features = features)


before<-seu@meta.data %>% group_by(orig.ident) %>% summarise(n=n()) %>% data.frame() %>% dplyr::rename(before_qc=n)
after<-seu_filtered@meta.data %>% group_by(orig.ident) %>% summarise(n=n()) %>% data.frame() %>% dplyr::rename(after_qc=n)
cell_counts<-before %>% left_join(after)
cell_counts$keep<-round(100*cell_counts$after_qc/cell_counts$before_qc,2)
colnames(cell_counts)[4]<-"keep(%)"


# Normalization and Dimensionality reduction
seu_clustered <- NormalizeData(seu_filtered) %>% 
            FindVariableFeatures() %>% 
            ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>%
            RunPCA(features = VariableFeatures(.)) %>% 
            RunHarmony(group.by.vars = "orig.ident") %>% 
            RunUMAP(reduction = 'harmony',dims = 1:30) %>%
            RunTSNE(reduction = 'harmony',dims = 1:30) 

elbow<-ElbowPlot(seu_clustered,reduction = "harmony", ndims = 30)
n.pcs=20
resolutions <- seq(0.2,0.8,by = 0.2)
seu_clustered <- seu_clustered %>%  
                 FindNeighbors(reduction = 'harmony', dims = 1:n.pcs) %>% 
                 FindClusters(resolution = resolutions)

marker_list<- list(Rods=c('Pde6a','Ppef2','Nr2e3','Rho','Tulp1'),
                  Cones=c('Gnat2','Opn1mw','Opn1sw','Pde6h'),
                  `Bipolar cells`=c('Camk2b','Grm6','Tmem215','Trpm1','Nyap2','Kcnma1'),
                  `Amacrine cells`=c('Gad1','C1ql2','Meg3','Nrxn1'),
                  Microglia=c('C1qa','Tmem119','Aif1','Laptm5','Tyrobp','Cx3cr1','Adgre1','Abca9'),
                  `Vascular cells`=c('Cd34','Cdh5','Rgs5'),
                  Astrocytes=c('Aqp4','Aldoc'),
                   `Muller glia`=c('Cp','Slc1a3'),
                   Immune=c('Ptprc'),
                   Schwann=c('Ttr','Mbp')
                  )

major_dot<-DotPlot(sample = seu_clustered, features = marker_list)

seu_annotated<- seu_clustered
seu_annotated$classI_anno<- seu_annotated$seurat_clusters
seu_annotated$classI_anno<- gsub('^0$|^10$|^12$','Rod',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^3$','Cone',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^7$|^11$','Bipolar',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^15$','Amacrine',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^17$','Immune',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^4$|^8$|^14$|^6$','Microglia',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^2$|^9$','Schwann',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^18$|^5$|^19$','Astrocytes',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^1$|^16$','Muller glia',seu_annotated$classI_anno)
seu_annotated$classI_anno<- gsub('^13$','Vascular cells',seu_annotated$classI_anno)

fig02<- CellDimPlot(seu_annotated,group.by = 'classI_anno',show_stat = F,label_insitu = T,
                    reduction = 'tsne',label = T,label.bg = 'white',label.fg = 'black')

FeaturePlot(seu_annotated,features = c('Rho','Tulp1'),cols = c('grey','red'),reduction = 'tsne')
FeaturePlot(seu_annotated,features = c('Cp','Slc1a3'),cols = c('grey','red'),reduction = 'tsne')
