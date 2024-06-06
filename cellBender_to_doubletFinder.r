

library(Seurat)
library(scCustomize)
library(DoubletFinder)

cb.sham.data <- Read_CellBender_h5_Mat(file_name = "")

sham$orig.ident <- "sham"
sham <- Add_Mito_Ribo_Seurat(sham, species = "mouse") #Store in "percent_mito", "percent_ribo"
sham <- PercentageFeatureSet(sham, "^Hb[ab]-", col.name = "percent_Hb")
Idents(sham) <- "orig.ident"

QC_fig.1 <- FeatureScatter(sham, feature1 = "percent_mito", 
                           feature2 = "nFeature_RNA") +
  FeatureScatter(opto4h, feature1 = "percent_mito", 
                 feature2 = "nFeature_RNA") +
  FeatureScatter(opto24h, feature1 = "percent_mito", 
                 feature2 = "nFeature_RNA")
QC_fig.2 <- FeatureScatter(sham, feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA") +
  FeatureScatter(opto4h, feature1 = "nCount_RNA", 
                 feature2 = "nFeature_RNA") +
  FeatureScatter(opto24h, feature1 = "nCount_RNA", 
                 feature2 = "nFeature_RNA")

sham <- subset(sham, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 50 & percent_Hb < 5)



#proceed doubletFinder sham
test.sham <- sham %>% 
  SCTransform(vst.flavor = "v2", verbose = FALSE) %>% 
  RunPCA(dims = 1:20) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters()

sweep.res.list_sham <- paramSweep(test.sham, PCs = 1:20, sct = T)
#unloadNamespace("KernSmooth")
#library(KernSmooth)
#Find pK
sweep.stats_sham <- summarizeSweep(sweep.res.list_sham, GT = FALSE)
bcmvn_sham <- find.pK(sweep.stats_sham)
ggplot(bcmvn_sham, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() # 0.16
pK <- bcmvn_sham %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))
annotations <- test.sham@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(test.sham@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
test.sham <- doubletFinder(test.sham, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
DimPlot(test.sham, reduction = 'umap', group.by = "DF.classificatiseurat_clustersons_0.25_0.16_291")
table(test.sham@meta.data$DF.classifications_0.25_0.16_291)
Idents(test.sham) <- "DF.classifications_0.25_0.16_291"
sham.sig <- subset(test.sham, idents = "Singlet")
DimPlot(sham.sig, reduction = 'umap', group.by = "seurat_clusters")

###Final filtered sham
sham.fil <- sham[,colnames(sham) %in% colnames(sham.sig)]

##export to rds
SaveSeuratRds(sham.fil, file = "sham_fil.rds")

