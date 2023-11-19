#setwd("your working directory")
#LOAD DATA
PBMSC<- Read10X(data.dir = "Blood/outs/filtered_feature_bc_matrix/")
BMMSC<- Read10X(data.dir = "Bone_Marrow/outs/filtered_feature_bc_matrix/")
ADMSC<-Read10X(data.dir = "Adipose/outs/filtered_feature_bc_matrix/")
PBMSC <-CreateSeuratObject(counts = PBMSC, project = "PBMSC", min.cells = 3, min.features = 200)
BMMSC <-CreateSeuratObject(counts = BMMSC, project = "BMMSC", min.cells = 3, min.features = 200)
ADMSC <-CreateSeuratObject(counts = ADMSC, project = "ADMSC", min.cells = 3, min.features = 200)
#add tissue variable
PBMSC$stim <-"PB-MSC"
BMMSC$stim <-"BM-MSC"
ADMSC$stim <-"AT-MSC"
#calculating mitochondria percent
mitochondrial.genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
PBMSC[["percent.mt"]] <- PercentageFeatureSet(PBMSC, features = mitochondrial.genes)
BMMSC[["percent.mt"]] <- PercentageFeatureSet(BMMSC, features = mitochondrial.genes)
ADMSC[["percent.mt"]] <- PercentageFeatureSet(ADMSC, features = mitochondrial.genes)

#calculating ribosome percent
PBMSC <- PercentageFeatureSet(PBMSC, "^RP[SL]", col.name = "percent_ribo")
BMMSC <- PercentageFeatureSet(BMMSC, "^RP[SL]", col.name = "percent_ribo")
ADMSC <- PercentageFeatureSet(ADMSC, "^RP[SL]", col.name = "percent_ribo")


#removing low quality cells
PBMSC<- subset(PBMSC, subset = nFeature_RNA >= 200 & nFeature_RNA <=6000  & percent.mt <=15 & nCount_RNA>=1000 & nCount_RNA<=40000)
BMMSC<- subset(BMMSC, subset = nFeature_RNA >= 200 & nFeature_RNA <=6000  & percent.mt <=15 & nCount_RNA>=1000 & nCount_RNA<=40000)
ADMSC<- subset(ADMSC, subset = nFeature_RNA >= 200 & nFeature_RNA <=6000  & percent.mt <=15 & nCount_RNA>=1000 & nCount_RNA<=40000)

#calculate cell cyle
library(biomaRt)
library(knitr)
convertHumanGeneList <- function(x) {
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                     version = 105)
  cow = useEnsembl("ensembl", dataset = "btaurus_gene_ensembl",
                   version = 105)
  
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                 values = x, mart = human, attributesL = c("external_gene_name"),
                 martL = cow, uniqueRows = T)
  
  humanx <- unique(genes[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

c.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
c.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)


##integrate data into one Seurat object
hypo.list<-c(PBMSC, BMMSC, ADMSC)
for (i in 1:length(hypo.list)) {
  hypo.list[[i]] <- NormalizeData(hypo.list[[i]], verbose = FALSE)
  hypo.list[[i]] <- CellCycleScoring(hypo.list[[i]],
                        s.features = c.s.genes, g2m.features = c.g2m.genes, set.ident = TRUE)
  Idents(hypo.list[[i]]) <- "orig.ident"
  hypo.list[[i]] <- FindVariableFeatures(hypo.list[[i]], selection.method = "vst", nfeatures = 5000, verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = hypo.list)
hypo.anchors <- FindIntegrationAnchors(object.list = hypo.list, anchor.features = features)
hypo.integrated <- IntegrateData(anchorset = hypo.anchors)

DefaultAssay(hypo.integrated) <- "integrated"
hypo.integrated <- ScaleData(hypo.integrated, verbose = T,vars.to.regress = c("stim","S.Score", "G2M.Score","percent.mt", "nFeature_RNA","percent_ribo"))
#perform PCA and UMAP analysis
hypo.integrated <- RunPCA(hypo.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(hypo.integrated,ndims = 50)
hypo.integrated <- RunUMAP(hypo.integrated, reduction = "pca", dims = 1:25)
hypo.integrated <- FindNeighbors(hypo.integrated, dims = 1:25)
hypo.integrated <- FindClusters(hypo.integrated, resolution = 0.2)

#finding markers and cluster anotation 
DefaultAssay(hypo.integrated) <- "RNA" 
hypo.markers <- FindAllMarkers(hypo.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
new.cluster.ids<-c("CD29++ MSC","CD29+ MSC","Macrophage/Monocyte","CD34+/CD105+ AT")
names(new.cluster.ids)<-levels(hypo.integrated)
hypo.integrated<-RenameIdents(hypo.integrated,new.cluster.ids)


#heatmap for top10 cluster markers
hypo.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(subset(hypo.integrated, downsample = 500),  features = top10$gene, size = 3, label = T,angle = 25) + scale_fill_viridis( na.value = "white")+guides(color="none")



