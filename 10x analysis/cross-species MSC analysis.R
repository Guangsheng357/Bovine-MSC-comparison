#get the 1:1 orthologs across different species using Ensembl Biomart
#read the alignment output from the filtered_feature_bc_matrix folder to create the gene expression dataframe,using cow for example
library(Matrix)
matrix_dir = "your alignment output directory"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

b<-as.data.frame(mat)#create the expression dataframe 
#filter the expression dataframe rows to only keep the orthologs, using cow for example
names<-read.table("cow.names",header = F)#cow.names is the ortholog gene list in cow
d<-b[names$V1,]#filter the rows using cow ortholog gene name vector
d<-na.omit(d)#remove the rows containing NA
write.csv(d,"cow_gene.csv",quote = F)#save the cow gene expression result

#for other species, repeat the steps above to get the expression data, make sure to convert the gene names into the same format for cross-species comparison

#read the expression csv into the seurat and perform integration
#read the expression data that only contains the ortholog
cow<-read.csv("cow.csv",header = T,row.names = 1)
mouse<-read.csv("mouse.csv",header = T,row.names = 1)
human_1<-read.csv("human1.csv",header = T,row.names = 1)
human_2<-read.csv("human2.csv",header = T,row.names = 1)
human_3<-read.csv("human3.csv",header=T,row.names = 1)

#create the seurat object
cow.cell <-CreateSeuratObject(counts = cow, project = "cow", min.cells = 3, min.features = 200)
mouse.cell <-CreateSeuratObject(counts = mouse, project = "mouse", min.cells = 3, min.features = 200)
h_1 <-CreateSeuratObject(counts = human_1, project = "human 1",  min.features = 200,min.cells = 3)
h_2 <-CreateSeuratObject(counts = human_2, project = "human 2", min.cells = 3, min.features = 200)
h_3 <-CreateSeuratObject(counts = human_3, project = "human 3", min.cells = 3, min.features = 200)

#add tissue variable
cow.cell$species <-"COW"
mouse.cell$species <-"MOUSE"
h_1$species <-"HUMAN"
h_2$species <-"HUMAN"
h_3$species <-"HUMAN"

#add sequence batch
cow.cell$batch <-"cow"
mouse.cell$batch <-"mouse"
h_1$batch <-"human_1"
h_2$batch <-"human_2"
h_3$batch <-"human_3"

#mitochondria percent
#since I have converted the all the gene names into human species, I can directly use MT feature to choose the mitochondrial genes
cow.cell <- PercentageFeatureSet(cow.cell, "^MT", col.name = "percent_mito")
mouse.cell <- PercentageFeatureSet(mouse.cell, "^MT", col.name = "percent_mito")
h_1 <- PercentageFeatureSet(h_1, "^MT", col.name = "percent_mito")
h_2 <- PercentageFeatureSet(h_2, "^MT", col.name = "percent_mito")
h_3 <- PercentageFeatureSet(h_3, "^MT", col.name = "percent_mito")

#remove low quality cells
mouse.cell<- subset(mouse.cell, subset = nFeature_RNA <=5000  & percent_mito <=10 & nCount_RNA>=200 )
cow.cell<- subset(cow.cell, subset = nFeature_RNA <=5000  & percent_mito <=10 & nCount_RNA>=200 )
h_1<- subset(h_1, subset = nFeature_RNA <=5000  & percent_mito <=10 & nCount_RNA>=200 )
h_2<- subset(h_2, subset = nFeature_RNA <=5000  & percent_mito <=10 & nCount_RNA>=200 )
h_3<- subset(h_3, subset = nFeature_RNA <=5000  & percent_mito <=10 & nCount_RNA>=200 )

#get the cell cycle gene list
c.s.genes <- cc.genes.updated.2019$s.genes
c.g2m.genes <- cc.genes.updated.2019$g2m.genes

##integrate data into one Seurat object
hypo.list<-c(cow.cell, mouse.cell, h_1,h_2,h_3)
for (i in 1:length(hypo.list)) {
  hypo.list[[i]] <- NormalizeData(hypo.list[[i]], verbose = FALSE)
  hypo.list[[i]] <- CellCycleScoring(hypo.list[[i]],
                                     s.features = c.s.genes, g2m.features = c.g2m.genes, set.ident = TRUE)
  Idents(hypo.list[[i]]) <- "orig.ident"
  hypo.list[[i]] <- FindVariableFeatures(hypo.list[[i]], selection.method = "vst", nfeatures = 5000, verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = hypo.list)
hypo.anchors <- FindIntegrationAnchors(object.list = hypo.list, anchor.features = features)
hypo.integrated <- IntegrateData(anchorset = hypo.anchors)#integrate the dataset

DefaultAssay(hypo.integrated) <- "integrated"#swith to integrated assay
hypo.integrated <- ScaleData(hypo.integrated, verbose = T,vars.to.regress = c("batch","S.Score", "G2M.Score","percent_mito"))#scale the data, including the batch effect
#PCA and UMAP analysis
hypo.integrated <- RunPCA(hypo.integrated, npcs = 50, verbose = TRUE)
ElbowPlot(hypo.integrated,ndims = 50)
hypo.integrated <- FindNeighbors(hypo.integrated, dims = 1:30)
hypo.integrated <- FindClusters(hypo.integrated, resolution = 0.2)
hypo.integrated <- RunUMAP(hypo.integrated, reduction = "pca", dims = 1:30)


DimPlot(hypo.integrated,group.by = "species",label = FALSE)+theme(plot.title = element_blank())#plot UMAP
#find the differential expression genes
DefaultAssay(hypo.integrated) <- "RNA"
de_obj<-hypo.integrated
Idents(de_obj)<-de_obj$species#set the identity to species
#pair-wise compariosn
bulk_MAST<- FindMarkers(de_obj, ident.1 = "species_1" , ident.2 = "species_2", verbose = T, test.use = "MAST", only.pos = FALSE, logfc.threshold = 0, latent.vars = c("nCount_RNA"))



