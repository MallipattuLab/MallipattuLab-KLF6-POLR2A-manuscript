##Comparing Control, OE, Control_AA, OE_AA
#File names are prepared for pre_mRNA alignment

library(dplyr)
library(Seurat)
library(Matrix)
library(cowplot)
library(rmarkdown)

#Create Seurat objects for each sample first
data_dir <- 'filtered_feature_bc_matrix_C/'
C_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_C/')
C_Data <- CreateSeuratObject(counts = C_Data,
                              min.cells = 3,
                              min.features = 200,
                              project = "Control")

data_dir <- 'filtered_feature_bc_matrix_SC11_3/'
C_AA_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_C_AA/')
C_AA_Data <- CreateSeuratObject(counts = C_AA_Data,
                                  min.cells = 3,
                                   min.features = 200,
                                  project = "Control_AA")

data_dir <- 'filtered_feature_bc_matrix_OE/'
OE_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_OE/')
OE_Data <- CreateSeuratObject(counts = OE_Data,
                              min.cells = 3,
                              min.features = 200,
                              project = "OE")

data_dir <- 'filtered_feature_bc_matrix_OE_AA/'
OE_AA_Data <- Read10X(data.dir = 'filtered_feature_bc_matrix_OE_AA/')
OE_AA_Data <- CreateSeuratObject(counts = OE_AA_Data,
                              min.cells = 3,
                              min.features = 200,
                              project = "OE_AA")



#Creates new merged seurat object
SC11.combined <- merge(C_Data, y = c(C_AA_Data, OE_Data, OE_AA_Data),   add.cell.ids = c("C_AA", "OE", "C", "OE_AA"), project = "Combined")
SC11 <- SC11.combined

# notice the cell names now have an added identifier
head(colnames(SC11))

#See original cell numbers
table(SC11$orig.ident)



# Get number of cells per cluster and per sample of origin
table(data.combined@meta.data[["seurat_clusters"]], data.combined@meta.data$orig.ident)

#Assess Mito%
SC11[["percent.mt"]] <- PercentageFeatureSet(SC11, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(SC11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
plot1 <- FeatureScatter(SC11, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SC11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#mitochondrial content less than 10 and other QC features (will change based on the graphs above (Data))
SC11<- subset(SC11, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)


#Identify Variable Features
obj.list <- SplitObject(SC11, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Takes a long time to run (~30 min)
data.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:30)

DefaultAssay(data.combined) <- "integrated"


# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)

#UMAP and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(data.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE)
plot_grid(p1,p2)


#show C_AA and KO side by side
DimPlot(data.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)



# Get number of cells per cluster and per sample of origin
table(data.combined@meta.data[["seurat_clusters"]], data.combined@meta.data$orig.ident)

#Get average expression of genes in each cluster
cluster.averages <- AverageExpression(data.combined, assays = "RNA", add.ident = "orig.ident")
write.csv(cluster.averages[["RNA"]], file = "CombinedClusterAveragespremRNA.csv")

#save file 
saveRDS(data.combined, file = "SC11_first.rds")




#load saved file
data.combined<- readRDS("SC11_first.rds")

#set default assay to RNA (as opposed to integrated)
DefaultAssay(data.combined) <- "RNA"




new.cluster.idss <- c("0", "1","2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13","14", "15", "16", "17", "18", "19", "20", "21", "22")
for(i in new.cluster.idss) {
  clus.markers <- FindMarkers(data.combined, ident.1 = i, ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
  head(clus.markers)
  str1 = i
  str2 = "deg.csv"
  print(i)
  result = paste(str1,str2, sep = ".")
  write.csv(clus.markers, file = result)
}
C0 <- read.csv("0.deg.csv")
C1 <- read.csv("1.deg.csv")
C2 <- read.csv("2.deg.csv")
C3 <- read.csv("3.deg.csv")
C4 <- read.csv("4.deg.csv")
C5 <- read.csv("5.deg.csv")
C6 <- read.csv("6.deg.csv")
C7 <- read.csv("7.deg.csv")
C8 <- read.csv("8.deg.csv")
C9 <- read.csv("9.deg.csv")
C10 <- read.csv("10.deg.csv")
C11 <- read.csv("11.deg.csv")
C12 <- read.csv("12.deg.csv")
C13 <- read.csv("13.deg.csv")
C14 <- read.csv("14.deg.csv")
C15 <- read.csv("15.deg.csv")
C16 <- read.csv("16.deg.csv")
C17 <- read.csv("17.deg.csv")
C18 <- read.csv("18.deg.csv")
C19 <- read.csv("19.deg.csv")
C20 <- read.csv("20.deg.csv")
C21 <- read.csv("21.deg.csv")
C22 <- read.csv("22.deg.csv")
#C23 <- read.csv("23.deg.csv")
#C24 <- read.csv("24.deg.csv")


C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2,  "C3" = C3,"C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12,  "C13" = C13,"C14" = C14,  "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22 )
#list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,"C12" = C12,  "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18,"C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22 )
write.xlsx(list_of_datasets, file = "Summary.DEGforclusters.Genes.xlsx")


#DEG between genotypes (choose ident.1 and ident.2)
for(i in new.cluster.idss) {
  a = paste(i, "Control", sep = "_")
  print(a)
  b = paste(i, "Control_AA", sep = "_")
  c = paste(i, "OE", sep = "_")
  d = paste(i, "OE_AA", sep = "_")
  MC.response <- FindMarkers(data.combined, ident.1 = a, ident.2 = b, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}



C0 <- read.csv("0comparison.csv")
C1 <- read.csv("1comparison.csv")
C2 <- read.csv("2comparison.csv")
C3 <- read.csv("3comparison.csv")
C4 <- read.csv("4comparison.csv")
C5 <- read.csv("5comparison.csv")
C6 <- read.csv("6comparison.csv")
C7 <- read.csv("7comparison.csv")
C8 <- read.csv("8comparison.csv")
C9 <- read.csv("9comparison.csv")
C10 <- read.csv("10comparison.csv")
C11 <- read.csv("11comparison.csv")
C12 <- read.csv("12comparison.csv")
C13 <- read.csv("13comparison.csv")
C14 <- read.csv("14comparison.csv")
C15 <- read.csv("15comparison.csv")
C16 <- read.csv("16comparison.csv")
C17 <- read.csv("17comparison.csv")
C18 <- read.csv("18comparison.csv")
C19 <- read.csv("19comparison.csv")
C20 <- read.csv("20comparison.csv")
C21 <- read.csv("21comparison.csv")
C22 <- read.csv("c22comparison.csv")
#C23 <- read.csv("c23comparison.csv")
#C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

#list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3,"C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11,"C12" = C12, "C13" = C13, "C14" = C14,"C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21)
#write.xlsx(list_of_datasets, file = "SummaryXXX vs XXX.DEG.Genes.xlsx")

# label cluster
data.combined <- RenameIdents(data.combined, `0` = "PT-S3-1", `1` = "PT-S1-S2-1", `2` = "PT-S1-S2-2", 
                                `3` = "PT-S3-2", `4` = "PT-S3-3", `5` = "DCT", `6` = "InjuredPT-A", `7` = "LH-AL-1",`8` = "EC-1", `9` = "CNT", 
                                `10` = "EC-2", `11` = "LH-AL-2", `12` = "InjuredPT-B", '13' = "CD/PC", '14'= "MC", '15' = "IC-B-1", 
                                '16' = "DCT/CNT", '17' = "Podocytes", '18' = "IC-A", '19' = "IC-B-2", '20' = "Macrophage", '21' = "LH-DL", '22' = "LH-AL-3")


# rearrange cluster
levels(data.combined)
names(Idents(data.combined)) <- levels(data.combined)
levels(x=data.combined) <- c("Podocytes", "MC", "EC-1", "EC-2", "PT-S1-S2-1", "PT-S1-S2-2", "PT-S3-1", "PT-S3-2", "PT-S3-3", "InjuredPT-A", "InjuredPT-B", "LH-DL", "LH-AL-1", "LH-AL-2", "LH-AL-3", "DCT", "DCT/CNT", "CNT", "CD/PC", "IC-A", "IC-B-1", "IC-B-2", "Macrophage")


