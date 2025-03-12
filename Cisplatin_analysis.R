setwd("/Klf6Project/AdditionalDataSet/")

# Get Dataset
# library(GEOquery)
# 
# getGEOSuppFiles("GSE197266")
# getGEOSuppFiles("GSE190887")

# z <- read.delim("GSE190887/GSM5733424_count.MM.txt", sep="\t")


library(Seurat)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(harmony)

CTRL <- read.delim("GSE197266/GSE197266_Ctrl_merge_as_matrix_RNA_data.txt.gz", sep="\t", header = TRUE , row.names = "Gene") #
AKI <- read.delim("GSE197266/GSE197266_AKI_merge_as_matrix_RNA_data.txt.gz", sep="\t", header = TRUE, row.names = "Gene") # 

SeuratObject_Control <- CreateSeuratObject(CTRL, project = "CTRL", assay = "RNA",
                       min.cells = 0, min.features = 0, names.field = 1,
                       names.delim = "_", meta.data = NULL)

SeuratObject_AKI <- CreateSeuratObject(AKI, project = "AKI", assay = "RNA",
                                           min.cells = 0, min.features = 0, names.field = 1,
                                           names.delim = "_", meta.data = NULL)
Klf6_combined.list <- lapply(X = c(SeuratObject_Control,
                                   SeuratObject_AKI), FUN = function(x) {
                                     x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
                                     VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
                                     
                                     x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 3)
                                     
                                     # run sctransform
                                     x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
                                   })


features <- SelectIntegrationFeatures(object.list = Klf6_combined.list)

Klf6.anchors <- FindIntegrationAnchors(object.list = Klf6_combined.list, anchor.features = features)

Klf6_Additional_data <- IntegrateData(anchorset = Klf6.anchors)

saveRDS(Klf6_Additional_data, file = "Klf6_Additional_data.rds")

Klf6_Additional_data <- readRDS("Klf6_Additional_data.rds")

Idents(Klf6_Additional_data) <- Klf6_Additional_data$orig.ident
Klf6_Additional_data <- subset(Klf6_Additional_data, idents =c("con1", "con2", "CP1", "CP2") )

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Klf6_Additional_data), 10)

all.genes <- rownames(Klf6_Additional_data)
Klf6_Additional_data <- ScaleData(Klf6_Additional_data, features = all.genes)

Klf6_Additional_data <- RunPCA(Klf6_Additional_data, features = VariableFeatures(object = Klf6_Additional_data))

print(Klf6_Additional_data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Klf6_Additional_data, dims = 1:2, reduction = "pca")

ElbowPlot(Klf6_Additional_data, ndims = 50)


Klf6_Additional_data$Sample <- Klf6_Additional_data[["orig.ident"]]

Klf6_Additional_data <- RunHarmony(Klf6_Additional_data, "Sample")
ElbowPlot(Klf6_Additional_data, ndims = 50, reduction = "harmony")

Klf6_Additional_data <- FindNeighbors(Klf6_Additional_data, dims = 1:50, reduction = "harmony")
Klf6_Additional_data <- FindClusters(Klf6_Additional_data, resolution = 0.8)

Klf6_Additional_data <- RunUMAP(Klf6_Additional_data, reduction = "harmony", dims = 1:50)

Klf6_Additional_data$group <- Klf6_Additional_data$Sample
Klf6_Additional_data$group[which(Klf6_Additional_data$group %in% c("con1", "con2"))] <- "control"
Klf6_Additional_data$group[which(Klf6_Additional_data$group %in% c("CP1", "CP2"))] <- "Cisplatin"

Klf6_Additional_data.markers <- FindAllMarkers(Klf6_Additional_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Klf6_Additional_data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

Klf6_Additional_data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top1
DoHeatmap(Klf6_Additional_data, features = top1$gene) + NoLegend()

new.cluster.ids <- c("PT",  #0
                     "Macrophage",  #1
                     "PT", #2
                     "PT", #3
                     "PT", #4 
                     "PT",  #5
                     "DCT_CNT",  #6 
                     "PT", #7
                     "PT",        #8
                     "Macrophage", #9
                     "EC", #10
                     "Pod", #11
                     "PT",    #12
                     "PT", #13
                     "Fib",   #14
                     "ICA",  #15
                     "CDPC", #16
                     "HLAL", #17
                     "Per", #18
                     "Tcell",  #19
                     "Macrophage",  #20
                     "EC2", #21
                     "EC",        #22
                     "PEC", #23
                     "Macrophage", #24
                     "Macrophage", #25
                     "ICB",    #26
                     "Macrophage" #27
)
names(new.cluster.ids) <- levels(Klf6_Additional_data)
Klf6_Additional_data <- RenameIdents(Klf6_Additional_data, new.cluster.ids)

DotPlot(Klf6_Additional_data, features = c("Slc34a1", "Lrp2")
        )+ RotatedAxis()

saveRDS(Klf6_Additional_data, file = "Klf6_Additional_data_anno.rds")

Klf6_Additional_data <- readRDS("Klf6_Additional_data_anno.rds")


PT_Subset <- subset(Klf6_Additional_data, idents =c("PT") )

all.genes <- rownames(PT_Subset)
PT_Subset <- ScaleData(PT_Subset, features = all.genes)

PT_Subset <- RunPCA(PT_Subset, features = VariableFeatures(object = PT_Subset))

print(PT_Subset[["pca"]], dims = 1:5, nfeatures = 5)

PT_Subset <- RunHarmony(PT_Subset, "Sample")

PT_Subset <- FindNeighbors(PT_Subset, dims = 1:50, reduction = "harmony")
PT_Subset <- FindClusters(PT_Subset, resolution = 0.8)

PT_Subset <- RunUMAP(PT_Subset, reduction = "harmony", dims = 1:50)

DimPlot(PT_Subset, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(PT_Subset, reduction = "umap", split.by = "group", label = TRUE) + NoLegend()

PT_Subset.markers <- FindAllMarkers(PT_Subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PT_Subset.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DotPlot(PT_Subset, features = c("Slc34a1", "Slc5a2", # Healthy PT S1
                                "S100a11",
                                "Vcam1", "Akap12","Havcr1",# Injury PT S1
                                "Fos",
                                "Slc22a6", #Healthy PT S2
                                "Krt20",  # New PT
                                "Slc5a10", #Healthy PT S3
                                "Atp11a",  # Injured PT S3
                                "Fxyd5",  # Maladaptive PT S1 
                                "Mki67", "Top2a", # Proliferative
                                "Klf6", "Polr2a"
                                ), assay = "SCT"
        )+ RotatedAxis()
ggplot2::ggsave("Klf6OE_Additional_Data_RNA_Dotplot.pdf", width = 10, height = 9)

new.cluster.ids <- c("Healthy_PTS1",  #0
                     "AcuteInjuryPT",  #1
                     "Healthy_PTS1", #2
                     "Healthy_PTS1", #3
                     "Healthy_PTS1", #4 
                     "Injury_PTS1",  #5
                     "PreInjuryPT",  #6 
                     "Healthy_PTS2", #7
                     "NewPT",        #8
                     "Healthy_PTS2", #9
                     "Injured_PTS3", #10
                     "Maladaptive_PTS1", #11
                     "Healthy_PTS1",    #12
                     "Maladaptive_PTS1", #13
                     "ProliferativePT"   #14
)
names(new.cluster.ids) <- levels(PT_Subset)
PT_Subset <- RenameIdents(PT_Subset, new.cluster.ids)

PT_Subset$ClusterName <- Idents(PT_Subset)

saveRDS(PT_Subset, file = "PT_Subset_annotated.rds")
PT_Subset <- readRDS("PT_Subset_annotated.rds")

Idents(PT_Subset) <- PT_Subset$ClusterName

newPTs <- Cells(PT_Subset)[which(PT_Subset$ClusterName == "NewPT")]
Controls <- Cells(PT_Subset)[which(PT_Subset$ClusterName %in% c("Healthy_PTS1","Healthy_PTS2"))]

PT_Subset$ClusterName2 <- levels(PT_Subset$ClusterName)[as.integer(PT_Subset$ClusterName)]
PT_Subset$ClusterName2[Controls] <- "Control"

Idents(PT_Subset) <- PT_Subset$ClusterName2
comparison <- FindMarkers(PT_Subset, 
                          ident.1 = "NewPT", 
                          ident.2 = "Control",
                         logfc.threshold = 0.25,
                         verbose = FALSE)
comparison$gene <- row.names(comparison)
openxlsx::write.xlsx(comparison, "DEGs_clusters_48_CP_newPTvsHealthyPT.xlsx")

DotPlot(PT_Subset, features = c("Polr2a", "Klf6"
                                            ), assay = "SCT"
)+ RotatedAxis()
ggplot2::ggsave("Klf6OE_Additional_Data_Polr2aKlf6.pdf", width = 6, height = 9)

Idents(PT_Subset) <- PT_Subset$Sample
PT_Subset_CP_Ctrl <- subset(PT_Subset, idents =c("con1", "con2", "CP1", "CP2") )

all.genes <- rownames(PT_Subset_CP_Ctrl)
PT_Subset_CP_Ctrl <- ScaleData(PT_Subset_CP_Ctrl, features = all.genes)

PT_Subset_CP_Ctrl <- RunPCA(PT_Subset_CP_Ctrl, features = VariableFeatures(object = PT_Subset_CP_Ctrl))

PT_Subset_CP_Ctrl <- RunHarmony(PT_Subset_CP_Ctrl, "Sample")
ElbowPlot(PT_Subset_CP_Ctrl, ndims = 50, reduction = "harmony")

PT_Subset_CP_Ctrl <- FindNeighbors(PT_Subset_CP_Ctrl, dims = 1:50, reduction = "harmony")
PT_Subset_CP_Ctrl <- FindClusters(PT_Subset_CP_Ctrl, resolution = 1)

PT_Subset_CP_Ctrl <- RunUMAP(PT_Subset_CP_Ctrl, reduction = "harmony", dims = 1:50)

PT_Subset_CP_Ctrl.markers <- FindAllMarkers(PT_Subset_CP_Ctrl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PT_Subset_CP_Ctrl.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


PT_Subset_CP_Ctrl.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top1

DotPlot(PT_Subset_CP_Ctrl, features = unique(top1$gene)
        , assay = "SCT"
)+ RotatedAxis()

DotPlot(PT_Subset_CP_Ctrl, features = c("Slc5a12","Slc34a1", "Slc5a2",  # Healthy PT S1
                                        "Slc13a3","Slc22a6",  #Healthy PT S2
                                        "Slc7a13","Slc5a10","Atp11a", #Healthy PT S3
                                "Ccn1","Klf6","Fos",
                                "Krt20","Cdkn1a",  "Sfn",
                                "Mt1", "Mt2","Cryab",
                                # "Vcam1", "Akap12","Havcr1",# Injury PT S1
                                # "Cd74",
                                
                                  # New PT
                                  # Injured PT S3
                                "Polr2a"
), assay = "SCT"
)+ RotatedAxis()

ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_RNA_Dotplot.pdf", width = 10, height = 9)

Idents(PT_Subset_CP_Ctrl) <- PT_Subset_CP_Ctrl$seurat_clusters
new.cluster.ids <- c("Healthy_PTS1_1",  #0
                     "Healthy_PTS1_2",  #1 
                     "PreInjuryPT",     #2
                     "Healthy_PTS1_2",  #3
                     "AcuteInjuryPT_1",   #4 
                     "AcuteInjuryPT_2",   #5
                     "Healthy_PTS1_1",  #6 
                     "Healthy_PTS1_2",  #7
                     "NewPT",           #8
                     "Healthy_PTS2",  #9
                     "Healthy_PTS1_3",  #10
                     "Healthy_PTS1_2",  #11
                     "AcuteInjuryPT_3",   #12
                     "Healthy_PTS3",    #13
                     "Healthy_PTS1_1",  #14
                     "Injury_PTS3",    #15
                     "Healthy_PTS1_1"   #16
)
names(new.cluster.ids) <- levels(PT_Subset_CP_Ctrl)
PT_Subset_CP_Ctrl <- RenameIdents(PT_Subset_CP_Ctrl, new.cluster.ids)
PT_Subset_CP_Ctrl$ClusterName <- Idents(PT_Subset_CP_Ctrl)
PT_Subset_CP_Ctrl$ClusterNameReodered <- factor(PT_Subset_CP_Ctrl$ClusterName, levels = c("Healthy_PTS1_1", "Healthy_PTS1_2","Healthy_PTS1_3",
                                                                                          "Healthy_PTS2", "Healthy_PTS3", 
                                                                                          "PreInjuryPT", 
                                                                                          "AcuteInjuryPT_1", "AcuteInjuryPT_2","AcuteInjuryPT_3",
                                                                                          "Injury_PTS3",
                                                                                          "NewPT"))

Idents(PT_Subset_CP_Ctrl) <- PT_Subset_CP_Ctrl$ClusterNameReodered

PT_Subset_CP_Ctrl.markers <- FindAllMarkers(PT_Subset_CP_Ctrl, only.pos = FALSE, 
                                            min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")

clusters <- levels(PT_Subset_CP_Ctrl$ClusterName)
Cluster_DEGs <- list()
for (cluster in clusters) {
  clusterDEGs <- PT_Subset_CP_Ctrl.markers[which(PT_Subset_CP_Ctrl.markers$cluster == cluster),]
  Cluster_DEGs[[cluster]] <- clusterDEGs
}
openxlsx::write.xlsx(Cluster_DEGs, "DEGs_clusters_48_CtrlCP.xlsx")

comparison <- FindMarkers(PT_Subset_CP_Ctrl, 
                          ident.1 = "NewPT", 
                          ident.2 = c("Healthy_PTS1_1", "Healthy_PTS1_2","Healthy_PTS1_3",
                                      "Healthy_PTS2", "Healthy_PTS3"),
                          logfc.threshold = 0.25,
                          assay = "SCT",
                          verbose = FALSE)
comparison$gene <- row.names(comparison)
openxlsx::write.xlsx(comparison, "DEGs_clusters_48_CtrlCP_newPTvsHealthyPT.xlsx")

DimPlot(PT_Subset_CP_Ctrl, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_UMAP.pdf", width = 10, height = 9)

DimPlot(PT_Subset_CP_Ctrl, reduction = "umap",  label = TRUE, pt.size = 0.5, split.by = "group")+ NoLegend()
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_SplitUMAP.pdf", width = 16, height = 9)

DotPlot(PT_Subset_CP_Ctrl, features = c("Polr2a", "Klf6"
                                            ), assay = "SCT"
 )+ RotatedAxis()
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_Polr2aKlf6.pdf", width = 6, height = 9)
PT_Subset_CP_Ctrl[["RNA"]] <- NULL

DefaultAssay(PT_Subset_CP_Ctrl) <- "SCT"
FeaturePlot(PT_Subset_CP_Ctrl, features = "Polr2a", split.by = "group", order = T, min.cutoff = 0.5)
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_Polr2a_SplitUMAPbyGroup.pdf", width = 16, height = 9)
FeaturePlot(PT_Subset_CP_Ctrl, features = "Klf6", split.by = "group", order = T, min.cutoff = 1)
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_Klf6_SplitUMAPbyGroup.pdf", width = 16, height = 9)

############### Stacked Histogram
library(ggplot2)

cluster_composition <- data.frame(cluster = PT_Subset_CP_Ctrl$ClusterNameReodered,
                                  sample = PT_Subset_CP_Ctrl$group
)

data <- as.data.frame(table(cluster_composition))
# levels(data[,2]) <- c("con", "con", "CP", "CP")
# for(Sample in sort(unique(PT_Subset_CP_Ctrl$Sample))){
#   data[which(data$sample %in% Sample),]$Freq <- data[which(data$sample %in% Sample),]$Freq/sum(data[which(data$sample %in% Sample),]$Freq)
# }

g <- ggplot(data, aes(fill=sample, y=Freq, x=cluster))+ geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
g
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_stacked_bar.pdf", width = 10, height = 9)

g2 <- ggplot(data, aes(fill=sample, y=Freq, x=cluster)) +
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
g2
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_stacked_bar_Percentage.pdf", width = 10, height = 9)

################ Violin Plot
Idents(PT_Subset_CP_Ctrl) <- PT_Subset_CP_Ctrl$ClusterName
PT_Subset_NewPT <- subset(PT_Subset_CP_Ctrl, idents =c("NewPT") )
Idents(PT_Subset_NewPT) <- PT_Subset_NewPT$group

DotPlot(PT_Subset_NewPT, features = c("Polr2a", "Klf6"
), assay = "SCT"
)+ RotatedAxis()
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_Polr2aKlf6_newPT_bySample.pdf", width = 6, height = 9)
VlnPlot(PT_Subset_NewPT, features =  c("Polr2a", "Klf6"), assay = "SCT")


PT_Subset_PreInjuryPT<- subset(PT_Subset_CP_Ctrl, idents =c("PreInjuryPT") )
Idents(PT_Subset_PreInjuryPT) <- PT_Subset_NewPT$group

DotPlot(PT_Subset_PreInjuryPT, features = c("Polr2a", "Klf6"
), assay = "SCT"
)+ RotatedAxis()
ggplot2::ggsave("Klf6OE_Additional_Data_CtrlCP_Polr2aKlf6_PreInjuryPT_bySample.pdf", width = 6, height = 9)
VlnPlot(PT_Subset_PreInjuryPT, features =  c("Polr2a", "Klf6"), assay = "SCT")
