library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
memory.limit(size = 9999999999)
options(future.globals.maxSize= 15912896000)
library(openxlsx)

#set working directory into project folder
setwd("/Klf6Project")

addArchRThreads(threads = 1) 
addArchRGenome("mm10")

#Get Input Fragment Files
inputFiles <-  c("SC521" = "/SC52_ATAC_KLF6/SC521ATAC/outs/fragments.tsv.gz",
                 "SC522" = "/SC52_ATAC_KLF6/SC522ATAC/outs/fragments.tsv.gz",
                 "SC523" = "/SC52_ATAC_KLF6/SC523ATAC/outs/fragments.tsv.gz",
                 "SC524" = "/SC52_ATAC_KLF6/SC524ATAC/outs/fragments.tsv.gz")
#Create Arrow Files (disabled here)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#ArchRProject
proj <- ArchRProject(ArrowFiles)

#Filter Cells
proj <- proj[proj$nFrags > 2500 & proj$TSSEnrichment > 8] #& !is.na(proj$Gex_nUMI))

#Doublet Filtration. Currently disabled just for tutorial. If you want to filter doublets uncomment below.
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)


#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "IterativeLSI"
)

#UMAP
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

#Batch Correction Harmony
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#UMAP
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

rd <- getReducedDims(proj, reducedDims = "Harmony",dimsToUse = 1:30)
sd <- colSds(rd)
ggplot(data.frame(dims = 1:30, stdev = sd[1:30]))+geom_point(mapping = aes_string(x = 'dims', y = 'stdev'))+labs(x = "LSI",y="Standard Deviation")+ylim(0.05,2)

#Add Clusters
proj <- addClusters(proj, reducedDims = "Harmony", name = "Clusters", resolution = 0.6, force = TRUE, maxClusters = 30, dimsToUse = 2:30)

saveArchRProject(ArchRProj = proj, outputDirectory = "Klf6_Multimome_analysis", load = FALSE)
proj <- loadArchRProject(path = "Klf6_Multimome_analysis", force = FALSE, showLogo = TRUE)

#load previously annotated Seurat RNA Object
seRNA <- readRDS("Klf6_Clustered.rds")
seRNA

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "ClusterName",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

unique(seRNA$ClusterName)

plotEmbedding(
  proj, 
  embedding = "UMAPHarmony",
  colorBy = "cellColData", 
  name = "predictedGroup_Un"
)

cPTs <- paste0(unique(seRNA$ClusterName)[c(7:13,15:18)], collapse = "|")
cPTs

cNPTs <- paste0(unique(seRNA$ClusterName)[c(1:6,14, 19:28)], collapse = "|")
cNPTs

clustPT <- rownames(cM)[grep(cPTs, preClust)]
clustPT

clustNPT <- rownames(cM)[grep(cNPTs, preClust)]
clustNPT

rnaPTs <- colnames(seRNA)[grep(cPTs, seRNA$ClusterName)]
head(rnaPTs)

rnaNPTs <- colnames(seRNA)[grep(cNPTs, seRNA$ClusterName)]
head(rnaNPTs)

groupList <- SimpleList(
  PTs = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% clustPT],
    RNA = rnaPTs
  ),
  NonPTs = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% clustNPT],
    RNA = rnaNPTs
  )    
)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupList = groupList,
  groupRNA = "ClusterName",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
  force = TRUE
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneIntegrationMatrix", 
  depthCol = "nFrags",
  name = "IterativeLSI_RNA",
  force = TRUE
)

proj <- addCombinedDims(proj, reducedDims = c("IterativeLSI", "IterativeLSI_RNA"), dimsToUse = c(2:30, 2:30), name =  "LSI_Combined")

#Batch Correction Harmony
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "Harmony_Combined",
  groupBy = "Sample",
  force = TRUE
)

#UMAP
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony_Combined", 
  name = "UMAPHarmony_Combined", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)


#Add Clusters
proj <- addClusters(proj, reducedDims = "Harmony_Combined", name = "Clusters", resolution = 2, force = TRUE, maxClusters = 40, dimsToUse = 2:30)

saveArchRProject(ArchRProj = proj, outputDirectory = "Klf6_Multimome_analysis", load = FALSE)
proj <- loadArchRProject(path = "Klf6_Multimome_analysis", force = FALSE, showLogo = TRUE)

cM <- confusionMatrix(proj$Clusters, proj$predictedGroup_Co)
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
cbind(colnames(cM)[apply(cM, 1, which.max)], rownames(cM)) #Assignments

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

remapClust <- c(
  "C12" = "PTS2",     
  "C9" = "PTS3"  ,    
  "C17" = "ICB-2"   ,   
  "C15" = "Injury"   ,  
  "C11" = "PTS1-1"    , 
  "C24" = "DCT"        ,
  "C27" = "Per"        ,
  "C10" = "PTS3"     ,
  "C22" = "MTAL"       ,
  "C18" = "PTS2"     ,
  "C28" = "PC"        ,
  "C14" = "DCT-CNT-2"  ,
  "C7" = "Injury"      ,
  "C3" = "EC"          ,
  "C26" = "ICB-1"      ,
  "C5" = "Per"         ,
  "C6" = "Fib"         ,
  "C20" = "DTL-ATL"    ,
  "C19" = "Pod"        ,
  "C16" = "PTS1-1"     ,
  "C1" = "Macrophage"  ,
  "C21" = "DTL-ATL"    ,
  "C4" = "vascular"    ,
  "C2" = "Tcell"  ,
  "C8" = "Injury"      ,
  "C13" = "PTS2"     ,
  "C23" = "MTAL"       ,
  "C25" = "CNT"        
)

proj$Clusters2 <- mapLabels(proj$Clusters, newLabels = remapClust, oldLabels = names(remapClust))
proj$Clusters2[which(proj$predictedGroup_Co %in% c("CNT","DCT-CNT-1", "DCT") )] <- proj$predictedGroup_Co[which(proj$predictedGroup_Co %in% c("CNT","DCT-CNT-1", "DCT") )]
proj$Clusters2[which(proj$predictedGroup_Co %in% "PTS1-2" )] <- proj$predictedGroup_Co[which(proj$predictedGroup_Co %in% "PTS1-2" )]

proj$Clusters2[which(proj$Clusters2 %in% "PTS2-2" )] <- "Injury"
proj$Clusters2[which(proj$Clusters2 %in% "PTS2-1" )] <- "PTS2"

#Rename clusters
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.xlsx(markerList@listData, "DEGs_default.xlsx")

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters2", force = TRUE)

# Adding modification and treatment labels
modificationClust <- c(
  "SC521" = "WT",
  "SC522" = "Klf6_OE",
  "SC523" = "Klf6_OE",
  "SC524" = "WT"
)
proj$modification <- mapLabels(proj$Sample, oldLabels = names(modificationClust), newLabels = modificationClust)
treatmentClust <-  c(
  "SC521" = "AA",
  "SC522" = "AA",
  "SC523" = "Control",
  "SC524" = "Control"
)
proj$Treatment <- mapLabels(proj$Sample, oldLabels = names(treatmentClust), newLabels = treatmentClust)


proj$Identity <- paste(proj$modification, proj$Treatment, sep = "_")

proj$modificationClusters2Combined <- paste(proj$Clusters2, proj$modification, sep = "_")
proj$treatmentClusters2Combined <- paste(proj$Clusters2, proj$Treatment, sep = "_")
proj$identityClusters2Combined <- paste(proj$Clusters2, proj$Identity, sep = "_")

########################################## Motif Enrichment ##################################################################

pathToMacs2 <- findMacs2()
addArchRThreads(threads = 16) 

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)

proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
write.xlsx(markersMotifs, "Motifs.xlsx")

saveArchRProject(ArchRProj = proj, outputDirectory = "Klf6_Multimome_analysis", load = FALSE)
proj <- loadArchRProject(path = "/gpfs/projects/MallipattuGroup/Jiakang/Klf6Project/Klf6_Multimome_analysis", force = FALSE, showLogo = TRUE)

########################################## Motif Enrichment ##################################################################

pathToMacs2 <- findMacs2()
addArchRThreads(threads = 16) 

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)


proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
# proj <- addMotifAnnotations(ArchRProj = proj, name = "Motif", force = TRUE, motifPWMs = pwmList)

proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

for (cluster in names(table(proj$Clusters2))) {
  markerTest <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters2",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = cluster,
  )
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1.5"
  )
  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(10)), ], aes(x = rank, y = mlog10Padj, label = TF),
      size = 5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) + 
    ggtitle(cluster)
  plotPDF(ggUp+ theme(text = element_text(size = 20)) ,name = paste0(cluster,"-Markers-Motifs-Enriched-UP"), width = 8, height = 11, ArchRProj = proj, addDOC = FALSE)
  
  
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC <= -1.5"
  )
  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(10)), ], aes(x = rank, y = mlog10Padj, label = TF),
      size = 5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))+ 
    ggtitle(cluster)
  plotPDF(ggDo+ theme(text = element_text(size = 20)) ,name = paste0(cluster,"-Markers-Motifs-Enriched-DOWN"), width = 8, height = 11, ArchRProj = proj, addDOC = FALSE)
}

############################## Peak2GeneLinkage ###########################################
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony_Combined",
  useMatrix = "GeneIntegrationMatrix"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1000,
  returnLoops = TRUE
)


proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony_Combined",
  useMatrix = "GeneIntegrationMatrix"
)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c(
    "PTS2",
    "PTS3",
    "Injury",
    "PTS1-1",
    "PTS1-2"
  )
)

markerGenes  <- c(
  "Polr2a"
)

markers <- getMarkers(markersPeaks, 
                      cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  useGroups = c(
    "PTS2",
    "PTS3",
    "Injury",
    "PTS1-1",
    "PTS1-2"
  ),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  groupBy = "Clusters2", 
  tileSize = 100,
  geneSymbol = markerGenes, 
  sizes = c(10, 3, 3, 4),
  features = getMarkers(markersPeaks, 
                        cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE),
  ylim = c(0,0.99),
  upstream = 30000,
  downstream = 10000,
  loops = getPeak2GeneLinks(proj)
)

grid::grid.newpage()
grid::grid.draw(p$Polr2a)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-30000-10000.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 9, height = 7)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  useGroups = c(
    "PTS2",
    "PTS3",
    "Injury",
    "PTS1-1",
    "PTS1-2"
  ),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  groupBy = "Clusters2", 
  tileSize = 40,
  ylim = c(0,0.99),
  geneSymbol = markerGenes, 
  sizes = c(10, 3, 3, 4),
  features = getMarkers(markersPeaks, 
                        cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE),
  upstream = 10000,
  downstream = 10000,
  loops = getPeak2GeneLinks(proj)
)

grid::grid.newpage()
grid::grid.draw(p$Polr2a)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-10000-10000.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 9, height = 7)

######################################## All PT Trajectory #########
All_PT_Subset <- subsetArchRProject(
  ArchRProj = proj,
  cells = proj$cellNames[which(proj$Clusters2 %in% c("PTS1-1", "PTS1-2", "PTS2", "PTS3", "Injury"))],
  outputDirectory = "All_PT_Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

All_PT_Subset <- addHarmony(
  ArchRProj = All_PT_Subset,
  reducedDims = "IterativeLSI_RNA",
  name = "Harmony_RNA",
  groupBy = "Sample",
  force = TRUE
)

#UMAP
All_PT_Subset <- addUMAP(
  ArchRProj = All_PT_Subset, 
  reducedDims = "Harmony_RNA", 
  name = "UMAPHarmony_RNA", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

All_PT_Subset <- addUMAP(
  ArchRProj = All_PT_Subset, 
  reducedDims = "IterativeLSI_RNA", 
  name = "UMAP_RNA", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

All_PT_Subset <- addClusters(All_PT_Subset, reducedDims = "Harmony_RNA", name = "reClusters", resolution = 2, force = TRUE, maxClusters = 20)
plotEmbedding(
  All_PT_Subset,
  embedding = "UMAPHarmony_RNA",
  colorBy = "cellColData",
  name = "reClusters",
  labelAsFactors = FALSE
)

trajectory <- c( "C9", "C20","C14","C16","C11","C10","C8","C6","C7","C3","C4","C5")

All_PT_Subset <- addTrajectory(
  ArchRProj = All_PT_Subset, 
  name = "Trajectory", 
  groupBy = "reClusters",
  trajectory = trajectory, 
  embedding = "UMAPHarmony_RNA",
  force = TRUE
)


trajMM  <- getTrajectory(ArchRProj = All_PT_Subset, name = "Trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = All_PT_Subset, name = "Trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
trajPM  <- getTrajectory(ArchRProj = All_PT_Subset, name = "Trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelTop = 50)
plotPDF(p1, name = "TrajectoryHeatMap_MotifMatrix", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), labelTop = 70)
plotPDF(p2, name = "TrajectoryHeatMap_GeneExpressionMatrix", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), labelTop = 50)
plotPDF(p3, name = "TrajectoryHeatMap_PeakMatrix", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)

p4 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), labelTop = 1000)
plotPDF(p4, name = "TrajectoryHeatMap_GeneExpressionMatrix_extended", width = 8, height = 250, ArchRProj = All_PT_Subset, addDOC = FALSE)

corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1 + ht2 , name = "TrajectoryHeatMap_Combined", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)

p1 <- plotTrajectory(All_PT_Subset, embedding = "UMAPHarmony_RNA", 
                     trajectory = "Trajectory", colorBy = "GeneIntegrationMatrix", name = "POLR2A", continuousSet = "blueYellow")
p1[[1]]
p2 <- plotTrajectory(All_PT_Subset, embedding = "UMAPHarmony_RNA", 
                     trajectory = "Trajectory", colorBy = "GeneIntegrationMatrix", name = "KLF6", continuousSet = "blueYellow")
p2[[1]]
plotPDF(p1, name = "Trajectory_Polr2a", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)
plotPDF(p2, name = "Trajectory_Klf6", width = 8, height = 12, ArchRProj = All_PT_Subset, addDOC = FALSE)
