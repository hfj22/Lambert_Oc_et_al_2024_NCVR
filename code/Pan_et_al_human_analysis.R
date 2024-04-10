
library(Seurat)
pan<-readRDS("pan_annotated.rds")
DimPlot(pan, reduction = "umap")
pan<-subset(pan, idents=c("VSMC 1", "VSMC 2", "VSMC 3", "Fibro-like VSMC"))
DefaultAssay(pan)<-"RNA"

pan.list <- SplitObject(pan, split.by = "experiment")
pan.list <- lapply(X = pan.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pan.list, nfeatures = 3000)
pan.list <- PrepSCTIntegration(object.list = pan.list, anchor.features = features)

pan.anchors <- FindIntegrationAnchors(object.list = pan.list, normalization.method = "SCT",
                                      anchor.features = features)
pan.combined.sct <- IntegrateData(anchorset = pan.anchors, normalization.method = "SCT")

pan.combined.sct <- RunPCA(pan.combined.sct, verbose = FALSE)
ElbowPlot(pan.combined.sct, ndims=50)

pct <- pan.combined.sct[["pca"]]@stdev / sum(pan.combined.sct[["pca"]]@stdev) * 100
pc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1]
pc

pan.combined.sct <- RunUMAP(pan.combined.sct, reduction = "pca", dims = 1:17)
DimPlot(pan.combined.sct, reduction = "umap", group.by = "experiment")

pan.combined.sct <- FindNeighbors(pan.combined.sct, reduction = "pca", dims = 1:17)
pan.combined.sct <- FindClusters(pan.combined.sct, resolution = 0.5)

DimPlot(pan.combined.sct, reduction = "umap", cols=c('pan1'='magenta','pan2'='blue','pan3'='purple'), group.by = "experiment",pt.size=0.1)

DefaultAssay(pan.combined.sct)<-"RNA"
pan.combined.sct<-NormalizeData(pan.combined.sct)
pan.combined.sct<-ScaleData(pan.combined.sct, features = rownames(pan.combined.sct))

FeaturePlot(pan.combined.sct, c("MYH11", "SPP1", "COL8A1", "TNFRSF11B", "VCAM1", "TIMP1","CCND1", "MKI67", "CRYAB"), ncol=3, pt.size=0.1)

pan.markers <- FindAllMarkers(pan.combined.sct, only.pos=TRUE)
pan.markers<-pan.markers[pan.markers$p_val_adj < 0.05,]
write.csv(pan.markers, file="pan markers.csv")
