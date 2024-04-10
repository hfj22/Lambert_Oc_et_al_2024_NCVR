library(Seurat)

pan<-readRDS("apoe.integrated3.rds")
DimPlot(pan, reduction = "umap")
DimPlot(pan, reduction = "umap", split.by="orig.ident")
DefaultAssay(pan)<-"RNA"

pan.list <- SplitObject(pan, split.by = "orig.ident")
pan.list <- pan.list [-1]
rm(pan)

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

pan.combined.sct <- RunUMAP(pan.combined.sct, reduction = "pca", dims = 1:20)
DimPlot(pan.combined.sct, reduction = "umap", split.by = "orig.ident")

current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
new.cluster.ids <- c("old-0","old-1","old-2","old-3","old-4","old-5","old-6","old-7","old-8","old-9","old-10","old-11","old-12","old-13","old-14")
pan.combined.sct@meta.data[["seurat_clusters"]] <- plyr::mapvalues(x = pan.combined.sct@meta.data[["seurat_clusters"]], from = current.cluster.ids, to = new.cluster.ids)
pan.combined.sct@meta.data[["seurat_clusters"]] -> pan.combined.sct@meta.data[["seurat_clusters_old"]]

pan.combined.sct <- FindNeighbors(pan.combined.sct, reduction = "pca", dims = 1:20)
pan.combined.sct <- FindClusters(pan.combined.sct, resolution = 0.3)

DimPlot(pan.combined.sct, reduction = "umap", label = TRUE)
DimPlot(pan.combined.sct, reduction = "umap", group.by = "orig.ident", pt.size=0.1)

DefaultAssay(pan.combined.sct)<-"RNA"
pan.combined.sct<-NormalizeData(pan.combined.sct)
pan.combined.sct<-ScaleData(pan.combined.sct, features = rownames(pan.combined.sct))

FeaturePlot(pan.combined.sct, c("Myh11", "Spp1", "Col8a1", "Tnfrsf11b", "Vcam1", "Timp1","Ccnd1", "Mki67", "Cryab"), ncol=3, pt.size=0.1)

pan.markers <- FindAllMarkers(pan.combined.sct, only.pos=TRUE)
pan.markers<-pan.markers[pan.markers$p_val_adj < 0.05,]
write.csv(pan.markers, file="pan markers.csv")

pan.markers.selected2 <- FindMarkers(pan.combined.sct, only.pos=TRUE, ident.1 = c("3","0"), ident.2="1")
pan.markers.selected2<-pan.markers.selected2[pan.markers.selected2$p_val_adj < 0.05,]
write.csv(pan.markers.selected2, file="pan markers selected2.csv")

pan.markers.selected2.all <- FindMarkers(pan.combined.sct, ident.1 = c("3","0"), ident.2="1", only.pos=FALSE, min.pct = -Inf, logfc.threshold = -Inf)
write.csv(pan.markers.selected2.all, file="pan markers selected2 all.csv")

