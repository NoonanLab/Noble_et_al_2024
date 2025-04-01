# Pipeline for mapping non-targeting control (NTC) guide NSCs to fetal human cortex

#------------------------------------------------------------------------------------

# Load in fetal human cortex Seurat
human_cortex_seurat <- readRDS(file = "./../human_chimp_organoid_project/Eze_et_al_2021_analysis/New_Eze_analysis/Eze_human_cortex_object.rds")

# Downsample cell types to 2000 cells each
human_cortex_seurat$cell_type <- Idents(human_cortex_seurat)
human_cortex_seurat$barcode <- rownames(human_cortex_seurat@meta.data)
human_cortex_seurat.ds.list <- human_cortex_seurat@meta.data %>% group_by(cell_type) %>% sample_n(size = 2000) %>% ungroup() %>% select(barcode)
human_cortex_seurat.ds <- human_cortex_seurat[,human_cortex_seurat.ds.list$barcode]
human_cortex_seurat.ds <- RunUMAP(human_cortex_seurat.ds, reduction = "mnn", dims = 1:50, return.model = TRUE)
human_cortex_seurat.ds.markers <- FindAllMarkers(human_cortex_seurat.ds, only.pos = TRUE, min.pct = 0.25)

# Load in NTC Seurat and preprocess
Singles_NTC <- readRDS(file = "./Single_NTC.rds")
Singles_NTC <- DietSeurat(Singles_NTC)

# Merge reference and query
cortex.NTC.merge <- merge(x = human_cortex_seurat.ds, y = Singles_NTC, add.cell.ids = c('ref', 'query'))
cortex.NTC.merge <- FindVariableFeatures(cortex.NTC.merge, nfeatures = 3000)
cortex.NTC.merge.varfeatures <- cortex.NTC.merge@assays$RNA@var.features
cortex.NTC.merge$ref_or_query <- str_split_fixed(colnames(cortex.NTC.merge), "_", 2)[,1]
cortex.NTC.split <- SplitObject(cortex.NTC.merge, split.by = "ref_or_query")

# Split into reference and query datasets and reintegrate
cortex.NTC.split[["ref"]] <- NormalizeData(cortex.NTC.split[["ref"]])
cortex.NTC.split[["ref"]]@assays$RNA@var.features <- cortex.NTC.merge.varfeatures
cortex.NTC.split[["ref"]] <- ScaleData(cortex.NTC.split[["ref"]])
cortex.NTC.split[["ref"]] <- RunPCA(cortex.NTC.split[["ref"]])

cortex.NTC.split[["query"]] <- NormalizeData(cortex.NTC.split[["query"]])
cortex.NTC.split[["query"]]@assays$RNA@var.features <- cortex.NTC.merge.varfeatures
cortex.NTC.split[["query"]] <- ScaleData(cortex.NTC.split[["query"]])
cortex.NTC.split[["query"]] <- RunPCA(cortex.NTC.split[["query"]])

# Find transfer anchors and transfer data
cortex.NTC.anchors <- FindTransferAnchors(reference = cortex.NTC.split[["ref"]], query = cortex.NTC.split[["query"]], dims = 1:30, reference.reduction = "pca")
NTC.predictions <- TransferData(anchorset = cortex.NTC.anchors, refdata = cortex.NTC.split[["ref"]]$cell_type, dims = 1:30)
cortex.NTC.split[["query"]]$cortex.predictions <- NTC.predictions$predicted.id

# Incorporate the MNN integrated UMAP back into the object and project data
colnames(human_cortex_seurat.ds) <- paste("ref", colnames(human_cortex_seurat.ds), sep = "_")
cortex.NTC.split[["ref"]][["umap"]] <- human_cortex_seurat.ds@reductions$umap
cortex.NTC.split[["query"]] <- MapQuery(anchorset = cortex.NTC.anchors, reference = cortex.NTC.split[["ref"]], query = cortex.NTC.split[["query"]], refdata = list(cell_type = "cell_type"), reference.reduction = "pca", reduction.model = "umap")
cortex.NTC.split[["ref"]][["ref.umap"]] <- cortex.NTC.split[["ref"]][["umap"]]
cortex.NTC.remerged <- merge(x = cortex.NTC.split[["ref"]], y = cortex.NTC.split[["query"]], add.cell.ids = c("ref", "query"))

# Save reference and query datasets
saveRDS(cortex.NTC.split, file = "./cortex.NTC.split.rds")
cortex.NTC.split <- readRDS(file = "./cortex.NTC.split.rds")

# Merge reference and query datasets maintaining dimensionality reductions
cortex.NTC.merge <- merge(cortex.NTC.split[["ref"]], y = cortex.NTC.split[["query"]], merge.dr = TRUE, add.cell.ids = c("", ""))
cortex.NTC.merge$ref_or_query <- str_split_fixed(colnames(cortex.NTC.merge), "_", 3)[,2]

# Create DimPlots from the cell types and the predictions
cortex.NTC.cell_type_palette <- paletteer::paletteer_c("grDevices::Dynamic", n = 17)
names(cortex.NTC.cell_type_palette) <- unique(names(table(cortex.NTC.merge$cell_type)))
cortex.NTC.DimPlots <- list()
cortex.NTC.DimPlots[["cell_type_ref"]] <- DimPlot(cortex.NTC.merge, group.by = "cell_type", order = TRUE, reduction = "ref.umap") + scale_color_manual(values = cortex.NTC.cell_type_palette, na.value = "#B8B8B8") + th_cell_type_dimplot_theme + labs(title = "Developing Human Cortex Cell Types") + theme(plot.title = element_text(size = 14, family = "Helvetica"), legend.text = element_text(size = 10))
cortex.NTC.DimPlots[["cell_type_pred_query"]] <- DimPlot(cortex.NTC.merge, group.by = "cortex.predictions", order = TRUE, reduction = "ref.umap") + scale_color_manual(values = cortex.NTC.cell_type_palette, na.value = "#B8B8B8") + th_cell_type_dimplot_theme + labs(title = "Predicted Cell Types for NTC NSCs") + theme(plot.title = element_text(size = 14, family = "Helvetica"), legend.position = "none")

# Create pie chart showing the percentage of cell type predictions
cortex.NTC.pie_chart <- as.data.frame(table(cortex.NTC.merge$cortex.predictions)) %>% mutate(Freq = Freq/sum(Freq)) %>% ggplot() + geom_col(aes(x = "", y = Freq, fill = Var1)) + scale_fill_manual(values = cortex.NTC.cell_type_palette) + coord_polar(theta = "y") + theme_minimal() + labs(fill = "Cell Type", title = "Proportion of Cell Type Predictions in NTCs") + theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 16), legend.text = element_text(size = 12), legend.title = element_text(size = 14), text = element_text(family = "Helvetica"))

# Create dot plot showing representation of cell type predictions by cluster
cortex.NTC.merge.metadata <- cortex.NTC.merge@meta.data %>% select(predicted.cell_type.score, predicted.cell_type) %>% filter(!is.na(predicted.cell_type))
rownames(cortex.NTC.merge.metadata) <- str_split_fixed(rownames(cortex.NTC.merge.metadata), "_", n = 3)[,3]
cortex.NTC.merge.metadata$seurat_clusters <- Singles_NTC$seurat_clusters
cortex.NTC.merge.metadata <- cortex.NTC.merge.metadata %>% group_by(predicted.cell_type, seurat_clusters) %>% summarize(avg_cell_type_score = mean(predicted.cell_type.score), num_cells_in_group = n())

# Hierarchical clustering of Louvain clusters based on average cell type prediction score
cortex.NTC.merge.cluster_pivot <- cortex.NTC.merge.metadata %>% select(predicted.cell_type, seurat_clusters, avg_cell_type_score) %>% pivot_wider(id_cols = seurat_clusters, names_from = predicted.cell_type, values_from = avg_cell_type_score)
cortex.NTC.merge.cluster_hclust <- hclust(dist(cortex.NTC.merge.cluster_pivot))
cortex.NTC.merge.metadata$seurat_clusters <- fct_relevel(cortex.NTC.merge.metadata$seurat_clusters, as.character(cortex.NTC.merge.cluster_hclust$order-1))

# Creating dot plot showing representation of cell type predictions across sample identity
cortex.NTC.merge.metadata %>% ggplot() + geom_point(aes(x = seurat_clusters, y = predicted.cell_type, fill = avg_cell_type_score, size = num_cells_in_group), pch = 21, color = "black") + scale_fill_viridis_c() + labs(x = "Cluster Assignment", y = "Predicted Cell Type", size = "# Cells", fill = "Average Prediction Score") + theme(panel.background = element_blank(), axis.ticks.x = element_blank(), text = element_text(family = "Helvetica"), panel.grid = element_line(color = "#F8F8F8")) + scale_y_discrete(expand = c(0,0.4)) + scale_size_continuous(range = c(3, 10))

# Write predictions and scores to TSV
cortex.NTC.merge.metadata <- cortex.NTC.merge@meta.data %>% select(predicted.cell_type.score, predicted.cell_type) %>% filter(!is.na(predicted.cell_type))
rownames(cortex.NTC.merge.metadata) <- str_split_fixed(rownames(cortex.NTC.merge.metadata), "_", n = 3)[,3]
cortex.NTC.merge.metadata$seurat_clusters <- Singles_NTC$seurat_clusters
cortex.NTC.merge.metadata$rep <- Singles_NTC$rep
write.table(cortex.NTC.merge@meta.data, sep = "\t", file = "./NTC_cortex_preds_metadata.tsv")
