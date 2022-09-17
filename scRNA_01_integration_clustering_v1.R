#integration of datasets from output of CellRanger and cell clusering
#05/30/19
#Cite Ni*,Lou*,Yao*, et al.
#ZIP1+ fibroblasts protect lung cancer against chemotherapy 
# via connexin-43 mediated intercellular Zn2+ transfer
#Created by Dekang Lv

#attach libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tibble)
library(tidyverse)
library(RColorBrewer)
library(glue)

####################################################
#Function
####################################################

# calculate density normalized to 1, independently for each facet variable
plot_density_split <- function(metadata_tbl, x_var, y_var, split_var, num_bins) {
  # ran into some issues with merging split geom_hex
  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    # geom_hex(aes(fill = stat(ndensity)), bins = num_bins) +
    stat_bin_2d(aes(fill = stat(ndensity)), bins = num_bins) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      strip.background = element_blank()
    ) +
    scale_fill_gradient2(low = "white", high = "darkred") +
    facet_wrap(vars(!!sym(split_var)))
}

# get table for density plot
get_density_diff_table <- function(metadata_tbl, x_var, y_var, split_var, num_bins){
  # generate a density plot split by stage
  density_plot = plot_density_split(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)
  # produce an object that can be rendered
  density_plot_tbl = ggplot_build(density_plot)
  # panel labels
  panels_tbl =
    tibble(
      PANEL = density_plot_tbl$layout$layout$PANEL,
      stage = density_plot_tbl$layout$layout$orig.ident
    )
  # merge panel contents and panel names
  density_tbl = density_plot_tbl$data[[1]]
  density_tbl = density_tbl %>% full_join(panels_tbl, by = "PANEL")
  return(density_tbl)
}

#density plot split by specified variable
#min_density = quantile(density_tbl$density, 0)
plot_density_diff <- function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", split_var, num_bins, group_pos, group_neg, interpolate = FALSE) {
  
  density_tbl = get_density_diff_table(metadata_tbl = metadata_tbl, 
                                       x_var = x_var, y_var = y_var, 
                                       split_var = split_var, num_bins = num_bins
  )
  min_density = quantile(density_tbl$density, 0)
  density_pos_tbl =
    density_tbl %>%
    filter(stage == group_pos) %>%
    select(x, y, cells_pos = count, density_pos = density)
  density_neg_tbl =
    density_tbl %>%
    filter(stage == group_neg) %>%
    select(x, y, cells_neg = count, density_neg = density)
  density_split_tbl = full_join(density_pos_tbl, density_neg_tbl, by = c("x", "y"))
  density_split_tbl[is.na(density_split_tbl)] = min_density
  density_split_tbl = density_split_tbl %>% mutate(density_diff = density_pos - density_neg)
  density_split_tbl = density_split_tbl %>% mutate(density_ratio = log(density_pos/density_neg))
  
  min_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.01)
  max_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.99)
  min_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.01)
  max_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.99)
  
  density_split_tbl =
    density_split_tbl %>%
    mutate(
      cells = cells_pos + cells_neg,
      log_density = log(density_pos + density_neg),
      density_ratio = if_else(density_ratio < min_density_ratio, min_density_ratio, density_ratio),
      density_ratio = if_else(density_ratio > max_density_ratio, max_density_ratio, density_ratio)
    ) %>%
    filter(cells > 0)
  
  ggplot(density_split_tbl, aes(x = x, y = y)) +
    # geom_tile(aes(fill = density_ratio)) +
    geom_raster(aes(fill = density_ratio), interpolate = interpolate) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    labs(title = glue("{group_pos} vs {group_neg}"), x = x_var, y = y_var) +
    scale_fill_gradient2(low = "#352a86", mid = "#41b899", high = "#f8fa0d")
}


####################################################
#Input Data
####################################################

set.min.cell <- 10
set.min.features <- 1000
#set.min.features <- 2000

set.nfeatures <- 1000
set.dims <- 20

set.resolution <- 0.2
#set.resolution <- 0.5

#directly load matrix data from the PBS output folder of CellRanger to dgCMatrix object
pbs <- Read10X(data.dir = "CellRanger output folder for PBS/filtered_feature_bc_matrix/")

#check the main quality of single cell transcriptomes
quantile(pbs@meta.data$nFeature_RNA)
quantile(pbs@meta.data$nCount_RNA)

#cellnames of the two dataset must have no overlap one for merging the two seurat object
#in the following codes
pbs@Dimnames[[2]] <- paste0(pbs@Dimnames[[2]],1)

#combine dgcMatrix with meta data to Seurat object
pbs <- CreateSeuratObject(counts = pbs, 
                          project = "tumor_pbs",
                          min.cells = set.min.cell,
                          min.features = set.min.features
)
pbs@meta.data$stim <- "PBS"

#normalize gene expression levels for each cell
pbs <- NormalizeData(pbs, verbose = FALSE)
#find variable genes in the dataset
pbs <- FindVariableFeatures(pbs, selection.method = "vst", nfeatures = set.nfeatures)

#load matrix data from the DOX output folder of CellRanger to Seurat object
dox <- Read10X(data.dir = "CellRanger output folder for DOX/filtered_feature_bc_matrix/")
quantile(dox@meta.data$nFeature_RNA)
quantile(dox@meta.data$nCount_RNA)

dox@Dimnames[[2]] <- paste0(dox@Dimnames[[2]],2)

dox <- CreateSeuratObject(counts = dox, 
                          project = "tumor_dox",
                          min.cells = set.min.cell,
                          min.features = set.min.features
)
dox@meta.data$stim <- "DOX"

dox <- NormalizeData(dox, verbose = FALSE)
dox <- FindVariableFeatures(dox, selection.method = "vst", nfeatures = set.nfeatures)

#integrate two datasets
#find cell anchors for dataset integration and integrate datasets
anchors <- FindIntegrationAnchors(object.list = list(pbs, dox), dims = 1:set.dims)
combined <- IntegrateData(anchorset = anchors, dims = 1:set.dims)

#DefaultAssay(combined) <- "integrated"
#scales and centers features in the dataset
combined <- ScaleData(combined, verbose = FALSE)

#run dimensionality reductions
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:set.dims)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:set.dims)

#computes the nearest neighbors for a given dataset
combined <- FindNeighbors(combined, reduction = "pca", k.param=20, dims = 1:set.dims)
#identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
combined <- FindClusters(combined, resolution = set.resolution)

#Save Output
saveRDS(combined,file = "combined.RDS")

#visualize clusters 
p1 <- DimPlot(combined, reduction = "umap", label = TRUE)
p2 <- DimPlot(combined, reduction = "umap", group.by = "stim")
p3 <- DimPlot(combined, reduction = "umap", split.by = "stim")

pg <- plot_grid(p1,p2,p3, ncol = 1)
#Save Output
ggsave(pg,filename = "clusters_umap.pdf")

freq <- as.data.frame(table(combined$stim,combined$seurat_clusters))
percentage=data.frame()
for(i in unique(freq$Var1)){
  freqi <- freq[freq$Var1 == i,]
  freqi$perc <- freqi[,3] / sum(freqi[,3]) * 100
  percentage <- rbind(percentage,freqi)
}
gp <- ggplot(percentage,aes(x = Var2,y = perc)) + 
  geom_bar(aes(fill = Var1),stat = "identity",position = "dodge") +
  scale_fill_manual(values =c("#7CAE00","#F8766D"))
#Save Output
ggsave(pg,filename = "percentage_clusters.pdf")

#visualize clusters density of DOX/PBS
new.dat <- cbind(combined@reductions$umap@cell.embeddings,combined@meta.data)
pdd <- plot_density_diff(new.dat, x_var = "UMAP_1", y_var = "UMAP_2", "orig.ident", 100, 
                         "tumor_dox", "tumor_pbs", interpolate = FALSE) 
#Save Output
ggsave(pdd,filename = "cell_density.pdf")
