#velocity plots using velocyto.R
#05/30/19
#Cite Ni*,Lou*,Yao*, et al.
#ZIP1+ fibroblasts protect lung cancer against chemotherapy 
# via connexin-43 mediated intercellular Zn2+ transfer
#Created by Dekang Lv

#attach libraries
library(Seurat)
library(velocyto.R)
library(plyr)

####################################################
#Input Data
####################################################

combined <- readRDS("combined.RDS")

#.loom files are output of velocyto
ldat1 <- list(dox1=read.loom.matrices("DOX.loom"),
              pbs1=read.loom.matrices("PBS.loom"))
#merge datasets
matrix.name <- names(ldat[1])
ldat1 <- lapply(matrix.name, function(x){
  dat.list <- lapply(ldat1, function(y){
    y[[x]]
  })
  dat.merged <- do.call(cbind, dat.list)
  dat.merged
})
names(ldat1) <- matrix.name

#change cell name
ldat2 <- lapply(ldat1, function(x) {
  colnames(x) <- ifelse(substr(colnames(x),1, 3)=="PBS",
                        paste0(substr(colnames(x),5, 20),"1"),
                        paste0(substr(colnames(x),5, 20),"2"))
  x
})

#samping to save time
sample.num <- 2000
new.dat <- cbind(combined@reductions$umap@cell.embeddings,combined@meta.data)
set.seed(666)
new.dat <- new.dat[sample(1:nrow(new.dat),sample.num),]
new.dat.pbs <- new.dat[new.dat$stim == "PBS",]
new.dat.dox <- new.dat[new.dat$stim == "DOX",]


#extract counts 
ldat3 <- lapply(ldat2, function(x) {
  x[, rownames(new.dat)]
})
emat1 <- ldat3$spliced
nmat1 <- ldat3$unspliced
#select gene
emat2 <- filter.genes.by.cluster.expression(emat1, cell_cluster1, min.max.cluster.average = .5)
nmat2 <- filter.genes.by.cluster.expression(nmat1, cell_cluster1, min.max.cluster.average = .5)
#length(intersect(rownames(emat2), rownames(nmat2)))

#refer state manifolds and save as rvel.qf1
fit.quantile  <-  0.05 
deltaT  <-  1 # default: 1
kCells  <-  10
rvel.qf1 <- gene.relative.velocity.estimates(emat2, nmat2, 
                                             deltaT = deltaT, kCells = kCells, 
                                             fit.quantile = fit.quantile)

#Save Output
pdf("velocity_all.pdf")
##set point color for each cell
cell_cluster1 <- as.character(new.dat$seurat_clusters)
names(cell_cluster1) <- rownames(new.dat)
colors1 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7","#BEBEBE")
names(colors1) <- sort(unique(new.dat$seurat_clusters))
cell.colors1  <-  plyr::mapvalues(cell_cluster1, names(colors1), colors1)
n  <-  100 
scale  <-  "sqrt" 
cell.alpha  <-  0.2 
cell.cex  <-  1 
arrow.scale  <-  1 
arrow.lwd  <-  1.5 
grid.n  <-  50
show.velocity.on.embedding.cor(emb1, rvel.qf1, n, scale=scale, 
                               cell.colors = ac(cell.colors1, alpha = cell.alpha),
                               cex = cell.cex, arrow.scale = arrow.scale, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 1,
                               grid.n = grid.n, arrow.lwd = arrow.lwd)
dev.off()


#extract counts 
ldat3 <- lapply(ldat2, function(x) {
  x[, rownames(new.dat.pbs)]
})
emat1 <- ldat3$spliced
nmat1 <- ldat3$unspliced
#select gene
emat2 <- filter.genes.by.cluster.expression(emat1, cell_cluster1, min.max.cluster.average = .5)
nmat2 <- filter.genes.by.cluster.expression(nmat1, cell_cluster1, min.max.cluster.average = .5)
#length(intersect(rownames(emat2), rownames(nmat2)))

#refer state manifolds and save as rvel.qf1
fit.quantile  <-  0.05
deltaT  <-  1 # default: 1
kCells  <-  10
rvel.qf1 <- gene.relative.velocity.estimates(emat2, nmat2, 
                                             deltaT = deltaT, kCells = kCells, 
                                             fit.quantile = fit.quantile)

#Save Output
pdf("velocity_pbs.pdf")
#set point color for each cell
cell_cluster1 <- as.character(new.dat.pbs$seurat_clusters)
names(cell_cluster1) <- rownames(new.dat.pbs)
colors1 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7","#BEBEBE")
names(colors1) <- sort(unique(new.dat.pbs$seurat_clusters))
cell.colors1  <-  plyr::mapvalues(cell_cluster1, names(colors1), colors1)

n  <-  100 
scale  <-  "sqrt" 
cell.alpha  <-  0.2 
cell.cex  <-  1 
arrow.scale  <-  1 
arrow.lwd  <-  1.5 
grid.n  <-  50
show.velocity.on.embedding.cor(emb1, rvel.qf1, n, scale=scale, 
                               cell.colors = ac(cell.colors1, alpha = cell.alpha),
                               cex = cell.cex, arrow.scale = arrow.scale, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 1,
                               grid.n = grid.n, arrow.lwd = arrow.lwd)
dev.off()


#extract counts 
ldat3 <- lapply(ldat2, function(x) {
  x[, rownames(new.dat.dox)]
})
emat1 <- ldat3$spliced
nmat1 <- ldat3$unspliced
#select gene
emat2 <- filter.genes.by.cluster.expression(emat1, cell_cluster1, min.max.cluster.average = .5)
nmat2 <- filter.genes.by.cluster.expression(nmat1, cell_cluster1, min.max.cluster.average = .5)
#length(intersect(rownames(emat2), rownames(nmat2)))

#refer state manifolds and save as rvel.qf1
fit.quantile  <-  0.05
deltaT  <-  1 # default: 1
kCells  <-  10
rvel.qf1 <- gene.relative.velocity.estimates(emat2, nmat2, 
                                             deltaT = deltaT, kCells = kCells, 
                                             fit.quantile = fit.quantile)

#Save Output
pdf("velocity_dox.pdf")
#set point color for each cell
cell_cluster1 <- as.character(new.dat.dox$seurat_clusters)
names(cell_cluster1) <- rownames(new.dat.dox)
colors1 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7","#BEBEBE")
names(colors1) <- sort(unique(new.dat.dox$seurat_clusters))
cell.colors1  <-  plyr::mapvalues(cell_cluster1, names(colors1), colors1)

n  <-  100 
scale  <-  "sqrt" 
cell.alpha  <-  0.2 
cell.cex  <-  1 
arrow.scale  <-  1 
arrow.lwd  <-  1.5 
grid.n  <-  50
show.velocity.on.embedding.cor(emb1, rvel.qf1, n, scale=scale, 
                               cell.colors = ac(cell.colors1, alpha = cell.alpha),
                               cex = cell.cex, arrow.scale = arrow.scale, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 1,
                               grid.n = grid.n, arrow.lwd = arrow.lwd)
dev.off()


