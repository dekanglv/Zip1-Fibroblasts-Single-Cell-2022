#identify markers and compute functional profiles of each gene cluster
#05/30/19
#Cite Ni*,Lou*,Yao*, et al.
#ZIP1+ fibroblasts protect lung cancer against chemotherapy 
# via connexin-43 mediated intercellular Zn2+ transfer
# Created by Dekang Lv

#attach libraries
library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)

####################################################
#Input Data
####################################################

combined <- readRDS("combined.RDS")

#calculate the average gene expression of each cell cluster
cluster_mean_expression=data.frame()
meta=combined@meta.data
matexpr=as.matrix(combined@assays$integrated@data)

for(clust in sort(unique(combined$seurat_clusters))){
  print(clust)
  cells=rownames(meta)[meta$seurat_clusters==clust]
  exprmat=matexpr[,cells]
  rowmean=rowMeans(exprmat)
  if(ncol(cluster_mean_expression)==0){
    cluster_mean_expression=as.data.frame(rowmean)
  }else{
    cluster_mean_expression=cbind(cluster_mean_expression,rowmean)  
  }
}
colnames(cluster_mean_expression)=sort(unique(combined$seurat_clusters))

#find the cluster with the largest average expression
maxcell=apply(cluster_mean_expression,1,function(x){
  colnames(cluster_mean_expression)[order(x,decreasing = T)[1]]
})
cluster_mean_expression$maxcell=maxcell

#Save Output
write.csv(cluster_mean_expression,file = "cluster_mean_expression.csv")


#method1 to select markers for each cluster
combined@active.ident = combined$seurat_clusters
#set only.pos=TRUE
diff.expr.genes <- Seurat::FindAllMarkers(combined,only.pos = T,logfc.threshold = 0.25) 
#filter markers with p_val_adj>0.05
diff.expr.genes = diff.expr.genes[diff.expr.genes$p_val_adj<=0.05,]
#diff.expr.genes =diff.expr.genes[diff.expr.genes$avg_logFC>=0.58 &diff.expr.genes$p_val_adj<0.05,]
#select markers have the highest mean expression level
diff.expr.genes$maxcell=maxcell[diff.expr.genes$gene]
diff.expr.genes = diff.expr.genes[as.character(diff.expr.genes$cluster)==diff.expr.genes$maxcell,]
#Save Output
write.csv(diff.expr.genes,file = "clustermarker.csv")

#method2 to select markers for each cluster
diff.expr.genes2 <- Seurat::FindAllMarkers(combined,only.pos = T,logfc.threshold = 0.25,test.use = "roc")
diff.expr.genes2 = diff.expr.genes2[diff.expr.genes2$myAUC>0.7,]
diff.expr.genes2=diff.expr.genes2[ order(diff.expr.genes2$myAUC,decreasing=T),]
diff.expr.genes2=diff.expr.genes2[ !duplicated(diff.expr.genes2$gene),]
diff.expr.genes2=diff.expr.genes2[order(diff.expr.genes2$cluster),]
#Save Output
write.csv(diff.expr.genes2,file = "clustermarker.csv")

#find cluster markers within cluster 0-7 from 12 clusters
conservemarkers_0_7<-data.frame()
for(i in 0:7){
  imarkers<-FindConservedMarkers(oldcombined, ident.1 = i, ident.2=setdiff(0:7,i),
                                 grouping.var = "stim", verbose = FALSE
                                 )
  imarkers<-cbind(gene=rownames(imarkers),cluster=i,imarkers)
  #print(nrow(imarkers))
  conservemarkers_0_7<-rbind(conservemarkers_0_7,imarkers)
}
write.csv(conservemarkers_0_7,file="clusterconservemarkers_0_7.csv",row.names = F)

#perform enrichment analysis with markers in conservemarkers_0_7
unique.marker<- unique(conservemarkers_0_7$gene)
eg = bitr(unique.marker, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
conservemarkers_0_7$SYMBOL<- conservemarkers_0_7$gene
conservemarkers_0_7.entriz<-merge(conservemarkers_0_7,eg,by="SYMBOL")
entrizlist<-list()
for(i in sort(unique(conservemarkers_0_7.entriz$cluster))){
  name=paste0("cluster",i)
  entrizlist[[name]]<-conservemarkers_0_7.entriz$ENTREZID[conservemarkers_0_7.entriz$cluster==i]
}
cc <- compareCluster(geneCluster = entrizlist, fun = "enrichKEGG", OrgDb="org.Mm.eg.db")
#Save Output
pdf("keggEnrichment.pdf")
dotplot(cc)
dev.off()

