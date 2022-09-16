#pathway activity inference, limma t test and GO enrichment analysis
#05/30/19
#Cite Ni*,Lou*,Yao*, et al.
#ZIP1+ fibroblasts protect lung cancer against chemotherapy 
# via connexin-43 mediated intercellular Zn2+ transfer
# Created by Dekang Lv

#attach libraries
library(Seurat)
library(progeny)
library(msigdbr)
library(gsva)
library(limma)
library(GO.db)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pheatmap)
library(ggplot2)
library(openxlsx)

####################################################
#Function
####################################################

#calculate PROGENy pathway scores from gene expression for each  cluster
progeny.cluster <- function(seurat.obj, set.sample=FALSE,set.num=100, progeny.org="Mouse", progeny.top=500){
  RNAexpr = as.matrix(seurat.obj@assays$RNA@data)
  if(set.sample == TRUE){
    cellall=c()
    for(i in sort(unique(seurat.obj$seurat_clusters))){
      celli=which(seurat.obj$seurat_clusters==i)
      if(length(celli)>set.num){
        set.seed(666)
        celli100=sort(sample(celli,set.num))
      }else{celli100=celli}
      cellall=c(cellall,celli100)
    }
    str(cellall)
    RNAexpr=RNAexpr[,cellall]
  }
  prog_RNAall=progeny(RNAexpr,scale=F,organism = progeny.org, top=progeny.top)
  group_list <- data.frame(cell=rownames(prog_RNAall),subcluster=paste0("C",seurat.obj@meta.data[rownames(prog_RNAall),]$seurat_clusters))
  meandf=data.frame()
  for(i in sort(unique(group_list$subcluster))) {
    ci=prog_RNAall[group_list$cell[group_list$subcluster==i],]
    meanv=apply(ci,2,mean)
    meandf=rbind(meandf,meanv)
  }
  colnames(meandf)=names(meanv)
  rownames(meandf)=sort(unique(group_list$subcluster))
  return(meandf)
}

#ID transformation from gene symbol to entriz ID
IDtrans<-function(symbol){
  library(org.Mm.eg.db)
  #keytypes(org.Hs.eg.db)
  symbol2entrezID.df=AnnotationDbi::select(org.Mm.eg.db, keys = symbol, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
  symbol2entrezID.df=symbol2entrezID.df[!is.na(symbol2entrezID.df$ENTREZID),]
  symbol2entrezID.df=unique(symbol2entrezID.df)
  return(symbol2entrezID.df)
}

#infer pathway activity using GSVA for cells from a seurat object
GSVA_Hallmark<-function(seurat.obj,set.sample=FALSE,set.num=100){
  
  expr=as.matrix(seurat.obj@assays$RNA@data)
  if(set.sample == TRUE){
    cellall=c()
    for(i in sort(unique(seurat.obj$seurat_clusters))){
      celli=which(seurat.obj$seurat_clusters==i)
      if(length(celli)>set.num){
        set.seed(666)
        celli100=sort(sample(celli,set.num))
      }else{celli100=celli}
      cellall=c(cellall,celli100)
    }
    str(cellall)
    expr=expr[,cellall]
  }
  
  Hallmark=msigdbr(species = "Mus musculus", category = "H")
  Hallmarkdf=Hallmark[,c("entrez_gene","gs_name")]
  Hallmarkdf$entrez_gene=as.character(Hallmarkdf$entrez_gene)
  Hallmarklist=unstack(Hallmarkdf)
  HallmarkVsEntrezID=lapply(Hallmarklist,function(x){sort(unique(x))})
  
  symbol2entrezID=IDtrans(rownames(expr))
  symbol2entrezID=symbol2entrezID[!duplicated(symbol2entrezID$SYMBOL),]
  symbol2entrezID=symbol2entrezID[!duplicated(symbol2entrezID$ENTREZID),]
  expr_filter=expr[symbol2entrezID$SYMBOL,]
  rownames(expr_filter)=symbol2entrezID$ENTREZID
  
  gsvadf<-as.data.frame(GSVA::gsva(as.matrix(expr_filter),HallmarkVsEntrezID,method="ssgsea"))
  #gsvadf$Functional.Pathway <- sub("HALLMARK_","",rownames(gsvadf))
  gsvadf
}

#enrichment analysis for a set of gene symbol
enrichment <- function(pathIDVsEntrezID,deg.symbol,low=10,high=500){
  
  deg.trans.df<- IDtrans(deg.symbol)
  deg.trans.vec <- deg.trans.df$SYMBOL
  names(deg.trans.vec)<- deg.trans.df$ENTREZID
  queryset<-deg.trans.df$ENTREZID
  AllKEGG.EntrezID<-unique(unlist(pathIDVsEntrezID))
  geneset<-intersect(queryset,AllKEGG.EntrezID)
  pathway.length<-sapply(pathIDVsEntrezID,length)
  pathIDVsEntrezID<-pathIDVsEntrezID[pathway.length>=low & pathway.length<=high]
  statistics=lapply(pathIDVsEntrezID,function(x){
    refset=x
    overlap <- intersect(geneset,refset)
    overlap.entrez <-paste(overlap,collapse = ",")
    k <- length(overlap)
    M <- length(refset)
    N <- length(AllKEGG.EntrezID)
    n <- length(geneset)
    pvalues <-  phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    fold<- (k/n)/(M/N)
    c(N,M,n,k,fold,pvalues)
  })
  statistics.df<-as.data.frame(do.call("rbind",statistics))
  colnames(statistics.df)<-c("#genes.in.background",
                             "#genes.in.the.pathway",
                             "#genes.in.the.query",
                             "#genes.overlapped",
                             "Enrichment.fold",
                             "P.value")
  statistics.df$BH.FDR<-p.adjust(statistics.df$P.value,method = "BH")
  overlap.ID=lapply(pathIDVsEntrezID,function(x){
    refset=x
    overlapEntrezID <- intersect(geneset,refset)
    overlapSymbol <- deg.trans.vec[overlapEntrezID]
    overlap.entrez <- paste(overlapEntrezID,collapse = ",")
    overlap.symbol <- paste(overlapSymbol,collapse = ",")
    c(overlap.entrez,overlap.symbol)
  })
  overlap.df<-as.data.frame(do.call("rbind",overlap.ID))
  colnames(overlap.df)<-c("EntrezID.overlapped",
                          "Genes.overlapped")
  res.df<-cbind(statistics.df,overlap.df)
  res.df
}

#GO enrichment analysis for a set of gene symbol
enrich_GO<-function(GOIDVsEntrezID,deg.symbol,GOTERM){
  enrichment<-enrichment(GOIDVsEntrezID,deg.symbol)
  GO2term.vec <-GOTERM$Functional.Pathway
  names(GO2term.vec) <- GOTERM$go_id
  GO2ontology.vec  <-GOTERM$Ontology
  names(GO2ontology.vec) <- GOTERM$go_id
  enrichment$Functional.Pathway <- GO2term.vec[rownames(enrichment)]
  enrichment$Ontology <- GO2ontology.vec[rownames(enrichment)]
  enrichment<-enrichment[order(enrichment$Ontology),]
  enrichment
}

####################################################
#Input Data
####################################################

combined <- readRDS("combined.RDS")

#calculate PROGENy pathway scores from gene expression
progeny.cluster.activity <- progeny.cluster(combined, progeny.top=500)
progeny.cluster.activity <- progeny.cluster(combined, set.sample=T,progeny.top=500)
#Save Output
pdf("progeny_activity_cluster.pdf")
pheatmap::pheatmap(t(progeny.cluster.activity),scale = "row",
                   clustering_method="ward.D2",
                   border_color="white")
dev.off()

#infer hallmark pathway activity
gsvadf <- GSVA_Hallmark(seurat.obj,set.sample=T,set.num=100)

#limma t test analysis of pathway activity scores for each cluster
group_list <- data.frame(cell=colnames(gsvadf),
                         subcluster=paste0("C",combined@meta.data[colnames(gsvadf),]$seurat_clusters)
)
lst <- list()
lst2 <- list()
for (i in sort(unique(group_list$subcluster))) {
  group_df <- group_list
  group_df$subcluster <- ifelse(group_df$subcluster==i,group_df$subcluster,"else")
  design <- model.matrix(~0+subcluster,group_df)
  rownames(design) <- group_df$cell
  contrast.matrix <- makeContrasts(contrasts = paste0("subcluster",i,"-subclusterelse"), levels = design)
  fit <- lmFit(gsvadf, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  toptable <- topTable(fit2, coef= 1 ,number=Inf,adjust.method = "BH", sort.by="P")
  toptable[["pathway"]] <- rownames(toptable)
  toptable[["subcluster"]] <- paste0(i,"VSotherclusters")
  lst[[i]] <- toptable[,c("pathway","t","P.Value","adj.P.Val","subcluster")]
  lst2[[i]] <- toptable[,"t",drop=F]
}

dt <- do.call("rbind",lst)
dt$pathway <- sub("HALLMARK_","",dt$pathway)
#Save Output
write.csv(dt,file = "ttest_hallmarks.csv")

#GO enrichment analysis of GEP top 100 gene
#input file is output of cNMF 
gep6=t(read.delim("cNMF.gene_spectra_score.k_6.dt_0_1.txt"))[-1,]
top100=c()
for(n in 1:ncol(gep6)){
  top100gene=rownames(gep6[order(gep6[,n],decreasing = T),])[1:100]
  top100=cbind(top100,top100gene)
}
colnames(top100)=paste0("GEP",1:ncol(gep6))
write.xlsx(top100,"GEPtop100gene.xlsx")

#prepare gene ontology term and gene id
goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
GOTERM.df <- toTable(GOTERM)
GOTERM.df <- GOTERM.df[, c("go_id", "Term", "Ontology")]
GOTERM.df <- unique(GOTERM.df)
colnames(GOTERM.df)[2]="Functional.Pathway"
go2gene <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                 keys = names(goterms), 
                                 column = "ENTREZID", 
                                 keytype = "GOALL", 
                                 multiVals = "list")
go2gene <- go2gene[sapply(go2gene,function(x){!is.na(x[1])})]
GOIDVsEntrezID <- lapply(go2gene,unique)

#perform GO enrichment for top 100 genes of each GEP
wb <- createWorkbook()
for(i in 1:ncol(top100)){
  deg=top100[,i]
  #enrichment analysis for all go terms
  enrich.GO<-enrich_GO(GOIDVsEntrezID,deg,GOTERM.df)
  significant.GO=enrich.GO[enrich.GO$BH.FDR<=0.05,]
  sheetname=paste0("GEP",i,"_GO_significant")
  addWorksheet(wb,sheetName = sheetname)
  writeData(wb,sheet = sheetname, x = significant.GO,rowNames =T)
}
#Save Output
saveWorkbook(wb, "GEP_GOenrichment.xlsx", overwrite = TRUE)

