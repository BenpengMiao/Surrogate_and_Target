# GO term enrichment

library(clusterProfiler)
library(enrichplot)
library("ggnewscale")
library(DOSE)
library(org.Mm.eg.db)
mm<-org.Mm.eg.db

#female
data<-read.table("Adult_female_Exposure_specific_DEG.txt",header = T)
k<-data[rowSums(data!=0)>0,];group<-c()
for(i in 1:ncol(k)){
  kk<-k[k[,i]!=0,]
  my.symbols<-as.character(rownames(kk))
  IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  Entrez<-unique(IDs[!is.na(IDs$ENTREZID),c(1,2)])
  group$x<-Entrez$ENTREZID
  names(group)<-colnames(k)[1:i]
}
xx <- compareCluster(group, fun="enrichGO",OrgDb=mm, pvalueCutoff=0.05,ont="BP")
xx <- pairwise_termsim(xx)
emapplot(xx,legend_n = 3,cex_line = 0.4,cex_label_category = 0.4)

result<-xx@compareClusterResult
for(x in 1:nrow(result)){
  symbols<- AnnotationDbi::select(mm, keys = strsplit(result[x,"geneID"], "/")[[1]], 
                                  columns = c("SYMBOL"), keytype = "ENTREZID")
  result[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
}
write.table(result,"Adult_female_Exposure_specific_DEG_GO.txt",sep="\t",quote=F,row.names = F)

#male
data<-read.table("Adult_male_Exposure_specific_DEG.txt",header = T)
k<-data[rowSums(data!=0)>0,];group<-c()
for(i in 1:ncol(k)){
  kk<-k[k[,i]!=0,]
  my.symbols<-as.character(rownames(kk))
  IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  Entrez<-unique(IDs[!is.na(IDs$ENTREZID),c(1,2)])
  group$x<-Entrez$ENTREZID
  names(group)<-colnames(k)[1:i]
}
xx <- compareCluster(group, fun="enrichGO",OrgDb=mm, pvalueCutoff=0.05,ont="BP")
xx <- pairwise_termsim(xx)
emapplot(xx,legend_n = 3,cex_line = 0.4,cex_label_category = 0.4)

result<-xx@compareClusterResult
for(x in 1:nrow(result)){
  symbols<- AnnotationDbi::select(mm, keys = strsplit(result[x,"geneID"], "/")[[1]], 
                                  columns = c("SYMBOL"), keytype = "ENTREZID")
  result[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
}
write.table(result,"Adult_male_Exposure_specific_DEG_GO.txt",sep="\t",quote=F,row.names = F)

