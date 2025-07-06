
# plot the blood and liver common DAR enriched TF binding motifs
library(ggplot2)
library(ggrepel)
data<-read.table("DAR_enriched_Common_motif_between_blood_and_liver.txt",header = T)
ggplot(data, aes(x=Blood, y=Liver,color = as.character(Share), shape=DAR))+
  geom_point(size=1.5)+theme_classic()+facet_grid(Sex~Type,scales = "free")+
  xlab("Percent of DARs with TF binding motif in Blood") + ylab("Percent of DARs with TF binding motif in Liver") +
  theme(legend.title = element_text(size=8), legend.text = element_text(size=8),axis.title=element_text(size=8),
        legend.position="bottom", legend.box="vertical")+
  geom_text_repel(aes(label = Class),size = 2.5)

#GO enrichment based on DARs with TF binding motif:
#show the DARs nearest genes enriched biology process from different groups
library(clusterProfiler)
library(enrichplot)
library("ggnewscale")
library(DOSE)
library(org.Mm.eg.db)
mm<-org.Mm.eg.db

files<-list.files(pattern = "peak.txt")
group<-c()
i<-1
data<-read.table(files[i],sep="\t",header = F)
my.symbols<-as.character(data$V1)
IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                           keytype = "SYMBOL")
Entrez<-unique(IDs[!is.na(IDs$ENTREZID),c(1,2)])

group$BPA10mg_Bl_M<-Entrez$ENTREZID

xx <- compareCluster(group, fun="enrichGO",OrgDb=mm, pvalueCutoff=0.05,ont="BP")
xx <- pairwise_termsim(xx)
emapplot(xx,legend_n = 3,cex_line = 0.4,cex_label_category = 0.4)

result<-xx@compareClusterResult
for(x in 1:nrow(result)){
  symbols<- AnnotationDbi::select(mm, keys = strsplit(result[x,"geneID"], "/")[[1]], 
                                  columns = c("SYMBOL"), keytype = "ENTREZID")
  result[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
}

write.table(result,"GO.txt",sep="\t",quote=F,row.names = F)


