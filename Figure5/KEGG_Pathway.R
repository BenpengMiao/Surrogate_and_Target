library(clusterProfiler)
library(org.Mm.eg.db)
mm<-org.Mm.eg.db

data<-readRDS("input_file") 
#Please get the data from here: 
data$new<-paste(as.character(data$lab),as.character(data$cond1),sep="_")
tp<-as.character(data$type)
tp<-tp[!(duplicated(tp))]
for(i in 1:length(tp)){
  dat<-data[data$type==tp[i],]
  con<-as.character(dat$new)
  con<-con[!(duplicated(con))]
  for(j in 1:length(con)){
    da<-dat[dat$new==con[j],]
    d<-da[da$sig=="TRUE",] 
    my.symbols<-d$gene
    IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
    Entrez<-unique(IDs[!is.na(IDs$ENTREZID),2])
    rP<-enrichKEGG(gene=Entrez,organism= 'mmu',pvalueCutoff = 0.05)
    Pathway<-rP@result;Pathway<-Pathway[Pathway$pvalue<0.05,]
    if(nrow(Pathway)>0){
      for(x in 1:nrow(Pathway)){
        symbols<- AnnotationDbi::select(mm, keys = strsplit(Pathway[x,"geneID"], "/")[[1]], columns = c("SYMBOL"), keytype = "ENTREZID")
        Pathway[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
      }
    }
    write.table(Pathway, file=paste0("KEGG_Pathway_", con[j],"_",tp[i],".txt"), 
                quote=F, sep="\t", row.names=F, col.names=T)
    
  }
}

#combine the results
files<-list.files(pattern="txt")
i<-1;data<-read.csv(files[i],sep="\t");data<-data[data$qvalue<0.1,]
data$Pathway<-paste(data$ID,data$Description,sep=",")
da<-data[,c("Pathway","pvalue")];colnames(da)<-c("Pathway",files[i])
d<-da
for(i in 2:length(files)){
  data<-read.csv(files[i],sep="\t");data<-data[data$qvalue<0.1,]
  data$Pathway<-paste(data$ID,data$Description,sep=",")
  da<-data[,c("Pathway","pvalue")];
  colnames(da)<-c("Pathway",files[i])
  d<-merge(d,da,by="Pathway",all=TRUE)
}
rownames(d)<-as.character(d$Pathway);d<-d[,2:ncol(d)];d[is.na(d)]<-1
write.table(d,"KEGG_Pathway_All_Qvalue_0.1.txt",sep="\t",quote=F)

#find the liver and blood shared pathway
data<-read.csv("KEGG_Pathway_All_Qvalue_0.1.txt",header = T,sep="\t")
n<-colnames(data);x<-c()
for(i in 1:length(n)){x[i]<-paste(strsplit(n[i],"_")[[1]][3:5],collapse = "_")}
x<-x[!(duplicated(x))];list<-c()
for(i in 1:length(x)){
  dat<-data[,grep(x[i],colnames(data))]
  da<-dat[rowSums(dat<1)>1,]
  list<-c(list,rownames(da))
}
list<-list[!(duplicated(list))];d<-data[list,]
write.table(d,"KEGG_Pathway_All_Qvalue_0.1_share_liver_blood.txt",sep="\t",quote=F)




