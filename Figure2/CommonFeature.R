
#1. number of DEG/DAR/shared by female and male in different epxousres

data<-read.table("common_feature.txt",header = T)
a<-as.data.frame(table(data[,c(4,7)]))
write.table(a,"Number_of_common_feature_female_male.txt",sep="\t",quote=F,row.names = F)

#2. The common DEGs between female and male were identified in at least 3 exposures.
da1<-data[data$feature=="DEG",]
da2<-data[data$feature=="DAR",]
da3<-data[data$feature=="DMR",]
a1<-as.data.frame(table(da1$name));b1<-as.data.frame(table(a1$Freq))
a2<-as.data.frame(table(da2$name));b2<-as.data.frame(table(a2$Freq))
a3<-as.data.frame(table(da3$name));b3<-as.data.frame(table(a3$Freq))

l1<-a1[a1$Freq>=3,]
d<-data[data$name %in% l1$Var1,];
d1<-d[,c(1,3,4)];d2<-d[,c(5,3,4)]
library(tidyr)
f1<-spread(d1,type,Female);f2<-spread(d2,type,Male);
colnames(f1)<-c("Gene",paste0("female_",colnames(f1)[2:6]))
colnames(f2)<-c("Gene",paste0("male_",colnames(f2)[2:6]))
f<-merge(f1,f2,by="Gene",all=TRUE)
rownames(f)<-as.character(f$Gene);f<-f[,2:ncol(f)]
f[is.na(f)]<-0; k<-f[,c(1,6,2,7,3,8,4,9,5,10)]
ann<-as.data.frame(rep(c("female","male"),5))
colnames(ann)<-"sex";rownames(ann)<-colnames(k)
ann$exposure<-c("BPA10mg","BPA10mg","BPA10ug","BPA10ug","PM2.5-JHU","PM2.5-JHU","TBT","TBT","As","As")
pheatmap(k,cluster_cols = F,breaks = myBreaks,color=myColor,
         gaps_col = c(2,4,6,8),annotation_col = ann)

#3. calculate the correlation
data<-read.table("common_feature.txt",header = T)
ft<-as.character(data$feature);ft<-ft[!(duplicated(ft))];f<-c()
for(j in ft){
  da<-data[data$feature==j,]
  group<-as.character(da$type);group<-group[!(duplicated(group))];corExp<-c();corNum<-c()
  for(i in 1:length(group)){
    d<-da[da$type==group[i],];a1<-as.numeric(d$Female);a2<-as.numeric(d$Male);corExp[i]<-cor(a1,a2);corNum[i]<-length(a1)}
  e<-cbind(corExp,corNum);e<-as.data.frame(e);e$exp<-group;e$type<-j;f<-rbind(f,e)}
write.table(f,"correlation_sex.txt",sep="\t",quote = F,row.names = F)

f$type<-factor(f$type,levels = c("DEG","DAR","DMR"))
ggplot(f, aes(x=log2(corNum), y=corExp,color = exp,size=corExp))+
  geom_hline(yintercept=0, linetype="dashed",col="blue")+geom_vline(xintercept = 2, linetype="dashed",col="blue")+
  geom_point(size=1.5,alpha=0.8)+theme_classic()+facet_wrap(~type,ncol = 3,nrow = 1)+
  xlab("log2 Number of features shared between female and male") + ylab("correaltion between sex") +
  theme(legend.title = element_text(size=8), legend.text = element_text(size=8),legend.position="bottom", legend.box="vertical")

#4. distribution of common DETGs/DARs/DMRs in female and male
data<-read.table("common_feature.txt",header = T)
ggplot(data, aes(x=Female, y=Male,color = type))+geom_point()+facet_grid(feature~type)+theme_classic()+
  geom_hline(yintercept=0, linetype="dashed",col="blue")+
  geom_vline(xintercept = 0, linetype="dashed",col="blue")
  







