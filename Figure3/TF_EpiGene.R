#1. number of TF for each exposures
#female
data<-read.table("DEG_Female_adult_TF_family.txt",header = T)
data$TF<-ifelse(data$Family=="zf-C2H2"|data$Family=="bHLH"|data$Family=="TF_bZIP"|data$Family=="Homeobox", data$Family,"Other")
a<-as.data.frame(table(data[,c(11,12)]))
a$TF<-factor(a$TF,levels = c("zf-C2H2","bHLH","TF_bZIP","Homeobox","Other"))
ggplot(a, aes(x=group, y=Freq,fill=TF)) + geom_bar(stat = "identity",position = "stack") + 
  theme_classic()+ scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8,vjust = 0.5), 
        strip.text.x = element_text(size=6), legend.title = element_text(size = 8))

#male
data<-read.table("DEG_Male_adult_TF_family.txt",header = T)
data$TF<-ifelse(data$Family=="zf-C2H2"|data$Family=="bHLH"|data$Family=="TF_bZIP"|data$Family=="Homeobox", data$Family,"Other")
a<-as.data.frame(table(data[,c(11,12)]))
a$TF<-factor(a$TF,levels = c("zf-C2H2","bHLH","TF_bZIP","Homeobox","Other"))
ggplot(a, aes(x=group, y=Freq,fill=TF)) + geom_bar(stat = "identity",position = "stack") + 
  theme_classic()+ scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8,vjust = 0.5), 
        strip.text.x = element_text(size=6), legend.title = element_text(size = 8))


#2. dot and verticle line show the TF C2H2
data<-read.table("DEG_TF_blood_C2H2.txt",header = T)
dat<-data[data$type=="Bartolomei_BPA10mg"|data$type=="Biswal_PM2.5"|data$type=="Walker_TBT"|data$type=="Zhibin_As",]
t<-as.character(dat$type);t<-t[!(duplicated(t))];f<-c()
for(i in 1:length(t)){
  a<-dat[dat$type==t[i],];b<-as.data.frame(table(a$gene));c<-as.character(b[b$Freq>1,]$Var1);a$share<-a$gene %in% c;f<-rbind(f,a)
}
ggplot(f, aes(x=gene, y=log2FoldChange,col=share)) +
  geom_point(size=0.7) +
  geom_hline(yintercept=0) +
  geom_linerange(aes(x=gene, ymax=log2FoldChange, ymin=0))+
  facet_grid(type~sex)

#3. TF Gene expression
data<-read.table("Gene_CPM_table_alltissue.txt",header = T)
list<-c("Bcl11a","Gfi1b","Tal1","Gata2")
sample<-colnames(data);s<-sample[grep("Bl",sample)];s<-s[grep("adt",s)];dat<-data[,s]
da<-dat[list,];da$Gene<-rownames(da);f<-gather(da,sample,CPM,1:(ncol(da)-1))
n<-as.character(f$sample);lab<-c();trt<-c();sex<-c();for(j in 1:length(n)){
  lab[j]<-strsplit(as.character(n[j]),"_")[[1]][1];trt[j]<-strsplit(as.character(n[j]),"_")[[1]][3]
  if(trt[j] != "Ctrl"){if((lab[j] == "BI") || (lab[j] == "MU")){trt[j]<-paste(lab[j],trt[j],sep=".")}}
  sex[j]<-strsplit(as.character(n[j]),"_")[[1]][6]};f$type<-trt;f$sex<-sex;

f$type<-factor(f$type,levels = c("Ctrl","As","BPA10mg","BPA10ug","DEHP","Pb","MU.PM2.5","BI.PM2.5","TBT","TCDD"))
ggplot(f, aes(type, CPM,fill=type))+scale_fill_brewer(palette = "Set3") +
  geom_boxplot(outlier.size = 2)+theme_classic()+facet_wrap(Gene~sex,scales = "free",ncol = 2,nrow = 4)+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1,vjust = 0.5),axis.text.y = element_text(size=10))

group<-c("As","BPA10mg","BPA10ug","DEHP","Pb","MU.PM2.5","BI.PM2.5","TBT","TCDD")
g1<-c();g2<-c();g3<-c();g4<-c()
for(l in list){for(i in c("M","F")){for(j in group){
  c1<-as.numeric(f[f$Gene==l & f$sex==i & f$type=="Ctrl",]$CPM);
  c2<-as.numeric(f[f$Gene==l & f$sex==i & f$type==j,]$CPM)
  a<-t.test(c1,c2)$p.value
  g1<-c(g1,l);g2<-c(g2,i);g3<-c(g3,j);g4<-c(g4,a)}}}
g<-as.data.frame(cbind(g1,g2,g3,g4));colnames(g)<-c("Gene","Sex","Type","TtestPvalue")
write.table(g,"T-test-Pvalue.txt",sep="\t",quote=F,row.names = F)

#4. EpiGene
data<-read.csv("DEG_adult_EpiGene_inf.txt",sep="\t")
rownames(data)<-as.character(data$MGI_symbol);d<-data[,2:15];
ann<-data[,c(18,19)];colnames(ann)<-c("Function","Modification")
x<-pheatmap(d);d<-d[x$tree_row$order,x$tree_col$order]
d<-d[,c(grep("Female",colnames(d)),grep("Male",colnames(d)))]
ann<-ann[x$tree_row$order,];d<-cbind(d,ann);a<-as.data.frame(table(ann$Function))
a<-a[order(a$Freq,decreasing = T),];name<-as.character(a$Var1)
k<-c()
for(i in 1:length(name)){
  x<-d[d$Function==name[i],];y<-pheatmap(x[,1:(ncol(x)-2)],cluster_cols = F);
  z<-x[y$tree_row$order,];k<-rbind(k,z)
}
k$Female_sum<-rowSums(k[,1:7]!=0);k$Male_sum<-rowSums(k[,8:14]!=0)
k$sum<-rowSums(k[,17:18]>0)


