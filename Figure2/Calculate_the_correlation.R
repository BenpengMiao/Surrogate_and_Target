# calculate the correlation
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