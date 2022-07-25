
load("prel_data_oriter_nov21.RData")
XX<-read.delim("distances_per_species.csv",header=FALSE,stringsAsFactors = FALSE,sep=",")
TX<-read.delim("bygenoriter_distNOMI.txt",header=FALSE,stringsAsFactors = FALSE,sep=" ")
ORI<-as.data.frame(ORI)
m<-match(XX$V1,rownames(ORI))
length(which(is.na(m)))
XX<-XX[-which(is.na(m)),]
XX<-cbind(XX,log2(1/XX$V3))
m<-m[-which(is.na(m))]

fx<-ORI$`pval>2`
fx[fx<0.01]<-0
fx[fx>=0.01]<-1

fx<-fx[m]


boxplot(XX$`log2(1/XX$V3)`~as.factor(fx))

t.test((XX$`log2(1/XX$V3)`~as.factor(fx)))





bygen<-NULL
ORITERNEW<-ORI[m,7]
names(ORITERNEW)<-rownames(ORI)[m]
#cor.test(XX$`log2(1/XX$V3)`[fx==0],ORITERNEW[fx==0])
#cor.test(XX$`log2(1/XX$V3)`[fx==1],ORITERNEW[fx==1])

#cor.test(ORITERNEW,XX$`log2(1/XX$V3)`)
modSpAll<-lm(ORITERNEW~XX$`log2(1/XX$V3)`)
summary(modSpAll)

w<-grep(x=names(ORITERNEW),pattern = "Staphylococcus")
modSpAllnoStaph<-lm(ORITERNEW[-w]~XX$`log2(1/XX$V3)`[-w])
summary(modSpAllnoStaph)
#cor.test(ORITERNEW[-w],XX$`log2(1/XX$V3)`[-w])




w<-grep(x=names(ORITERNEW),pattern = "Escherichia|Bordetella|Neisseria|Campylobacter")
cor.test(ORITERNEW[w],XX$`log2(1/XX$V3)`[w])




head(ORITERNEW)
nomi<-unlist(lapply(strsplit(x=names(ORITERNEW),split = "_"),function(x)x[1]))
un<-unique(nomi)
for(i in 1:length(un)){
  w<-grep(XX$V1,pattern = un[i])
  bygen<-rbind(bygen,c(mean(XX$`log2(1/XX$V3)`[w]),mean(ORITERNEW[w])))
  }
rownames(bygen)<-un
head(bygen)

cor.test(bygen[,1],bygen[,2])



m<-match(rownames(bygen),TX$V1)
TX<-TX[m,]
mod<-lm(bygen[,2]~bygen[,1])
summary(mod)

colnames(bygen)<-c("PS","oriter")
bygen<-as.data.frame(bygen)
modProteo<-lm(bygen$oriter[which(TX$V4==0)]~bygen$PS[which(TX$V4==0)])
summary(modProteo)

modAll<-lm(bygen$oriter~bygen$PS)
summary(modAll)


modGram<-lm(bygen$oriter[which(TX$V4==1)]~bygen$PS[which(TX$V4==1)])
summary(modGram)


R2<-c((summary(modSpAll))$adj.r.squared,(summary(modSpAllnoStaph))$adj.r.squared,(summary(modAll))$adj.r.squared,(summary(modProteo))$adj.r.squared,(summary(modGram))$adj.r.squared)
barplot(R2)







bygen<-bygen[-which(rownames(bygen)=="Staphylococcus"),]
modAllNoStaph<-lm(bygen$ori~bygen$PS)
summary(modAllNoStaph)


modGramNoStaph<-lm(bygen$ori[which(TX$V4==1)]~bygen$PS[which(TX$V4==1)])
summary(modGramNoStaph)

R2<-c((summary(modSpAll))$adj.r.squared,(summary(modSpAllnoStaph))$adj.r.squared,(summary(modAll))$adj.r.squared,(summary(modAllNoStaph))$adj.r.squared,(summary(modProteo))$adj.r.squared,(summary(modGram))$adj.r.squared,(summary(modGramNoStaph))$adj.r.squared)
barplot(R2)




# m<-match(tree_filt_root_d$tip.label,names(ORITERNEW))
# tree_filt_root_dF<- drop.tip(phy = tree_filt_root_d,
#                              tip = tree_filt_root_d$tip.label[which(is.na(m))])
# 
# 
# X<-XX$`log2(1/XX$V3)`
# names(X)<-names(ORITERNEW)
# 
# cm<-contMap(tree=tree_filt_root_dF,x = X)
# library(RColorBrewer)
# objp<-setMap(cm,invert=TRUE)
# plot(objp)
# 
# X<-ORITERNEW
# X[X>2]<-2
# cm<-contMap(tree=tree_filt_root_dF,x = X,plot=FALSE)
# 
# objp<-setMap(cm,invert=TRUE)
# plot(objp)



