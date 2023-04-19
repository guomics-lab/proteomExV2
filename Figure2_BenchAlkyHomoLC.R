# unique(sapply(sapply(list.files("."),strsplit,"_"),"[[",2))[1:3]
# ReductionAlkylations####
peptideFiles<-list.files(pattern="peptide.tsv","1_ReductionAlkylations/",full.names = T)
proteinFiles<-list.files(pattern="protein.tsv","1_ReductionAlkylations/",full.names = T)

condFullNames<-unique(sapply(sapply(list.files("1_ReductionAlkylations/"
                                               ,pattern = "tsv"),strsplit,"-"),"[[",1))
condNames<-unique(sapply(sapply(condFullNames,strsplit,"_"),"[[",1))
names(condFullNames)<-sapply(sapply(condFullNames,strsplit,"_"),"[[",1)
condCol<-hcl.colors(length(condNames), palette = "Dynamic", alpha = 0.9)
names(condCol)<-condNames
# "#DB9D85E6" "#ABB065E6" "#5CBD92E6" "#4CB9CCE6" "#ACA4E2E6" "#E093C3E6"
peptideMats<-list()
for(i in 1:length(peptideFiles)){
  peptideMats[[i]]<-read.delim(peptideFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}

proteinMats<-list()
for(i in 1:length(proteinFiles)){
  proteinMats[[i]]<-read.delim(proteinFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
names(peptideMats)<-names(proteinMats)<-condFullNames

df_ReducAkyl<-as.data.frame(matrix(NA,nrow=length(condFullNames),ncol=4))
row.names(df_ReducAkyl)<-condFullNames
colnames(df_ReducAkyl)<-c("NumPeptide","NumProtein","CysModPerc","CysPeptidePerc")
df_ReducAkyl$condName<-names(condFullNames)
df_ReducAkyl$Rep<-c(rep(c(1,1,2,2,3,3),4),rep(1,7))

df_ReducAkyl$NumPeptide<-as.vector(sapply(peptideMats,nrow))
df_ReducAkyl$NumProtein<-as.vector(sapply(proteinMats,nrow))

cysMod<-function(pepMat){
  cysModNum<-length(which(grepl("C\\(57.0214\\)",pepMat$Assigned.Modifications)))
  cysPepNum<-sum(grepl("C",row.names(pepMat)))
  return(list(round(cysModNum/cysPepNum*100,2),round(cysPepNum/nrow(pepMat)*100,2)))
}
df_ReducAkyl$CysModPerc<-unlist(sapply(peptideMats,cysMod)[1,])
df_ReducAkyl$CysPeptidePerc<-unlist(sapply(peptideMats,cysMod)[2,])
write.csv(df_ReducAkyl,"1_ReductionAlkylations.csv")

df_ReducAkyl$condName<-as.factor(df_ReducAkyl$condName)
par(mfrow=c(1,4),mar=c(8,4,2,2))
stripchart(NumPeptide~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,34000))
stripchart(NumProtein~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,3400))
stripchart(CysModPerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(50,100))
stripchart(CysPeptidePerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,15))

#noMLQC
df_ReducAkylNoML<-df_ReducAkyl[!grepl("QC",row.names(df_ReducAkyl)),]
df_ReducAkylNoML$condName<-as.vector(df_ReducAkylNoML$condName)
df_ReducAkylNoML$condName<-as.factor(df_ReducAkylNoML$condName)
df_ReducAkylNoML$condName<-factor(df_ReducAkylNoML$condName,
                                  levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1"))

Mean_NumPeptide<-aggregate(df_ReducAkylNoML$NumPeptide, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(df_ReducAkylNoML$NumProtein, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysModPerc<-aggregate(df_ReducAkylNoML$CysModPerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysPeptidePerc<-aggregate(df_ReducAkylNoML$CysPeptidePerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_NumPeptide)<-names(Mean_NumProtein)<-names(Mean_CysModPerc)<-names(Mean_CysPeptidePerc)<-c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")



df_ReducAkylNoML4<-df_ReducAkylNoML3<-df_ReducAkylNoML2<-df_ReducAkylNoML1<-df_ReducAkylNoML
df_ReducAkylNoML1$condName<-factor(df_ReducAkylNoML1$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_NumPeptide,decreasing = T)])
df_ReducAkylNoML2$condName<-factor(df_ReducAkylNoML2$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_NumProtein,decreasing = T)])
df_ReducAkylNoML3$condName<-factor(df_ReducAkylNoML3$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_CysModPerc,decreasing = T)])
df_ReducAkylNoML4$condName<-factor(df_ReducAkylNoML4$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_CysPeptidePerc,decreasing = T)])

pValueCal<-function(vectorIn,factorIn){
  # vectorIn=df_ReducAkylNoML3$NumPeptide
  # factorIn=df_ReducAkylNoML3$condName
  
  levelNames<-levels(factorIn)
  comb_level = combn(levelNames,2)
  combNames<-paste(t(comb_level)[,1],t(comb_level)[,2],sep="-")
  pvalues<-c()
  for(j in 1:ncol(comb_level)){
    vector1<-vectorIn[factorIn==comb_level[1,j]]
    vector2<-vectorIn[factorIn==comb_level[2,j]]
    pvalues[j]<-t.test(vector1,vector2)$p.value
  }
  names(pvalues)<-combNames
  return(pvalues)
}

pValueCalSort<-function(vectorIn,factorIn){
  # vectorIn=df_ReducAkylNoML3$CysModPerc
  # factorIn=df_ReducAkylNoML3$condName
  levelNames<-levels(factorIn)
  combNames<-paste(levelNames[1:(length(levelNames)-1)],levelNames[2:length(levelNames)])
  pvalues<-c()
  for(j in 1:(length(levelNames)-1)){
    vector1<-vectorIn[factorIn==levelNames[j]]
    vector2<-vectorIn[factorIn==levelNames[j+1]]
    pvalues[j]<-t.test(vector1,vector2)$p.value
  }
  names(pvalues)<-combNames
  return(pvalues)
}

pvalueStars<-function(pvalues){
  stars<-rep(NA,length(pvalues))
  stars[pvalues<=0.05&pvalues>0.01]<-"*"
  stars[pvalues<=0.01&pvalues>0.001]<-"**"
  stars[pvalues<=0.001&pvalues>0.0001]<-"***"
  stars[pvalues<=0.0001]<-"****"
  return(stars)
}

{
pdf("1_ReductionAlkylations_pvalue.pdf",width=8.27,height=2.71)
par(mfrow=c(1,4),mar=c(8,4,1,1))
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==1]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(94,100),col=condCol[levels(df_ReducAkylNoML3$condName)],ylab="CysModPerc")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_CysModPerc,decreasing = T),y1=sort(Mean_CysModPerc,decreasing = T))
sem3<-aggregate(df_ReducAkylNoML[,3], list(df_ReducAkylNoML3$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)-sem3)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysModPerc,decreasing = T)+sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==1]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(94,100),col=condCol[levels(df_ReducAkylNoML3$condName)],ylab="CysModPerc")
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==2]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML3$condName)])
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==3]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML3$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML3$CysModPerc,df_ReducAkylNoML3$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.01*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],y1=1.01*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.011*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])

stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==1]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,15),col=condCol[levels(df_ReducAkylNoML4$condName)],ylab="CysPeptidePerc")
sem4<-aggregate(df_ReducAkylNoML[,4], list(df_ReducAkylNoML4$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysPeptidePerc,decreasing = T)-sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)+sem4)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysPeptidePerc,decreasing = T)-sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)-sem4)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysPeptidePerc,decreasing = T)+sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)+sem4)
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==1]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,15),col=condCol[levels(df_ReducAkylNoML4$condName)],ylab="CysPeptidePerc")
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==2]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML4$condName)])
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==3]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML4$condName)])
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_CysPeptidePerc,decreasing = T),y1=sort(Mean_CysPeptidePerc,decreasing = T))
pvalues<-pValueCalSort(df_ReducAkylNoML4$CysPeptidePerc,df_ReducAkylNoML4$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])



stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==1]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==1],vertical=T,method = "jitter",las=2,cex=0,pch=15,ylim=c(0,30000),col=condCol[levels(df_ReducAkylNoML1$condName)]
           ,ylab="NumPeptide")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_NumPeptide,decreasing = T),y1=sort(Mean_NumPeptide,decreasing = T))
sem1<-aggregate(df_ReducAkylNoML[,1], list(df_ReducAkylNoML1$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_NumPeptide,decreasing = T)-sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)+sem1)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumPeptide,decreasing = T)-sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)-sem1)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumPeptide,decreasing = T)+sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)+sem1)
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==1]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,30000),col=condCol[levels(df_ReducAkylNoML1$condName)]
           ,ylab="NumPeptide")
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==2]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML1$condName)])
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==3]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML1$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML1$NumPeptide,df_ReducAkylNoML1$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])


stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==1]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,3500),col=condCol[levels(df_ReducAkylNoML2$condName)],ylab="NumProtein")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_NumProtein,decreasing = T),y1=sort(Mean_NumProtein,decreasing = T))
sem2<-aggregate(df_ReducAkylNoML[,2], list(df_ReducAkylNoML2$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_NumProtein,decreasing = T)-sem2,
         y1=sort(Mean_NumProtein,decreasing = T)+sem2)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumProtein,decreasing = T)-sem2,
         y1=sort(Mean_NumProtein,decreasing = T)-sem2)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumProtein,decreasing = T)+sem2,
         y1=sort(Mean_NumProtein,decreasing = T)+sem2)
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==1]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,3500),col=condCol[levels(df_ReducAkylNoML2$condName)],ylab="NumProtein")
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==2]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML2$condName)])
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==3]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML2$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML2$NumProtein,df_ReducAkylNoML2$condName)
stars<-pvalueStars(pvalues)
if(sum(is.na(stars))!=length(stars)){
  arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
  text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])
}
dev.off()
}

# Homogenization####
preMat_filename = '2_Homogenization_report.pr_matrix.tsv'
proMat_filename = '2_Homogenization_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[13],"_"))[c(2,5)],collapse = "_"))
}
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
pepAll<-unique(preMatIn$Stripped.Sequence)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll
for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}
Homo_pepNum<-apply(df_PepID,2,sum)

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

Homo_proNum<-apply(df_proID,2,sum)
Homo_Num<-as.data.frame(cbind(Homo_pepNum,Homo_proNum))
Homo_Num$GelName<-factor(sapply(sapply(row.names(Homo_Num),strsplit,"_"),"[[",1))
Homo_col<-hcl.colors(4, palette = "Warm", alpha = 0.9)
names(Homo_col)<-as.vector(levels(Homo_Num$GelName))


Mean_NumPeptide<-aggregate(Homo_Num$Homo_pepNum, list(Homo_Num$GelName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(Homo_Num$Homo_proNum, list(Homo_Num$GelName), FUN=mean)[,2]

{
pdf("2_Homogenization_pvalue.pdf",width=5,height=2.71)
par(mfrow=c(1,2),mar=c(4,4,1,1))
stripchart(Homo_Num$Homo_pepNum~Homo_Num$GelName,col=Homo_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,40000),cex=0)
segments(x0=seq(1,4)-0.1,x1=seq(1,4)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
sem1<-aggregate(Homo_Num[,1], list(Homo_Num$GelName), FUN=sd)[,2]
segments(x0=seq(1,4),x1=seq(1,4),y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide+sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide-sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumPeptide+sem1,
         y1=Mean_NumPeptide+sem1)
stripchart(Homo_Num$Homo_pepNum~Homo_Num$GelName,col=Homo_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,40000),add=T,cex=0.7)
pvalues<-pValueCal(Homo_Num$Homo_pepNum,Homo_Num$GelName)
stars<-pvalueStars(pvalues)
combSel<-combn(4,2)[,!is.na(stars)]
arrows(x0=combSel[1,],x1=combSel[2,],y0=1.5*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,y1=1.5*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,code=3,angle=90,length = 0)
text(x=(combSel[1,]+combSel[2,])/2,y=1.51*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,labels=stars[!is.na(stars)])



stripchart(Homo_Num$Homo_proNum~Homo_Num$GelName,col=Homo_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,5000),cex=0)
segments(x0=seq(1,4)-0.1,x1=seq(1,4)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
sem1<-aggregate(Homo_Num[,2], list(Homo_Num$GelName), FUN=sd)[,2]
segments(x0=seq(1,4),x1=seq(1,4),y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein+sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein-sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumProtein+sem1,
         y1=Mean_NumProtein+sem1)
stripchart(Homo_Num$Homo_proNum~Homo_Num$GelName,col=Homo_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,5000),add=T,cex=0.7)
pvalues<-pValueCal(Homo_Num$Homo_proNum,Homo_Num$GelName)
stars<-pvalueStars(pvalues)
combSel<-combn(4,2)[,!is.na(stars)]
arrows(x0=combSel[1,],x1=combSel[2,],y0=1.25*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,y1=1.25*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,code=3,angle=90,length = 0)
text(x=(combSel[1,]+combSel[2,])/2,y=1.26*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,labels=stars[!is.na(stars)])
dev.off()
}
# LC####
preMat_filename = '3_LC_gradients_report.pr_matrix.tsv'
proMat_filename = '3_LC_gradients_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[13],"_"))[c(3,5)],collapse = "_"))
}
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
pepAll<-unique(preMatIn$Stripped.Sequence)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll
for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}
LC_pepNum<-apply(df_PepID,2,sum)

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

LC_proNum<-apply(df_proID,2,sum)
LC_Num<-as.data.frame(cbind(LC_pepNum,LC_proNum))
LC_Num$GradLength<-factor(sapply(sapply(row.names(LC_Num),strsplit,"_"),"[[",1))
LC_col<-hcl.colors(4, palette = "Cold", alpha = 0.9)
names(LC_col)<-as.vector(levels(LC_Num$GradLength))

LC_Num$MeanPeak<-c(2.46,2.51,2.53,3.47,3.55,3.4)

Mean_NumPeptide<-aggregate(LC_Num$LC_pepNum, list(LC_Num$GradLength), FUN=mean)[,2]
Mean_NumProtein<-aggregate(LC_Num$LC_proNum, list(LC_Num$GradLength), FUN=mean)[,2]
Mean_MeanPeak<-aggregate(LC_Num$MeanPeak, list(LC_Num$GradLength), FUN=mean)[,2]

sem1<-aggregate(LC_Num[,1], list(LC_Num$GradLength), FUN=sd)[,2]
sem2<-aggregate(LC_Num[,2], list(LC_Num$GradLength), FUN=sd)[,2]
sem3<-aggregate(LC_Num[,4], list(LC_Num$GradLength), FUN=sd)[,2]


{
  pdf("3_LC_gradients_pvalue.pdf",width=5,height=2.71)
  par(mfrow=c(1,3),mar=c(4,4,1,1))
  stripchart(LC_Num$LC_pepNum~LC_Num$GradLength,col=LC_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,50000),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumPeptide-sem1,
           y1=Mean_NumPeptide+sem1)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide-sem1,
           y1=Mean_NumPeptide-sem1)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide+sem1,
           y1=Mean_NumPeptide+sem1)
  stripchart(LC_Num$LC_pepNum~LC_Num$GradLength,col=LC_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,50000),add=T,cex=0.7,at=c(1,2),xlim=c(0,3))
  pvalues<-pValueCal(LC_Num$LC_pepNum,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.2*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,y1=1.2*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.21*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,labels=stars[!is.na(stars)])
  
  
  stripchart(LC_Num$LC_proNum~LC_Num$GradLength,col=LC_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,6000),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumProtein-sem2,
           y1=Mean_NumProtein+sem2)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein-sem2,
           y1=Mean_NumProtein-sem2)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein+sem2,
           y1=Mean_NumProtein+sem2)
  stripchart(LC_Num$LC_proNum~LC_Num$GradLength,col=LC_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,6000),add=T,cex=0.7,at=c(1,2))
  pvalues<-pValueCal(LC_Num$LC_proNum,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.1*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,y1=1.1*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.11*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,labels=stars[!is.na(stars)])
  
  
  stripchart(LC_Num$MeanPeak~LC_Num$GradLength,col=LC_col,ylab="MeanPeakFWHM",vertical=T,method='jitter',pch=19,las=2,ylim=c(0,4),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MeanPeak,y1=Mean_MeanPeak)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MeanPeak-sem3,
           y1=Mean_MeanPeak+sem3)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MeanPeak-sem3,
           y1=Mean_MeanPeak-sem3)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MeanPeak+sem3,
           y1=Mean_MeanPeak+sem3)
  stripchart(LC_Num$MeanPeak~LC_Num$GradLength,col=LC_col,ylab="MeanPeakFWHM",vertical=T,method='jitter',pch=19,las=2,ylim=c(0,4),add=T,cex=0.7,at=c(1,2))
  pvalues<-pValueCal(LC_Num$MeanPeak,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.24*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,y1=1.24*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.25*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,labels=stars[!is.na(stars)])
  dev.off()
}
