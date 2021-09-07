library(KEGGgraph)
library(org.Hs.eg.db)
library(stringr)
library(MASS)
library(SKAT)
library(kernlab)
library(KRLS)
library(dplyr)

gene<-read.table('curatedExpressionLiver/expression.txt',header = TRUE)
outcome<-read.table('curatedPhenotype/phenotype.txt',header = TRUE)
indi<-read.table('curatedPhenotype/individuals.txt',header = TRUE)
feature<-read.table('curatedExpressionLiver/features.txt',fill = TRUE)


#combine gene expression data with response variable
length(which(indi$expression=='Yes')) # 466,we use all gene information
patientid<-colnames(gene)
outcomeid<-colnames(outcome)
comid<-patientid[patientid%in%outcomeid]
y<-outcome[15,]
index<-0
for (i in 1:466) {
  index<-c(index,which(outcomeid%in%comid[i]))
}

index<-index[-1]
comy<-y[index]

geney<-gene[,-1]
geney<-rbind(geney,comy)
rownames(geney)<-c(as.character(gene[,1]),'y')


#import gene pathway, match gene expression data
a<-system.file("extdata/hsa00982.xml",package = "KEGGgraph")
aGraph <- parseKGML2Graph(a, genesOnly=FALSE)
plot(aGraph)
aGraph

# To get the node gene names from pathway
outs <- sapply(edges(aGraph), length) > 0
ins <- sapply(inEdges(aGraph), length) > 0
ios <- outs | ins
ioGeneID <- translateKEGGID2GeneID(names(ios))
nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
nodesNames<-nodesNames[-which(is.na(nodesNames))]


comgene<-feature$V1[which(feature$V3%in%nodesNames)]
comgene<-as.character(comgene)

pathgeney<-geney[c(which(rownames(geney)%in%comgene),nrow(geney)),]
rowname<-rownames(pathgeney)
pathgeney<-as.data.frame(sapply(pathgeney, as.numeric))
rownames(pathgeney)<-rowname


#now gene expreesion data matches pathway, but one gene may mathes more than one id, we need to take average of those id
uniqcomgene<-unique(feature$V3[which(feature$V3%in%nodesNames)])
uniqcomgene<-as.character(uniqcomgene)
newobs <- NULL
for(i in uniqcomgene){
  temp =as.character(feature$V1[which(feature$V3==i)])
  temp.data = data.frame(pathgeney[which(rownames(pathgeney)%in%temp), ])
  tempobs <-colMeans(temp.data)
  newobs <- rbind(newobs, tempobs)
}

newobs<-rbind(newobs,pathgeney[nrow(pathgeney),])
rownames(newobs) <-c(uniqcomgene,'y')


#remove individual whose response variable is NA

remy<-which(is.na(comy))
remgeney<-newobs[,-remy]


#get pathway subgraph
## translate the KEGG IDs into Gene Symbol
if(require(org.Hs.eg.db)) {
  ioGeneID <- translateKEGGID2GeneID(names(ios))
  nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
} else {
  nodesNames <- names(ios)
}
names(nodesNames) <- names(ios)

comname<-nodesNames[nodesNames%in%uniqcomgene]

comname<-as.data.frame(comname)
rownames(comname)
mapkGsub=subGraph(rownames(comname),aGraph) 
plot(mapkGsub)
usub=as(mapkGsub, "matrix")
mapkGsub

Gsub=usub
for (i in 1:dim(usub)[1])
{
  for (j in 1:dim(usub)[1])
  {if (usub[j,i]==1) Gsub[i,j]=1}
}



#put the subset in order that is the same as the order of diffusion kernel
comname<-as.character(nodesNames[nodesNames%in%uniqcomgene])
orderdat<-remgeney
for (i in 1:(nrow(orderdat)-1)) {
  for (j in 1:(nrow(orderdat)-1)) {
    if(uniqcomgene[j]%in%comname[i])
      orderdat[i,]<-remgeney[j,]
  }
}

rownames(orderdat)<-c(comname,"y")


#remove individual whose gene expression contains NA
remindi<-NA
for (i in 1:nrow(orderdat)) {
  for (j in 1:ncol(orderdat)) {
    if(is.na(orderdat[i,j]))
      remindi<-c(remindi,j)
  }
}

remindi<-remindi[-1]
remindi<-unique(remindi)
finaldata<-orderdat[,-remindi]

write.csv(finaldata,file = '00982CYP2E1.csv')
write.csv(Gsub,file = 'Gsub00982.csv')
####################################################################
finaldata<-read.csv("result/hsa00982/00982CYP2E1/00982CYP2E1.csv",header = TRUE)
rownames(finaldata)<-finaldata[,1]
finaldata<-finaldata[,-1]

Gsub<-read.csv("result/hsa00982/00982Adjacency.csv")
rownames(Gsub)<-Gsub[,1]
Gsub<-Gsub[,-1]
Gsub<-as.matrix(Gsub)
####################################################################

D<-diag(colSums(Gsub))
L<-D-Gsub
H=-L

decom<-eigen(L)
R<-decom$vectors
Lamda<-decom$values

input<-as.matrix(finaldata)
remy0<-which(input[nrow(input),]==0)
input<-input[,-remy0]
input<-t(input)

Q1<-quantile(log(input[,ncol(input)]),1/4)
Q3<-quantile(log(input[,ncol(input)]),3/4)
iqr<-IQR(log(input[,ncol(input)]))
remyIQR<-which(log(input[,ncol(input)])<Q1-2*iqr | log(input[,ncol(input)])>Q3+2*iqr)
input<-input[-remyIQR,]

#SKAT test
a<-input[,-ncol(input)]
a<-scale(a)

logy<-log(input[,ncol(input)])
obj1<-SKAT_Null_Model(logy~1, out_type = 'C',Adjustment = FALSE)

#Gaussian kernel or linear kernel
lker=gausskernel(a,sigma =2*dim(a)[2])
lker=a%*%t(a)
SKATlker<-SKAT.linear.Other(obj1$res,a,obj1$X1,kernel = lker,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
SKATlker$p.value


#Diffuison kernel combined with Gaussian kernel or linear kernel
#tune<-seq(-2,2,0.1)
tune<-seq(-1,1,0.1)
difp<-NULL

for (i in 1:length(tune)) {
  dif=R%*%diag(exp(-tune[j]*Lamda))%*%t(R)
  difker<-Kernelsim(a,dif)
  #difker<-a%*%dif%*%t(a)
  SKATdifker<-SKAT.linear.Other(obj1$res,a,obj1$X1,kernel = difker,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
  difp<-rbind(difp,SKATdifker$p.value)
}

pvalue<-cbind(tune,difp)
difp[11]
min(difp)
tune[which(difp==min(difp))]

write.csv(pvalue,file = 'g00982CYP2E1P.csv')


plotdata <- data.frame('Bandwidth' = tune, "P_Value" = as.vector(difp))
p_min <- min(plotdata$P_Value)
plabel <- (plotdata$P_Value == p_min) + 1
plotdata[['label']] <- as.factor(plabel)
library(ggplot2)
myplot <- ggplot(plotdata, aes(x = Bandwidth, y = P_Value, color = label)) + geom_point()
show(myplot)


ccp<-sum((1/length(difp))*tan((0.5-difp)*pi))
p_cauchy<-1-pcauchy(ccp)
p_cauchy

1/2-atan(ccp)/pi

