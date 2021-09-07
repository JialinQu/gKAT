library(ttutils)
library(Matrix)
library(SPAtest)
library(stringr)
library(MASS)
library(SKAT)
library(kernlab)
library(KRLS)
library("igraph")
library('EQL')

source('gKAT.R')

#define number of features and number of samples
np<-100
nsamp<-400

#create adjacency matrix
A<-matrix(1,nrow = np,ncol = np)
I=sample(1:np,100*np,replace = TRUE)
J=sample(1:np,100*np,replace = TRUE)
for (i in 1:np) {
  A[I[i],J[i]]<-0
}

diag(A)<-0

#plot network based on adjacency matrix
net=graph.adjacency(A, mode="undirected", weighted=NULL)
plot(net)

#calculate variance-covariance matrix and negative Laplacian matrix
sigma=similarity(net, method = "jaccard", mode="all")
D=diag(colSums(A))
L=D-A
H=-L

decom<-eigen(L)
R<-decom$vectors
Lamda<-decom$values

gp1<-NULL
gp11<-NULL
gp12<-NULL


dp1<-NULL
dp11<-NULL
dp12<-NULL


tune<-seq(-1,1,0.1)

runnum<-1000
for (i in 1:runnum) {
  print(i)
  
  X<-mvrnorm(nsamp,mu=rep(0,np),Sigma = sigma)
  
  # beta1<-matrix(c(rep(c(0,0.035),25),rep(0,50)))
  # #beta1<-matrix(rep(c(0,0.035),25))
  # z1=X%*%beta1
  # y1=rnorm(nsamp,mean = z1,sd=1)
  
  z1<-matrix(0,nrow=nsamp,ncol=1)
  y1<-matrix(0,nsamp,1)
  
  # for (i in 1:nsamp) {
  #   z1[i]<-(0.15*X[i,1]+0.15*X[i,3]+0.15*X[i,5]+0.15*X[i,7]+0.15*X[i,9] +0.6*X[i,1]*X[i,3]+0.6*X[i,5]*X[i,7]+ 0.6*X[i,9]*X[i,11]+0.6*X[i,13]*X[i,15]+0.6*X[i,17]*X[i,19])
  #   y1[i]<-rnorm(1,mean = z1[i],sd=1)
  # }
  
  
  for (i in 1:nsamp) {
    sum1=0
    for (j in 1:16) {
      sum1=sum1+0.01*(X[i,j])^2+0.015*(X[i,j])^3+exp(-X[i,j]^2/10)
    }
    z1[i]<-sum1
    y1[i]<-rnorm(1,mean = z1[i],sd=1)
  }
  
  
  
  #SKAT test
  obj1<-SKAT_Null_Model(y1~1, out_type = 'C',Adjustment = FALSE)
  
  #multiple kernels
  gker1=gausskernel(X,sigma =2*dim(X)[2])
  gker2=X%*%t(X)
  SKATgker11<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = gker1,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
  gp11<-rbind(gp11,SKATgker11$p.value)
  SKATgker12<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = gker2,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
  gp12<-rbind(gp12,SKATgker12$p.value)
  gp1<-1/2*tan((0.5-gp11)*pi)+1/2*tan((0.5-gp12)*pi)
  
  print(length(which(1-pcauchy(gp1)<0.05)))
  
  
  #diffusion kernel combined with multiple kernels
  difp11<-NULL
  difp12<-NULL
  for (i in 1:length(tune)) {
    dif=R%*%diag(exp(-tune[i]*Lamda))%*%t(R)
    difker1<-Kernelsim(X,dif)
    difker2<-X%*%dif%*%t(X)
    SKATdifker11<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = difker1,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
    difp11<-rbind(difp11,SKATdifker11$p.value)
    SKATdifker12<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = difker2,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
    difp12<-rbind(difp12,SKATdifker12$p.value)
  }
  
  
  ccp11<-sum((1/length(difp11))*tan((0.5-difp11)*pi))
  p_cauchy11<-1-pcauchy(ccp11)
  dp11<-rbind(dp11,p_cauchy11)
  
  ccp12<-sum((1/length(difp12))*tan((0.5-difp12)*pi))
  p_cauchy12<-1-pcauchy(ccp12)
  dp12<-rbind(dp12,p_cauchy12)
  
  dp1<-1/2*tan((0.5-dp11)*pi)+1/2*tan((0.5-dp12)*pi)
  
  print(length(which(1-pcauchy(dp1)<0.05)))
}



write.csv(gp1,file = 'gene50_400n_c1.csv')
write.csv(dp1,file = 'gene50_400n_dc1.csv')



c(length(which(1-pcauchy(gp1)<0.05))/runnum, length(which(1-pcauchy(dp1)<0.05))/runnum)
