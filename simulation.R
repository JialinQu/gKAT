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
np<-50
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

dp1<-NULL

tune<-seq(-1,1,0.1)
#tune<-seq(0,1,0.1)

runnum<-1000
for (i in 1:runnum) {
  print(i)
  
  X<-mvrnorm(nsamp,mu=rep(0,np),Sigma = sigma, tol = 0)
  
  # beta1<-matrix(c(rep(c(0,0.035),25),rep(0,50)))
  # beta1<-matrix(rep(c(0,0.035),25))
  # z1=X%*%beta1
  # y1=rnorm(nsamp,mean = z1,sd=1)
  
  #z1<-matrix(0,nrow=nsamp,ncol=1)
  #y1=rnorm(nsamp,mean = z1,sd=1)
  
  z1<-matrix(0,nrow=nsamp,ncol=1)
  y1<-matrix(0,nrow=nsamp,ncol=1)
  
  # for (k in 1:nsamp) {
  #   z1[k]<-(0.15*X[k,1]+0.15*X[k,3]+0.15*X[k,5]+0.15*X[k,7]+0.15*X[k,9] +0.6*X[k,1]*X[k,3]+0.6*X[k,5]*X[k,7]+ 0.6*X[k,9]*X[k,11]+0.6*X[k,13]*X[k,15]+0.6*X[k,17]*X[k,19])
  #   y1[k]<-rnorm(1,mean = z1[k],sd=1)
  # }
  
  
  for (k in 1:nsamp) {
    sum1=0
    for (s in 1:16) {
      sum1=sum1+0.01*(X[k,s])^2+0.015*(X[k,s])^3+exp(-X[k,s]^2/10)
    }
    z1[k]<-sum1
    y1[k]<-rnorm(1,mean = z1[k],sd=1)
  }
  
  
  #SKAT test
  obj1<-SKAT_Null_Model(y1~1, out_type = 'C',Adjustment = FALSE)
  
  #Gaussian kernel or linear kernel
  #gker=gausskernel(X,sigma =2*dim(X)[2])
  gker=X%*%t(X)
  SKATgker1<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = gker,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
  gp1<-rbind(gp1,SKATgker1$p.value)
  
  print(length(which(gp1<0.05)))
  
  #diffusion kernel combined with Gaussian kernel or linear kernel
  difp1<-NULL
  for (j in 1:length(tune)) {
    #dif=diag(1,nrow(H),ncol(H))+tune[i]*H+(tune[i]^2)/factorial(2)*H%*%H
    dif=R%*%diag(exp(-tune[j]*Lamda))%*%t(R)
    #difker<-Kernelsim(X,dif)
    difker<-X%*%dif%*%t(X)
    SKATdifker1<-SKAT.linear.Other(obj1$res,X,obj1$X1,kernel = difker,weights = NULL,obj1$s2,method = 'davies',obj1$res.out,obj1$n.Resampling)
    difp1<-rbind(difp1,SKATdifker1$p.value)
  }
  
  
  ccp1<-sum((1/length(difp1))*tan((0.5-difp1)*pi))
  p_cauchy1<-1-pcauchy(ccp1)
  #p_cauchy1<-min(difp1)
  dp1<-rbind(dp1,p_cauchy1)
  
  print(length(which(dp1<0.05)))
  
}


write.csv(gp1,file = 'gene50_400n_l1C.csv')
write.csv(dp1,file = 'gene50_400n_dl1C.csv')


print(length(which(gp1<0.05))/runnum)

print(length(which(dp1<0.05))/runnum)
