SKAT.logistic.Other = function(res, Z, X1, kernel , weights = NULL, pi_1, method,res.out,n.Resampling){
  
  n = nrow(Z) 
  m = ncol(Z)   
  
  # If m >> p and ( linear or linear.weight) kernel than call 
  # Linear function
  
  if (class(kernel) == "matrix") {
    K = kernel
  } 
  
  Q = t(res)%*%K%*%res/2
  Q.res = NULL
  
  D  = diag(pi_1)   
  gg = X1%*%solve(t(X1)%*%(X1 * pi_1))%*%t(X1 * pi_1)  ### Just a holder... not all that useful by itself
  P0 = D-(gg * pi_1)      ### This is the P0 or P in Zhang and Lin
  # P0 = D-D%*%gg  
  
  if(method == "davies"){
    P0_half = Get_Matrix_Square.1(P0)
    #print(dim(P0_half))
    W1 = P0_half %*% K %*% t(P0_half)
  } else {
    #W    = P0%*%K
    W = K * pi_1 - (X1 *pi_1) %*%solve(t(X1)%*%(X1 * pi_1))%*% ( t(X1 * pi_1) %*% K) 
  }
  
  if( method == "davies" ){
    out<-Get_Davies_PVal(Q, W1, Q.res)    
  } 
  
  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling, Test.Type = method, Q = Q, param=out$param )  
  return(re)
  
}

##############################################################
Get_Matrix_Square.1<-function(A){
  
  out<-eigen(A,symmetric=TRUE)
  ID1<-which(out$values > 0)
  if(length(ID1)== 0){
    stop("Error to obtain matrix square!")
  }
  out1<-t(out$vectors[,ID1]) * sqrt(out$values[ID1])
  return(out1)
}
#############################################################
Get_Davies_PVal<-function(Q, W, Q.resampling = NULL){
  
  K<-W/2
  
  Q.all<-c(Q,Q.resampling)
  
  re<-Get_PValue(K,Q.all)
  param<-list()
  param$liu_pval<-re$p.val.liu[1]
  param$Is_Converged<-re$is_converge[1]
  
  
  p.value.resampling = NULL
  if(length(Q.resampling) > 0){
    p.value.resampling<-re$p.value[-1]
    param$liu_pval.resampling<-re$p.val.liu[-1]
    param$Is_Converged.resampling<-re$is_converge[-1]
    
  }
  
  
  re<-list(p.value = re$p.value[1], param=param,p.value.resampling = p.value.resampling
           , pval.zero.msg=re$pval.zero.msg )  
  return(re)
}

###############################################################
Get_PValue<-function(K,Q){
  
  lambda<-Get_Lambda(K)
  re<-Get_PValue.Lambda(lambda,Q)
  return(re)
}

##############################################################
#only focus on positive eigenvalue
Get_Lambda<-function(K){
  
  out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)
  #print(out.s$values)
  
  #out.s1<-eigen(K,symmetric=TRUE)
  #print(out.s1$values)
  
  lambda1<-out.s$values
  IDX1<-which(lambda1 >= 0)
  
  # eigenvalue bigger than sum(eigenvalues)/1000
  IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
  #cat("Lambda:", lambda1, "\n")
  #K1<<-K
  
  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda<-lambda1[IDX2]
  return(lambda)
  
}
################################################################
Get_PValue.Lambda<-function(lambda,Q){
  
  #print(lambda)
  n1<-length(Q)
  
  p.val<-rep(0,n1)
  p.val.liu<-rep(0,n1)
  is_converge<-rep(0,n1)
  p.val.liu<-Get_Liu_PVal.MOD.Lambda(Q, lambda)
  
  for(i in 1:n1){
    out<-SKAT_davies(Q[i],lambda,acc=10^(-6))
    
    p.val[i]<-out$Qq
    #p.val.liu[i]<-SKAT_liu(Q[i],lambda)
    
    is_converge[i]<-1
    
    # check convergence
    if(length(lambda) == 1){
      p.val[i]<-p.val.liu[i]
    } else if(out$ifault != 0){
      is_converge[i]<-0
    }
    
    # check p-value
    if(p.val[i] > 1 || p.val[i] <= 0 ){
      is_converge[i]<-0
      p.val[i]<-p.val.liu[i]
    }
  }
  
  p.val.msg = NULL
  p.val.log=NULL
  #cat(p.val[1])
  if(p.val[1] == 0){
    
    param<-Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
    p.val.log<-Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p=TRUE)[1]
    
  }
  
  return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, p.val.log=p.val.log, pval.zero.msg=p.val.msg))
  
}
########################################################
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  
  return(p.value)
  
}
##############################################################
Get_Liu_PVal.MOD.Lambda.Zero<-function(Q, muQ, muX, sigmaQ, sigmaX, l, d){
  
  
  Q.Norm<-(Q - muQ)/sigmaQ
  Q.Norm1<-Q.Norm * sigmaX + muX
  
  temp<-c(0.05,10^-10, 10^-20,10^-30,10^-40,10^-50, 10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  #qchisq(temp, df=1000000000,lower.tail=FALSE)	
  out<-qchisq(temp,df = l,ncp=d, lower.tail=FALSE)
  #cat(c(Q.Norm1,l,d, out))
  #cat("\n")
  IDX<-max(which(out < Q.Norm1))
  
  pval.msg<-sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)
  
}

##################################################################
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}

#######################################
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation
  
  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }
  
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}








#########################################################################
#continuous y 
############################################################################
SKAT.linear.Other = function(res,Z,X1, kernel, weights = NULL, s2, method,res.out,n.Resampling){
  
  
  n<-nrow(Z)
  m = ncol(Z) 
  if (is.matrix(kernel)) {
    
    K = kernel
    
  }
  
  
  Q = t(res)%*%K%*%res/(2*s2)
  
  Q.res = NULL
  if(n.Resampling > 0){
    Q.res<-rep(0,n.Resampling)
    for(i in 1:n.Resampling){
      
      Q.res[i] = t(res.out[,i])%*%K%*%res.out[,i]/(2*s2)
    }
  }
  
  
  W = K - X1%*%solve( t(X1)%*%X1)%*%( t(X1) %*% K)	# W = P0 K
  
  if(method == "davies"){
    # P0_half = P0
    W1 = W - (W %*% X1) %*%solve( t(X1)%*%X1)%*% t(X1)
  } 
  
  
  if( method == "liu" ){
    
    out<-Get_Liu_PVal(Q, W,Q.res)    
    pval.zero.msg=NULL
    
  } else if( method == "liu.mod" ){
    
    out<-Get_Liu_PVal.MOD(Q, W,Q.res)   
    pval.zero.msg = NULL 
    
  } else if( method == "davies" ){
    
    out<-Get_Davies_PVal(Q, W1,Q.res)  
    pval.zero.msg = out$pval.zero.msg  
    
  } else {
    stop("Invalid Method!")
  }
  
  
  
  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling
           , Test.Type = method, Q = Q, param=out$param, pval.zero.msg=pval.zero.msg )  
  
  return(re)
  
}


#####################################
#similarity kernel function
Kernelsim=function(x,sim){
  K=matrix(,nrow(x),nrow(x))
  C<-sum(diag(sim))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x)) {
      K[i,j]=exp(-0.5*t((x[i,]-x[j,]))%*%sim%*%(x[i,]-x[j,])/C)}}
  K
}