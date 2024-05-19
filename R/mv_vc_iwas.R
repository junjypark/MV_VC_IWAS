library(CompQuadForm)
library(MASS)

vc_sim=function(stat, lambda, n.sim=1e6){
  n.lambda = length(lambda)
  result = rep(0, n.sim)
  for (i in 1:n.lambda){
    result=result+lambda[i]*rchisq(n.sim, df=1)
  }
  pv=(1+sum(result>stat))/(n.sim+1)
  return(pv)
}

#corr   : p times 1 vector (from GWAS)
#LD     : p times p matrix of correlations between variants
#A1     : p times q1 matrix of coefficients predicting imaging data from genotypes
#A2     : p times q2 matrix of coefficients predicting imaging data from genotypes
#method : Either "davies" or "Liu". Davies method is used as a default.

mv_vc_iwas=function(corr,R,A1,A2 = NULL, egger = T, method = "liu", n.sim = 1e6){
  if (is.null(A2)){
    if (egger){ 
      A2=cbind(rep(1,nrow(A1)),A2) 
      U=crossprod(A1, corr)-crossprod(A1, R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, corr)
      V=crossprod(A1, R)%*%A1+crossprod(A1,R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, R)%*%A1
    } else{
      U=crossprod(A1, corr)
      V=crossprod(A1, R)%*%A1
    }    
  } else{
    if (egger){ A2=cbind(1,A2) }
    U=crossprod(A1, corr)-crossprod(A1, R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, corr)
    V=crossprod(A1, R)%*%A1+crossprod(A1,R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, R)%*%A1
  }
  

  
  # if (is.null(A2)){
  #   U=crossprod(A1, corr)
  #   V=crossprod(A1, R)%*%A1
  # } else{
  #   A2=cbind(1,A2)
  #   U=crossprod(A1, corr)-crossprod(A1, R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, corr)
  #   V=crossprod(A1, R)%*%A1-crossprod(A1,R)%*%A2%*%ginv(crossprod(A2,R)%*%A2)%*%crossprod(A2, R)%*%A1
  # }
  
  V=(V+t(V))/2 #Stability purposes
  
  if (nrow(V)==1){
    pval=1-pchisq(as.numeric(U^2/V), df = 1)
  } else if (nrow(V)<5){
    stat=sum(U^2)
    wts=eigen(V)$values
    pval=vc_sim(stat, wts, n.sim = n.sim)
  } else{
    stat=sum(U^2)
    wts=eigen(V)$values
    
    if (method=="davies"){ pval=davies(stat, lambda=wts)$Qq }
    if (method=="liu"){ pval=liu(stat, lambda=wts) }    
  }

  return(pval)
}
