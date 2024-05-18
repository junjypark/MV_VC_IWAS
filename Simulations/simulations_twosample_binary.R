library(bindata)
library(foreach)

#Number of simulations
n.sim=5000

#Number of samples
#assuming an equal sample size (n) in 2 datasets
n1=2000
n2=2000
n=n1+n2

#Number of SNPs
LD=readRDS("LD_chr16_gene1029.rds")
p=nrow(LD)

#Generate genotype data
G_s1=rmvbin(n1, margprob=rep(0.2,p), sigma = LD)
G_s2=rmvbin(n2, margprob=rep(0.2,p), sigma = LD)
gsum_s1=c(scale(apply(G_s1,1,mean)))
gsum_s2=c(scale(apply(G_s2,1,mean)))

#Sample correlations of SNPs
R=cor(G_s2)

#Number of imaging features in each modality
q=10

#tau is the signal: 
# --tau=0 is the null hypothesis.
# --tau value greater than 0 is the alternative hypothesis
tau=0

cl <- makeCluster(5)
registerDoParallel(cl)

results = foreach(sim=1:n.sim, .combine="c",.packages = c("MASS","CompQuadForm"))%dopar%{
  #Unmeasured confounding variable
  U_s1=rnorm(n1)
  U_s2=rnorm(n2)
  
  # Generate imaging features
  ## genetic effects to imaging features
  alpha1=matrix(rnorm(p*q,0,0.3),p)
  alpha2=matrix(rnorm(p*q,0,0.3),p)
  
  ## confounding effects
  gamma1=rnorm(q,0,1)
  gamma2=rnorm(q,0,1)
  M1_s1=G_s1%*%alpha1+U_s1%*%t(gamma1)+matrix(rnorm(n1*q,0,3),n1)
  M2_s1=G_s1%*%alpha2+U_s1%*%t(gamma2)+matrix(rnorm(n1*q,0,3),n1)
  M1_s2=G_s2%*%alpha1+U_s2%*%t(gamma1)+matrix(rnorm(n2*q,0,3),n2)
  M2_s2=G_s2%*%alpha2+U_s2%*%t(gamma2)+matrix(rnorm(n2*q,0,3),n2)

  #generate phenotype for sample 2
  pred_s2=c(-2+M1_s2%*%rnorm(q,0,tau)+M2_s2%*%rnorm(q, 0, 0.05)+gsum_s2)+U_s2
  prob_s2=exp(pred_s2)/(1+exp(pred_s2))
  y_s2=rbinom(n2, 1, prob_s2)
  
  #Get marginal GWAS summary statistics using sample 2
  #Assuming a logistic regression model
  gwas=rep(NA,p)
  for (j in 1:p){ gwas[j]=summary(glm(y_s2~G_s2[,j],family=binomial))$coef[-1,3] }
  
  #Get predictive models for imaging features using genotypes
  # Use Sample 1
  A1=lm(M1_s1~G_s1)$coef[-1,]
  A2=lm(M2_s1~G_s1)$coef[-1,]
  
  #p-value of the proposed method
  mv_vc_iwas(gwas, R, alpha1, alpha2)
}

#Type 1 error rate (tau=0) or power (tau>0)
mean(results<0.05,na.rm=T)

stopCluster(cl)
