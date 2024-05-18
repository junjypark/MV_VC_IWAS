library(bindata)
library(foreach)
library(doParallel)


n.sim = 5000  #Number of simulations
n = 2000      #Number of samples

LD = readRDS("LD_chr16_gene1029.rds")
p = nrow(LD)  #Number of SNPs

#Generate genotype data
G = rmvbin(n, margprob=rep(0.2,p), sigma = LD)
gsum = scale(rowMeans(G))

#Sample correlations of SNPs
R = cor(G)

#Number of imaging features in each modality
q = 10

#tau is the signal: 
# --tau=0 is the null hypothesis.
# --tau value greater than 0 is the alternative hypothesis
tau=0.1

cl = makeCluster(5)
registerDoParallel(cl)

results = foreach(sim=1:n.sim, .combine="rbind", .packages = c("MASS","CompQuadForm"))%dopar%{
  #Unmeasured confounding variable
  U=rnorm(n) 
  
  # Generate imaging features
  ## genetic effects to imaging features
  alpha1 = matrix(rnorm(p*q,0,0.3),p)
  alpha2 = matrix(rnorm(p*q,0,0.3),p)
  
  ## confounding effects
  gamma1 = rnorm(q,0,1)
  gamma2 = rnorm(q,0,1)
  
  M1 = G%*%alpha1 + U%*%t(gamma1) + matrix(rnorm(n*q,0,5),n)
  M2 = G%*%alpha2 + U%*%t(gamma2) + matrix(rnorm(n*q,0,5),n)

  #generate phenotype.  
  pred = c(-2 + M1%*%rnorm(q,0,tau) + M2%*%rnorm(q, 0, 0.05) + gsum + U)
  y = pred + rnorm(n, 0, 4)

  #Get GWAS summary statistics
  #Assuming a logistic regression model
  gwas = vector()
  for (j in 1:p){ gwas[j]=summary(lm(y~G[,j]))$coef[-1,3] }

  #Get predictive models for imaging features using genotypes
  A1 = lm(M1 ~ G)$coef[-1,]
  A2 = lm(M2 ~ G)$coef[-1,]

  #p-value of the proposed method (Method 1)
  pvalue_m1 = mv_vc_iwas(gwas, R, A1, A2)

  #p-value of the proposed method (Method 2)
  pvec_m2 = vector()
  for (j in 1:q){ pvec_m2[j] = mv_vc_iwas(gwas, R, A1[,j], A2) }
  pvalue_m2 = min(pvec_m2)
  
  #p-value of the proposed method (Method 3)
  pvalue_m3 = mv_vc_iwas(gwas, R, A1)
  
  #p-value of the proposed method (Method 4)
  pvec_m4 = vector()
  for (j in 1:q){ pvec_m4[j] = mv_vc_iwas(gwas, R, A1[,j]) }
  pvalue_m4 = min(pvec_m4)
  
  c(pvalue_m1, pvalue_m2, pvalue_m3, pvalue_m4)
}

stopCluster(cl)
