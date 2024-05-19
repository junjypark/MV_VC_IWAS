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
