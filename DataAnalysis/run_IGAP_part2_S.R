args=(commandArgs(TRUE))

function_call=as.character(args[1])      # ='power_one_changepoint'
if(as.character(args[2])=='ALL'){
  job_num=as.character(args[2])          # = a number between 1 and 6
}else{
  job_num=as.numeric(args[2])
}
ncores=as.numeric(args[3])               # 80 (niagara number of cores per node)
total_nodes_to_use=as.numeric(args[4])   # always equal to 1
compute_server=as.character(args[5])     # name of compute server

chr = job_num

library(foreach)
library(doParallel)

source("/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/files/method.R")

ncores = 1#Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

hg=readRDS(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/hg19info/hg19_chr",chr,".rds"))
pvecs = pvecs2 = pvecs3 = pvecs4 = rep(NA, nrow(hg))
for (gene in 1:nrow(hg)){
  pvecs[gene]=NA
  tryCatch({
    #Load coefficients
    setwd(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/chr",chr))
    coef=readRDS(paste0("gene",gene,".rds"))
    
    A1=as.matrix(cbind(coef$S))
    A2=as.matrix(cbind(coef$F, coef$D))
    
    #Load LD
    setwd(paste0("/gpfs/fs1/home/j/junpark/junpark/imgGenet/chr",chr))
    LD=readRDS(paste0("combine_gene",gene,".rds"))$LD
    index=match(rownames(A1), rownames(LD))
    LD=LD[index,index]
    
    #Load GWAS summary statistics
    setwd(paste0("/gpfs/fs0/project/j/junpark/tianyu47/IGAP_0.5/chr",chr))
    summ=readRDS(paste0("IGAP_chr",chr,"_gene", gene,".rds"))
    index=match(rownames(A1), summ$MarkerName)
    summ=summ[index,]
    
    tstats=summ$Beta/summ$SE
    pvecs[gene]=mv_vc_iwas(tstats,LD,A1,A2)
    pvecs2[gene]=mv_vc_iwas(tstats,LD,A1,A2=NULL)
    pvecs3[gene]=mv_vc_iwas(tstats,LD,A1,A2=NULL, egger = F)
    pvecs4[gene]=mv_vc_iwas(tstats,LD,A1,A2, egger=F)
    
  }, error = function(e) {print("error")})
}

hg$pval=pvecs
hg$pval2=pvecs2
hg$pval3=pvecs3
hg$pval4=pvecs4


setwd("/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/IGAP_results")
saveRDS(hg, paste0("chr",chr,"_results_S.rds"))


