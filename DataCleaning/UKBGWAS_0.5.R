#######################################
# cluster setup
#######################################
library(foreach)
library(doParallel)
library(data.table)
library(genio)

args=(commandArgs(TRUE))
set = as.numeric(args[1])

gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)
cumgeneset <- cumsum(ceiling(gene_n_v/160))
geneset <- ceiling(gene_n_v/160)

chr <- sum(!(set <= cumgeneset)) + 1
if (chr == 1) {
  gene_nums <- 1:160+(set-1)*160
  if (max(gene_nums) > gene_n_v[chr]) {
    gene_nums <- min(gene_nums):gene_n_v[chr]
  }
} else {
  gene_nums <- 1:160+(set-cumgeneset[chr-1]-1)*160
  if (max(gene_nums) > gene_n_v[chr]) {
    gene_nums <- min(gene_nums):gene_n_v[chr]
  }
}

print(chr)
print(gene_nums)

ncores = 80 #Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

UKB_summarydata <- function(chr, IDP, geneindex, MAFdata) {
  tryCatch({
    
    # read SNPlist
    IDP <- sprintf("%04d", IDP)
    readname <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5/chr", 
                      chr, "/Intersect3_0.5_chr", chr, "_IDP", IDP, "_gene", geneindex, ".txt", sep = "")
    SNPlist <- unlist(fread(readname, header = F))
    
    # summary data
    IDPname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/", "IDP", IDP, "_chr", chr, ".txt")
    IDP_summarydata <- fread(IDPname, header = T)
    
    # outdata
    betas <- rep(NA, length(SNPlist))
    se <- rep(NA, length(SNPlist))
    p <- rep(NA, length(SNPlist))
    MAF <- rep(NA, length(SNPlist))
    for (i in 1:length(SNPlist)) {
      sumindex <- which(IDP_summarydata$rsid == SNPlist[i])
      MAFindex <- which(MAFdata$ID == SNPlist[i])
      
      betas[i] <- IDP_summarydata$beta[sumindex]
      se[i] <- IDP_summarydata$se[sumindex]
      p[i] <- IDP_summarydata$`pval(-log10)`[sumindex]
      MAF[i] <- MAFdata$ALT_FREQS[MAFindex]
    }
    
    # final data
    finaldata <- data.frame(SNPlist, betas, se, p, MAF)
    colnames(finaldata) <- c("rsID", "beta", "se", "pval(-log10)", "MAF")
    finaloutname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5/chr", chr, "/UKBGWAS_0.5_chr", chr, 
                           "_IDP", IDP, "_gene", geneindex, ".txt")
    write.table(finaldata, finaloutname, row.names = F, col.names = T)
    
  }, error = function(e){print(geneindex)})
}



# MAF
MAFname <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKB_MAF_bychr/", "UKBB_ImgSampleGenotype_MAF_Chr", chr, ".afreq", sep = "")
MAFdata <- fread(MAFname, header = T)

all_MRI_ids <- readRDS("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/all_MRI_ids.rds")

#################################################
# results
#################################################

result=foreach(geneindex=gene_nums, .combine="rbind", .packages = c("data.table"))%dopar%{
  for (IDP in all_MRI_ids) {
    UKB_summarydata(chr, IDP, geneindex, MAFdata)
  }
  print(geneindex)
}