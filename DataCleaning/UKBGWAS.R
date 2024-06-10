
library(foreach)
library(doParallel)
library(data.table)
library(genio)
library(dplyr)
library(plink2R)

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
getDoParWorkers()

GWASfunc <- function(chr, index) {
  tryCatch({
    a_name <- paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/genedata1M/chr", chr, "/chr", 
                     chr, "_gene", index)
    a <- read_plink(a_name)
    b <- readRDS("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/behavioral.rds")
    sort_b <- b %>%
      filter(eid %in% a$fam[, 2]) %>%
      arrange(match(eid, a$fam[, 2]))
    
    GWAS_result <- as.data.frame(matrix(NA, ncol(a$bed), 5))
    colnames(GWAS_result) <- c("rsID", "beta", "SE", "z", "p")
    
    for (snp in 1:ncol(a$bed)) {
      dat <- sort_b[, 1:14]
      dat <- cbind(dat, a$bed[, snp])
      names(dat)[15] <- colnames(a$bed)[snp]
      dat <- dat[, -2]
      dat$sex <- as.factor(dat$sex)
      fit <- glm(AD ~ ., data = dat, family = binomial(link = "logit"))
      GWAS_result[snp, 1] <- colnames(a$bed)[snp]
      tryCatch({
        GWAS_result[snp, 2:5] <- summary(fit)$coef[14, ]
      }, error = function(e){print("error")})
    }
    
    outnames <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKB_ADGWAS/UKB_rawADGWAS/chr", chr, "/UKB_GWAS_chr", chr, "_gene", index, ".rds")
    saveRDS(GWAS_result, outnames)
    
  }, error = function(e){print("error")})
}

#################################################
# results
#################################################

result=foreach(geneindex=gene_nums, .combine="rbind", .packages = c("data.table"))%dopar%{
  GWASfunc(chr, geneindex)
  print(geneindex)
}

