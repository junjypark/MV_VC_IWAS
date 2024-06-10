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


data_intersect <- function(chr, IDP, geneindex, IGAP_snp) {
  tryCatch({
    
    # read UKB IDP
    familyID2 <- readRDS("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/familyID2.rds")
    UKBB_filename <- paste("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/genedata1M/chr", chr, "/chr", chr, "_gene", geneindex, sep = "")
    UKBB_gene <- read_plink(UKBB_filename, verbose = F)
    UKBB_gene$fam <- UKBB_gene$fam[which(UKBB_gene$fam$fam %in% familyID2), ]
    UKBB_gene$X <- UKBB_gene$X[, which(colnames(UKBB_gene$X) %in% familyID2)]
    UKBB_snp <- unlist(UKBB_gene$bim[, 2])
    
    # read 1000G clump
    IDP <- sprintf("%04d", IDP)
    file1000G_name <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/clump_0.5_SNPlist/chr", chr, "/SNPlist_IDP", IDP, "_chr", chr, "_0.5.txt", sep = "")
    dft_1000G <- fread(file1000G_name, header = F)
    snp_1000G <- unlist(dft_1000G[, 1])
    
    
    # intersection
    int1_snp <- Reduce(intersect, list(UKBB_snp, snp_1000G))
    int2_snp <- Reduce(intersect, list(snp_1000G, IGAP_snp))
    full_int_snp <- Reduce(intersect, list(int1_snp, int2_snp))
    full_int_snp <- as.data.frame(full_int_snp)
    
    # save data
    if (nrow(full_int_snp) > 0) {
      outname <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5/chr", chr, "/Intersect3_0.5_chr", chr, "_IDP", IDP, "_gene", geneindex, ".txt", sep = "")
      write.table(full_int_snp, outname, col.names = F, row.names = F)
    }
    
  }, error = function(e){print(geneindex)})
}



# IGAP
IGAPfilename <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_chr", chr, ".txt", sep = "")
IGAP <- fread(IGAPfilename)
IGAP_snp <- unlist(IGAP[, 4])


all_MRI_ids <- readRDS("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB_IDP_chr/all_MRI_ids.rds")

#################################################
# simulation
#################################################

result=foreach(geneindex=gene_nums, .combine="rbind", .packages = c("data.table", "genio"))%dopar%{
  for (IDP in all_MRI_ids) {
    data_intersect(chr, IDP, geneindex, IGAP_snp)
  }
  print(geneindex)
}  

