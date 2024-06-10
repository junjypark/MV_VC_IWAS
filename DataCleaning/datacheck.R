library(data.table)
library(foreach)
library(doParallel)

gene_n_v <- c(2568, 1648, 1359, 1018, 1172, 1282, 1235, 919, 1028, 990, 1523, 1282, 550, 843, 954, 1083, 1489, 383, 1806, 723, 330, 583)

checkSNP_datas <- function(chr, geneindex) {
      
    # SNPlist
    SNPlistname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/intersection_snplist_0.5/chr", chr, "/SNPslist_0.5_chr", chr, "_gene", geneindex, ".txt")
    SNPlist <- unlist(fread(SNPlistname, header = F))
      
    # UKB GWAS
    chr_path <- paste("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5/chr", chr, "/", sep = "")
    all_files <- list.files(chr_path)
    matching_files <- grep(paste0("gene", geneindex, ".txt"), all_files, value = TRUE)
    UKBGWASchecks <- rep(NA, length(matching_files))
    for (i in 1:length(matching_files)) {
      UKBfilename <- paste0(chr_path, matching_files[i])
      UKBGWAS <- fread(UKBfilename, header = T)
      UKBGWASchecks[i] <- all(UKBGWAS$rsID %in% SNPlist)
    }
    UKBGWAScheckall <- NA
    UKBGWAScheckall <- all(UKBGWASchecks)
      
    # 1000G LD
    LDname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/1000G/R_LD_0.5/chr", chr, "/LD_chr", chr, "_gene", geneindex, ".rds")
    LD <- readRDS(LDname)
    LDcheck <- all(colnames(LD) %in% SNPlist)
      
    # IGAP data
    IGAPname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/IGAP/IGAP_0.5_bygene/chr", chr, "/IGAP_chr", chr, "_gene", geneindex, ".rds")
    IGAP <- readRDS(IGAPname)
    IGAPcheck <- all(IGAP$MarkerName %in% SNPlist)
    
    result <- list(geneindex, UKBGWAScheckall, LDcheck, IGAPcheck)
    return(result)
      
}


ncores = 40 #Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

for (chr in 5:22) {
  checkresult=foreach(geneindex=1:gene_n_v[2], .combine="rbind", .packages = c("data.table"), .errorhandling = 'remove')%dopar%{
    checkSNP_datas(chr, geneindex)
  }
  colnames(checkresult) <- c("geneindex", "UKBGWAS", "LD", "IGAP")
  
  outname <- paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/datacheck/chr", chr, ".rds")
  saveRDS(checkresult, outname)
  
  table(unlist(checkresult[,2]))
  table(unlist(checkresult[,3]))
  table(unlist(checkresult[,4]))
}








