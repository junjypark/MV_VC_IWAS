
# read in IDP id and category data
MRIIDP <- read.csv("/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/IDP_category/IDP_description.csv", header = T)

# for each Chr and each gene
for (chr in 1:22) {
  
  # Stage 2 final pvalue read in
  setwd("/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/UKB_results")
  chr_S_results <- readRDS(paste0("chr", chr, "_results_S.rds"))
  chr_D_results <- readRDS(paste0("chr", chr, "_results_D.rds"))
  chr_F_results <- readRDS(paste0("chr", chr, "_results_F.rds"))
  chr_S_results <- chr_S_results[, -c(6:8)]
  chr_D_results <- chr_D_results[, -c(6:8)]
  chr_F_results <- chr_F_results[, -c(6:8)]
  
  # chr gene length read in
  hg19_chr <- readRDS(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/hg19info/hg19_chr", chr, ".rds"))
  chr_IDP_results <- vector("list", length = nrow(hg19_chr))
  
  
  for (gene in 1:nrow(hg19_chr)) {
    tryCatch({
      
      # Stage 1 IDPs being included in Stage 2 analysis read in
      stage1_coef <- readRDS(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/chr", chr, 
                                    "/gene", gene, ".rds"))
      
      # Structural MRI IDP ID and category
      gene_S_colns <- as.numeric(colnames(stage1_coef$S))
      gene_S_colns <- rbind(MRIIDP$Category.name[which(MRIIDP$Pheno %in% gene_S_colns)], gene_S_colns)
      gene_S_colns <- rbind(MRIIDP$IDP.description[which(MRIIDP$Pheno %in% gene_S_colns)], gene_S_colns)
      
      # Diffusion MRI IDP ID and category
      gene_D_colns <- as.numeric(colnames(stage1_coef$D))
      gene_D_colns <- rbind(MRIIDP$Category.name[which(MRIIDP$Pheno %in% gene_D_colns)], gene_D_colns)
      gene_D_colns <- rbind(MRIIDP$IDP.description[which(MRIIDP$Pheno %in% gene_D_colns)], gene_D_colns)
      
      # Functional MRI IDP ID and category
      gene_F_colns <- as.numeric(colnames(stage1_coef$F))
      
      # Save MRI IDP ID, category, MRI significance results
      gene_IDP <- list(chr = chr, gene = gene, genename = chr_S_results[gene, 4],
                       S = gene_S_colns, D = gene_D_colns, F = gene_F_colns, 
                       sig_S = chr_S_results[gene, 5], sig_D = chr_D_results[gene, 5], sig_F = chr_F_results[gene, 5])
      chr_IDP_results[[gene]] <- gene_IDP
    }, error = function(e) {
      print(paste0("gene", gene, " failed"))
    })
  }  
  
  saveRDS(chr_IDP_results, file = paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/IDP_category/chr", chr, "_IDP_results.rds"))
}
