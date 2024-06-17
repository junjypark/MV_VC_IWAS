setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results")


for (chr in 1:22) {
  chr_result <- readRDS(paste0("chr", chr, "_IDP_results.rds"))
  chr_hg <- readRDS(paste0("hg19_chr", chr, ".rds"))
  # structural
  chr_S <- as.data.frame(matrix(NA, nrow = length(chr_result), ncol = 10))
  colnames(chr_S) <- c("# regional and tissue volume", "# cortical area", "# cortical thickness", 
                       "# cortical grey-white contrast", "# regional and tissue intensity", 
                       "regional and tissue volume", "cortical area", "cortical thickness", 
                       "cortical grey-white contrast", "regional and tissue intensity")
  # diffusion
  chr_D <- as.data.frame(matrix(NA, nrow = length(chr_result), ncol = 16))
  colnames(chr_D) <- c("# white matter hyperintensity volume", "# regional T2*", 
                       "# WM tract FA", "# WM tract MO", "# WM tract diffusivity", 
                       "# WM tract ICVF", "# WM tract OD", "# WM tract ISOVF", 
                       "white matter hyperintensity volume", "regional T2*", 
                       "WM tract FA", "WM tract MO", "WM tract diffusivity", 
                       "WM tract ICVF", "WM tract OD", "WM tract ISOVF")
  IDP_S_description <- list()
  IDP_S_gene <- c()
  IDP_D_description <- list()
  IDP_D_gene <- c()
  
  for (i in 1:length(chr_result)) {
    
    if (length(chr_result[[i]]) == 0) {
      chr_S[i, ] <- rep(NA, 10)
      chr_D[i, ] <- rep(NA, 16)
    } else {
      
      S_result <- t(as.data.frame(chr_result[[i]]$S))
      colnames(S_result) <- c("description", "category", "col_num")
      S_result <- as.data.frame(S_result)
      # structural
      chr_S[i, 1] <- length(which(S_result$category == "regional and tissue volume"))
      chr_S[i, 2] <- length(which(S_result$category == "cortical area"))
      chr_S[i, 3] <- length(which(S_result$category == "cortical thickness"))
      chr_S[i, 4] <- length(which(S_result$category == "cortical grey-white contrast"))
      chr_S[i, 5] <- length(which(S_result$category == "regional and tissue intensity"))
      
      if (!is.na(chr_result[[i]]$sig_S )) {
        if (chr_result[[i]]$sig_S < 2.153873e-06) {
          # whether regional and tissue volume is significant
          if (chr_S[i, 1] > 0) {
            chr_S[i, 6] <- T
          } else {
            chr_S[i, 6] <- F
          }
          # whether cortical area is significant
          if (chr_S[i, 2] > 0) {
            chr_S[i, 7] <- T
          } else {
            chr_S[i, 7] <- F
          }
          # whether cortical thickness is significant
          if (chr_S[i, 3] > 0) {
            chr_S[i, 8] <- T
          } else {
            chr_S[i, 8] <- F
          }
          # whether cortical grey-white contrast is significant
          if (chr_S[i, 4] > 0) {
            chr_S[i, 9] <- T
          } else {
            chr_S[i, 9] <- F
          }
          # whether regional and tissue intensity is significant
          if (chr_S[i, 5] > 0) {
            chr_S[i, 10] <- T
          } else {
            chr_S[i, 10] <- F
          }
          
          IDP_S_description <- append(IDP_S_description, list(S_result$description))
          IDP_S_gene <- c(IDP_S_gene, chr_hg[i, 4])
          
        } else {
          chr_S[i, 6:10] <- rep(F, 5)
        }
      } else {
        chr_S[i, 6:10] <- rep(F, 5)
      }
      
      D_result <- t(as.data.frame(chr_result[[i]]$D))
      colnames(D_result) <- c("description", "category", "col_num")
      D_result <- as.data.frame(D_result)
      # diffusion
      chr_D[i, 1] <- length(which(D_result$category == "white matter hyperintensity volume"))
      chr_D[i, 2] <- length(which(D_result$category == "regional T2*"))
      chr_D[i, 3] <- length(which(D_result$category == "WM tract FA"))
      chr_D[i, 4] <- length(which(D_result$category == "WM tract MO"))
      chr_D[i, 5] <- length(which(D_result$category == "WM tract diffusivity"))
      chr_D[i, 6] <- length(which(D_result$category == "WM tract ICVF"))
      chr_D[i, 7] <- length(which(D_result$category == "WM tract OD"))
      chr_D[i, 8] <- length(which(D_result$category == "WM tract ISOVF"))
      
      # significant results
      if (!is.na(chr_result[[i]]$sig_D)) {
        if (chr_result[[i]]$sig_D < 2.283522e-06) {
          # whether white matter hyperintensity volume is significant
          if (chr_D[i, 1] > 0) {
            chr_D[i, 9] <- T
          } else {
            chr_D[i, 9] <- F
          }
          # whether regional T2* is significant
          if (chr_D[i, 2] > 0) {
            chr_D[i, 10] <- T
          } else {
            chr_D[i, 10] <- F
          }
          # whether WM tract FA is significant
          if (chr_D[i, 3] > 0) {
            chr_D[i, 11] <- T
          } else {
            chr_D[i, 11] <- F
          }
          # whether WM tract MO is significant
          if (chr_D[i, 4] > 0) {
            chr_D[i, 12] <- T
          } else {
            chr_D[i, 12] <- F
          }
          # whether WM tract diffusivity is significant
          if (chr_D[i, 5] > 0) {
            chr_D[i, 13] <- T
          } else {
            chr_D[i, 13] <- F
          }
          # whether WM tract ICVF is significant
          if (chr_D[i, 6] > 0) {
            chr_D[i, 14] <- T
          } else {
            chr_D[i, 14] <- F
          }
          # whether WM tract OD is significant
          if (chr_D[i, 7] > 0) {
            chr_D[i, 15] <- T
          } else {
            chr_D[i, 15] <- F
          }
          # whether WM tract ISOVF is significant
          if (chr_D[i, 8] > 0) {
            chr_D[i, 16] <- T
          } else {
            chr_D[i, 16] <- F
          }
          # append the description vector to the list of description
          IDP_D_description <- append(IDP_D_description, list(D_result$description))
          IDP_D_gene <- c(IDP_D_gene, chr_hg[i, 4])
          
        } else {
          chr_D[i, 9:16] <- rep(F, 8)
        }
      } else {
        chr_D[i, 9:16] <- rep(F, 8)
      }
    }
  }
  gene <- chr_hg[, 4]
  chr_S <- cbind(gene, chr_S)
  chr_D <- cbind(gene, chr_D)
  description_S <- list(IDP_S_description, IDP_S_gene)
  description_D <- list(IDP_D_description, IDP_D_gene)
  setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results")
  saveRDS(chr_S, file = paste0("chr", chr, "_S.rds"))
  saveRDS(chr_D, file = paste0("chr", chr, "_D.rds"))
  if (chr %in% c(8, 19)) {
    saveRDS(description_S, file = paste0("chr", chr, "_S_description.rds"))
    saveRDS(description_D, file = paste0("chr", chr, "_D_description.rds"))
  }
}


############################################################################################################
# IDP results
############################################################################################################


MRIIDP <- read.csv("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP_description.csv", header = T)
MRIIDP <- as.data.frame(MRIIDP)

# chr8
chr <- 8
chr_S <- readRDS(paste0("chr", chr, "_S.rds"))
chr_S <- as.data.frame(chr_S)

chr_S$gene[which(chr_S[, 7] == T)]
chr_S$gene[which(chr_S[, 8] == T)]
chr_S$gene[which(chr_S[, 9] == T)]
chr_S$gene[which(chr_S[, 10] == T)]
chr_S$gene[which(chr_S[, 11] == T)]
chr_D <- readRDS(paste0("chr", chr, "_D.rds"))
chr_D <- as.data.frame(chr_D)
chr_S$gene[which(chr_D[, 10] == T)]
chr_S$gene[which(chr_D[, 11] == T)]
chr_S$gene[which(chr_D[, 12] == T)]
chr_S$gene[which(chr_D[, 13] == T)]
chr_S$gene[which(chr_D[, 14] == T)]
chr_S$gene[which(chr_D[, 15] == T)]
chr_S$gene[which(chr_D[, 16] == T)]
chr_S$gene[which(chr_D[, 17] == T)]

# chr8 IDP description
des_S <- readRDS(paste0("chr", chr, "_S_description.rds"))
length(unique(unlist(des_S[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_S[[1]])))])
# sMRI IDP description by gene
des_S_gene <- rep(NA, length(des_S[[1]]))
for (i in 1:length(des_S_gene)) {
  des_S_gene[i] <- paste(unlist(des_S[[1]][i]), collapse = ";", sep = "")
}
des_S_gene <- cbind(des_S[[2]], des_S_gene)
colnames(des_S_gene) <- c("gene", "description")
write.csv(des_S_gene, file = "chr8_sMRI_IDP_description.csv", row.names = F)
# dMRI IDP description by gene
des_D <- readRDS(paste0("chr", chr, "_D_description.rds"))
length(unique(unlist(des_D[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_D[[1]])))])
des_D_gene <- rep(NA, length(des_D[[1]]))
for (i in 1:length(des_D_gene)) {
  des_D_gene[i] <- paste(unlist(des_D[[1]][i]), collapse = ";", sep = "")
}
des_D_gene <- cbind(des_D[[2]], des_D_gene)
colnames(des_D_gene) <- c("gene", "description")
write.csv(des_D_gene, file = "chr8_dMRI_IDP_description.csv", row.names = F)

# chr19
chr <- 19
chr_S <- readRDS(paste0("chr", chr, "_S.rds"))
chr_S <- as.data.frame(chr_S)
chr_S$gene[which(chr_S[, 7] == T)]
chr_S$gene[which(chr_S[, 8] == T)]
chr_S$gene[which(chr_S[, 9] == T)]
chr_S$gene[which(chr_S[, 10] == T)]
chr_S$gene[which(chr_S[, 11] == T)]

des_S <- readRDS(paste0("chr", chr, "_S_description.rds"))
length(unique(unlist(des_S[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_S[[1]])))])

chr_D <- readRDS(paste0("chr", chr, "_D.rds"))
chr_D <- as.data.frame(chr_D)
chr_S$gene[which(chr_D[, 10] == T)]
chr_S$gene[which(chr_D[, 11] == T)]
chr_S$gene[which(chr_D[, 12] == T)]
chr_S$gene[which(chr_D[, 13] == T)]
chr_S$gene[which(chr_D[, 14] == T)]
chr_S$gene[which(chr_D[, 15] == T)]
chr_S$gene[which(chr_D[, 16] == T)]
chr_S$gene[which(chr_D[, 17] == T)]

des_D <- readRDS(paste0("chr", chr, "_D_description.rds"))
length(unique(unlist(des_D[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_D[[1]])))])

# chr19 IDP description
des_S <- readRDS(paste0("chr", chr, "_S_description.rds"))
length(unique(unlist(des_S[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_S[[1]])))])
# sMRI IDP description by gene
des_S_gene <- rep(NA, length(des_S[[1]]))
for (i in 1:length(des_S_gene)) {
  des_S_gene[i] <- paste(unlist(des_S[[1]][i]), collapse = ";", sep = "")
}
des_S_gene <- cbind(des_S[[2]], des_S_gene)
colnames(des_S_gene) <- c("gene", "description")
write.csv(des_S_gene, file = "chr19_sMRI_IDP_description.csv", row.names = F)
# dMRI IDP description by gene
des_D <- readRDS(paste0("chr", chr, "_D_description.rds"))
length(unique(unlist(des_D[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_D[[1]])))])
des_D_gene <- rep(NA, length(des_D[[1]]))
for (i in 1:length(des_D_gene)) {
  des_D_gene[i] <- paste(unlist(des_D[[1]][i]), collapse = ";", sep = "")
}
des_D_gene <- cbind(des_D[[2]], des_D_gene)
colnames(des_D_gene) <- c("gene", "description")
write.csv(des_D_gene, file = "chr19_dMRI_IDP_description.csv", row.names = F)