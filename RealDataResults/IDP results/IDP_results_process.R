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

  for (i in 1:length(chr_result)) {
    
    if (length(chr_result[[i]]) == 0) {
      chr_S[i, ] <- rep(NA, 10)
      chr_D[i, ] <- rep(NA, 16)
    } else {
      
      # structural
      chr_S[i, 1] <- length(which(chr_result[[i]]$S== "regional and tissue volume"))
      chr_S[i, 2] <- length(which(chr_result[[i]]$S== "cortical area"))
      chr_S[i, 3] <- length(which(chr_result[[i]]$S== "cortical thickness"))
      chr_S[i, 4] <- length(which(chr_result[[i]]$S== "cortical grey-white contrast"))
      chr_S[i, 5] <- length(which(chr_result[[i]]$S== "regional and tissue intensity"))
      
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
        }else {
          chr_S[i, 6:10] <- rep(F, 5)
        }
      } else {
        chr_S[i, 6:10] <- rep(F, 5)
      }
      
      # diffusion
      chr_D[i, 1] <- length(which(chr_result[[i]]$D== "white matter hyperintensity volume"))
      chr_D[i, 2] <- length(which(chr_result[[i]]$D== "regional T2*"))
      chr_D[i, 3] <- length(which(chr_result[[i]]$D== "WM tract FA"))
      chr_D[i, 4] <- length(which(chr_result[[i]]$D== "WM tract MO"))
      chr_D[i, 5] <- length(which(chr_result[[i]]$D== "WM tract diffusivity"))
      chr_D[i, 6] <- length(which(chr_result[[i]]$D== "WM tract ICVF"))
      chr_D[i, 7] <- length(which(chr_result[[i]]$D== "WM tract OD"))
      chr_D[i, 8] <- length(which(chr_result[[i]]$D== "WM tract ISOVF"))
      
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
  setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results")
  saveRDS(chr_S, file = paste0("chr", chr, "_S.rds"))
  saveRDS(chr_D, file = paste0("chr", chr, "_D.rds"))
}


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

# chr19
chr <- 19
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

