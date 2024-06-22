############################################################################################################
# IGAP
############################################################################################################
setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results")

# read in all diffusion data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_D.rds")
  if (chr == 1) {
    all_d <- readRDS(filename)[[1]]
  } else {
    all_d <- c(all_d, readRDS(filename)[[1]])
  }
}
all_d <- unlist(all_d)
# remove NA in data_d
all_d <- all_d[!is.na(all_d)]
length(all_d)

# diffusion by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_D.rds")
  data_d <- readRDS(filename)[[1]]
  
  data_t <- vector("list", length(data_d))
  for (gene in 1:length(data_d)) {
    tryCatch({
      data_t[[gene]] <- data_d[[gene]] < (0.05/270127)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_D_t.rds"))
}

# read in all structural data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_S.rds")
  if (chr == 1) {
    all_s <- readRDS(filename)[[1]]
  } else {
    all_s <- c(all_s, readRDS(filename)[[1]])
  }
}

all_s <- unlist(all_s)
# remove NA in data_s
all_s <- all_s[!is.na(all_s)]
length(all_s)

# structural by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_S.rds")
  data_s <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_s))
  for (gene in 1:length(data_s)) {
    tryCatch({
      data_t[[gene]] <- data_s[[gene]] < (0.05/400995)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_S_t.rds"))
}

# read in all functional data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_F.rds")
  if (chr == 1) {
    all_f <- readRDS(filename)[[1]]
  } else {
    all_f <- c(all_f, readRDS(filename)[[1]])
  }
}

all_f <- unlist(all_f)
# remove NA in data_s
all_f <- all_f[!is.na(all_f)]
length(all_f)

# functional by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_F.rds")
  data_f <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_f))
  for (gene in 1:length(data_f)) {
    tryCatch({
      data_t[[gene]] <- data_f[[gene]] < (0.05/27885)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_F_t.rds"))
}

IDP_description <- read.csv("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP_description.csv")
# diffusion IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_D_t.rds")
  data_D_t <- readRDS(filename)
  result_D <- data.frame()
  for (gene in 1:length(data_D_t)) {
    tryCatch({
      if (!is.na(unlist(data_D_t[[gene]]))) {
        if (length(which(unlist(data_D_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_D_t[[gene]]) == T)), 
                                        paste(names(data_D_t[[gene]])[which(unlist(data_D_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_D))
}

# structural IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_S_t.rds")
  data_S_t <- readRDS(filename)
  result_S <- data.frame()
  for (gene in 1:length(data_S_t)) {
    tryCatch({
      if (!is.na(unlist(data_S_t[[gene]]))) {
        if (length(which(unlist(data_S_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_S_t[[gene]]) == T)), 
                                        paste(names(data_S_t[[gene]])[which(unlist(data_S_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_S))
}

# functional IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_F_t.rds")
  data_F_t <- readRDS(filename)
  result_F <- data.frame()
  for (gene in 1:length(data_F_t)) {
    tryCatch({
      if (!is.na(unlist(data_F_t[[gene]]))) {
        if (length(which(unlist(data_F_t[[gene]]) == T)) > 0) {
          result_F <- rbind(result_F, c(chr_hg[gene, 4], length(which(unlist(data_F_t[[gene]]) == T)), 
                                        paste(names(data_F_t[[gene]])[which(unlist(data_F_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_F))
}


# read in all diffusion data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_D.rds")
  if (chr == 1) {
    all_d <- readRDS(filename)[[1]]
  } else {
    all_d <- c(all_d, readRDS(filename)[[1]])
  }
}
all_d <- unlist(all_d)
# remove NA in data_d
all_d <- all_d[!is.na(all_d)]
length(all_d)

# diffusion by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_D.rds")
  data_d <- readRDS(filename)[[1]]
  
  data_t <- vector("list", length(data_d))
  for (gene in 1:length(data_d)) {
    tryCatch({
      data_t[[gene]] <- data_d[[gene]] < (0.05/270127)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_D_t.rds"))
}

# read in all structural data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_S.rds")
  if (chr == 1) {
    all_s <- readRDS(filename)[[1]]
  } else {
    all_s <- c(all_s, readRDS(filename)[[1]])
  }
}

all_s <- unlist(all_s)
# remove NA in data_s
all_s <- all_s[!is.na(all_s)]
length(all_s)

# structural by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_S.rds")
  data_s <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_s))
  for (gene in 1:length(data_s)) {
    tryCatch({
      data_t[[gene]] <- data_s[[gene]] < (0.05/400995)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_S_t.rds"))
}

# read in all functional data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_F.rds")
  if (chr == 1) {
    all_f <- readRDS(filename)[[1]]
  } else {
    all_f <- c(all_f, readRDS(filename)[[1]])
  }
}

all_f <- unlist(all_f)
# remove NA in data_s
all_f <- all_f[!is.na(all_f)]
length(all_f)

# functional by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_results_F.rds")
  data_f <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_f))
  for (gene in 1:length(data_f)) {
    tryCatch({
      data_t[[gene]] <- data_f[[gene]] < (0.05/27885)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_results_F_t.rds"))
}

IDP_description <- read.csv("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP_description.csv")
# diffusion IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_D_t.rds")
  data_D_t <- readRDS(filename)
  result_D <- data.frame()
  for (gene in 1:length(data_D_t)) {
    tryCatch({
      if (!is.na(unlist(data_D_t[[gene]]))) {
        if (length(which(unlist(data_D_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_D_t[[gene]]) == T)), 
                                        paste(names(data_D_t[[gene]])[which(unlist(data_D_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_D))
}

# structural IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_S_t.rds")
  data_S_t <- readRDS(filename)
  result_S <- data.frame()
  for (gene in 1:length(data_S_t)) {
    tryCatch({
      if (!is.na(unlist(data_S_t[[gene]]))) {
        if (length(which(unlist(data_S_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_S_t[[gene]]) == T)), 
                                        paste(names(data_S_t[[gene]])[which(unlist(data_S_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_S))
}

# functional IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_results_F_t.rds")
  data_F_t <- readRDS(filename)
  result_F <- data.frame()
  for (gene in 1:length(data_F_t)) {
    tryCatch({
      if (!is.na(unlist(data_F_t[[gene]]))) {
        if (length(which(unlist(data_F_t[[gene]]) == T)) > 0) {
          result_F <- rbind(result_F, c(chr_hg[gene, 4], length(which(unlist(data_F_t[[gene]]) == T)), 
                                        paste(names(data_F_t[[gene]])[which(unlist(data_F_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_F))
}


############################################################################################################
# UKB
############################################################################################################

setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results")

# read in all diffusion data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_D.rds")
  if (chr == 1) {
    all_d <- readRDS(filename)[[1]]
  } else {
    all_d <- c(all_d, readRDS(filename)[[1]])
  }
}
all_d <- unlist(all_d)
# remove NA in data_d
all_d <- all_d[!is.na(all_d)]
length(all_d)

# diffusion by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_D.rds")
  data_d <- readRDS(filename)[[1]]
  
  data_t <- vector("list", length(data_d))
  for (gene in 1:length(data_d)) {
    tryCatch({
      data_t[[gene]] <- data_d[[gene]] < (0.05/235065)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_D_t.rds"))
}

# read in all structural data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_S.rds")
  if (chr == 1) {
    all_s <- readRDS(filename)[[1]]
  } else {
    all_s <- c(all_s, readRDS(filename)[[1]])
  }
}

all_s <- unlist(all_s)
# remove NA in data_s
all_s <- all_s[!is.na(all_s)]
length(all_s)

# structural by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_S.rds")
  data_s <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_s))
  for (gene in 1:length(data_s)) {
    tryCatch({
      data_t[[gene]] <- data_s[[gene]] < (0.05/355606)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_S_t.rds"))
}

# read in all functional data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_F.rds")
  if (chr == 1) {
    all_f <- readRDS(filename)[[1]]
  } else {
    all_f <- c(all_f, readRDS(filename)[[1]])
  }
}

all_f <- unlist(all_f)
# remove NA in data_s
all_f <- all_f[!is.na(all_f)]
length(all_f)

# functional by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_F.rds")
  data_f <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_f))
  for (gene in 1:length(data_f)) {
    tryCatch({
      data_t[[gene]] <- data_f[[gene]] < (0.05/24464)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_F_t.rds"))
}

IDP_description <- read.csv("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP_description.csv")
# diffusion IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_D_t.rds")
  data_D_t <- readRDS(filename)
  result_D <- data.frame()
  for (gene in 1:length(data_D_t)) {
    tryCatch({
      if (!is.na(unlist(data_D_t[[gene]]))) {
        if (length(which(unlist(data_D_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_D_t[[gene]]) == T)), 
                                        paste(names(data_D_t[[gene]])[which(unlist(data_D_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_D))
}

# structural IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_S_t.rds")
  data_S_t <- readRDS(filename)
  result_S <- data.frame()
  for (gene in 1:length(data_S_t)) {
    tryCatch({
      if (!is.na(unlist(data_S_t[[gene]]))) {
        if (length(which(unlist(data_S_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_S_t[[gene]]) == T)), 
                                        paste(names(data_S_t[[gene]])[which(unlist(data_S_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_S))
}

# functional IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_F_t.rds")
  data_F_t <- readRDS(filename)
  result_F <- data.frame()
  for (gene in 1:length(data_F_t)) {
    tryCatch({
      if (!is.na(unlist(data_F_t[[gene]]))) {
        if (length(which(unlist(data_F_t[[gene]]) == T)) > 0) {
          result_F <- rbind(result_F, c(chr_hg[gene, 4], length(which(unlist(data_F_t[[gene]]) == T)), 
                                        paste(names(data_F_t[[gene]])[which(unlist(data_F_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_F))
}


# read in all diffusion data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_D.rds")
  if (chr == 1) {
    all_d <- readRDS(filename)[[1]]
  } else {
    all_d <- c(all_d, readRDS(filename)[[1]])
  }
}
all_d <- unlist(all_d)
# remove NA in data_d
all_d <- all_d[!is.na(all_d)]
length(all_d)

# diffusion by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_D.rds")
  data_d <- readRDS(filename)[[1]]
  
  data_t <- vector("list", length(data_d))
  for (gene in 1:length(data_d)) {
    tryCatch({
      data_t[[gene]] <- data_d[[gene]] < (0.05/235065)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_D_t.rds"))
}

# read in all structural data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_S.rds")
  if (chr == 1) {
    all_s <- readRDS(filename)[[1]]
  } else {
    all_s <- c(all_s, readRDS(filename)[[1]])
  }
}

all_s <- unlist(all_s)
# remove NA in data_s
all_s <- all_s[!is.na(all_s)]
length(all_s)

# structural by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_S.rds")
  data_s <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_s))
  for (gene in 1:length(data_s)) {
    tryCatch({
      data_t[[gene]] <- data_s[[gene]] < (0.05/400995)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_S_t.rds"))
}

# read in all functional data
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_F.rds")
  if (chr == 1) {
    all_f <- readRDS(filename)[[1]]
  } else {
    all_f <- c(all_f, readRDS(filename)[[1]])
  }
}

all_f <- unlist(all_f)
# remove NA in data_s
all_f <- all_f[!is.na(all_f)]
length(all_f)

# functional by chr
for (chr in 1:22) {
  filename <- paste0("chr", chr, "_UKBresults_F.rds")
  data_f <- readRDS(filename)[[1]]
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/hg19_chr", chr, ".rds"))
  data_t <- vector("list", length(data_f))
  for (gene in 1:length(data_f)) {
    tryCatch({
      data_t[[gene]] <- data_f[[gene]] < (0.05/27885)
    }, error = function(e) {print(gene)})
  }
  saveRDS(data_t, file = paste0("chr", chr, "_UKBresults_F_t.rds"))
}

IDP_description <- read.csv("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP_description.csv")
# diffusion IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_D_t.rds")
  data_D_t <- readRDS(filename)
  result_D <- data.frame()
  for (gene in 1:length(data_D_t)) {
    tryCatch({
      if (!is.na(unlist(data_D_t[[gene]]))) {
        if (length(which(unlist(data_D_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_D_t[[gene]]) == T)), 
                                        paste(names(data_D_t[[gene]])[which(unlist(data_D_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_D))
}

# structural IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_S_t.rds")
  data_S_t <- readRDS(filename)
  result_S <- data.frame()
  for (gene in 1:length(data_S_t)) {
    tryCatch({
      if (!is.na(unlist(data_S_t[[gene]]))) {
        if (length(which(unlist(data_S_t[[gene]]) == T)) > 0) {
          result_D <- rbind(result_D, c(chr_hg[gene, 4], length(which(unlist(data_S_t[[gene]]) == T)), 
                                        paste(names(data_S_t[[gene]])[which(unlist(data_S_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_S))
}

# functional IDP
for (chr in 1:22) {
  chr_hg <- readRDS(paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results//hg19_chr", chr, ".rds"))
  filename <- paste0("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/MVIWAS/MVIWAS results/chr", chr, "_UKBresults_F_t.rds")
  data_F_t <- readRDS(filename)
  result_F <- data.frame()
  for (gene in 1:length(data_F_t)) {
    tryCatch({
      if (!is.na(unlist(data_F_t[[gene]]))) {
        if (length(which(unlist(data_F_t[[gene]]) == T)) > 0) {
          result_F <- rbind(result_F, c(chr_hg[gene, 4], length(which(unlist(data_F_t[[gene]]) == T)), 
                                        paste(names(data_F_t[[gene]])[which(unlist(data_F_t[[gene]]) == T)], collapse = "; ")))
        }
      }
    }, error = function(e) {})
  }
  print(dim(result_F))
}


