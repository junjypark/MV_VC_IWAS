setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB")


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
        if (abs(chr_result[[i]]$sig_S) < 2.153873e-06) {
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
        if (abs(chr_result[[i]]$sig_D) < 2.283522e-06) {
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
  setwd("/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB")
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

############################################################################################################
# chr8
############################################################################################################

############################################################################################################
# structural
chr <- 8
chr_S <- readRDS(paste0("chr", chr, "_S.rds"))
chr_S <- as.data.frame(chr_S)
chr_S$gene[which(chr_S[, 7] == T)]
chr_S$gene[which(chr_S[, 8] == T)]
chr_S$gene[which(chr_S[, 9] == T)]
chr_S$gene[which(chr_S[, 10] == T)]
chr_S$gene[which(chr_S[, 11] == T)]

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
# write.csv(des_S_gene, file = "chr8_sMRI_IDP_description.csv", row.names = F)

# chr8 region
chr_hg <- readRDS(paste0("hg19_chr", chr, ".rds"))
chr_result <- readRDS(paste0("chr", chr, "_IDP_results.rds"))
# obtain $sig_S from each element of the list
chr_sig_S <- lapply(chr_result, function(x) x$sig_S)
# replace NULLs with NAs
chr_sig_S <- lapply(chr_sig_S, function(x) if (is.null(x)) NA else x)
chr_sig_S <- unlist(chr_sig_S)
chr_hg <- cbind(chr_hg, chr_sig_S)
# remove NAs in chr_hg[, 5]
which(chr_hg[, 5] < 0.05/19662)
chr_result_sig <- chr_result[which(chr_hg[, 5] < 0.05/19662)]
MRI_S_region <- c()
MRI_S_category <- c()
for (i in 1:length(chr_result_sig)) {
  if (is.null(MRI_S_region)) {
    MRI_S_region <- MRIIDP$region[as.numeric(chr_result_sig[[i]]$S[3, ])]
    MRI_S_category  <- MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$S[3, ])]
  } else {
    MRI_S_region <- c(MRI_S_region, MRIIDP$region[as.numeric(chr_result_sig[[i]]$S[3, ])])
    MRI_S_category <- c(MRI_S_category, MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$S[3, ])])
  }
}
# remove leading empty space in MRI_S_region
MRI_S_region <- gsub("^\\s+", "", MRI_S_region)
# make MRI_S_region all lower case
MRI_S_region <- tolower(MRI_S_region)
# count the frequency of each region
MRI_S_region_results <- as.data.frame(table(MRI_S_region, MRI_S_category))
# remove rows with freq 0 in MRI_S_region_results
MRI_S_region_results <- MRI_S_region_results[which(MRI_S_region_results$Freq != 0), ]
write.csv(MRI_S_region_results, file = "chr8_sMRI_region_results.csv", row.names = F)

############################################################################################################
# diffusion
chr_D <- readRDS(paste0("chr", chr, "_D.rds"))
chr_D <- as.data.frame(chr_D)
chr_D$gene[which(chr_D[, 10] == T)]
chr_D$gene[which(chr_D[, 11] == T)]
chr_D$gene[which(chr_D[, 12] == T)]
chr_D$gene[which(chr_D[, 13] == T)]
chr_D$gene[which(chr_D[, 14] == T)]
chr_D$gene[which(chr_D[, 15] == T)]
chr_D$gene[which(chr_D[, 16] == T)]
chr_D$gene[which(chr_D[, 17] == T)]

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

# chr8 region
chr_hg <- readRDS(paste0("hg19_chr", chr, ".rds"))
chr_result <- readRDS(paste0("chr", chr, "_IDP_results.rds"))
# obtain $sig_S from each element of the list
chr_sig_D <- lapply(chr_result, function(x) x$sig_D)
# replace NULLs with NAs
chr_sig_D <- lapply(chr_sig_D, function(x) if (is.null(x)) NA else x)
chr_sig_D <- unlist(chr_sig_D)
chr_hg <- cbind(chr_hg, chr_sig_D)
# remove NAs in chr_hg[, 5]
which(abs(chr_hg[, 5]) < 0.05/20937)
chr_result_sig <- chr_result[which(abs(chr_hg[, 5]) < 0.05/20937)]
MRI_D_region <- c()
MRI_D_category <- c()
for (i in 1:length(chr_result_sig)) {
  if (is.null(MRI_D_region)) {
    MRI_D_region <- MRIIDP$region[as.numeric(chr_result_sig[[i]]$D[3, ])]
    MRI_D_category <- MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$D[3, ])]
  } else {
    MRI_D_region <- c(MRI_D_region, MRIIDP$region[as.numeric(chr_result_sig[[i]]$D[3, ])])
    MRI_D_category <- c(MRI_D_category, MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$D[3, ])])
  }
}
# remove leading empty space in MRI_D_region
MRI_D_region <- gsub("^\\s+", "", MRI_D_region)
# make MRI_D_region all lower case
MRI_D_region <- tolower(MRI_D_region)
# count the frequency of each region by category
MRI_D_region_results <- as.data.frame(table(MRI_D_region, MRI_D_category))
# remove rows with freq 0 in MRI_D_region_results
MRI_D_region_results <- MRI_D_region_results[which(MRI_D_region_results$Freq != 0), ]
write.csv(MRI_D_region_results, file = "chr8_dMRI_region_results.csv", row.names = F)

############################################################################################################
# chr19
############################################################################################################

############################################################################################################
# structural
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

# chr19 region
chr_hg <- readRDS(paste0("hg19_chr", chr, ".rds"))
chr_result <- readRDS(paste0("chr", chr, "_IDP_results.rds"))
# obtain $sig_S from each element of the list
chr_sig_S <- lapply(chr_result, function(x) x$sig_S)
# replace NULLs with NAs
chr_sig_S <- lapply(chr_sig_S, function(x) if (is.null(x)) NA else x)
chr_sig_S <- unlist(chr_sig_S)
chr_hg <- cbind(chr_hg, chr_sig_S)
# remove NAs in chr_hg[, 5]
which(chr_hg[, 5] < 0.05/19662)
chr_result_sig <- chr_result[which(chr_hg[, 5] < 0.05/19662)]
MRI_S_region <- c()
MRI_S_category <- c()
for (i in 1:length(chr_result_sig)) {
  if (is.null(MRI_S_region)) {
    MRI_S_region <- MRIIDP$region[as.numeric(chr_result_sig[[i]]$S[3, ])]
    MRI_S_category  <- MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$S[3, ])]
  } else {
    MRI_S_region <- c(MRI_S_region, MRIIDP$region[as.numeric(chr_result_sig[[i]]$S[3, ])])
    MRI_S_category <- c(MRI_S_category, MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$S[3, ])])
  }
}
# remove leading empty space in MRI_S_region
MRI_S_region <- gsub("^\\s+", "", MRI_S_region)
# make MRI_S_region all lower case
MRI_S_region <- tolower(MRI_S_region)
# count the frequency of each region
MRI_S_region_results <- as.data.frame(table(MRI_S_region, MRI_S_category))
# remove rows with freq 0 in MRI_S_region_results
MRI_S_region_results <- MRI_S_region_results[which(MRI_S_region_results$Freq != 0), ]
write.csv(MRI_S_region_results, file = "chr19_sMRI_region_results.csv", row.names = F)


############################################################################################################
# diffusion
chr_D <- readRDS(paste0("chr", chr, "_D.rds"))
chr_D <- as.data.frame(chr_D)
chr_D$gene[which(chr_D[, 10] == T)]
chr_D$gene[which(chr_D[, 11] == T)]
chr_D$gene[which(chr_D[, 12] == T)]
chr_D$gene[which(chr_D[, 13] == T)]
chr_D$gene[which(chr_D[, 14] == T)]
chr_D$gene[which(chr_D[, 15] == T)]
chr_D$gene[which(chr_D[, 16] == T)]
chr_D$gene[which(chr_D[, 17] == T)]

des_D <- readRDS(paste0("chr", chr, "_D_description.rds"))
length(unique(unlist(des_D[[1]])))
table(MRIIDP$Category.name[which(MRIIDP$IDP.description %in% unique(unlist(des_D[[1]])))])

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

# chr19 region
chr_hg <- readRDS(paste0("hg19_chr", chr, ".rds"))
chr_result <- readRDS(paste0("chr", chr, "_IDP_results.rds"))
# obtain $sig_S from each element of the list
chr_sig_D <- lapply(chr_result, function(x) x$sig_D)
# replace NULLs with NAs
chr_sig_D <- lapply(chr_sig_D, function(x) if (is.null(x)) NA else x)
chr_sig_D <- unlist(chr_sig_D)
chr_hg <- cbind(chr_hg, chr_sig_D)
# remove NAs in chr_hg[, 5]
which(chr_hg[, 5] < 0.05/20937)
chr_result_sig <- chr_result[which(chr_hg[, 5] < 0.05/20937)]
MRI_D_region <- c()
MRI_D_category <- c()
for (i in 1:length(chr_result_sig)) {
  if (is.null(MRI_D_region)) {
    MRI_D_region <- MRIIDP$region[as.numeric(chr_result_sig[[i]]$D[3, ])]
    MRI_D_category <- MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$D[3, ])]
  } else {
    MRI_D_region <- c(MRI_D_region, MRIIDP$region[as.numeric(chr_result_sig[[i]]$D[3, ])])
    MRI_D_category <- c(MRI_D_category, MRIIDP$Category.name[as.numeric(chr_result_sig[[i]]$D[3, ])])
  }
}
# remove leading empty space in MRI_D_region
MRI_D_region <- gsub("^\\s+", "", MRI_D_region)
# make MRI_D_region all lower case
MRI_D_region <- tolower(MRI_D_region)
# count the frequency of each region by category
MRI_D_region_results <- as.data.frame(table(MRI_D_region, MRI_D_category))
# remove rows with freq 0 in MRI_D_region_results
MRI_D_region_results <- MRI_D_region_results[which(MRI_D_region_results$Freq != 0), ]
write.csv(MRI_D_region_results, file = "chr19_dMRI_region_results.csv", row.names = F)

# paste(MRI_S_region_results$MRI_S_region[which(MRI_S_region_results$MRI_S_category == "regional and tissue volume")], collapse = ", ")


############################################################################################################
# chr8 sMRI
############################################################################################################
# Load the data
file_path <- "/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB/chr8_sMRI_region_results.csv"
df <- read.csv(file_path)
colnames(df) <- c("IDP_region", "IDP_category", "Freq")

# Summarize the number of unique IDP regions by IDP category
summary_data <- df %>%
  group_by(IDP_category) %>%
  summarise(unique_regions = n_distinct(IDP_region))
print(summary_data)
sum(summary_data$unique_regions)

# check regions shared by different IDP categories
df %>%
  group_by(IDP_region) %>%
  summarise(unique_categories = n_distinct(IDP_category)) %>%
  filter(unique_categories > 1)

############################################################################################################
# chr8 dMRI
############################################################################################################
# Load the data
file_path <- "/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB/chr8_dMRI_region_results.csv"
df <- read.csv(file_path)
colnames(df) <- c("IDP_region", "IDP_category", "Freq")
length(unique(df$IDP_region))

# Summarize the number of unique IDP regions by IDP category
summary_data <- df %>%
  group_by(IDP_category) %>%
  summarise(unique_regions = n_distinct(IDP_region))
print(summary_data)
sum(summary_data$unique_regions)

df %>%
  group_by(IDP_region) %>%
  summarise(unique_categories = n_distinct(IDP_category)) %>%
  filter(unique_categories > 1)

############################################################################################################
# chr19 sMRI
############################################################################################################
# Load the data
file_path <- "/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB/chr19_sMRI_region_results.csv"
df <- read.csv(file_path)
colnames(df) <- c("IDP_region", "IDP_category", "Freq")
length(unique(df$IDP_region))

# Summarize the number of unique IDP regions by IDP category
summary_data <- df %>%
  group_by(IDP_category) %>%
  summarise(unique_regions = n_distinct(IDP_region))
print(summary_data)
sum(summary_data$unique_regions)

# check regions shared by different IDP categories
sharedf <- df %>%
  group_by(IDP_region) %>%
  summarise(unique_categories = n_distinct(IDP_category)) %>%
  filter(unique_categories > 1)

############################################################################################################
# chr19 sMRI
############################################################################################################
# Load the data
file_path <- "/Users/tianyuan/Documents/GitHub/MV_VC_IWAS/RealDataResults/IDP results/UKB/chr19_dMRI_region_results.csv"
df <- read.csv(file_path)
colnames(df) <- c("IDP_region", "IDP_category", "Freq")
length(unique(df$IDP_region))

# Summarize the number of unique IDP regions by IDP category
summary_data <- df %>%
  group_by(IDP_category) %>%
  summarise(unique_regions = n_distinct(IDP_region))
print(summary_data)
sum(summary_data$unique_regions)

# check regions shared by different IDP categories
sharedf <- df %>%
  group_by(IDP_region) %>%
  summarise(unique_categories = n_distinct(IDP_category)) %>%
  filter(unique_categories > 1)