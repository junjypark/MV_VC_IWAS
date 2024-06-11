############################################################################################################################################
setwd("/Users/tianyuan/Desktop/UoT/Research/UKBB/May 21 results/UKB_results")

library(qqman)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(ggforce)

################################################
# UKB diffusion data
################################################
# read in files chr1_results_D.rds to chr22_results_D.rds
# combine them into one file
# save the combined file as all_results_D.rds
for (i in 1:22) {
  file_name <- paste0("chr", i, "_results_D.rds")
  data <- readRDS(file_name)
  if (i == 1) {
    UKB_d <- data
  } else {
    UKB_d <- rbind(UKB_d, data)
  }
}
# remove column 5-7
UKB_d <- UKB_d[, -c(5, 7:8)]
colnames(UKB_d) <- c("CHR", "POS1", "POS2", "Gene", "P")
# take mean of POS1 and POS2 and append to data
UKB_d$BP <- (UKB_d$POS1 + UKB_d$POS2) / 2
# remove POS1 and POS2
UKB_d <- UKB_d[, -c(2:3)]
# remove row with NAs in P
UKB_d <- UKB_d[!is.na(UKB_d$P), ]
UKB_d$CHR <- factor(UKB_d$CHR, levels = 1:22)
min(UKB_d$P[which(UKB_d$P != 0)])
which(UKB_d$P == 0)
# check values that UKB_d$P == 0 and replace with a random value between 1e-15 to 1e-20
UKB_d$P[which(UKB_d$P == 0)] <- runif(length(which(UKB_d$P == 0)), min = 0, max = 1)* (1e-15 - 1e-20) + 1e-20
head(UKB_d)

################################################
# UKB structural data
################################################
# read in files chr1_results_S.rds to chr22_results_S.rds
# combine them into one file
# save the combined file as all_results_S.rds
for (i in 1:22) {
  file_name <- paste0("chr", i, "_results_S.rds")
  data <- readRDS(file_name)
  if (i == 1) {
    UKB_s <- data
  } else {
    UKB_s <- rbind(UKB_s, data)
  }
}
# remove pvalues in column 6-8
UKB_s <- UKB_s[, -c(5, 7:8)]
colnames(UKB_s) <- c("CHR", "POS1", "POS2", "Gene", "P")
# take mean of POS1 and POS2 and append to data
UKB_s$BP <- (UKB_s$POS1 + UKB_s$POS2) / 2
# remove POS1 and POS2
UKB_s <- UKB_s[, -c(2:3)]
# remove row with NAs in P
UKB_s <- UKB_s[!is.na(UKB_s$P), ]
UKB_s$CHR <- factor(UKB_s$CHR, levels = 1:22)
min(UKB_s$P[which(UKB_s$P != 0)])
which(UKB_s$P == 0)
# check values that UKB_s$P == 0 and replace with a random value between 1e-15 to 1e-20
UKB_s$P[which(UKB_s$P == 0)] <- runif(length(which(UKB_s$P == 0)), min = 0, max = 1)* (1e-15 - 1e-20) + 1e-20
head(UKB_s)
nrow(UKB_s)


################################################
# UKB functional data
################################################
# read in files chr1_results_F.rds to chr22_results_F.rds
# combine them into one file
# save the combined file as all_results_F.rds
for (i in 1:22) {
  file_name <- paste0("chr", i, "_results_F.rds")
  data <- readRDS(file_name)
  if (i == 1) {
    UKB_f <- data
  } else {
    UKB_f <- rbind(UKB_f, data)
  }
}
# remove pvalues in column 6-8
UKB_f <- UKB_f[, -c(5, 7:8)]
colnames(UKB_f) <- c("CHR", "POS1", "POS2", "Gene", "P")
# take mean of POS1 and POS2 and append to data
UKB_f$BP <- (UKB_f$POS1 + UKB_f$POS2) / 2
# remove POS1 and POS2
UKB_f <- UKB_f[, -c(2:3)]
# remove row with NAs in P
UKB_f <- UKB_f[!is.na(UKB_f$P), ]
UKB_f$CHR <- factor(UKB_f$CHR, levels = 1:22)
min(UKB_f$P[which(UKB_f$P != 0)])
which(UKB_f$P == 0)
# check values that UKB_f$P == 0 and replace with a random value between 1e-15 to 1e-20
UKB_f$P[which(UKB_f$P == 0)] <- runif(length(which(UKB_f$P == 0)), min = 0, max = 1)* (1e-15 - 1e-20) + 1e-20
head(UKB_f)
nrow(UKB_f)


################################################
# UKB diffusion figure
################################################

# Prepare the dataset for plotting
don_d <- UKB_d %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(UKB_d, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each Gene
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(CHR %in% c(8, 19), "yes", "no")) %>%
  mutate( is_annotate=ifelse((-log10(P)) > (-log10(0.05/nrow(UKB_d))), "yes", "no"))

# Prepare X axis
axisdf_d <- don_d %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
p_d <- ggplot(don_d, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=CHR), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf_d$CHR, breaks= axisdf_d$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +     # remove space between plot area and x axis
  
  # modality threshold
  geom_hline(yintercept = -log10(0.05/nrow(UKB_d)), linetype = "dashed", color = "blue") +
  
  # global threshold
  geom_hline(yintercept = -log10(0.05/(nrow(UKB_d) + nrow(UKB_s) + nrow(UKB_f))), linetype = "dashed", color = "red") +
  
  # # Add highlighted points
  geom_point(data=subset(don_d, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don_d, is_annotate=="yes"), aes(label=Gene), size = 3, 
                    max.overlaps = 80, segment.alpha = 0.5) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels 45 degrees
    axis.text.y = element_text(size = 10),  # Increase y-axis text size
    plot.margin = unit(rep(0.5, 4), "cm")
  ) +
  xlab(NULL) # Remove x-axis label

ggsave("manhattan_plot_UKB_d_NoAdj.png", plot = p_d, width = 12, height = 5, dpi = 300)


################################################
# UKB structural figure
################################################

# Prepare the dataset for plotting
don_s <- UKB_s %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(UKB_s, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each Gene
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(CHR %in% c(8, 19), "yes", "no")) %>%
  mutate( is_annotate=ifelse((-log10(P)) > (-log10(0.05/nrow(UKB_s))), "yes", "no"))

# Prepare X axis
axisdf_s <- don_s %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
p_s <- ggplot(don_s, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=CHR), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf_s$CHR, breaks= axisdf_s$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +     # remove space between plot area and x axis
  
  # modality threshold
  geom_hline(yintercept = -log10(0.05/nrow(UKB_s)), linetype = "dashed", color = "blue") +
  
  # global threshold
  geom_hline(yintercept = -log10(0.05/(nrow(UKB_d) + nrow(UKB_s) + nrow(UKB_f))), linetype = "dashed", color = "red") +
  
  # # Add highlighted points
  geom_point(data=subset(don_s, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don_s, is_annotate=="yes"), aes(label=Gene), size = 3, 
                    max.overlaps = 80, segment.alpha = 0.5) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels 45 degrees
    axis.text.y = element_text(size = 10),  # Increase y-axis text size
    plot.margin = unit(rep(0.5, 4), "cm")
  ) +
  xlab(NULL) # Remove x-axis label

ggsave("manhattan_plot_UKB_s_NoAdj.png", plot = p_s, width = 12, height = 5, dpi = 300)


################################################
# UKB functional figure
################################################

# Prepare the dataset for plotting
don_f <- UKB_f %>%
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(UKB_f, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each Gene
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_annotate=ifelse((-log10(P) )> (-log10(0.05/nrow(UKB_f))), "yes", "no"))

# Prepare X axis
axisdf_f <- don_f %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
p_f <- ggplot(don_f, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=CHR), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf_f$CHR, breaks= axisdf_f$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +     # remove space between plot area and x axis
  
  # modality threshold
  geom_hline(yintercept = -log10(0.05/nrow(UKB_f)), linetype = "dashed", color = "blue") +
  
  # global threshold
  geom_hline(yintercept = -log10(0.05/(nrow(UKB_d) + nrow(UKB_s) + nrow(UKB_f))), linetype = "dashed", color = "red") +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don_f, is_annotate=="yes"), aes(label=Gene), size = 3, 
                    max.overlaps = 80, segment.alpha = 0.5) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels 45 degrees
    axis.text.y = element_text(size = 10),  # Increase y-axis text size
    plot.margin = unit(rep(0.5, 4), "cm")
  ) +
  xlab(NULL) # Remove x-axis label

ggsave("manhattan_plot_UKB_f_NoAdj.png", plot = p_f, width = 12, height = 5, dpi = 300)


################################################
# IGAP chr8 overlap
################################################

chr8_d <- don_d$Gene[which(don_d$is_annotate == "yes" & don_d$CHR == "8")]
chr8_s <- don_s$Gene[which(don_s$is_annotate == "yes" & don_s$CHR == "8")]
overlap_8 <- intersect(chr8_d, chr8_s)

# Define unique and overlap genes
unique_diffusion <- setdiff(chr8_d, chr8_s)
unique_structural <- setdiff(chr8_s, chr8_d)
overlap_genes <- intersect(chr8_d, chr8_s)

# Create basic Venn plot with transparent fills and adjusted circle positions
p <- ggplot() +
  geom_circle(aes(x0 = -1.5, y0 = 0, r = 2.5), fill = "lightblue", color = "black", size = 0.5, alpha = 0.5) +
  geom_circle(aes(x0 = 1.5, y0 = 0, r = 2.5), fill = "lightyellow", color = "black", size = 0.5, alpha = 0.5) +
  annotate("text", x = -1.5, y = 0, label = paste(unique_diffusion, collapse = "\n"), size = 3, hjust = 1) +
  annotate("text", x = 1.5, y = 0, label = paste(unique_structural, collapse = "\n"), size = 3, hjust = 0) +
  annotate("text", x = 0, y = 0, label = paste(overlap_genes, collapse = "\n"), size = 3) +
  annotate("text", x = -2, y = 1.5, label = "Diffusion", size = 5, fontface = "bold") +  # Add title for Diffusion
  annotate("text", x = 2, y = 1.5, label = "Structural", size = 5, fontface = "bold") +  # Add title for Structural
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA))

# Print the plot
print(p)

ggsave("venn_plot_UKB8_NoAdj.png", plot = p, width = 6, height = 4, dpi = 300)

################################################
# IGAP chr19 overlap
################################################

chr19_d <- don_d$Gene[which(don_d$is_annotate == "yes" & don_d$CHR == "19")]
chr19_s <- don_s$Gene[which(don_s$is_annotate == "yes" & don_s$CHR == "19")]
overlap_19 <- intersect(chr19_d, chr19_s)

# Define unique and overlap genes
unique_diffusion <- setdiff(chr19_d, chr19_s)
unique_structural <- setdiff(chr19_s, chr19_d)
overlap_genes <- intersect(chr19_d, chr19_s)

# Split the overlap genes into two columns
overlap_genes_col1 <- overlap_genes[1:ceiling(length(overlap_genes)/2)]
overlap_genes_col2 <- overlap_genes[(ceiling(length(overlap_genes)/2) + 1):length(overlap_genes)]

# Create basic Venn plot with transparent fills and adjusted circle positions
p <- ggplot() +
  geom_circle(aes(x0 = -1.5, y0 = 0, r = 2.5), fill = "lightblue", color = "black", size = 0.5, alpha = 0.5) +
  geom_circle(aes(x0 = 1.5, y0 = 0, r = 2.5), fill = "lightyellow", color = "black", size = 0.5, alpha = 0.5) +
  annotate("text", x = -2, y = 0, label = paste(unique_diffusion, collapse = "\n"), size = 3, hjust = 1) +
  annotate("text", x = 2, y = 0, label = paste(unique_structural, collapse = "\n"), size = 3, hjust = 0) +
  annotate("text", x = -0.25, y = seq(1.5, -1.5, length.out = length(overlap_genes_col1)), 
           label = overlap_genes_col1, size = 3, hjust = 1) +
  annotate("text", x = 0.95, y = seq(1.5, -1.5, length.out = length(overlap_genes_col2)), 
           label = overlap_genes_col2, size = 3, hjust = 1) +
  annotate("text", x = -2, y = 2, label = "Diffusion", size = 5, fontface = "bold") +  # Add title for Diffusion
  annotate("text", x = 2, y = 2, label = "Structural", size = 5, fontface = "bold") +  # Add title for Structural
  theme_void()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA))

# Print the plot
print(p)

ggsave("venn_plot_UKB19_NoAdj.png", plot = p, width = 6, height = 4, dpi = 300)