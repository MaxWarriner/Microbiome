library(phyloseq)
library(BiocManager)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(rlang)
library(MicrobiotaProcess)
library(vegan)
library(dplyr)
library(ALDEx2)
library(microbiomeMarker)
library(ggsci)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(parallel)
library(doParallel)
library(pROC)

ps <- read_rds('categorized_data.RDS')
sam <- data.frame(ps@sam_data)

# ## Part I : Alpha Diversity (Chao1 and Shannon indices) -----------------


create_a_diversity_plot <- function(ps, variable) {
  # Convert sample data to data frame
  sam_data <- as(sample_data(ps), "data.frame")
  
  # First remove NAs
  keep_samples <- !is.na(sam_data[[variable]])
  ps_filtered <- prune_samples(keep_samples, ps)
  
  # Get updated sample data after NA removal
  sam_data_filtered <- as(sample_data(ps_filtered), "data.frame")
  
  # Count occurrences of each value in the variable
  value_counts <- table(sam_data_filtered[[variable]])
  
  # Identify values that appear at least 5 times
  valid_values <- names(value_counts)[value_counts >= 5]
  
  # Filter samples to only include those with valid values
  final_keep_samples <- sam_data_filtered[[variable]] %in% valid_values
  ps_final <- prune_samples(final_keep_samples, ps_filtered)
  
  # Get alpha diversity metrics
  alpha_div <- estimate_richness(ps_final, measures = c("Chao1", "Shannon"))
  alpha_div$sample_id <- rownames(alpha_div)
  alpha_div <- merge(alpha_div, sam_data_filtered, by.x = "sample_id", by.y = "row.names")
  
  # Create the plot
  a_diversity_factor <- plot_richness(ps_final, x = variable, color = variable, 
                                      measures = c("Chao1", "Shannon"))
  a_diversity_factor$layers[[2]] = NULL 
  
  # Add statistical annotations
  a_diversity_factor <- a_diversity_factor + 
    geom_boxplot() + 
    geom_jitter(alpha = 0.25) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 18, face = "bold"), 
          legend.position = "none", 
          axis.text.x = element_text(size = 14)) + 
    labs(x = gsub("_", " ", variable)) + 
    # Add p-values for each metric
    stat_compare_means(
      method = "kruskal.test",
      label = "p.format",
      label.x = 1,    # Left side
      label.y = Inf,  # Top
      hjust = 0,      # Left-align text
      vjust = 1.5,    # Pushes it down slightly from top edge
      size = 10
    )
  # scale_x_discrete(limits = c("never", "sometimes", "always"))
  
  print(a_diversity_factor)
  
  return(a_diversity_factor)
}

create_a_diversity_plot(ps, "Kitchen_Material")
create_a_diversity_plot(ps, "House_Floor_Material")


# ## Part II : Beta Diversity -----------------


bray <- phyloseq::distance(ps, method = "bray")
jaccard <- phyloseq::distance(ps, method = "jaccard")

jaccard_pcoa <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")
bray_pcoa <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")

create_pcoa_plot <- function(variable, jaccard_dist, bray_dist, jaccard_pcoa, bray_pcoa, sam) {
  
  permanova_jaccard <- vegan::adonis2(jaccard_dist ~ sam[[variable]])
  permanova_bray <- vegan::adonis2(bray_dist ~ sam[[variable]])
  
  pval_jaccard <- permanova_jaccard$`Pr(>F)`[1]
  pval_bray <- permanova_bray$`Pr(>F)`[1]
  
  p_text_jaccard <- ifelse(pval_jaccard < 0.001, "p < 0.001",
                           paste("p =", format(round(pval_jaccard, 3), nsmall = 3)))
  p_text_bray <- ifelse(pval_bray < 0.001, "p < 0.001",
                        paste("p =", format(round(pval_bray, 3), nsmall = 3)))
  
  pcoa_jaccard_plot <- ggordpoint(obj = jaccard_pcoa, biplot = FALSE, speciesannot = TRUE,
                                  factorNames = c(variable), ellipse = TRUE, linesize = 1.5,
                                  ellipse_linewd = 1, ellipse_lty = 2) +
    ggtitle(paste(gsub("_", " ", variable), "(Jaccard)")) +
    guides(color=guide_legend(title=gsub("_", " ", variable), override.aes = list(size = 4))) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
    annotate("text", x = Inf, y = Inf, label = p_text_jaccard,
             hjust = 1.1, vjust = 1.5, size = 12, fontface = "bold")
  
  pcoa_bray_plot <- ggordpoint(obj = bray_pcoa, biplot = FALSE, speciesannot = TRUE,
                               factorNames = c(variable), ellipse = TRUE, linesize = 1.5,
                               ellipse_linewd = 1, ellipse_lty = 2) +
    ggtitle(paste(gsub("_", " ", variable), "(Bray-Curtis)")) +
    guides(color=guide_legend(title=gsub("_", " ", variable), override.aes = list(size = 4))) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
    annotate("text", x = Inf, y = Inf, label = p_text_bray,
             hjust = 1.1, vjust = 1.5, size = 12, fontface = "bold")
  
  combined_plot <- pcoa_jaccard_plot + pcoa_bray_plot
  
  ggsave(combined_plot,
         filename = paste(variable, "_pcoa_combined.png", sep = ""),
         device = "png",
         height = 6, width = 14, units = "in")
  
  return(combined_plot)
}


# House Floor Material
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis")
ps <- readRDS("categorized_data.RDS")
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/Machine Learning/Heatmaps")
ps <- subset_samples(ps, House_Floor_Material %in% c("Cement", "Dust"))
sam <- ps@sam_data
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/beta diversity PCOA plots")

bray <- phyloseq::distance(ps, method = "bray")
jaccard <- phyloseq::distance(ps, method = "jaccard")

jaccard_pcoa <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")
bray_pcoa <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")

create_pcoa_plot("House_Floor_Material", jaccard, bray, jaccard_pcoa, bray_pcoa, sam)


# Kitchen Floor Material
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis")
ps <- readRDS("categorized_data.RDS")
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/Machine Learning/Heatmaps")
ps <- subset_samples(ps, Kitchen_Material %in% c("Cement", "Dust"))
sam <- ps@sam_data
setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/beta diversity PCOA plots")

bray <- phyloseq::distance(ps, method = "bray")
jaccard <- phyloseq::distance(ps, method = "jaccard")

jaccard_pcoa <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")
bray_pcoa <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")

create_pcoa_plot("Kitchen_Material", jaccard, bray, jaccard_pcoa, bray_pcoa, sam)




