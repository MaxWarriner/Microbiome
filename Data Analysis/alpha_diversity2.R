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

setwd("C:/Users/12697/Documents/Microbiome/Data Analysis")
ps <- read_rds('categorized_data.RDS')
sam <- data.frame(ps@sam_data)

create_a_diversity_plot <- function(ps, variable, alpha_level = 0.05, save_dir = NULL) {
  # Convert sample data to data frame
  sam_data <- as(sample_data(ps), "data.frame")
  
  # Remove NAs
  keep_samples <- !is.na(sam_data[[variable]])
  ps_filtered <- prune_samples(keep_samples, ps)
  sam_data_filtered <- as(sample_data(ps_filtered), "data.frame")
  
  # Keep groups with at least 5 samples
  value_counts <- table(sam_data_filtered[[variable]])
  valid_values <- names(value_counts)[value_counts >= 5]
  final_keep_samples <- sam_data_filtered[[variable]] %in% valid_values
  ps_final <- prune_samples(final_keep_samples, ps_filtered)
  sam_data_final <- as(sample_data(ps_final), "data.frame")
  
  # Compute diversity metrics
  alpha_div <- estimate_richness(ps_final, measures = c("Chao1", "Shannon"))
  sam_data_final$Chao1 <- alpha_div$Chao1
  sam_data_final$Shannon <- alpha_div$Shannon

  # Run significance tests
  chao_formula <- as.formula(paste("Chao1 ~ ", variable, sep = ""))
  shannon_formula <- as.formula(paste("Shannon ~ ", variable, sep = ""))
  
  chao_p <- kruskal.test(chao_formula, data = sam_data_final)$p.value
  shan_p <- kruskal.test(shannon_formula, data = sam_data_final)$p.value
  
  # Check if either metric is significant
  if (chao_p < alpha_level | shan_p < alpha_level) {
    message(sprintf("Significant results for '%s' (Chao1 p=%.4f, Shannon p=%.4f)", variable, chao_p, shan_p))
    
    # Create plot
    p <- plot_richness(ps_final, x = variable, color = variable, 
                       measures = c("Chao1", "Shannon"))
    p$layers[[2]] <- NULL
    
    p <- p + 
      geom_boxplot() + 
      geom_jitter(alpha = 0.25) + 
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18, face = "bold"),
            legend.position = "none",
            axis.text.x = element_text(size = 14)) +
      labs(x = gsub("_", " ", variable)) +
      stat_compare_means(
        method = "kruskal.test",
        label = "p.format",
        label.x = 1,
        label.y = Inf,
        hjust = 0,
        vjust = 1.5,
        size = 10
      )
    
      ggsave(filename = paste0(variable, "_diversity_plot.png"),
             plot = p, width = 8, height = 6, dpi = 300)

    
    return(p)
  } else {
    message(sprintf("No significant difference for '%s' (Chao1 p=%.4f, Shannon p=%.4f) — skipped.", variable, chao_p, shan_p))
    return(NULL)
  }
}

setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/alpha diversity boxplots/Significant Plots")
for(factor in colnames(sam)){
  try(create_a_diversity_plot(ps, factor))
}




