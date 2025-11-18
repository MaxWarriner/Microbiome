library(tidyverse)
library(phyloseq)
library(vegan)
library(stats)
library(MicrobiotaProcess)
library(patchwork)


#Microbiome Analysis

setwd("C:/Users/12697/Documents/Microbiome/Data Analysis")

ps <- readRDS('categorized_data.RDS')

sam <- data.frame(ps@sam_data) |>
  rename("Deworming Pill in the Last Year" = Dewormingin1yr)

sample_data(ps) <- sam

bdiv <- tibble(variable = colnames(sam), 
               jaccard_p = rep(NA, 128), 
               bray_p = rep(NA, 128))

bray <- phyloseq::distance(ps, method = "bray")
jaccard <- phyloseq::distance(ps, method = "jaccard")

for (i in 1:128) {
  try({
    var <- sam[[i]]
    keep <- !is.na(var)
    
    # Skip variables with too few valid samples or constant values
    if (sum(keep) < 3 || length(unique(var[keep])) < 2) {
      bdiv$bray_p[i] <- NA
      bdiv$jaccard_p[i] <- NA
    } else {
      # Clean variable and corresponding distance matrices
      var_clean <- var[keep]
      bray_sub <- as.dist(as.matrix(bray)[keep, keep])
      jaccard_sub <- as.dist(as.matrix(jaccard)[keep, keep])
      
      # Make a clean data frame for adonis2
      df <- data.frame(var = var_clean)
      
      # Bray–Curtis PERMANOVA
      permanova_bray <- vegan::adonis2(bray_sub ~ var, data = df, na.action = na.omit)
      bdiv$bray_p[i] <- permanova_bray$`Pr(>F)`[1]
      
      # Jaccard PERMANOVA
      permanova_jaccard <- vegan::adonis2(jaccard_sub ~ var, data = df, na.action = na.omit)
      bdiv$jaccard_p[i] <- permanova_jaccard$`Pr(>F)`[1]
    }
  })
}

sig <- bdiv |>
  filter(bray_p <= 0.05 & jaccard_p <= 0.05) |>
  pull(variable)

bdiv$adjusted_bray <- p.adjust(bdiv$bray_p, method = "BH")
bdiv$adjusted_jaccard <- p.adjust(bdiv$jaccard_p, method = "BH")


#Create PCOA plots for significant stuff

create_pcoa_plot <- function(variable, jaccard_dist, bray_dist, jaccard_pcoa, bray_pcoa, sam) {
  
  permanova_jaccard <- vegan::adonis2(jaccard_dist ~ sam[[variable]])
  permanova_bray <- vegan::adonis2(bray_dist ~ sam[[variable]])
  
  pval_jaccard <- permanova_jaccard$`Pr(>F)`[1]
  pval_bray <- permanova_bray$`Pr(>F)`[1]
  
  r2_jaccard <- permanova_jaccard$R2[1]
  r2_bray <- permanova_bray$R2[1]
  
  p_text_jaccard <- ifelse(pval_jaccard < 0.001, "p < 0.001",
                           paste0("p = ", format(round(pval_jaccard, 3), nsmall = 3)))
  p_text_bray <- ifelse(pval_bray < 0.001, "p < 0.001",
                        paste0("p = ", format(round(pval_bray, 3), nsmall = 3)))
  
  r2_text_jaccard <- paste0("R² = ", format(round(r2_jaccard, 3), nsmall = 3))
  r2_text_bray <- paste0("R² = ", format(round(r2_bray, 3), nsmall = 3))
  
  pcoa_jaccard_plot <- ggordpoint(
    obj = jaccard_pcoa,
    biplot = FALSE,
    speciesannot = TRUE,
    factorNames = c(variable),
    ellipse = TRUE,
    linesize = 1.5,
    ellipse_linewd = 1,
    ellipse_lty = 2
   ) +
  #   ggtitle(paste(gsub(pattern = "_", replacement = " ", variable), " (Jaccard)")) +
    ggtitle("Deworming Pill in the Last Year (Jaccard)") + 
    guides(color = guide_legend(
      title = gsub("_", " ", variable),
      override.aes = list(size = 4)
    )) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
    annotate("text", x = Inf, y = Inf,
             label = paste(p_text_jaccard, r2_text_jaccard, sep = " • "),
             hjust = 1.1, vjust = 1.5, size = 10, fontface = "bold")
  
  pcoa_bray_plot <- ggordpoint(
    obj = bray_pcoa,
    biplot = FALSE,
    speciesannot = TRUE,
    factorNames = c(variable),
    ellipse = TRUE,
    linesize = 1.5,
    ellipse_linewd = 1,
    ellipse_lty = 2
  ) +
    # ggtitle(paste(gsub(pattern = "_", replacement = " ", variable), " (Bray-Curtis)")) +
    ggtitle("Deworming Pill in the Last Year (Bray-Curtis)") + 
    guides(color = guide_legend(
      title = gsub("_", " ", variable),
      override.aes = list(size = 4)
    )) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
    annotate("text", x = Inf, y = Inf,
             label = paste(p_text_bray, r2_text_bray, sep = " • "),
             hjust = 1.1, vjust = 1.5, size = 10, fontface = "bold")
  
  combined_plot <- pcoa_jaccard_plot + pcoa_bray_plot
  
  ggsave(
    combined_plot,
    filename = paste(variable, "_pcoa_combined.png", sep = ""),
    device = "png",
    height = 6,
    width = 14,
    units = "in"
  )
  
  combined_plot
}

setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/beta diversity PCOA plots/Significant Combined")

prune_na_samples <- function(ps, variable) {
  # Convert sample data to a data frame
  sam <- as.data.frame(sample_data(ps))
  
  # Identify samples without NA for the variable
  keep <- !is.na(sam[[variable]])
  
  # Print how many samples were dropped
  message(sum(!keep), " samples removed due to NA in '", variable, "'")
  
  # Prune phyloseq object
  ps_pruned <- prune_samples(keep, ps)
  
  return(ps_pruned)
}



for (i in 1:5){
  setwd("C:/Users/12697/Documents/Microbiome/Data Analysis")
  ps <- readRDS('categorized_data.RDS')
  
  ps <- prune_na_samples(ps, sig[4])
  
  bray <- phyloseq::distance(ps, method = "bray")
  jaccard <- phyloseq::distance(ps, method = "jaccard")
  
  jaccard_pcoa <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")
  bray_pcoa <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")
  
  sam <- data.frame(ps@sam_data)
  setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/beta diversity PCOA plots/Significant Combined")
  try({create_pcoa_plot(sig[i], jaccard, bray, jaccard_pcoa, bray_pcoa, sam)})
}

try({create_pcoa_plot(sig[4], jaccard, bray, jaccard_pcoa, bray_pcoa, sam)})


