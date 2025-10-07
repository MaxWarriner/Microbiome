
# Libraries & Data --------------------------------------------------------

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



# Alpha Diversity ---------------------------------------------------------

div <- estimate_richness(ps)

sam$Shannon <- div$Shannon
sam$chao1 <- div$Chao1

# Socioeconomic Index
shannon_SES <- ggplot(data = sam, aes(x = as.factor(SES_index))) + 
  geom_boxplot(aes(y = Shannon)) + 
  theme_bw() + 
  xlab('Socioeconomic Index') + 
  ylab('Shannon Diversity')

chao1_SES <- ggplot(data = sam, aes(x = as.factor(SES_index))) + 
  geom_boxplot(aes(y = chao1)) +
  theme_bw() + 
  xlab('Socioeconomic Index') +
  ylab('Chao1 Diversity')

shannon_SES + chao1_SES

summary(lm(Shannon ~ SES_index, data = sam))
summary(lm(chao1 ~ SES_index, data = sam))

# Hygiene Index
shannon_hygiene <- ggplot(data = sam, aes(x = as.factor(hygiene_index))) + 
  geom_boxplot(aes(y = Shannon)) + 
  theme_bw() + 
  xlab('Hygiene Index') + 
  ylab('Shannon Diversity')

chao1_hygiene <- ggplot(data = sam, aes(x = as.factor(hygiene_index))) + 
  geom_boxplot(aes(y = chao1)) +
  theme_bw() + 
  xlab('Hygiene Index') +
  ylab('Chao1 Diversity')

shannon_hygiene + chao1_hygiene

summary(lm(Shannon ~ hygiene_index, data = sam)) # significant model with Shannon index
summary(lm(chao1 ~ hygiene_index, data = sam))

# Household Infrastructure Quality Index
shannon_HIQ <- ggplot(data = sam, aes(x = as.factor(HIQ_index))) + 
  geom_boxplot(aes(y = Shannon)) + 
  theme_bw() + 
  xlab('Household Infrastructure Index') + 
  ylab('Shannon Diversity')

chao1_HIQ <- ggplot(data = sam, aes(x = as.factor(HIQ_index))) + 
  geom_boxplot(aes(y = chao1)) +
  theme_bw() + 
  xlab('Household Infrastructure Index') +
  ylab('Chao1 Diversity')

shannon_HIQ + chao1_HIQ

summary(lm(Shannon ~ HIQ_index, data = sam))
summary(lm(chao1 ~ HIQ_index, data = sam))


#Household Sanitation Index

shannon_HSI <- ggplot(data = sam, aes(x = as.factor(HSI_index))) + 
  geom_boxplot(aes(y = Shannon)) + 
  theme_bw() + 
  xlab('Household Sanitation Index') + 
  ylab('Shannon Diversity')

chao1_HSI <- ggplot(data = sam, aes(x = as.factor(HSI_index))) + 
  geom_boxplot(aes(y = chao1)) +
  theme_bw() + 
  xlab('Household Sanitation Index') +
  ylab('Chao1 Diversity')

shannon_HSI + chao1_HSI

summary(lm(Shannon ~ HSI_index, data = sam))
summary(lm(chao1 ~ HSI_index, data = sam))


# Beta Diversity ----------------------------------------------------------

pcoares <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")

physeq_dist <- phyloseq::distance(ps, method = "bray")

create_pcoa_plot <- function(variable, pcoares, physeq_dist, ps, sam) {
  
  permanova <- vegan::adonis2(physeq_dist ~ sam[[variable]])
  p_value <- permanova$`Pr(>F)`[1]
  
  # Format p-value for display
  p_text <- ifelse(p_value < 0.001, "p < 0.001", 
                   paste("p =", format(round(p_value, 3), nsmall = 3)))
  
  # Create PCoA plot
  pcoaplot <- ggordpoint(obj = pcoares, biplot = FALSE, speciesannot = TRUE,
                         factorNames = c(variable), ellipse = TRUE, linesize = 1.5, 
                         ellipse_linewd = 1, ellipse_lty = 2) +
    ggtitle(gsub("_", " ", variable)) + 
    guides(color=guide_legend(title=gsub("_", " ", variable), 
                              override.aes = list(size = 4))) +
    theme( 
      legend.title = element_blank(), 
      legend.text = element_text(size = 20)) +
    # Add p-value annotation
    annotate("text", x = Inf, y = Inf, 
             label = p_text, 
             hjust = 1.1, vjust = 1.5, 
             size = 12, fontface = "bold")
  
  print(pcoaplot)
  
  # ggsave(pcoaplot,
  #        filename = paste(variable, "_pcoa_bray.png", sep = ""),
  #        device = "png",
  #        height = 6, width = 8, units = "in")
  
  return(pcoaplot)
}

create_pcoa_plot("SES_group", pcoares, physeq_dist, ps, sam)

create_pcoa_plot("hygiene_group", pcoares, physeq_dist, ps, sam)

create_pcoa_plot("HSI_group", pcoares, physeq_dist, ps, sam)

create_pcoa_plot("HIQ_group", pcoares, physeq_dist, ps, sam)

#nothing significant


# Differential Abundance --------------------------------------------------

SVs <- as.data.frame(otu_table(ps))         # Extract OTU table
SVs <- data.frame(t(SVs))
taxonomy <- as.data.frame(phyloseq::tax_table(ps))    # Extract taxonomy

create_volcano_plot <- function(ps_obj, tax_col, condition_col, metadata, SVs, taxonomy, comparison_title = NULL) {
  
  # Create comparison title if not provided
  if(is.null(comparison_title)) {
    comparison_title <- paste(top_values, collapse = " vs ")
  }
  
  # Combine SVs with taxonomy to match the specified level
  SVs_with_taxonomy <- as.data.frame(SVs) %>%
    rownames_to_column(var = "OTU_name")
  
  taxonomy <- as.data.frame(taxonomy) %>%
    rownames_to_column(var = "OTU_name") |>
    dplyr::select("OTU_name", tax_col)
  
  SVs_with_taxonomy <- SVs_with_taxonomy |>
    left_join(taxonomy, by = "OTU_name")
  
  # Aggregate the SVs data to the specified taxonomic level
  tax_level_SVs <- SVs_with_taxonomy %>%
    dplyr::select(-OTU_name) %>%
    group_by(!!sym(tax_col))
  
  tax_level_SVs <- tax_level_SVs |>
    summarise(across(everything(), ~sum(.x, na.rm = TRUE)))
  
  tax_level_SVs <- tax_level_SVs |>
    filter(!!sym(tax_col) != "none") %>%
    column_to_rownames(var = tax_col)
  
  # Convert to matrix for ALDEx2 analysis
  tax_level_SVs <- as.matrix(tax_level_SVs)
  
  condition <- as.character(as.factor(metadata[[condition_col]]))
  
  # Run ALDEx2 analysis
  aldex_data <- aldex(tax_level_SVs, conditions = condition, mc.samples = 138, test = "t", effect = TRUE)
  
  # Combine t-test and effect size results
  results <- data.frame(aldex_data)
  
  # Merge results with taxonomy
  results <- results %>%
    rownames_to_column(var = tax_col)
  
  # Create the volcano plot (before multiple testing correction)
  p <- results %>%
    mutate(Significant = if_else(wi.eBH < 0.05, TRUE, FALSE)) %>%
    mutate(Taxon = as.character(!!sym(tax_col))) %>%
    mutate(TaxonToPrint = if_else(wi.eBH < 0.05, paste(Taxon, "(", round(wi.eBH,3) , ")", sep = ""), "")) |>
    ggplot(aes(x = diff.btw, y = -log10(wi.eBH), color = Significant, label = TaxonToPrint)) +
    geom_text_repel(size = 4, nudge_y = 0.05, max.overlaps = Inf) +
    geom_point(alpha = 0.6, shape = 16) +
    theme_minimal() +  
    xlab("log2(fold change)") +
    ylab("-log10(P-value)") +
    theme(legend.position = "none") +
    ggtitle(paste("Taxonomic Level:", tax_col, "\nComparison:", comparison_title))
  
  return(p)
}

create_combined_volcano <- function(ps, factor, metadata, SVs, taxonomy){
  
  # Get comparison title (top 2 values)
  var_values <- metadata[[factor]]
  non_na_values <- var_values[!is.na(var_values)]
  value_counts <- table(non_na_values)
  top_values <- names(sort(value_counts, decreasing = TRUE))[1:2]
  comparison_title <- paste(top_values, collapse = " vs ")
  
  # Create plots with consistent comparison title
  volcano_plot_phylum <- create_volcano_plot(ps, "Phylum", factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_family <- create_volcano_plot(ps, "family", factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_genus <- create_volcano_plot(ps, "Genus", factor, metadata, SVs, taxonomy, comparison_title)
  
  # Combine plots with a single main title
  combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus) +
    plot_annotation(
      title = paste("Differential Abundance Analysis for", factor),
      subtitle = paste("Comparison:", comparison_title),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  # ggsave(combined_volcano_plot, 
  #        filename = paste(factor, "_volcano.pdf", sep = ""),
  #        device = "pdf",
  #        height = 8, width = 12, units = "in")
  
  return(combined_volcano_plot)
}


create_combined_volcano(ps, "SES_group", sam, SVs, taxonomy)
