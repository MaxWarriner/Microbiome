
#Socioeconomic, Environmental, and Demographic Factor Analysis on Gut Microbiome

#Modified Code from Christine Bi by Max Warriner

# ## Part 0: data/phyloseq preprocessing -----------------------------------

#Import the phyloseq object 
#First you need to set working directory to source file location in the 'Session' tab above

library(BiocManager)
library(phyloseq)
library(tidyverse)
library(dplyr)

ps <- readRDS("categorized_data.RDS") # Load in using whatever file name you have

#Rarefy Data
set.seed(13)
ps <- rarefy_even_depth(ps, 1e+05)

#Changing the taxonomic table labels (From Rank 1~Rank 7 to Kingdom~Species)
# Access the existing taxonomy table from the phyloseq object
tax_table_ps <- tax_table(ps)[,1:7]

# Assign the modified taxonomy table back to the phyloseq object
tax_table(ps) <- tax_table_ps

#define list of factors to be tested
factors <- colnames(ps@sam_data)

# ## Part I : Alpha Diversity (Chao1 and Shannon indices) -----------------

library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(rlang)


create_a_diversity_plot <- function(ps,variable){ #function to create an alpha diversity boxplot
  
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
  
  #create the plot
  a_diversity_factor<- plot_richness(ps_final, x=variable, color=variable, measures=c("Chao1", "Shannon"))
  a_diversity_factor$layers[[2]] = NULL 
  a_diversity_factor <- a_diversity_factor  + geom_boxplot() + theme_bw()
  
  print(a_diversity_factor) #sneak-peek
  
  #save the file as a pdf with big dimensions
  ggsave(a_diversity_factor, 
         filename = paste(variable,"_boxplot.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 12, units = "in")
  
}

#run through all the factors in the dataset
for (var in factors){
  try(create_a_diversity_plot(ps,var)) #try() will continue running the loop even if there's an error(if certain variables can't work with the plot)
}


# ## Part II (a) Beta diversity (Jaccard and Bray Distance) ----------------------


library(MicrobiotaProcess)
# Define the function to create PCoA plots
create_pcoa_plot <- function(variable) {
  
  # Convert sample data to data frame
  sam_data <- as(sample_data(ps), "data.frame")
  
  # First remove NAs
  keep_samples <- !is.na(sam_data[[variable]])
  ps_filtered <- prune_samples(keep_samples, ps)
  
  # Get updated sample data after NA removal
  sam_data_filtered <- as(sample_data(ps_filtered), "data.frame")
  
  # Count occurrences of each value in the variable
  value_counts <- table(sam_data_filtered[[variable]])
  
  # Identify values that appear at least 3 times (PCOA plot can't draw ellipse with < 3 points)
  valid_values <- names(value_counts)[value_counts >= 3]
  
  # Filter samples to only include those with valid values
  final_keep_samples <- sam_data_filtered[[variable]] %in% valid_values
  ps_final <- prune_samples(final_keep_samples, ps_filtered)
  
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps_final, distmethod = "bray", method = "hellinger")
  
  # Create PCoA plot
  pcoaplot <- ggordpoint(obj = pcoares, biplot = FALSE, speciesannot = TRUE,
                         factorNames = c(variable), ellipse = TRUE, linesize = 1.5)
  print(pcoaplot)
  
  ggsave(pcoaplot, 
         filename = paste(variable, "_pcoa_bray.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(pcoaplot)
}

for (var in factors){
  try(pcoa_factor <- create_pcoa_plot(var))
}

# ##(b) Beta diversity using Jaccard distance


library(MicrobiotaProcess)
# Define the function to create PCoA plots
create_pcoa_plot <- function(variable) {
  
  # Convert sample data to data frame
  sam_data <- as(sample_data(ps), "data.frame")
  
  # First remove NAs
  keep_samples <- !is.na(sam_data[[variable]])
  ps_filtered <- prune_samples(keep_samples, ps)
  
  # Get updated sample data after NA removal
  sam_data_filtered <- as(sample_data(ps_filtered), "data.frame")
  
  # Count occurrences of each value in the variable
  value_counts <- table(sam_data_filtered[[variable]])
  
  # Identify values that appear at least 3 times (PCOA plot can't draw ellipse with < 3 points)
  valid_values <- names(value_counts)[value_counts >= 3]
  
  # Filter samples to only include those with valid values
  final_keep_samples <- sam_data_filtered[[variable]] %in% valid_values
  ps_final <- prune_samples(final_keep_samples, ps_filtered)
  
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps_final, distmethod = "jaccard", method = "hellinger")
  
  # Create PCoA plot
  pcoaplot <- ggordpoint(obj = pcoares, biplot = FALSE, speciesannot = TRUE,
                         factorNames = c(variable), ellipse = TRUE)
  
  ggsave(pcoaplot, 
         filename = paste(variable, "_pcoa_jaccard.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(pcoaplot)
}


for (var in factors){
  try(pcoa_factor <- create_pcoa_plot(var))
}


##Calculations of significance:

library(vegan)
library(dplyr)

calculate_permanova <- function(ps) {
  # Calculate both distance matrices
  bray_dist <- get_dist(ps, distmethod = "bray")
  jaccard_dist <- get_dist(ps, distmethod = "jaccard")
  
  # Get sample data
  sampleda <- data.frame(sample_data(ps), check.names = FALSE)
  
  # Initialize results data frame
  results <- data.frame()
  
  # Set seed for reproducibility
  set.seed(1024)
  
  # Process each distance matrix
  for (dist_name in c("bray", "jaccard")) {
    # Get the appropriate distance matrix
    dist_matrix <- if (dist_name == "bray") bray_dist else jaccard_dist
    
    # Ensure sample data ordering matches distance matrix
    sampleda_ordered <- sampleda[match(colnames(as.matrix(dist_matrix)), rownames(sampleda)), , drop = FALSE]
    
    # Test each variable
    for (variable in names(sampleda_ordered)) {
      # Get the variable values, remove NAs, and count frequencies
      var_values <- sampleda_ordered[[variable]]
      non_na_values <- var_values[!is.na(var_values)]
      value_counts <- table(non_na_values)
      
      # Skip if there are fewer than 2 unique non-NA values
      if (length(value_counts) < 2) {
        message("Skipping ", variable, " with ", dist_name, " distance - fewer than 2 non-NA categories")
        next
      }
      
      # Get top 2 most frequent non-NA values
      top_values <- names(sort(value_counts, decreasing = TRUE))[1:2]
      
      # Subset data to only include samples with top 2 non-NA values
      subset_idx <- var_values %in% top_values & !is.na(var_values)
      dist_subset <- as.dist(as.matrix(dist_matrix)[subset_idx, subset_idx])
      sampleda_subset <- sampleda_ordered[subset_idx, , drop = FALSE]
      
      # Ensure we have at least 2 samples per group
      if (min(table(sampleda_subset[[variable]])) < 2) {
        message("Skipping ", variable, " with ", dist_name, " distance - insufficient samples in one of the top groups")
        next
      }
      
      # Record which values are being compared
      comparison_note <- paste("Comparing:", paste(top_values, collapse = " vs "))
      
      # Try running PERMANOVA
      test_result <- tryCatch({
        formula <- as.formula(paste("dist_subset ~", variable))
        adores <- adonis2(formula, data = sampleda_subset, permutations = 999)
        
        data.frame(
          variable = variable,
          distance_method = dist_name,
          comparison = comparison_note,
          F_statistic = adores[1, "F"],
          p_value = adores[1, "Pr(>F)"],
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message("Skipped ", variable, " with ", dist_name, " distance (Error: ", e$message, ")")
        return(NULL)
      })
      
      # Add to results if successful
      if (!is.null(test_result)) {
        results <- bind_rows(results, test_result)
      }
    }
  }
  
  # Add significance stars
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(
        significance = case_when(
          p_value <= 0.01 ~ "***",
          p_value > 0.01 & p_value <= 0.05 ~ "**",
          p_value > 0.05 & p_value <= 0.1 ~ "*",
          p_value > 0.05 ~ ""
        )
      )
  }
  
  return(results)
}

permanova_results <- calculate_permanova(ps)

permanova_results <- permanova_results |>
  arrange("p_value")

write.csv(permanova_results, "PERMANOVA_results.csv")


# ## Part III Overall microbial abundance "at the Phylum, Genus and Family Level --------

# Define the function to create bar plots with titles
create_phylum_barplot <- function(variable) {
  
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
  
  phylumtaxa <- get_taxadf(obj = ps_final, taxlevel = 2)
  barplot <- ggbartax(obj = phylumtaxa, facetNames = variable, plotgroup = TRUE, topn = 5) +
    xlab(NULL) +
    ylab("relative abundance (%)") +
    labs(title = variable) +  # Add variable name as the title
    scale_fill_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(31))) +
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 2))
  
  print(barplot)
  
  ggsave(barplot, 
         filename = paste(variable, "_phylum_barplot.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(barplot)
}

for (variable in factors){
  try(create_phylum_barplot(variable))
}




# Define the function to create bar plots at the family level
create_family_barplot <- function(variable) {
  
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
  
  familytaxa <- get_taxadf(obj = ps_final, taxlevel = 5)
  barplot <- ggbartax(obj = familytaxa, facetNames = variable, plotgroup = TRUE, topn = 10) +
    xlab(NULL) +
    ylab("Relative Abundance (%)") +
    labs(title = variable) +  # Add variable name as the title
    scale_fill_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(31))) +
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 2))
  
  ggsave(barplot, 
         filename = paste(variable, "_family_barplot.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(barplot)
}

for (variable in factors){
  try(create_family_barplot(variable))
}


##Overall microbial abundance "at the Genus level" across different conditions (in one panel)


# Define the function to create bar plots at the genus level
create_genus_barplot <- function(variable) {
  
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
  
  genustaxa <- get_taxadf(obj = ps_final, taxlevel = 6)
  barplot <- ggbartax(obj = genustaxa, facetNames = variable, plotgroup = TRUE, topn = 15) +
    xlab(NULL) +
    ylab("Relative Abundance (%)") +
    labs(title = variable, fill = NULL) +  # Add variable name as the title and remove the legend title
    scale_fill_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(31))) +
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 2)) +
    theme_classic() +
    theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
  
  ggsave(barplot, 
         filename = paste(variable, "_genus_barplot.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(barplot)
}

for (variable in factors){
  try(create_genus_barplot(variable))
}



# ## Part IV: Differential Abundance Analysis using ALDEx2 and Volc --------


#BEFORE multiple testing correction
#Volcano plots of differentially abundant taxa according to ALEDX2 results between positive and negative at the phylum, family, and genus levels (BEFORE multiple testing correction)


library(ALDEx2)
# Function to create volcano plots for taxonomic levels

metadata <- as.data.frame(sample_data(ps))  # Extract metadata
SVs <- as.data.frame(otu_table(ps))         # Extract OTU table
SVs <- data.frame(t(SVs))
taxonomy <- as.data.frame(phyloseq::tax_table(ps))    # Extract taxonomy

create_volcano_plot <- function(ps_obj, taxlevel, condition_col, metadata, SVs, taxonomy) {
  
  # Ensure that taxonomy contains the specified level
  tax_col <- colnames(taxonomy)[taxlevel]
  taxonomy <- taxonomy %>%
    rownames_to_column(var = "OTU_name") %>%
    dplyr::select(OTU_name, !!sym(tax_col)) %>%
    distinct()
  
  # Combine SVs with taxonomy to match the specified level
  SVs_with_taxonomy <- as.data.frame(SVs) %>%
    rownames_to_column(var = "OTU_name")
  
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
    mutate(Significant = if_else(we.ep < 0.05, TRUE, FALSE)) %>%
    mutate(Taxon = as.character(!!sym(tax_col))) %>%
    mutate(TaxonToPrint = if_else(we.ep < 0.05, paste(Taxon, "(", round(we.ep,3) , ")", sep = ""), "")) |>
    ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant, label = TaxonToPrint)) +
    geom_text_repel(size = 4, nudge_y = 0.05, max.overlaps = Inf) +  # Increase max.overlaps
    geom_point(alpha = 0.6, shape = 16) +
    theme_minimal() +  
    xlab("log2(fold change)") +
    ylab("-log10(P-value)") +
    theme(legend.position = "none") +
    ggtitle(paste("Factor:", condition_col,"-- Taxonomic Level:", tax_col))
  
  
  return(p)
}

create_combined_volcano <- function(ps, factor, metadata, SVs, taxonomy){
  
  
  volcano_plot_phylum <- create_volcano_plot(ps, 2, factor, metadata, SVs, taxonomy)
  volcano_plot_family <- create_volcano_plot(ps, 5, factor, metadata, SVs, taxonomy)
  volcano_plot_genus <- create_volcano_plot(ps, 6, factor, metadata, SVs, taxonomy)
  combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus)
  
  ggsave(combined_volcano_plot, 
         filename = paste(factor, "_volcano.pdf", sep = ""),
         device = "pdf",
         height = 8, width = 12, units = "in")
  
}

create_combined_volcano(ps, "Sex", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Modeofdelivery", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Everhadvaccinated", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Childsfingernailtrimmed", metadata, SVs, taxonomy)
create_combined_volcano(ps, "BCGscar", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Arechildsfingernailsdirty", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Didyourparentsorhealthprofessionalsgaveyouotherantibiotics", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill", metadata, SVs, taxonomy)





## AFTER multiple testing correction




# ## Part V: Differential Abundance Using Lefse -----------------------------


library(microbiomeMarker)

#Adjust cutoff values if needed
#Note: If it gives error like "Error in names(marker) <- `*vtmp*` : attempt to set an attribute on NULL", it means no marker was identified, and you will have to increase your wilcoxon/kw/pvalue cutoffs and/or decrease your lda_cutoff.

create_lefse_edgar <- function(ps, factor){
  
  # Run LefSe
  mm_lefse <- run_lefse(
    ps,
    wilcoxon_cutoff = 0.05,
    group = factor,
    kw_cutoff = 0.05,
    multigrp_strat = TRUE,
    lda_cutoff = 2
  )
  
  # Run edgeR
  mm_edger <- run_edger(
    ps,
    group = factor,
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
  )
  
  # Create plots
  lefse_bar_plot <- plot_ef_bar(mm_lefse)
  lefse_dot_plot <- plot_ef_dot(mm_lefse)
  edger_bar_plot <- plot_ef_bar(mm_edger)
  #lefse_cladogram <- (plot_cladogram(mm_lefse, color = c(Negative = "darkgreen", Positive = "red")) +
  #theme(plot.margin = margin(0, 0, 0, 0))
  
  combined_plot <- (lefse_bar_plot | edger_bar_plot) / lefse_dot_plot
  
  print(combined_plot)
  
  ggsave(combined_plot, 
         filename = paste("lefse_edgar_plot_",factor, ".pdf", sep = ""),
         device = "pdf",
         height = 12, width = 16, units = "in")
  
}

create_lefse_edgar(ps, "Age")
create_lefse_edgar(ps, "Sex")
create_lefse_edgar(ps, "Modeofdelivery")
create_lefse_edgar(ps, "Everhadvaccinated")
create_lefse_edgar(ps, "Childsfingernailtrimmed")
create_lefse_edgar(ps, "BCGscar")
create_lefse_edgar(ps, "Arechildsfingernailsdirty")
create_lefse_edgar(ps, "Didyourparentsorhealthprofessionalsgaveyouotherantibiotics")
create_lefse_edgar(ps, "Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill")
