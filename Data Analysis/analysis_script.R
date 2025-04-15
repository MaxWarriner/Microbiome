
#Socioeconomic, Environmental, and Demographic Factor Analysis on Gut Microbiome

#Modified Code from Christine Bi by Max Warriner: 
#Script is designed to be able to run through a large amount of variables in a phyloseq object with efficiency

# ## Part 0: data/phyloseq preprocessing -----------------------------------

#Import the phyloseq object 
#First you need to set working directory to source file location in the 'Session' tab above

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
  
  print(a_diversity_factor)
  
  return(a_diversity_factor)
}

# #run through all the factors in the dataset
# for (var in factors){
#   try(create_a_diversity_plot(ps,var)) #try() will continue running the loop even if there's an error(if certain variables can't work with the plot)
# }

# Demographic Factors
a_diversity_age <- create_a_diversity_plot(ps, "Age")
a_diversity_weight <- create_a_diversity_plot(ps, "Weight")
a_diversity_vaccine <- create_a_diversity_plot(ps, "Vaccinated")
combined_demographic <-  a_diversity_age + a_diversity_weight + a_diversity_vaccine + 
  plot_annotation(
    title = "Demographic Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_demographic,
       filename = "demographic_alpha_diversity",
       device = "jpg",
       height = 12, width = 28, units = "in")

# Lifestyle Factors
a_diversity_raw_veg <- create_a_diversity_plot(ps, "Frequency_of_Eating_Raw.Undercooked_Vegetables")
a_diversity_washing_hands_method_eating <- create_a_diversity_plot(ps, "Method_of_Washing_Hands_Before_Eating")
a_diversity_use_school_latrine <- create_a_diversity_plot(ps, "Frequency_of_Using_School_Latrine")
a_diversity_bathing_river <- create_a_diversity_plot(ps, "Frequency_of_Bathing_in_River")
combined_lifestyle <- a_diversity_raw_veg + a_diversity_washing_hands_method_eating + a_diversity_use_school_latrine + 
   a_diversity_bathing_river +
  plot_annotation(
    title = "Lifestyle Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_lifestyle,
       filename = "lifestyle_alpha_diversity",
       device = "jpg",
       height = 18, width = 32, units = "in")

# Clinical Factors
a_diversity_fever <- create_a_diversity_plot(ps, "Fever_in_Last_Two_Weeks")
a_diversity_diarrhea <- create_a_diversity_plot(ps, "Diarrhea_in_Last_Two_Weeks")
a_diversity_wheezing <- create_a_diversity_plot(ps, "How_Many_Times_Wheezing_or_Whistling")
combined_clinical <- a_diversity_fever + a_diversity_diarrhea + a_diversity_wheezing + 
  plot_annotation(
    title = "Clinical Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))

ggsave(combined_clinical,
       filename = "clinical_alpha_diversity",
       device = "jpg",
       height = 12, width = 32, units = "in")


# Household Factors
a_diversity_kitchen_material <- create_a_diversity_plot(ps, "Kitchen_Material")
a_diversity_house_material <- create_a_diversity_plot(ps, "House_Floor_Material")
a_diversity_younger_siblings <- create_a_diversity_plot(ps, "Number_of_Siblings_Younger_than_12")
a_diversity_pet <- create_a_diversity_plot(ps, "Pet_in_House.Compound")
a_diversity_wood <- create_a_diversity_plot(ps, "Wood_as_Cooking_Fuel")
combined_household <- a_diversity_kitchen_material + a_diversity_house_material +
  a_diversity_younger_siblings + a_diversity_pet + a_diversity_wood + 
  plot_annotation(
    title = "Household Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_household,
       filename = "household_alpha_diversity",
       device = "jpg",
       height = 16, width = 32, units = "in")


# ## Part II (a) Beta diversity (Jaccard and Bray Distance) ----------------------

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
  
  # Get final sample data
  sam_data_final <- as(sample_data(ps_final), "data.frame")
  
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps_final, distmethod = "bray", method = "hellinger")
  
  # Perform PERMANOVA to get p-value
  physeq_dist <- phyloseq::distance(ps_final, method = "bray")
  permanova <- vegan::adonis2(physeq_dist ~ sam_data_final[[variable]])
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
    theme(title = element_text(size=4,face="bold"), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 32)) +
    # Add p-value annotation
    annotate("text", x = Inf, y = Inf, 
             label = p_text, 
             hjust = 1.1, vjust = 1.5, 
             size = 14, fontface = "bold")
  
  print(pcoaplot)
  
  # ggsave(pcoaplot, 
  #        filename = paste(variable, "_pcoa_bray.pdf", sep = ""),
  #        device = "pdf",
  #        height = 6, width = 8, units = "in")
  
  return(pcoaplot)
}
# for (var in factors){
#   try(pcoa_factor <- create_pcoa_plot(var))
# }

# ##(b) Beta diversity using Jaccard distance

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


# for (var in factors){
#   try(pcoa_factor <- create_pcoa_plot(var))
# }


# Demographic Factors
beta_diversity_age <- create_pcoa_plot("Age")
beta_diversity_weight <- create_pcoa_plot("Weight")
beta_diversity_vaccine <- create_pcoa_plot("Vaccinated")
beta_diversity_bcg <- create_pcoa_plot("BCG_Scar")
combined_demographic <- beta_diversity_age  + beta_diversity_weight  + beta_diversity_vaccine + beta_diversity_bcg + 
  plot_annotation(
    title = "Demographic Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_demographic,
       filename = "demographic_beta_diversity",
       device = "jpg",
       height = 24, width = 24, units = "in")

# Lifestyle Factors
beta_diversity_trimming <- create_pcoa_plot("Child.s_Fingernails_Trimmed")
beta_diversity_footwear <- create_pcoa_plot("Prefer_Sandals_or_Barefoot_in_House")
beta_diversity_raw_veg <- create_pcoa_plot("Frequency_of_Eating_Raw.Undercooked_Vegetables")
beta_diversity_washing_hands_method_eating <- create_pcoa_plot("Method_of_Washing_Hands_Before_Eating")
beta_diversity_use_school_latrine <- create_pcoa_plot("Frequency_of_Using_School_Latrine")
beta_diversity_clothes_river <- create_pcoa_plot("Frequency_of_Clothes_Washing_in_River")
combined_lifestyle <- beta_diversity_trimming + beta_diversity_washing_hands_method_eating +
  beta_diversity_footwear + beta_diversity_raw_veg + beta_diversity_use_school_latrine + beta_diversity_clothes_river  + 
  plot_annotation(
    title = "Lifestyle Factors") &
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5))
ggsave(combined_lifestyle,
       filename = "lifestyle_beta_diversity",
       device = "jpg",
       height = 32, width = 48, units = "in")

# Clinical Factors
beta_diversity_fever <- create_pcoa_plot("Fever_in_Last_Two_Weeks")
beta_diversity_wheezing <- create_pcoa_plot("How_Many_Times_Wheezing_or_Whistling")
combined_clinical <- beta_diversity_fever + beta_diversity_wheezing + 
  plot_annotation(
    title = "Clinical Factors") &
  theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5))

ggsave(combined_clinical,
       filename = "clinical_beta_diversity",
       device = "jpg",
       height = 12, width = 32, units = "in")


# Environmental Factors
beta_diversity_rural_urban <- create_pcoa_plot("Urban.Rural")
combined_environmental <- beta_diversity_rural_urban + 
  plot_annotation(
    title = "Environmental Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))

ggsave(combined_environmental,
       filename = "environmental_beta_diversity",
       device = "jpg",
       height = 12, width = 12, units = "in")

# Household Factors
beta_diversity_kitchen_material <- create_pcoa_plot("Kitchen_Material")
beta_diversity_house_material <- create_pcoa_plot("House_Floor_Material")
beta_diversity_younger_siblings <- create_pcoa_plot("Number_of_Siblings_Younger_than_12")
beta_diversity_animals <- create_pcoa_plot("No_Animals_in_House.Compound")
beta_diversity_pet <- create_pcoa_plot("Pet_in_House.Compound")
beta_diversity_wood <- create_pcoa_plot("Wood_as_Cooking_Fuel")
combined_household <- beta_diversity_kitchen_material + beta_diversity_house_material +
  beta_diversity_younger_siblings + beta_diversity_animals + beta_diversity_pet + beta_diversity_wood + 
  plot_annotation(
    title = "Household Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_household,
       filename = "household_alpha_diversity",
       device = "jpg",
       height = 16, width = 38, units = "in")






##Calculations of significance:

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
  
  # ggsave(barplot, 
  #        filename = paste(variable, "_phylum_barplot.pdf", sep = ""),
  #        device = "pdf",
  #        height = 6, width = 8, units = "in")
  
  return(barplot)
}

# for (variable in factors){
#   try(create_phylum_barplot(variable))
# }


# Demographic Factors
phylum_age <- create_phylum_barplot( "Age")
phylum_weight <- create_phylum_barplot( "Weight")
phylum_vaccine <- create_phylum_barplot( "Vaccinated")
combined_demographic <-  phylum_age + phylum_weight + phylum_vaccine + 
  plot_annotation(
    title = "Demographic Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_demographic,
       filename = "demographic_phylum",
       device = "jpg",
       height = 12, width = 28, units = "in")

# Lifestyle Factors
phylum_raw_veg <- create_phylum_barplot( "Frequency_of_Eating_Raw.Undercooked_Vegetables")
phylum_washing_hands_method_eating <- create_phylum_barplot( "Method_of_Washing_Hands_Before_Eating")
phylum_use_school_latrine <- create_phylum_barplot( "Frequency_of_Using_School_Latrine")
phylum_bathing_river <- create_phylum_barplot( "Frequency_of_Bathing_in_River")
combined_lifestyle <- phylum_raw_veg + phylum_washing_hands_method_eating + phylum_use_school_latrine + 
  phylum_bathing_river +
  plot_annotation(
    title = "Lifestyle Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_lifestyle,
       filename = "lifestyle_phylum",
       device = "jpg",
       height = 18, width = 32, units = "in")

# Clinical Factors
phylum_fever <- create_phylum_barplot( "Fever_in_Last_Two_Weeks")
phylum_diarrhea <- create_phylum_barplot( "Diarrhea_in_Last_Two_Weeks")
phylum_wheezing <- create_phylum_barplot( "How_Many_Times_Wheezing_or_Whistling")
combined_clinical <- phylum_fever + phylum_diarrhea + phylum_wheezing + 
  plot_annotation(
    title = "Clinical Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))

ggsave(combined_clinical,
       filename = "clinical_phylum",
       device = "jpg",
       height = 12, width = 32, units = "in")


# Household Factors
phylum_kitchen_material <- create_phylum_barplot( "Kitchen_Material")
phylum_house_material <- create_phylum_barplot( "House_Floor_Material")
phylum_younger_siblings <- create_phylum_barplot( "Number_of_Siblings_Younger_than_12")
phylum_pet <- create_phylum_barplot( "Pet_in_House.Compound")
phylum_wood <- create_phylum_barplot( "Wood_as_Cooking_Fuel")
combined_household <- phylum_kitchen_material + phylum_house_material +
  phylum_younger_siblings + phylum_pet + phylum_wood + 
  plot_annotation(
    title = "Household Factors") &
  theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5))
ggsave(combined_household,
       filename = "household_phylum",
       device = "jpg",
       height = 16, width = 32, units = "in")


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

# Function to create volcano plots for taxonomic levels

metadata <- as.data.frame(sample_data(ps))  # Extract metadata
SVs <- as.data.frame(otu_table(ps))         # Extract OTU table
SVs <- data.frame(t(SVs))
taxonomy <- as.data.frame(phyloseq::tax_table(ps))    # Extract taxonomy

# Before Multiple Testing Correction

create_volcano_plot <- function(ps_obj, taxlevel, condition_col, metadata, SVs, taxonomy, comparison_title = NULL) {
  
  # Get the variable values, remove NAs, and count frequencies
  var_values <- metadata[[condition_col]]
  non_na_values <- var_values[!is.na(var_values)]
  value_counts <- table(non_na_values)
  
  # Get top 2 most frequent non-NA values
  top_values <- names(sort(value_counts, decreasing = TRUE))[1:2]
  
  # Create comparison title if not provided
  if(is.null(comparison_title)) {
    comparison_title <- paste(top_values, collapse = " vs ")
  }
  
  # Subset data to only include samples with top 2 non-NA values
  subset_idx <- var_values %in% top_values & !is.na(var_values)
  metadata <- metadata[subset_idx, , drop = FALSE]
  
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
  
  row_names_sample <- rownames(metadata)
  matching_cols <- c("OTU_name", colnames(SVs_with_taxonomy)[ncol(SVs_with_taxonomy)] , paste0("X", gsub("-", ".", row_names_sample)))
  
  # Filter df1
  SVs_with_taxonomy <- SVs_with_taxonomy[, colnames(SVs_with_taxonomy) %in% matching_cols, drop = FALSE]
  
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
    mutate(Significant = if_else(wi.ep < 0.05, TRUE, FALSE)) %>%
    mutate(Taxon = as.character(!!sym(tax_col))) %>%
    mutate(TaxonToPrint = if_else(wi.ep < 0.05, paste(Taxon, "(", round(wi.ep,3) , ")", sep = ""), "")) |>
    ggplot(aes(x = diff.btw, y = -log10(wi.ep), color = Significant, label = TaxonToPrint)) +
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
  volcano_plot_phylum <- create_volcano_plot(ps, 2, factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_family <- create_volcano_plot(ps, 5, factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_genus <- create_volcano_plot(ps, 6, factor, metadata, SVs, taxonomy, comparison_title)
  
  # Combine plots with a single main title
  combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus) +
    plot_annotation(
      title = paste("Differential Abundance Analysis for", factor),
      subtitle = paste("Comparison:", comparison_title),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  ggsave(combined_volcano_plot, 
         filename = paste(factor, "_volcano.pdf", sep = ""),
         device = "pdf",
         height = 8, width = 12, units = "in")
  
  return(combined_volcano_plot)
}

for (variable in factors){
  try(create_combined_volcano(ps, variable, metadata, SVs, taxonomy))
}





#After Multiple Testing Correction

#Volcano plots of differentially abundant taxa according to ALEDX2 results between positive and negative at the phylum, family, and genus levels (BEFORE multiple testing correction)


create_volcano_plot <- function(ps_obj, taxlevel, condition_col, metadata, SVs, taxonomy, comparison_title = NULL) {
  
  # Get the variable values, remove NAs, and count frequencies
  var_values <- metadata[[condition_col]]
  non_na_values <- var_values[!is.na(var_values)]
  value_counts <- table(non_na_values)
  
  # Get top 2 most frequent non-NA values
  top_values <- names(sort(value_counts, decreasing = TRUE))[1:2]
  
  # Create comparison title if not provided
  if(is.null(comparison_title)) {
    comparison_title <- paste(top_values, collapse = " vs ")
  }
  
  # Subset data to only include samples with top 2 non-NA values
  subset_idx <- var_values %in% top_values & !is.na(var_values)
  metadata <- metadata[subset_idx, , drop = FALSE]
  
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
  
  row_names_sample <- rownames(metadata)
  matching_cols <- c("OTU_name", colnames(SVs_with_taxonomy)[ncol(SVs_with_taxonomy)] , paste0("X", gsub("-", ".", row_names_sample)))
  
  # Filter df1
  SVs_with_taxonomy <- SVs_with_taxonomy[, colnames(SVs_with_taxonomy) %in% matching_cols, drop = FALSE]
  
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
  volcano_plot_phylum <- create_volcano_plot(ps, 2, factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_family <- create_volcano_plot(ps, 5, factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_genus <- create_volcano_plot(ps, 6, factor, metadata, SVs, taxonomy, comparison_title)
  
  # Combine plots with a single main title
  combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus) +
    plot_annotation(
      title = paste("Differential Abundance Analysis for", factor),
      subtitle = paste("Comparison:", comparison_title),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  ggsave(combined_volcano_plot, 
         filename = paste(factor, "_volcano.pdf", sep = ""),
         device = "pdf",
         height = 8, width = 12, units = "in")
  
  return(combined_volcano_plot)
}

for (variable in factors){
  try(create_combined_volcano(ps, variable, metadata, SVs, taxonomy))
}




# ## Part V: Differential Abundance Using Lefse and Edger  -----------------------------

#Adjust cutoff values if needed
#Note: If it gives error like "Error in names(marker) <- `*vtmp*` : attempt to set an attribute on NULL", it means no marker was identified, and you will have to increase your wilcoxon/kw/pvalue cutoffs and/or decrease your lda_cutoff.

create_lefse_edgar <- function(ps, factor){
  
  # Extract sample data as a data frame
  sample_df <- data.frame(sample_data(ps))
  
  # Get the variable of interest
  var_vector <- sample_df[[factor]]
  
  # Exclude NAs and get top two most frequent categories
  var_vector_no_na <- var_vector[!is.na(var_vector)]
  if(length(unique(var_vector_no_na)) < 2) {
    message(paste("Skipping", factor, "- fewer than 2 categories after removing NAs"))
    return(NULL)
  }
  
  top_categories <- names(sort(table(var_vector_no_na), decreasing = TRUE))[1:2]
  
  # Get sample names to keep
  keep_samples <- rownames(sample_df)[!is.na(var_vector) & var_vector %in% top_categories]
  
  # Subset the phyloseq object using sample names
  ps_subset <- prune_samples(sample_names(ps) %in% keep_samples, ps)
  
  # Verify we have exactly two categories
  if(length(unique(sample_data(ps_subset)[[factor]])) != 2) {
    message(paste("Skipping", factor, "- doesn't have exactly 2 categories after subsetting"))
    return(NULL)
  }
  
  # Initialize plot objects as NULL
  lefse_bar_plot <- NULL
  lefse_dot_plot <- NULL
  edger_bar_plot <- NULL
  
  # Try to run LefSe
  mm_lefse <- tryCatch(
    run_lefse(
      ps_subset,
      wilcoxon_cutoff = 0.01,
      group = factor,
      kw_cutoff = 0.01,
      multigrp_strat = TRUE,
      lda_cutoff = 2
    ),
    error = function(e) {
      message(paste("LefSe failed for", factor, ":", e$message))
      return(NULL)
    }
  )
  
  # Create LefSe plots if successful
  if(!is.null(mm_lefse)) {
    lefse_bar_plot <- tryCatch(
      plot_ef_bar(mm_lefse),
      error = function(e) {
        message(paste("Failed to create LefSe bar plot:", e$message))
        NULL
      }
    )
    
    lefse_dot_plot <- tryCatch(
      plot_ef_dot(mm_lefse),
      error = function(e) {
        message(paste("Failed to create LefSe dot plot:", e$message))
        NULL
      }
    )
  }
  
  # Try to run edgeR
  mm_edger <- tryCatch(
    run_edger(
      ps_subset,
      group = factor,
      pvalue_cutoff = 0.05,
      p_adjust = "fdr"
    ),
    error = function(e) {
      message(paste("edgeR failed for", factor, ":", e$message))
      return(NULL)
    }
  )
  
  # Create edgeR plot if successful
  if(!is.null(mm_edger)) {
    edger_bar_plot <- tryCatch(
      plot_ef_bar(mm_edger),
      error = function(e) {
        message(paste("Failed to create edgeR bar plot:", e$message))
        NULL
      }
    )
  }
  
  # Create placeholder plots if any are missing
  if(is.null(lefse_bar_plot)) {
    lefse_bar_plot <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "LefSe analysis failed", size = 6) + 
      theme_void()
  }
  
  if(is.null(edger_bar_plot)) {
    edger_bar_plot <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "edgeR analysis failed", size = 6) + 
      theme_void()
  }
  
  if(is.null(lefse_dot_plot)) {
    lefse_dot_plot <- ggplot() + 
      annotate("text", x = 1, y = 1, label = "LefSe dot plot unavailable", size = 6) + 
      theme_void()
  }
  
  # Create combined plot layout
  combined_plot <- (lefse_bar_plot | edger_bar_plot) / lefse_dot_plot +
    plot_annotation(
      title = paste("LefSe/Edger Plots for:", factor),
      subtitle = paste("Top two categories:", paste(top_categories, collapse = " vs ")),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  print(combined_plot)
  
  # Save plot
  tryCatch({
    ggsave(combined_plot, 
           filename = paste("lefse_edgar_plot_", factor, ".pdf", sep = ""),
           device = "pdf",
           height = 12, width = 16, units = "in")
  }, error = function(e) {
    message(paste("Failed to save plot for", factor, ":", e$message))
  })
  
  return(combined_plot)
}

# Run for each factor
for (variable in factors) {
  try(create_lefse_edgar(ps, variable))
}

