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
sam <- data.frame(ps@sam_data) |>
  rename("Frequency_of_Eating_Raw_or_Undercooked_Vegetables" = Frequency_of_Eating_Raw.Undercooked_Vegetables)
sample_data(ps) <- sam

ps@sam_data$Frequency_of_Using_School_Latrine <- factor(ps@sam_data$Frequency_of_Using_School_Latrine, levels = c("never", "sometimes", "always"))

create_a_diversity_plot <- function(ps, variable, alpha_level = 0.05) {
  library(phyloseq)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  
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
  
  # Function to calculate effect size (eta squared)
  kruskal_eta_sq <- function(formula, data) {
    test <- kruskal.test(formula, data = data)
    H <- test$statistic
    k <- length(unique(data[[all.vars(formula)[2]]]))
    n <- nrow(data)
    eta_sq <- (H - k + 1) / (n - k)
    list(statistic = H, p = test$p.value, eta_sq = eta_sq)
  }
  
  # Run significance tests
  chao_res <- kruskal_eta_sq(as.formula(paste("Chao1 ~", variable)), sam_data_final)
  shan_res <- kruskal_eta_sq(as.formula(paste("Shannon ~", variable)), sam_data_final)
  
  if (chao_res$p < alpha_level | shan_res$p < alpha_level) {
    message(sprintf("Significant results for '%s' (Chao1 p=%.4f, η²=%.3f; Shannon p=%.4f, η²=%.3f)", 
                    variable, chao_res$p, chao_res$eta_sq, shan_res$p, shan_res$eta_sq))
    
    # Reshape data into long format
    long_data <- sam_data_final %>%
      dplyr::select(all_of(variable), Chao1, Shannon) %>%
      pivot_longer(cols = c("Chao1", "Shannon"), names_to = "Measure", values_to = "Value")
    
    # Add p-values and eta squared for each measure
    annotations <- data.frame(
      Measure = c("Chao1", "Shannon"),
      label = c(
        sprintf("p=%.4f, η²=%.3f", chao_res$p, chao_res$eta_sq),
        sprintf("p=%.4f, η²=%.3f", shan_res$p, shan_res$eta_sq)
      ),
      y = c(max(long_data$Value[long_data$Measure=="Chao1"]) * 1.05,
            max(long_data$Value[long_data$Measure=="Shannon"]) * 1.05)
    )
    
    # Plot
    p <- ggplot(long_data, aes_string(x = variable, y = "Value", color = variable)) +
      geom_boxplot() +
      geom_jitter(alpha = 0.25, width = 0.1) +
      facet_wrap(~Measure, scales = "free_y") +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18, face = "bold"),
            legend.position = "none",
            axis.text.x = element_text(size = 14)) +
      labs(x = "Weight (kg)") 
      # labs(x = gsub("_", " ", variable))
      # scale_x_discrete(breaks = c("never", "sometimes", "always"), labels = c("never", "sometimes", "always"))
    
    # Automatically center based on numeric positions of factor levels
    annotations$x <- sapply(annotations$Measure, function(m) {
      groups <- unique(long_data[[variable]][long_data$Measure == m])
      mean(as.numeric(factor(groups, levels = unique(long_data[[variable]]))))
    })
    
    # Add the text
    p <- p + geom_text(
      data = annotations,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 10
    )
    
    # Save plot if directory specified
      ggsave(filename = paste0(variable, "_diversity_plot.png"),
             plot = p, width = 12, height = 6, dpi = 300)
    
    return(p)
  } else {
    message(sprintf("No significant difference for '%s' (Chao1 p=%.4f, Shannon p=%.4f) — skipped.", 
                    variable, chao_res$p, shan_res$p))
    return(NULL)
  }
}

setwd("C:/Users/12697/Documents/Microbiome/Data Analysis/Figures/alpha diversity boxplots/Significant Plots")
for(factor in colnames(sam)){
  try(create_a_diversity_plot(ps, factor))
}


create_a_diversity_plot(ps, "Frequency_of_Using_School_Latrine")

create_a_diversity_plot(ps, "Frequency_of_Eating_Raw_or_Undercooked_Vegetables")

create_a_diversity_plot(ps, "Weight (kg)")


#Create tables for alpha diversity

sam <- ps@sam_data
alpha_div <- estimate_richness(ps, measures = c("Chao1", "Shannon"))
sam$Chao1 <- alpha_div$Chao1
sam$Shannon <- alpha_div$Shannon

sam <- data.frame(sam) |>
  filter(Frequency_of_Using_Soap_After_Using_Toilet %in% c("sometimes", "always"))

library(effectsize)

kruskal_eta_sq <- function(formula, data) {
  test <- kruskal.test(formula, data = data)
  H <- test$statistic
  k <- length(unique(data[[all.vars(formula)[2]]]))
  n <- nrow(data)
  eta_sq <- (H - k + 1) / (n - k)
  list(statistic = H, p = test$p.value, eta_sq = eta_sq)
}

latrine_eta_shannon <- kruskal_eta_sq(Shannon ~ Frequency_of_Using_Soap_After_Using_Toilet, sam)
latrine_eta_Chao1 <- kruskal_eta_sq(Chao1 ~ Frequency_of_Using_Soap_After_Using_Toilet, sam)


latrine_table <- tibble("Diversity Measure" = c("Shannon", "Chao1"), 
                        W = c(latrine_eta_shannon$statistic, latrine_eta_Chao1$statistic), 
                        "p-value" = c(latrine_eta_shannon$p, latrine_eta_Chao1$p), 
                        "Effect Size" = c(latrine_eta_shannon$eta_sq, latrine_eta_Chao1$eta_sq))

latrine_table$Interpretation <- interpret_eta_squared(latrine_table$`Effect Size`)

