
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
library(Maaslin2)

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

combined_SES <- shannon_SES + chao1_SES

ggsave(plot = combined_SES, filename = 'SES_alpha_boxplot.png', width = 7, height = 4)

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

combined_hygiene <- shannon_hygiene + chao1_hygiene

ggsave(plot = combined_hygiene, filename = 'hygiene_alpha_boxplot.png', width = 7, height = 4)

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

combined_HIQ <- shannon_HIQ + chao1_HIQ 

ggsave(plot = combined_HIQ, filename = 'HIQ_alpha_boxplot.png', width = 7, height = 4)

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

combined_HSI <- shannon_HSI + chao1_HSI

ggsave(plot = combined_HSI, filename = 'HSI_alpha_boxplot.png', width = 7, height = 4)

summary(lm(Shannon ~ HSI_index, data = sam))
summary(lm(chao1 ~ HSI_index, data = sam))


# Beta Diversity ----------------------------------------------------------

pcoares <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")

physeq_dist <- phyloseq::distance(ps, method = "jaccard")

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

  ggsave(pcoaplot,
         filename = paste(variable, "_pcoa_jaccard.png", sep = ""),
         device = "png",
         height = 6, width = 8, units = "in")
  
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

#before correction
create_volcano_plot <- function(ps_obj, tax_col, condition_col, metadata, SVs, taxonomy, comparison_title = NULL) {
  
  taxonomy <- as.data.frame(phyloseq::tax_table(ps))    # Extract taxonomy
  
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

#after correction
create_volcano_plot <- function(ps_obj, tax_col, condition_col, metadata, SVs, taxonomy, comparison_title = NULL) {
  
  taxonomy <- as.data.frame(phyloseq::tax_table(ps))    # Extract taxonomy
  
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
  volcano_plot_family <- create_volcano_plot(ps, "Family", factor, metadata, SVs, taxonomy, comparison_title)
  volcano_plot_genus <- create_volcano_plot(ps, "Genus", factor, metadata, SVs, taxonomy, comparison_title)
  
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


create_combined_volcano(ps, "SES_group", sam, SVs, taxonomy)
create_combined_volcano(ps, "hygiene_group", sam, SVs, taxonomy)
create_combined_volcano(ps, "HSI_group", sam, SVs, taxonomy)
create_combined_volcano(ps, "HIQ_group", sam, SVs, taxonomy)





# Machine Learning --------------------------------------------------------

cv_predict_clr_xgb <- function(
    ps_obj,
    outcome_var,
    meta_cols,
    min_prevalence = 0.05,
    nfolds = 10,
    seed = 100,
    n_cores = 4
) {
  set.seed(1313)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  get_tax_matrix <- function(ps, taxlevel, prefix) {
    ps_tax <- tax_glom(ps, taxlevel)
    otumat <- as.data.frame(otu_table(ps_tax))
    if (!taxa_are_rows(ps_tax)) otumat <- t(otumat)
    taxmat <- as.data.frame(tax_table(ps_tax))
    tax_names <- taxmat[[taxlevel]]
    tax_names[is.na(tax_names) | tax_names == ""] <- "Unassigned"
    rownames(otumat) <- paste0(prefix, make.unique(tax_names))
    present_counts <- rowSums(otumat > 0)
    keep <- present_counts >= (min_prevalence * ncol(otumat))
    otumat <- otumat[keep, , drop = FALSE]
    return(otumat)
  }
  
  fam_tab <- get_tax_matrix(ps_obj, "Family", "F__")
  gen_tab <- get_tax_matrix(ps_obj, "Genus", "G__")
  all_feat_tab <- rbind(fam_tab, gen_tab)
  all_feat_tab[all_feat_tab == 0] <- 1e-6
  clr_mat <- compositions::clr(all_feat_tab)
  clr_mat <- t(clr_mat)
  
  meta <- as.data.frame(sample_data(ps_obj))
  clr_samples <- rownames(clr_mat)
  meta <- meta[clr_samples, , drop = FALSE]
  y <- as.factor(meta[[outcome_var]])
  keep <- !is.na(y)
  clr_mat <- clr_mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  y <- y[keep]
  levels(y) <- make.names(levels(y))
  
  meta_cols <- intersect(meta_cols, colnames(meta))
  meta_features <- meta[, meta_cols, drop = FALSE]
  meta_features[] <- lapply(meta_features, function(x) as.numeric(as.character(x)))
  meta_features <- meta_features[, sapply(meta_features, function(x) all(x %in% c(0, 1, NA))), drop = FALSE]
  
  X_full <- cbind(clr_mat, meta_features)
  nzv <- caret::nearZeroVar(X_full)
  if (length(nzv) > 0) X_full <- X_full[, -nzv, drop = FALSE]
  
  # SIMPLIFIED hyperparameter grid
  xgb_grid <- expand.grid(
    nrounds = c(100, 300, 600),
    max_depth = c(3, 5, 7, 9),
    eta = c(0.01, 0.05, 0.1),
    gamma = c(0, 0.5),
    colsample_bytree = c(0.7, 0.9),
    min_child_weight = c(1, 3),
    subsample = c(0.7, 0.9)
  )
  
  fitControl <- caret::trainControl(
    method = "cv",
    number = nfolds,
    classProbs = TRUE,
    savePredictions = "final",
    allowParallel = TRUE, 
    verboseIter = T
  )
  
  xgb_fit <- caret::train(
    x = X_full,
    y = y,
    method = "xgbTree",
    tuneGrid = xgb_grid,
    trControl = fitControl,
    verbose = FALSE
  )
  
  overall_acc <- max(xgb_fit$results$Accuracy)
  best_params <- xgb_fit$bestTune
  
  varimp <- caret::varImp(xgb_fit, scale = FALSE)$importance
  varimp_df <- tibble(Feature = rownames(varimp), Importance = varimp[,1]) %>%
    arrange(desc(Importance)) %>%
    mutate(Level = dplyr::case_when(
      grepl("^F__", Feature) ~ "Family",
      grepl("^G__", Feature) ~ "Genus",
      Feature %in% meta_cols ~ "Metadata",
      TRUE ~ "Other"
    ))
  
  preds <- xgb_fit$pred$pred[order(xgb_fit$pred$rowIndex)]
  y_true <- xgb_fit$pred$obs[order(xgb_fit$pred$rowIndex)]
  overall_cm <- caret::confusionMatrix(preds, y_true)
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(list(
    overall_accuracy = overall_acc,
    best_params = best_params,
    confusion_matrix = overall_cm,
    feature_importance = varimp_df,
    predictions = preds,
    y_true = y_true,
    xgb_fit = xgb_fit
  ))
}

## ROC, Variable Importance, and Heatmap plots

plot_roc_curve_gg <- function(model_results, positive_class = NULL, factor) {
  
  y_true <- model_results$y_true
  if (is.null(positive_class)) {
    positive_class <- levels(y_true)[1]
  }
  pred_prob <- model_results$xgb_fit$pred %>%
    filter(obs %in% levels(y_true)) %>%
    arrange(rowIndex) %>%
    pull(!!as.name(positive_class))
  
  roc_obj <- roc(y_true, pred_prob, levels = rev(levels(y_true)), direction = "<")
  auc_value <- auc(roc_obj)
  
  roc <- ggroc(roc_obj, legacy.axes = TRUE, colour = "darkgreen", size = 1.3) +
    geom_abline(linetype = "dashed", color = "gray50") +
    labs(
      title = paste(factor, ": ROC Curve (AUC = ", round(auc_value, 3), ")", sep = ""),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(filename = paste(factor, "ROC.png", sep = ""), plot = roc, width = 7, height = 5)
  
}

plot_top_importance <- function(model_results, n_top = 10, bar_color = "#2c7bb6", factor) {
  
  varimp <- model_results$feature_importance %>% 
    arrange(desc(Importance)) %>%
    head(n_top)
  
  plot <- ggplot(varimp, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = bar_color) +
    coord_flip() +
    labs(
      title = paste(factor, ": Top", n_top, "Feature Importance"),
      x = "",
      y = "Importance"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(filename = paste(factor, "top_importance.png", sep = ""), plot = plot, width = 10, height = 5)
  
}

plot_top_feature_heatmap_clr <- function(
    ps_obj,
    model_results,
    n_top = 10,
    metadata_vars,
    outcome_var,
    min_prevalence = 0.05
) {
  
  # 1. Top N features
  top_features <- head(model_results$feature_importance$Feature, n_top)
  
  # 2. Get genus & family tables
  get_tax_matrix <- function(ps, taxlevel, prefix) {
    ps_tax <- tax_glom(ps, taxlevel)
    otumat <- as.data.frame(otu_table(ps_tax))
    if (!taxa_are_rows(ps_tax)) otumat <- t(otumat)
    taxmat <- as.data.frame(tax_table(ps_tax))
    tax_names <- taxmat[[taxlevel]]
    tax_names[is.na(tax_names) | tax_names == ""] <- "Unassigned"
    rownames(otumat) <- paste0(prefix, make.unique(tax_names))
    present_counts <- rowSums(otumat > 0)
    keep <- present_counts >= (min_prevalence * ncol(otumat))
    otumat <- otumat[keep, , drop=FALSE]
    return(otumat)
  }
  
  fam_tab <- get_tax_matrix(ps_obj, "Family", "F__")
  gen_tab <- get_tax_matrix(ps_obj, "Genus", "G__")
  all_feat_tab <- rbind(fam_tab, gen_tab)
  
  # 3. CLR transformation
  all_feat_tab[all_feat_tab == 0] <- 1e-6
  clr_mat <- compositions::clr(all_feat_tab)
  clr_mat <- t(clr_mat)  # samples as rows, features as columns
  
  # 4. Get metadata (as numeric)
  meta <- as.data.frame(sample_data(ps_obj))
  clr_samples <- rownames(clr_mat)
  meta <- meta[clr_samples, , drop = FALSE]
  metadata_vars <- intersect(metadata_vars, colnames(meta))
  metadata_features <- meta[, metadata_vars, drop = FALSE]
  metadata_features[] <- lapply(metadata_features, function(x) as.numeric(as.character(x)))
  
  # 5. Combine CLR and metadata
  X_full <- cbind(clr_mat, metadata_features)
  # Ensure all columns are numeric
  for (i in seq_len(ncol(X_full))) {
    if (!is.numeric(X_full[, i])) {
      X_full[, i] <- as.numeric(as.character(X_full[, i]))
    }
  }
  
  # 6. Subset to top features, as rows (order preserved)
  top_feats_present <- top_features[top_features %in% colnames(X_full)]
  if (length(top_feats_present) == 0) stop("No top features found in input matrix.")
  heatmap_mat <- t(X_full[, top_feats_present, drop = FALSE])
  
  # Remove rows (features) that are all NA or all zero
  keep_rows <- apply(heatmap_mat, 1, function(x) any(!is.na(x)) && any(x != 0))
  heatmap_mat <- heatmap_mat[keep_rows, , drop = FALSE]
  if (nrow(heatmap_mat) == 0) stop("No valid features to plot after removing all-NA/zero rows.")
  
  # 7. Get annotation for sample outcome
  sample_anno <- as(sample_data(ps_obj), "data.frame")[, outcome_var, drop = FALSE]
  sample_anno <- sample_anno[colnames(heatmap_mat), , drop = FALSE]
  sample_anno[[outcome_var]] <- as.factor(as.character(sample_anno[[outcome_var]]))
  rownames(sample_anno) <- colnames(heatmap_mat)
  
  # 8. Order columns by group
  group_order <- order(sample_anno[[outcome_var]])
  heatmap_mat <- heatmap_mat[, group_order, drop = FALSE]
  sample_anno <- sample_anno[group_order, , drop = FALSE]
  
  # 9. Z-score (standardize) each feature (row)
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
  heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0
  
  # 10. Plot heatmap (no clustering of samples)
  heatmap <- pheatmap::pheatmap(
    mat = heatmap_mat_scaled,
    annotation_col = sample_anno,
    main = paste("Top", n_top, "Features Heatmap"),
    clustering_method = "complete",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    fontsize_row = 10,
    fontsize_col = 7,
    scale = "none",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  
  ggsave(filename = paste(outcome_var, "_heatmap.png", sep = ""), plot = heatmap, height = 6, width = 18)
  
  return(heatmap)
  
}

# Hygiene
hygiene_results <- cv_predict_clr_xgb(ps, "hygiene_group", meta_cols = c("Age", "Sex"))
hygiene_results$confusion_matrix
head(hygiene_results$feature_importance, 50)

plot_top_feature_heatmap_clr(ps_obj = ps, model_results = hygiene_results,
                             n_top = 10, metadata_vars = c("Sex", "Age"), 
                             outcome_var = "hygiene_group", min_prevalence = 0.05)

plot_roc_curve_gg(hygiene_results, factor = "Hygiene_Group")

plot_top_importance(hygiene_results, n_top = 10, factor = "hygiene_group")

# SES
SES_results <- cv_predict_clr_xgb(ps, "SES_group", meta_cols = c("Age", "Sex"))
SES_results$confusion_matrix
head(SES_results$feature_importance, 50)

plot_top_feature_heatmap_clr(ps_obj = ps, model_results = SES_results,
                             n_top = 10, metadata_vars = c("Sex", "Age"), 
                             outcome_var = "SES_group", min_prevalence = 0.05)

plot_roc_curve_gg(SES_results, factor = "SES_Group")

plot_top_importance(SES_results, n_top = 10, factor = "SES_group")

# HSI
HSI_results <- cv_predict_clr_xgb(ps, "HSI_group", meta_cols = c("Age", "Sex"))
HSI_results$confusion_matrix
head(HSI_results$feature_importance, 50)

plot_top_feature_heatmap_clr(ps_obj = ps, model_results = HSI_results,
                             n_top = 10, metadata_vars = c("Sex", "Age"), 
                             outcome_var = "HSI_group", min_prevalence = 0.05)

plot_roc_curve_gg(HSI_results, factor = "HSI_Group")

plot_top_importance(HSI_results, n_top = 10, factor = "HSI_group")



# HIQ
HIQ_results <- cv_predict_clr_xgb(ps, "HIQ_group", meta_cols = c("Age", "Sex"))
HIQ_results$confusion_matrix
head(HIQ_results$feature_importance, 50)

plot_top_feature_heatmap_clr(ps_obj = ps, model_results = HIQ_results,
                             n_top = 10, metadata_vars = c("Sex", "Age"), 
                             outcome_var = "HIQ_group", min_prevalence = 0.05)

plot_roc_curve_gg(HIQ_results, factor = "HIQ_Group")

plot_top_importance(HIQ_results, n_top = 10, factor = "HIQ_results")



# Maaslin2 ----------------------------------------------------------------


# Extract sample data, OTU table, taxonomy
metadata <- as.data.frame(sample_data(ps))
SVs <- as.data.frame(otu_table(ps))


run_maaslin2_and_plot <- function(
    ps_obj, 
    variable, 
    title_name = variable, 
    taxlevel = "Genus",
    qval_sig_max = 0.1,
    min_prevalence = 0.05,
    min_samples = 5,               # Minimum number of samples with abundance > 0
    min_abundance_mean = 0.0005,    # Minimum mean relative abundance across all samples
    min_abundance_max  = 0.001,     # Minimum *max* relative abundance in any sample
    adjust_vars = c("Age", "sex"),
    remove_unassigned = TRUE
) {
  # Prepare feature table (taxa as rows, samples as columns)
  otumat <- as.data.frame(otu_table(ps_obj))
  if (!taxa_are_rows(ps_obj)) otumat <- data.frame(t(otumat))
  
  # Aggregate at desired taxonomic level
  taxmat <- as.data.frame(tax_table(ps_obj))
  if (!(taxlevel %in% colnames(taxmat))) stop("taxlevel not in taxonomy table!")
  
  # Assign "Unassigned" to missing genus
  otumat$Taxon <- taxmat[rownames(otumat), taxlevel]
  otumat$Taxon[is.na(otumat$Taxon) | otumat$Taxon == "" ] <- "Unassigned"
  
  # Aggregate counts at genus level
  feature_table <- otumat %>%
    group_by(Taxon) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    column_to_rownames("Taxon")
  
  # Convert to relative abundance (to apply abundance filters)
  feature_table_rel <- sweep(feature_table, 2, colSums(feature_table), "/")
  
  # Filtering step: prevalence, mean abundance, and max abundance
  present_counts <- rowSums(feature_table_rel > 0)
  mean_abund <- rowMeans(feature_table_rel)
  max_abund <- apply(feature_table_rel, 1, max)
  
  keep_taxa <- (present_counts > (min_prevalence * ncol(feature_table_rel))) &
    (present_counts >= min_samples) &
    (mean_abund >= min_abundance_mean | max_abund >= min_abundance_max)
  
  feature_table <- feature_table[keep_taxa, ]
  feature_table_rel <- feature_table_rel[keep_taxa, ]
  
  # Prepare metadata
  meta <- as(sample_data(ps_obj), "data.frame")
  meta <- meta[ , , drop = FALSE]
  
  # Match samples between metadata and feature table
  colnames(feature_table) <- gsub("X", "", colnames(feature_table))
  colnames(feature_table) <- gsub("\\.", "-", colnames(feature_table))
  common_samples <- intersect(colnames(feature_table), rownames(meta))
  feature_table <- feature_table[ , common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  
  # Compose fixed effects
  fixed_effects <- unique(c(variable, adjust_vars))
  fixed_effects <- fixed_effects[fixed_effects %in% colnames(meta)]
  
  # Run MaAsLin2 
  fit_data <- Maaslin2(
    input_data = as.data.frame(t(feature_table)),
    input_metadata = meta,
    output = tempdir(),   
    fixed_effects = fixed_effects,
    normalization = "TSS",
    min_prevalence = min_prevalence
  )
  
  # Results
  results <- fit_data$results
  results_var <- results[results$metadata == variable, ]
  results_var <- results_var %>%
    dplyr::mutate(sig = ifelse(qval < qval_sig_max, "Significant", "NS"))
  
  # Only keep one result per genus
  results_var_unique <- results_var %>%
    group_by(feature) %>%
    slice_min(order_by = qval, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Prevalence and abundance calculation (add to output)
  prevalence_df <- data.frame(
    feature = rownames(feature_table_rel),
    prevalence = rowSums(feature_table_rel > 0),
    prevalence_frac = rowSums(feature_table_rel > 0) / ncol(feature_table_rel),
    mean_abund = rowMeans(feature_table_rel),
    max_abund = apply(feature_table_rel, 1, max)
  )
  results_var_unique <- left_join(results_var_unique, prevalence_df, by = "feature")
  
  # Keep Only Genus Names That are in the Original Genus Column
  genus_list <- unique(taxmat[[taxlevel]])
  genus_list <- genus_list[!is.na(genus_list) & genus_list != ""]
  if (remove_unassigned) {
    genus_list <- genus_list[genus_list != "Unassigned"]
  }
  results_var_unique <- results_var_unique %>% filter(feature %in% genus_list)
  
  # Plot: Label rare significant genera
  results_var_unique <- results_var_unique %>%
    mutate(
      label_flag = ifelse(
        sig == "Significant" & prevalence < min_samples, 
        paste0(feature, " (rare)"), 
        ifelse(sig == "Significant", feature, "")
      )
    )
  
  p <- ggplot(results_var_unique, aes(x = coef, y = -log10(qval), color = sig)) +
    geom_point() +
    geom_text_repel(
      aes(label = label_flag),
      size = 3, 
      max.overlaps = 12, 
      color = "red",
      force = 2,
      box.padding = 0.5
    ) +
    scale_color_manual(values = c("Significant" = "red", "NS" = "black")) +
    labs(
      x = "Effect Size Coefficient",
      y = "-log10(FDR)",
      title = title_name
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(filename = paste(variable, "maslin_volcano.png", sep = ""), plot = p, width = 8, height = 6)
  
  print(p)
  n_rare_sig <- sum(results_var_unique$sig == "Significant" & results_var_unique$prevalence < min_samples)
  if (n_rare_sig > 0) {
    message(sprintf("Warning: %d significant genera are present in fewer than %d samples (labeled as (rare)). Interpret with caution!", n_rare_sig, min_samples))
  }
  return(list(plot = p, table = results_var_unique, all_results = results))
}

plot_all_significant_boxplots <- function(
    ps_obj, 
    maaslin2_table, 
    variable = "FoodSecure_vs_FoodInsecure", 
    variable_name = variable,
    taxlevel = "Genus"
) {
  
  sig_genera <- maaslin2_table %>%
    filter(sig == "Significant", metadata == variable) %>%
    pull(feature) %>%
    unique()
  
  # Helper function for a single genus
  plot_single <- function(genus) {
    qval <- maaslin2_table %>%
      filter(feature == genus, metadata == variable) %>%
      arrange(qval) %>%
      slice(1) %>%
      pull(qval)
    qval_label <- if(length(qval) == 0) "NA" else formatC(qval, digits = 3, format = "f")
    auto_title <- paste0(genus, " (q = ", qval_label, ")")
    
    otumat <- as.data.frame(otu_table(ps_obj))
    if (!taxa_are_rows(ps_obj)) otumat <- t(otumat)
    taxmat <- as.data.frame(tax_table(ps_obj))
    otumat$Taxon <- taxmat[rownames(otumat), taxlevel]
    feature_table <- otumat %>%
      group_by(Taxon) %>%
      summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
      column_to_rownames("Taxon")
    feature_table_rel <- sweep(feature_table, 2, colSums(feature_table), "/")
    if (!genus %in% rownames(feature_table_rel)) return(NULL)
    genus_abund <- as.numeric(feature_table_rel[genus, ])
    sample_names <- colnames(feature_table_rel)
    meta <- as(sample_data(ps_obj), "data.frame")
    meta <- meta[sample_names, , drop = FALSE]
    df <- data.frame(
      Sample = sample_names,
      Abundance = genus_abund,
      Group = as.factor(meta[[variable]])
    )
    ggplot(df, aes(x = Group, y = Abundance, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
      labs(
        title = auto_title,
        x = variable_name,
        y = "Relative Abundance"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # Loop and collect
  plot_list <- lapply(sig_genera, plot_single)
  names(plot_list) <- sig_genera
  plot_list <- plot_list[!sapply(plot_list, is.null)]  # Remove any nulls
  print(paste("Created", length(plot_list), "boxplots for significant genera."))
  plot_list
  
  if(plot_list > 0){
    ggsave(filename = paste(variable, "maaslin_boxplots.png", sep = ""), plot = plot_list, width = 12, height = 5)
  }
  
}

#hygiene
hygiene_maaslin2 <- run_maaslin2_and_plot(ps, "hygiene_group", qval_sig_max = 0.2)
hygiene_boxplots <- plot_all_significant_boxplots(
  ps_obj = ps,
  maaslin2_table = hygiene_maaslin2$table,
  variable = "hygiene_group",
  variable_name = "hygiene_group"
)
print(hygiene_boxplots[[1]])

#SES
SES_maaslin2 <- run_maaslin2_and_plot(ps, "SES_group", qval_sig_max = 0.2)
SES_boxplots <- plot_all_significant_boxplots(
  ps_obj = ps,
  maaslin2_table = SES_maaslin2$table,
  variable = "SES_group",
  variable_name = "SES_group"
)
print(SES_boxplots[[1]])

#HSI
HSI_maaslin2 <- run_maaslin2_and_plot(ps, "HSI_group", qval_sig_max = 0.2)
HSI_boxplots <- plot_all_significant_boxplots(
  ps_obj = ps,
  maaslin2_table = HSI_maaslin2$table,
  variable = "HSI_group",
  variable_name = "HSI_group"
)
print(HSI_boxplots[[1]])


#HIQ
HIQ_maaslin2 <- run_maaslin2_and_plot(ps, "HIQ_group", qval_sig_max = 0.2)
HIQ_boxplots <- plot_all_significant_boxplots(
  ps_obj = ps,
  maaslin2_table = HIQ_maaslin2$table,
  variable = "HIQ_group",
  variable_name = "HIQ_group"
)
print(HIQ_boxplots[[1]])

