#Socioeconomic Factor Analysis

#Modified Code from Christine Bi by Max Warriner

# ## Part 0: data/phyloseq preprocessing -----------------------------------



#Import the phyloseq object 
#Assume your R code, biom file, tree file, and new columns (csv) file are all in the same folder
#First you need to set working directory to source file location in the 'Session' tab above

library(BiocManager)
library(phyloseq)
library(tidyverse)
library(dplyr)

ps <- readRDS("categorized_data.RDS")
ps

#Rarefy Data
set.seed(13)
ps <- rarefy_even_depth(ps, 1e+05)

#Changing the taxonomic table labels (From Rank 1~Rank 7 to Kingdom~Species)
# Access the existing taxonomy table from the phyloseq object
tax_table_ps <- tax_table(ps)[,1:7]

# Assign the modified taxonomy table back to the phyloseq object
tax_table(ps) <- tax_table_ps

sam_data <- sam_data(ps)

convert_binary_to_yesno <- function(df, skip_non_numeric = TRUE) {
  for (col_name in names(df)) {
    col_data <- df[[col_name]]
    
    # Skip if not numeric (when skip_non_numeric = TRUE)
    if (skip_non_numeric && !is.numeric(col_data)) next
    
    # Get unique non-NA values
    unique_vals <- unique(na.omit(col_data))
    
    # Check if binary (0/1)
    if (all(unique_vals %in% c(0, 1)) && length(unique_vals) <= 2) {
      df[[col_name]] <- factor(col_data,
                               levels = c(0, 1),
                               labels = c("no", "yes"))
    }
  }
  return(df)
}

sam_data <- convert_binary_to_yesno(sam_data)

sam_data(ps) <- sam_data

# ## Part I : Alpha Diversity (Chao1 and Shannon indices) -----------------

library(ggrepel)
library(patchwork)
library(RColorBrewer)

create_a_diversity_plot <- function(ps, factor){

a_diversity_factor<- plot_richness(ps, x=factor, color=factor, measures=c("Chao1", "Shannon"))
a_diversity_factor$layers[[2]] = NULL 
a_diversity_factor <- a_diversity_factor  + geom_boxplot() + theme_bw()

print(a_diversity_factor)

ggsave(a_diversity_factor, 
       filename = paste(factor,"_boxplot.pdf", sep = ""),
       device = "pdf",
       height = 6, width = 5, units = "in")

}

create_a_diversity_plot(ps, "Household.Number")
create_a_diversity_plot(ps, "Siblings.Younger.than.12")
create_a_diversity_plot(ps, "HeardALnamebefore")
create_a_diversity_plot(ps, "HeardTTnamebefore")
create_a_diversity_plot(ps, "HeardHWnamebefore")
create_a_diversity_plot(ps, "HeardInWormnamebefore")
create_a_diversity_plot(ps, "HeardMalanamebefore")
create_a_diversity_plot(ps, "HeardTBnamebefore")
create_a_diversity_plot(ps, "HeardSChnamebefore")
create_a_diversity_plot(ps, "HPtold")
create_a_diversity_plot(ps, "Teachertold")
create_a_diversity_plot(ps, "Mediatold")
create_a_diversity_plot(ps, "Doyouknowhowintestinalwomstransmitted")
create_a_diversity_plot(ps, "Yourlivingaddress")
create_a_diversity_plot(ps, "Familyoccupation")
create_a_diversity_plot(ps, "Maternaleducationalstatus")


# ## Part II (a) Beta diversity using Bray distance ----------------------

library(MicrobiotaProcess)
# Define the function to create PCoA plots
create_pcoa_plot <- function(variable) {
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")
  
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

# Create plots for each variable
create_pcoa_plot("Household.Number")
create_pcoa_plot("Siblings.Younger.than.12")
create_pcoa_plot("HeardALnamebefore")
create_pcoa_plot("HeardTTnamebefore")
create_pcoa_plot("HeardHWnamebefore")
create_pcoa_plot("HeardInWormnamebefore")
create_pcoa_plot("HeardMalanamebefore")
create_pcoa_plot("HeardTBnamebefore")
create_pcoa_plot("HeardSChnamebefore")
create_pcoa_plot("HPtold")
create_pcoa_plot("Teachertold")
create_pcoa_plot("Mediatold")
create_pcoa_plot("Doyouknowhowintestinalwomstransmitted")
create_pcoa_plot("Yourlivingaddress")
create_pcoa_plot("Familyoccupation")
create_pcoa_plot("Maternaleducationalstatus")

##Calculations of Significance:

library(vegan)

distme <- get_dist(ps, distmethod ="bray")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
set.seed(1024)

factors <- c("Household.Number", "Siblings.Younger.than.12", "HeardALnamebefore", "HeardTTnamebefore", "HeardHWnamebefore", "HeardInWormnamebefore", "HeardMalanamebefore", "HeardTBnamebefore", "HeardSChnamebefore", "HPtold", "Teachertold", "Mediatold", "Doyouknowhowintestinalwomstransmitted", "Yourlivingaddress", "Familyoccupation", "Maternaleducationalstatus")

adores_summary <- data.frame(factor = factors, 
                             p = rep(NA, length(factors)))

run_adores <- function(factors){
  
  for (x in 1:length(factors)){
    adores_summary[x,"p"] <- adonis2(distme ~ factors[x], data = sampleda, permutations = 999)
  }
  
}
run_adores(factors)


adores <- adonis2(distme ~ Household.Number, data=sampleda, permutations=999) #p-value = 
summary <- as.data.frame(summary(adores))

adores_sex <- adonis2(distme ~ Sex, data=sampleda, permutations=999) #p-value = 
adores_sex

adores_weight <- adonis2(distme ~ weight, data=sampleda, permutations=999) #p-value = 
adores_weight

adores_height <- adonis2(distme ~ Height, data=sampleda, permutations=999) #p-value = 
adores_height

adores_delivery <- adonis2(distme ~ Modeofdelivery, data=sampleda, permutations=999) #p-value = 
adores_delivery

adores_vaccine <- adonis2(distme ~ Everhadvaccinated, data=sampleda, permutations=999) #p-value = 
adores_vaccine

adores_BCG <- adonis2(distme ~ BCGscar, data=sampleda, permutations=999) #p-value = 
adores_BCG

adores_trimmed <- adonis2(distme ~ Childsfingernailtrimmed, data=sampleda, permutations=999) #p-value = 
adores_trimmed

adores_dirtynails <- adonis2(distme ~ Arechildsfingernailsdirty, data=sampleda, permutations=999) #p-value = 
adores_dirtynails

adores_nailsoften <- adonis2(distme ~ Howoftendoyoutrimyourfingernails, data=sampleda, permutations=999) #p-value = 
adores_nailsoften

adores_antibiotics <- adonis2(distme ~ Didyourparentsorhealthprofessionalsgaveyouotherantibiotics, data=sampleda, permutations=999) #p-value = 
adores_antibiotics

adores_deworming <- adonis2(distme ~ Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill, data=sampleda, permutations=999) #p-value = 
adores_deworming


# ##(b) Beta diversity using Jaccard distance


library(MicrobiotaProcess)
# Define the function to create PCoA plots
create_pcoa_plot <- function(variable) {
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps, distmethod = "jaccard", method = "hellinger")
  
  # Create PCoA plot
  pcoaplot <- ggordpoint(obj = pcoares, biplot = FALSE, speciesannot = TRUE,
                         factorNames = c(variable), ellipse = TRUE)
  
  ggsave(pcoaplot, 
         filename = paste(variable, "_pcoa_jaccard.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(pcoaplot)
}

# Create plots for each variable
pcoa_Factor_age <- create_pcoa_plot("Age")
print(pcoa_Factor_age)

pcoa_Factor_sex <- create_pcoa_plot("Sex")
print(pcoa_Factor_sex)

pcoa_Factor_weight <- create_pcoa_plot("weight")
print(pcoa_Factor_weight)

pcoa_Factor_height <- create_pcoa_plot("Height")
print(pcoa_Factor_height)

pcoa_Factor_delivery <- create_pcoa_plot("Modeofdelivery")
print(pcoa_Factor_delivery)

pcoa_Factor_vaccine <- create_pcoa_plot("Everhadvaccinated")
print(pcoa_Factor_vaccine)

pcoa_Factor_BCG <- create_pcoa_plot("BCGscar")
print(pcoa_Factor_BCG)

pcoa_Factor_trimmed <- create_pcoa_plot("Childsfingernailtrimmed")
print(pcoa_Factor_trimmed)

pcoa_Factor_nailsdirty <- create_pcoa_plot("Arechildsfingernailsdirty")
print(pcoa_Factor_nailsdirty)

pcoa_Factor_nailsoften <- create_pcoa_plot("Howoftendoyoutrimyourfingernails")
print(pcoa_Factor_nailsoften)

pcoa_Factor_antibiotics <- create_pcoa_plot("Didyourparentsorhealthprofessionalsgaveyouotherantibiotics")
print(pcoa_Factor_antibiotics)

pcoa_Factor_deworming <- create_pcoa_plot("Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill")
print(pcoa_Factor_deworming)



##Calculations of significance:

#PERMANOVA: Jaccard's distance

library(vegan)
distme <- get_dist(ps, distmethod ="jaccard")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
set.seed(1024)

adores_age <- adonis2(distme ~ Age, data=sampleda, permutations=999) #p-value = 0.49
adores_age

adores_sex <- adonis2(distme ~ Sex, data=sampleda, permutations=999) #p-value = 0.272
adores_sex

adores_weight <- adonis2(distme ~ weight, data=sampleda, permutations=999) #p-value = 0.094 
adores_weight

adores_height <- adonis2(distme ~ Height, data=sampleda, permutations=999) #p-value = 0.466
adores_height

adores_delivery <- adonis2(distme ~ Modeofdelivery, data=sampleda, permutations=999) #p-value = 0.673
adores_delivery

adores_vaccine <- adonis2(distme ~ Everhadvaccinated, data=sampleda, permutations=999) #p-value = 0.526
adores_vaccine

adores_BCG <- adonis2(distme ~ BCGscar, data=sampleda, permutations=999) #p-value = 0.601
adores_BCG

adores_trimmed <- adonis2(distme ~ Childsfingernailtrimmed, data=sampleda, permutations=999) #p-value = 0.474
adores_trimmed

adores_dirtynails <- adonis2(distme ~ Arechildsfingernailsdirty, data=sampleda, permutations=999) #p-value = 0.653
adores_dirtynails

adores_nailsoften <- adonis2(distme ~ Howoftendoyoutrimyourfingernails, data=sampleda, permutations=999) #p-value = 0.689
adores_nailsoften

adores_antibiotics <- adonis2(distme ~ Didyourparentsorhealthprofessionalsgaveyouotherantibiotics, data=sampleda, permutations=999) #p-value = 0.879
adores_antibiotics

adores_deworming <- adonis2(distme ~ Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill, data=sampleda, permutations=999) #p-value = 0.213
adores_deworming


# ## Part III ##Overall microbial abundance "at the Phylum, Genus and Family Level --------

# Define the function to create bar plots with titles
create_phylum_barplot <- function(variable) {
  phylumtaxa <- get_taxadf(obj = ps, taxlevel = 2)
  barplot <- ggbartax(obj = phylumtaxa, facetNames = variable, plotgroup = TRUE, topn = 5) +
    xlab(NULL) +
    ylab("relative abundance (%)") +
    labs(title = variable) +  # Add variable name as the title
    scale_fill_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(31))) +
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 2))
  
  ggsave(barplot, 
         filename = paste(variable, "_phylum_barplot.pdf", sep = ""),
         device = "pdf",
         height = 6, width = 8, units = "in")
  
  return(barplot)
}

# Create plots for each variable
phylum_Factor_age <- create_phylum_barplot("Age")
print(phylum_Factor_age)

phylum_Factor_sex <- create_phylum_barplot("Sex")
print(phylum_Factor_sex)

phylum_Factor_weight <- create_phylum_barplot("weight")
print(phylum_Factor_weight)

phylum_Factor_height <- create_phylum_barplot("Height")
print(phylum_Factor_height)

phylum_Factor_delivery <- create_phylum_barplot("Modeofdelivery")
print(phylum_Factor_delivery)

phylum_Factor_vaccine <- create_phylum_barplot("Everhadvaccinated")
print(phylum_Factor_vaccine)

phylum_Factor_BCG <- create_phylum_barplot("BCGscar")
print(phylum_Factor_BCG)

phylum_Factor_trimmed <- create_phylum_barplot("Childsfingernailtrimmed")
print(phylum_Factor_trimmed)

phylum_Factor_nailsdirty <- create_phylum_barplot("Arechildsfingernailsdirty")
print(phylum_Factor_nailsdirty)

phylum_Factor_nailsoften <- create_phylum_barplot("Howoftendoyoutrimyourfingernails")
print(phylum_Factor_nailsoften)

phylum_Factor_antibiotics <- create_phylum_barplot("Didyourparentsorhealthprofessionalsgaveyouotherantibiotics")
print(phylum_Factor_antibiotics)

phylum_Factor_deworming <- create_phylum_barplot("Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill")
print(phylum_Factor_deworming)



#Microbial Abundance by Phylum Across Different Conditions

#Stacked bar plots showing the relative abundance of microbial phyla in all samples, grouped by a) FactorA and b) FactorB. Each bar represents the average relative abundance of the top five phyla, with "Others" representing the remaining phyla.

##Overall microbial abundance "at the Family level" across different conditions (in one panel)


# Define the function to create bar plots at the family level
create_family_barplot <- function(variable) {
  familytaxa <- get_taxadf(obj = ps, taxlevel = 5)
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

# Create plots for each variable 
family_Factor_age <- create_family_barplot("Age")
print(family_Factor_age)

family_Factor_sex <- create_family_barplot("Sex")
print(family_Factor_sex)

family_Factor_weight <- create_family_barplot("weight")
print(family_Factor_weight)

family_Factor_height <- create_family_barplot("Height")
print(family_Factor_height)

family_Factor_delivery <- create_family_barplot("Modeofdelivery")
print(family_Factor_delivery)

family_Factor_vaccine <- create_family_barplot("Everhadvaccinated")
print(family_Factor_vaccine)

family_Factor_BCG <- create_family_barplot("BCGscar")
print(family_Factor_BCG)

family_Factor_trimmed <- create_family_barplot("Childsfingernailtrimmed")
print(family_Factor_trimmed)

family_Factor_nailsdirty <- create_family_barplot("Arechildsfingernailsdirty")
print(family_Factor_nailsdirty)

family_Factor_nailsoften <- create_family_barplot("Howoftendoyoutrimyourfingernails")
print(family_Factor_nailsoften)

family_Factor_antibiotics <- create_family_barplot("Didyourparentsorhealthprofessionalsgaveyouotherantibiotics")
print(family_Factor_antibiotics)

family_Factor_deworming <- create_family_barplot("Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill")
print(family_Factor_deworming)




#Microbial Abundance by Family Across Different Conditions

#Stacked bar plots showing the relative abundance of microbial families in all samples, grouped by a) FactorA and b) FactorB. Each bar represents the average relative abundance of the top ten families, with "Others" representing the remaining families. 

##Overall microbial abundance "at the Genus level" across different conditions (in one panel)


# Define the function to create bar plots at the genus level
create_genus_barplot <- function(variable) {
  genustaxa <- get_taxadf(obj = ps, taxlevel = 6)
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

# Create plots for each variable 
genus_Factor_age <- create_genus_barplot("Age")
print(genus_Factor_age)

genus_Factor_sex <- create_genus_barplot("Sex")
print(genus_Factor_sex)

genus_Factor_weight <- create_genus_barplot("weight")
print(genus_Factor_weight)

genus_Factor_height <- create_genus_barplot("Height")
print(genus_Factor_height)

genus_Factor_delivery <- create_genus_barplot("Modeofdelivery")
print(genus_Factor_delivery)

genus_Factor_vaccine <- create_genus_barplot("Everhadvaccinated")
print(genus_Factor_vaccine)

genus_Factor_BCG <- create_genus_barplot("BCGscar")
print(genus_Factor_BCG)

genus_Factor_trimmed <- create_genus_barplot("Childsfingernailtrimmed")
print(genus_Factor_trimmed)

genus_Factor_nailsdirty <- create_genus_barplot("Arechildsfingernailsdirty")
print(genus_Factor_nailsdirty)

genus_Factor_nailsoften <- create_genus_barplot("Howoftendoyoutrimyourfingernails")
print(genus_Factor_nailsoften)

genus_Factor_antibiotics <- create_genus_barplot("Didyourparentsorhealthprofessionalsgaveyouotherantibiotics")
print(genus_Factor_antibiotics)

genus_Factor_deworming <- create_genus_barplot("Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill")
print(genus_Factor_deworming)



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


# Function to create volcano plots for taxonomic levels
create_volcano_plot <- function(ps_obj, taxlevel, condition_col, metadata, SVs, taxonomy) {
  
  # Ensure that taxonomy contains the specified level
  tax_col <- colnames(taxonomy)[taxlevel]
  taxonomy <- taxonomy %>%
    rownames_to_column(var = "OTU_name") %>%
    dplyr::select(OTU_name, !!sym(tax_col)) %>%
    distinct()
  
  # Combine SVs with taxonomy to match the specified level
  SVs_with_taxonomy <- as.data.frame(SVs) %>%
    rownames_to_column(var = "OTU_name") %>%
    left_join(taxonomy, by = "OTU_name")
  
  # Aggregate the SVs data to the specified taxonomic level
  tax_level_SVs <- SVs_with_taxonomy %>%
    dplyr::select(-OTU_name) %>%
    group_by(!!sym(tax_col)) %>%
    summarise(across(everything(), ~sum(.x, na.rm = TRUE)))%>%
    filter(!!sym(tax_col) != "none") %>%
    column_to_rownames(var = tax_col)
  
  # Convert to matrix for ALDEx2 analysis
  tax_level_SVs <- as.matrix(tax_level_SVs)
  
  condition <- as.character(as.factor(metadata[[condition_col]]))
  
  # Run ALDEx2 analysis
  aldex_data <- aldex(tax_level_SVs, conditions = condition, mc.samples = 128, test = "t", effect = TRUE)
  
  # Combine t-test and effect size results
  results <- data.frame(aldex_data)
  
  # Merge results with taxonomy
  results <- results %>%
    rownames_to_column(var = tax_col)
  
  # Create the volcano plot (after multiple testing correction)
  p <- results %>%
    mutate(Significant = if_else(wi.eBH < 0.1, TRUE, FALSE)) %>%
    mutate(Taxon = as.character(!!sym(tax_col))) %>%
    mutate(TaxonToPrint = if_else(wi.eBH < 0.05, paste(Taxon, "(", round(wi.eBH,3) , ")", sep = ""), "")) |>
    ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant, label = TaxonToPrint)) +
    geom_text_repel(size = 2, nudge_y = 0.05, max.overlaps = Inf) +  # Increase max.overlaps
    geom_point(alpha = 0.6, shape = 16) +
    theme_minimal() +  
    xlab("log2(fold change)") +
    ylab("-log10(P-value)") +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "red")) +
    ggtitle(paste("Taxonomic Level:", tax_col))
  
  return(p)
}


create_combined_volcano(ps, "Sex", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Modeofdelivery", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Everhadvaccinated", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Childsfingernailtrimmed", metadata, SVs, taxonomy)
create_combined_volcano(ps, "BCGscar", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Arechildsfingernailsdirty", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Didyourparentsorhealthprofessionalsgaveyouotherantibiotics", metadata, SVs, taxonomy)
create_combined_volcano(ps, "Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill", metadata, SVs, taxonomy)



# ## PartV: Differential Abundance Using Lefse -----------------------------


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
