---
title: "16s microbiome analysis pipeline"
output: html_document
date: "2025-02-06"
---
#If there's any questions, feel free to ask me in lab or email me: ybi@colgate.edu! --Christine

##Part 0: data/phyloseq preprocessing

#Import the phyloseq object 
#Assume your R code, biom file, tree file, and new columns (csv) file are all in the same folder
#First you need to set working directory to source file location in the 'Session' tab above

library(BiocManager)
library(phyloseq)
library(tidyverse)
ps <- readRDS("categorized_data.RDS")
ps

#Changing the taxonomic table labels (From Rank 1~Rank 7 to Kingdom~Species)
# Access the existing taxonomy table from the phyloseq object
tax_table_ps <- tax_table(ps)

# Define the new column names for the taxonomy levels
new_tax_labels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Assign the modified taxonomy table back to the phyloseq object
tax_table(ps) <- tax_table_ps
##Microbiome analysis:

## Part I
## Alpha Diversity (Chao1 and Shannon indices)

library(ggrepel)
library(patchwork)
library(RColorBrewer)

#Age
a_diversity_factor_age<- plot_richness(ps, x="Age", color="Age", measures=c("Chao1", "Shannon"))
a_diversity_factor_age$layers[[2]] = NULL 
a_diversity_factor_age <- a_diversity_factor_age  + geom_boxplot() + theme_bw()
print(a_diversity_factor_age)
ggsave(a_diversity_factor_age, 
       filename = "age_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#Sex
a_diversity_factor_sex <- plot_richness(ps, x="Sex", color="Sex", measures=c("Chao1", "Shannon"))
a_diversity_factor_sex$layers[[2]] = NULL 
a_diversity_factor_sex <- a_diversity_factor_sex  + geom_boxplot() + theme_bw()
print(a_diversity_factor_sex)
ggsave(a_diversity_factor_sex, 
       filename = "sex_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#Weight
a_diversity_factor_weight <- plot_richness(ps, x="weight", color="weight", measures=c("Chao1", "Shannon"))
a_diversity_factor_weight$layers[[2]] = NULL 
a_diversity_factor_weight <- a_diversity_factor_weight  + geom_boxplot() + theme_bw()
print(a_diversity_factor_weight)
ggsave(a_diversity_factor_sex, 
       filename = "weight_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#Height
a_diversity_factor_height <- plot_richness(ps, x="Height", color="Height", measures=c("Chao1", "Shannon"))
a_diversity_factor_height$layers[[2]] = NULL 
a_diversity_factor_height <- a_diversity_factor_height  + geom_boxplot() + theme_bw()
print(a_diversity_factor_height)
ggsave(a_diversity_factor_height, 
       filename = "height_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#Mode of delivery
a_diversity_factor_Modeofdelivery <- plot_richness(ps, x="Modeofdelivery", color="Modeofdelivery", measures=c("Chao1", "Shannon"))
a_diversity_factor_Modeofdelivery$layers[[2]] = NULL 
a_diversity_factor_Modeofdelivery <- a_diversity_factor_Modeofdelivery  + geom_boxplot() + theme_bw()
print(a_diversity_factor_Modeofdelivery)
ggsave(a_diversity_factor_Modeofdelivery, 
       filename = "delivery_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#Vaccination
a_diversity_factor_vaccine <- plot_richness(ps, x="Everhadvaccinated", color="Everhadvaccinated", measures=c("Chao1", "Shannon"))
a_diversity_factor_vaccine$layers[[2]] = NULL 
a_diversity_factor_vaccine <- a_diversity_factor_vaccine  + geom_boxplot() + theme_bw()
print(a_diversity_factor_vaccine)
ggsave(a_diversity_factor_vaccine, 
       filename = "vaccine_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#BCG scar
a_diversity_factor_BCG <- plot_richness(ps, x="BCGscar", color="BCGscar", measures=c("Chao1", "Shannon"))
a_diversity_factor_BCG$layers[[2]] = NULL 
a_diversity_factor_BCG <- a_diversity_factor_BCG  + geom_boxplot() + theme_bw()
print(a_diversity_factor_BCG)
ggsave(a_diversity_factor_BCG, 
       filename = "BCG_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#fingernails trimmed
a_diversity_factor_trimmed <- plot_richness(ps, x="Childsfingernailtrimmed", color="Childsfingernailtrimmed", measures=c("Chao1", "Shannon"))
a_diversity_factor_trimmed$layers[[2]] = NULL 
a_diversity_factor_trimmed <- a_diversity_factor_trimmed  + geom_boxplot() + theme_bw()
print(a_diversity_factor_trimmed)

#fingernails dirty
a_diversity_factor_dirtynails <- plot_richness(ps, x="Arechildsfingernailsdirty", color="Arechildsfingernailsdirty", measures=c("Chao1", "Shannon"))
a_diversity_factor_dirtynails$layers[[2]] = NULL 
a_diversity_factor_dirtynails <- a_diversity_factor_dirtynails  + geom_boxplot() + theme_bw()
print(a_diversity_factor_dirtynails)
ggsave(a_diversity_factor_dirtynails, 
       filename = "dirtynails_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#How often do you trim your fingernails
a_diversity_factor_nailsoften <- plot_richness(ps, x="Howoftendoyoutrimyourfingernails", color="Howoftendoyoutrimyourfingernails", measures=c("Chao1", "Shannon"))
a_diversity_factor_nailsoften$layers[[2]] = NULL 
a_diversity_factor_nailsoften <- a_diversity_factor_nailsoften  + geom_boxplot() + theme_bw()
print(a_diversity_factor_nailsoften)
ggsave(a_diversity_factor_nailsoften, 
       filename = "nailsoften_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#antibiotic
a_diversity_factor_antibiotic <- plot_richness(ps, x="Didyourparentsorhealthprofessionalsgaveyouotherantibiotics", color="Didyourparentsorhealthprofessionalsgaveyouotherantibiotics", measures=c("Chao1", "Shannon"))
a_diversity_factor_antibiotic$layers[[2]] = NULL 
a_diversity_factor_antibiotic <- a_diversity_factor_antibiotic  + geom_boxplot() + theme_bw()
print(a_diversity_factor_antibiotic)
ggsave(a_diversity_factor_antibiotic, 
       filename = "antibiotic_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

#deworming
a_diversity_factor_deworming <- plot_richness(ps, x="Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill", color="Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill", measures=c("Chao1", "Shannon"))
a_diversity_factor_deworming$layers[[2]] = NULL 
a_diversity_factor_deworming <- a_diversity_factor_deworming  + geom_boxplot() + theme_bw()
print(a_diversity_factor_deworming)
ggsave(a_diversity_factor_deworming, 
       filename = "deworming_boxplot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")


##Part II
##(a) Beta diversity using Bray distance

library(MicrobiotaProcess)
# Define the function to create PCoA plots
create_pcoa_plot <- function(variable) {
  # Perform PCoA
  pcoares <- get_pcoa(obj = ps, distmethod = "bray", method = "hellinger")
  
  # Create PCoA plot
  pcoaplot <- ggordpoint(obj = pcoares, biplot = FALSE, speciesannot = TRUE,
                         factorNames = c(variable), ellipse = TRUE, linesize = 1.5)
  
  ggsave(pcoaplot, 
         filename = paste(variable, "_pcoa_bray.pdf", sep = ""),
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

##Calculations of Significance:

library(vegan)

distme <- get_dist(ps, distmethod ="bray")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
set.seed(1024)

adores_age <- adonis2(distme ~ Age, data=sampleda, permutations=999) #p-value = 0.655
adores_age

adores_sex <- adonis2(distme ~ Sex, data=sampleda, permutations=999) #p-value = 0.34
adores_sex

adores_weight <- adonis2(distme ~ weight, data=sampleda, permutations=999) #p-value = 0.102
adores_weight

adores_height <- adonis2(distme ~ Height, data=sampleda, permutations=999) #p-value = 0.494
adores_height

adores_delivery <- adonis2(distme ~ Modeofdelivery, data=sampleda, permutations=999) #p-value = 0.683
adores_delivery

adores_vaccine <- adonis2(distme ~ Everhadvaccinated, data=sampleda, permutations=999) #p-value = 0.567
adores_vaccine

adores_BCG <- adonis2(distme ~ BCGscar, data=sampleda, permutations=999) #p-value = 0.546
adores_BCG

adores_trimmed <- adonis2(distme ~ Childsfingernailtrimmed, data=sampleda, permutations=999) #p-value = 0.503
adores_trimmed

adores_dirtynails <- adonis2(distme ~ Arechildsfingernailsdirty, data=sampleda, permutations=999) #p-value = 0.641
adores_dirtynails

adores_nailsoften <- adonis2(distme ~ Howoftendoyoutrimyourfingernails, data=sampleda, permutations=999) #p-value = 0.708
adores_nailsoften

adores_antibiotics <- adonis2(distme ~ Didyourparentsorhealthprofessionalsgaveyouotherantibiotics, data=sampleda, permutations=999) #p-value = 0.842
adores_antibiotics

adores_deworming <- adonis2(distme ~ Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill, data=sampleda, permutations=999) #p-value = 0.214
adores_deworming


##(b) Beta diversity using Jaccard distance

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





##Part III
##Overall microbial abundance "at the Phylum level" across different conditions (in one panel)


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



#Microbial Abundance by Genus Across Different Conditions

#Stacked bar plots showing the relative abundance of microbial genera in all samples, grouped by a) FactorA and b) FactorB. Each bar represents the average relative abundance of the top fifteen genera, with "Others" representing the remaining genera. 

##Part IV: Differential Abundance Analysis using ALDEx2 and Volcano plots

##FactorA, BEFORE multiple testing correction
#Volcano plots of differentially abundant taxa according to ALEDX2 results between "FactorA" positive and "FactorA" negative at the phylum, family, and genus levels (BEFORE multiple testing correction)


library(ALDEx2)
# Function to create volcano plots for taxonomic levels
create_volcano_plot <- function(ps_obj, taxlevel, condition_col) {
  metadata <- as.data.frame(sample_data(ps_obj))  # Extract metadata
  SVs <- as.data.frame(otu_table(ps_obj))         # Extract OTU table
  taxonomy <- as.data.frame(phyloseq::tax_table(ps_obj))    # Extract taxonomy
  
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
  
  # Create the volcano plot (before multiple testing correction)
  p <- results %>%
    mutate(Significant = if_else(we.ep < 0.05, TRUE, FALSE)) %>%
    mutate(Taxon = as.character(!!sym(tax_col))) %>%
    mutate(TaxonToPrint = if_else(we.ep < 0.05, Taxon, "")) %>%
    ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant, label = TaxonToPrint)) +
    geom_text_repel(size = 2, nudge_y = 0.05, max.overlaps = Inf) +  # Increase max.overlaps
    geom_point(alpha = 0.6, shape = 16) +
    theme_minimal() +  
    xlab("log2(fold change)") +
    ylab("-log10(P-value)") +
    theme(legend.position = "none") +
    ggtitle(paste("Taxonomic Level:", tax_col))
  
  return(p)
}

# Create volcano plots for each taxonomic level
volcano_plot_phylum_age <- create_volcano_plot(ps, 2, "Age")
volcano_plot_family_age <- create_volcano_plot(ps, 5, "Age")
volcano_plot_genus_age <- create_volcano_plot(ps, 6, "Age")
combined_volcano_plot_age <- (volcano_plot_phylum_age | volcano_plot_family_age) / (volcano_plot_genus_age)
# Display the combined plot
print(combined_volcano_plot_age)
ggsave(combined_volcano_plot_age, 
       filename = "age_volcano.pdf",
       device = "pdf",
       height = 6, width = 8, units = "in")



Differentially Abundant Taxa Between FactorA Positive and Negative Groups (Before multiple testing correction)

Volcano plots displaying differentially abundant taxa according to ALDEx2 results at various taxonomic levels: A) Phylum, B) Family, and C) Genus. Taxa with P-values (we.ep) less than 0.05, before multiple testing correction, are highlighted in red. 

##FactorB, BEFORE multiple testing correction
#Volcano plots of differentially abundant taxa according to ALEDX2 results between "FactorB" positive and "FactorB" negative at the phylum, family, and genus levels (BEFORE multiple testing correction)

```{r}
#After running the chunk above, run this.
# Create volcano plots for each taxonomic level
volcano_plot_phylum <- create_volcano_plot(ps, 2, "FactorB")
volcano_plot_family <- create_volcano_plot(ps, 5, "FactorB")
volcano_plot_genus <- create_volcano_plot(ps, 6, "FactorB")

# Combine all plots into a panel
combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus)

# Display the combined plot
print(combined_volcano_plot)
```

Differentially Abundant Taxa Between FactorB Positive and Negative Groups (Before multiple testing correction)

Volcano plots displaying differentially abundant taxa according to ALDEx2 results at various taxonomic levels: A) Phylum, B) Family, and C) Genus. Taxa with P-values (we.ep) less than 0.05, before multiple testing correction, are highlighted in red.


##FactorA, AFTER multiple testing correction
```{r}
# Function to create volcano plots for taxonomic levels
create_volcano_plot <- function(ps_obj, taxlevel, condition_col) {
  metadata <- as.data.frame(sample_data(ps_obj))  # Extract metadata
  SVs <- as.data.frame(otu_table(ps_obj))         # Extract OTU table
  taxonomy <- as.data.frame(phyloseq::tax_table(ps_obj))    # Extract taxonomy
  
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
    mutate(TaxonToPrint = if_else(wi.eBH < 0.1, Taxon, "")) %>%
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

# Create volcano plots for each taxonomic level
volcano_plot_phylum <- create_volcano_plot(ps, 2, "FactorA")
volcano_plot_family <- create_volcano_plot(ps, 5, "FactorA")
volcano_plot_genus <- create_volcano_plot(ps, 6, "FactorA")

# Combine all plots into a panel
combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus)

# Display the combined plot
print(combined_volcano_plot)
```

Differentially Abundant Taxa Between FactorA Positive and Negative Groups (After multiple testing correction)

Volcano plots displaying differentially abundant taxa according to ALDEx2 results at various taxonomic levels: A) Phylum, B) Family, and C) Genus. Taxa with P-values (we.eBH) less than 0.1, after BH-correction, are highlighted in red. 

##FactorB, AFTER multiple testing correction

```{r}
#After running the chunk above, run this.
# Create volcano plots for each taxonomic level
volcano_plot_phylum <- create_volcano_plot(ps, 2, "FactorB")
volcano_plot_family <- create_volcano_plot(ps, 5, "FactorB")
volcano_plot_genus <- create_volcano_plot(ps, 6, "FactorB")

# Combine all plots into a panel
combined_volcano_plot <- (volcano_plot_phylum | volcano_plot_family) / (volcano_plot_genus)

# Display the combined plot
print(combined_volcano_plot)
```


Differentially Abundant Taxa Between FactorB Positive and Negative Groups (After multiple testing correction)

Volcano plots displaying differentially abundant taxa according to ALDEx2 results at various taxonomic levels: A) Phylum, B) Family, and C) Genus. Taxa with P-values (we.eBH) less than 0.1, after BH-correction, are highlighted in red. 


##PartV: Differential Abundance Using Lefse

#FactorA
```{r}
library(microbiomeMarker)

#Adjust cutoff values if needed
#Note: If it gives error like "Error in names(marker) <- `*vtmp*` : attempt to set an attribute on NULL", it means no marker was identified, and you will have to increase your wilcoxon/kw/pvalue cutoffs and/or decrease your lda_cutoff.

# Run LefSe
mm_lefse <- run_lefse(
    ps,
    wilcoxon_cutoff = 0.05,
    group = "FactorA",
    kw_cutoff = 0.05,
    multigrp_strat = TRUE,
    lda_cutoff = 2
)

# Run edgeR
mm_edger <- run_edger(
    ps,
    group = "FactorA",
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
)

# Create plots
lefse_bar_plot <- plot_ef_bar(mm_lefse)
lefse_dot_plot <- plot_ef_dot(mm_lefse)
edger_bar_plot <- plot_ef_bar(mm_edger)
#lefse_cladogram <- (plot_cladogram(mm_lefse, color = c(Negative = "darkgreen", Positive = "red")) +
    #theme(plot.margin = margin(0, 0, 0, 0))

lefse_bar_plot
lefse_dot_plot
edger_bar_plot
#lefse_cladogram 
#This cladogram gives a lot of warnings and I don't know how to suppress them, so I commented it out. The information should be repetitive as the other three but just presented in a different format. 
```

#FactorB
```{r}
library(microbiomeMarker)

#Adjust cutoff values if needed
#Note: If it gives error like "Error in names(marker) <- `*vtmp*` : attempt to set an attribute on NULL", it means no marker was identified, and you will have to increase your wilcoxon/kw/pvalue cutoffs and/or decrease your lda_cutoff.

# Run LefSe
mm_lefse <- run_lefse(
    ps,
    wilcoxon_cutoff = 0.05,
    group = "FactorB",
    kw_cutoff = 0.05,
    multigrp_strat = TRUE,
    lda_cutoff = 2
)

# Run edgeR
mm_edger <- run_edger(
    ps,
    group = "FactorB",
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
)

# Create plots
lefse_bar_plot <- plot_ef_bar(mm_lefse)
lefse_dot_plot <- plot_ef_dot(mm_lefse)
edger_bar_plot <- plot_ef_bar(mm_edger)

lefse_bar_plot
lefse_dot_plot
edger_bar_plot
```

#Additional code: Correlation Heatmap (if wanted) (takes a long time to run, also not related to FactorA/FactorB, but looks pretty)
#Correlation between abundant taxa and MUAC,WEIGHTkg,HEIGHTcm (numerical variables)
```{r}
library(microViz)
library(pheatmap)

#First convert the variables to numeric and reassign to ps
metadata <- sample_data(ps)
metadata$MUAC <- as.numeric(metadata$MUAC)
metadata$WEIGHTkg <- as.numeric(metadata$WEIGHTkg)
metadata$HEIGHTcm <- as.numeric(metadata$HEIGHTcm)
sample_data(ps) <- metadata

#Tax fix
ps <- ps %>% tax_fix(unknowns = c("none"))

# Building correlations for phylum
top_phylum <- ps %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  group_by(Phylum) %>%
  summarise(Abundance = mean(Abundance)) %>%
  top_n(25, Abundance) %>%
  pull(Phylum)
ps_top_phylum <- subset_taxa(ps, Phylum %in% top_phylum)
correlations_df_phylum <- ps_top_phylum %>% 
  tax_model(
    trans = "clr",
    rank = "Phylum", variables = list("MUAC","WEIGHTkg", "HEIGHTcm"), 
    type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE
  ) %>% 
  tax_models2stats(rank = "Phylum")

# Building correlations for family
top_family <- ps %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  group_by(Family) %>%
  summarise(Abundance = mean(Abundance)) %>%
  top_n(25, Abundance) %>%
  pull(Family)
ps_top_family <- subset_taxa(ps, Family %in% top_family)
correlations_df_family <- ps_top_family %>% 
  tax_model(
    trans = "clr",
    rank = "Family", variables = list("MUAC","WEIGHTkg", "HEIGHTcm"), 
    type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE
  ) %>% 
  tax_models2stats(rank = "Family")

# Building correlations for genus
top_genus <- ps %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) %>%
  top_n(25, Abundance) %>%
  pull(Genus)
ps_top_genus <- subset_taxa(ps, Genus %in% top_genus)
correlations_df_genus <- ps_top_genus %>% 
  tax_model(
    trans = "clr",
    rank = "Family", variables = list("MUAC","WEIGHTkg", "HEIGHTcm"), 
    type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE
  ) %>% 
  tax_models2stats(rank = "Genus")


generate_correlation_heatmap <- function(correlations_df) {
  
# Get nice looking ordering of correlation estimates using hclust
  taxa_hclust <- correlations_df %>% 
    dplyr::select(term, taxon, estimate) %>% 
    tidyr::pivot_wider(
      id_cols = taxon, names_from = term, values_from = estimate
    ) %>% 
    tibble::column_to_rownames("taxon") %>% 
    as.matrix() %>% 
    stats::dist(method = "euclidean") %>% 
    hclust(method = "ward.D2") 
  
  taxa_order <- taxa_hclust$labels[taxa_hclust$order]
  
# Create the correlation heatmap
  heatmap_plot <- correlations_df %>% 
    mutate(p.FDR = p.adjust(p.value, method = "fdr")) %>% 
    ggplot(aes(x = term, y = taxon)) +
    geom_raster(aes(fill = estimate)) +
    geom_point(
      data = function(x) filter(x, p.value < 0.05),
      shape = "asterisk"
    ) +
    geom_point(
      data = function(x) filter(x, p.FDR < 0.05),
      shape = "circle", size = 3
    ) +
    scale_y_discrete(limits = taxa_order) +
    scale_fill_distiller(palette = "BrBG", limits = c(-0.45, 0.45)) + 
    theme_minimal() +
    labs(
      x = NULL, y = NULL, fill = "Spearman's\nRank\nCorrelation",
      caption = paste(
        "Asterisk indicates p < 0.05, not FDR adjusted",
        "Filled circle indicates FDR corrected p < 0.05", sep = "\n"
      ))
  
  return(heatmap_plot)
}

# Generate heatmaps for each taxonomic level
heatmap_phylum <- generate_correlation_heatmap(correlations_df_phylum)
heatmap_family <- generate_correlation_heatmap(correlations_df_family)
heatmap_genus <- generate_correlation_heatmap(correlations_df_genus)

# Display each plot
heatmap_phylum
heatmap_family
heatmap_genus
```
