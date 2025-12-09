README for the Microbiome Analysis Code

Poster: PDF of poster presented in the spring of 2025

ps_clean: The original phyloseq object without being cleaned

Categorized data vs. Continuous Data:
- Cleaned datasets
- Versions of the same dataset where variables are either left continous or turned into categories for analysis/classification

Statistics Tables: 
- Tables showing characteristics of the study population
- Shows the breakdown of each variable used in analysis and N in each group

Data Cleaning and Organization:
- dat.work (1-3) is just some data cleaning on the original dataset
- Things such as: grouping for tables, renaming variables, changing 1's and 0's to yes's and no's, and creating proxy variables

Data Analysis: 
- The bulk of the coding work is found in this folder
- proxy_analyis: running the main analysis pipeline on the proxy variables
- proxy: creating the proxy variables from other variables
- PERMANOVA_results: Spreadsheet of the PERMANOVA results for each categorical variable for Beta Diversity
- beta_diversity: Updated beta diversity script from the main pipeline (analysis_script)
- analysis_script: The main pipeline that contains every type of analysis found in this study
- Figures: folder containing all of the figures in the study

analysis_scipt: The functions are set to make figures and save them to your computer (change the working directory calls to match your file paths). The functions can either be looped through for a large number of variables, or you can just run it once for an individual variable, so run or comment out certain sections as needed. 

