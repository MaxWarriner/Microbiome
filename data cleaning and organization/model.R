library(phyloseq)
library(glmnet)
library(caret)
library(dplyr)

ps <- readRDS("categorized_data.rds") # read in ps file

#Rarefy Data
set.seed(13)
ps <- rarefy_even_depth(ps, 1e+05)

a_diversity <- estimate_richness(ps, measures = c("Shannon", "Chao1")) # calculate alpha diversity of each sample

sample_data <- ps@sam_data #extract sample data

sample_data$sample <- rownames(sample_data) #create column with sample

a_diversity$sample <- rownames(sample_data) #create same column in diversity object to be able to merge by

sample_data <- data.frame(sample_data) #convert sample data to data frame for merging

a_diversity <- merge(sample_data, a_diversity, by = c("sample")) #merge the data together

a_diversity <- a_diversity |> # remove sample column now that it's not needed
  dplyr::select(-sample)

a_diversity <- mutate_if(a_diversity, is.character, as.factor) #change all character columns to factor

a_diversity <- a_diversity %>%
  dplyr::select_if(~ !(is.factor(.) & nlevels(.) == 1)) # remove all columns where factor only has one level



a_diversity <- a_diversity |>
  mutate(Chao_plus = Chao1/mean(a_diversity$Chao1)*100,
         Shannon_plus = Shannon/mean(a_diversity$Shannon)*100,
         Chao_Shannon = (Chao_plus + Shannon_plus)/2)

a_diversity <- a_diversity %>%
  dplyr::select_if(~ length(unique(.)) > 1) # remove all columns where numeric has only one level

remove_na_columns <- function(df) {
  # Check which columns have any NA values
  na_columns <- sapply(df, function(x) any(is.na(x)))
  
  # Keep only columns without NAs
  df_clean <- df[, !na_columns, drop = FALSE]
  
  return(df_clean)
}

a_diversity <- remove_na_columns(a_diversity)

model <- lm(Chao_Shannon ~ ., data = a_diversity[, c(1:68, 74)])


removeAliased <- function(mod) {
  # Check if there are any aliased terms
  aliased_info <- alias(mod)
  if (is.null(aliased_info$Complete)) {
    return(mod)  # Return original model if no aliasing
  }
  
  # Retrieve aliased terms
  rn <- rownames(aliased_info$Complete)
  
  # Handle factor levels in coefficient names
  factor_vars <- mod$model[, sapply(mod$model, is.factor), drop = FALSE]
  regex.base <- if (ncol(factor_vars) > 0) {
    unique(unlist(lapply(factor_vars, levels)))
  } else {
    character(0)
  }
  
  # Clean coefficient names to match formula terms
  aliased <- if (length(regex.base) > 0) {
    gsub(paste0(regex.base, "$", collapse = "|"), "", 
         gsub(paste0(regex.base, ":", collapse = "|"), ":", rn))
  } else {
    rn
  }
  
  # Skip if no aliased terms left after cleaning
  if (length(aliased) == 0) {
    return(mod)
  }
  
  # Update formula (safely)
  uF <- reformulate(".", response = ".", 
                    drop.terms = terms(reformulate(aliased)), 
                    intercept = attr(mod$terms, "intercept"))
  
  update(mod, uF)
}

model <- removeAliased(model) # remove colinear variables from data


#Perform Multivariable Linear Regression
summary(model)

#perform stepwise regression
library(MASS)
model = stepAIC(model, direction = "both", trace = FALSE) # Stepwise regression model with shannon measure
model_summary <- summary(model)
model_coeffs <- data.frame(model_summary$coefficients)

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1)) #function that rounds all the numbers in the dataframe
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

model_coeffs <- round_df(model_coeffs, 4) #round each data frame to four decimal places

model_coeffs <- model_coeffs[-1,]

library(pwr)
r2  = summary(model)$r.squared          # R-squared for our linear model
my.f2 = r2 / (1-r2)                         # Effect Size

print(r2)
print(my.f2)
