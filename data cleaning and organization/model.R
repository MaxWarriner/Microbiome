library(phyloseq)
library(glmnet)
library(caret)
library(dplyr)

ps <- readRDS("continuous_data.rds") # read in ps file

#Rarefy Data
set.seed(13)
ps <- rarefy_even_depth(ps, 1e+05)

a_diversity <- estimate_richness(ps, measures = c("Shannon", "Chao1")) # calculate alpha diversity of each sample

sample_data <- ps@sam_data #extract sample data

sample_data$sample <- rownames(sample_data) #create column with sample

a_diversity$sample <- rownames(sample_data) #create same column in diversity object to be able to merge by

sample_data <- data.frame(sample_data) #convert sample data to data frame for merging

a_diversity <- merge(sample_data, a_diversity, by = c("sample")) #merge the data together

a_diversity <- a_diversity |> # remove unnecessary columns for the model
  dplyr::select(-sample, -se.chao1)

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

a_diversity <- a_diversity |>
  dplyr::select(-Chao1, -Chao_plus, -Shannon, -Shannon_plus)


a_diversity <- a_diversity[,c(1:16, 33:46, 54:78, 88:93, 116)]

model <- lm(Chao_Shannon ~ ., data = a_diversity, na.action = na.omit)

summary(model)

top10_rsquared <- function(model) {
  # Check if the input is a linear model
  if (!inherits(model, "lm")) {
    stop("Input must be a linear model (lm) object")
  }
  
  # Get the model data and formula
  model_data <- model.frame(model)
  response_var <- as.character(formula(model)[[2]])
  
  # Get all predictor variables
  predictors <- setdiff(colnames(model_data), response_var)
  
  # If there are 10 or fewer predictors, return them all
  if (length(predictors) <= 10) {
    message("Model has 10 or fewer predictors - returning all variables")
    return(predictors)
  }
  
  # Calculate individual R-squared contributions
  rsq_contributions <- sapply(predictors, function(var) {
    # Create a formula with just this predictor
    f <- as.formula(paste(response_var, "~", var))
    # Fit a simple model
    m <- lm(f, data = model_data)
    # Return R-squared
    summary(m)$r.squared
  })
  
  # Sort by R-squared in descending order
  sorted_vars <- names(sort(rsq_contributions, decreasing = TRUE))
  
  # Return the top 10
  return(sorted_vars[1:10])
}


top10_rsquared(model)

best_model <- lm(Chao_Shannon ~ Familyoccupation + Yourlatrineconnectedto + Doyouuseschoollatrine + Ifwithsoaphowoftenyouuseit +
                   Ifwithsoaphowoftenyouuseit_A + Doyoueatraworundercookedvegitables + Underweight + Isyourkitcheninyourhouseorsepatated + 
                   Howoftendoyoubathinriver + Doyoudeficateintheopenfield, data = a_diversity)

#Perform Multivariable Linear Regression
summary(best_model)

library(mlbench)     # For PimaIndiansDiabetes2 dataset
library(dplyr)       # For data manipulation (dplyr) 
library(broom)       # For making model summary tidy
library(visreg)      # For plotting logodds and probability 
library(rcompanion)  # To calculate pseudo R-squared
library(MASS)        # For stepwise model selection
library(ROCR)        # To find the probability threshold for best accuracy
library(car)         # For multicollinearity function vif()

nagelkerke(best_model)


#Normality of Residuals
hist(best_model$residuals, freq=FALSE, xlab = "Residuals", main="",cex.axis=1.5,cex.lab=1.5)
curve(dnorm(x,mean(best_model$residuals), sd(best_model$residuals)), -1.6, 2.1, add=TRUE,col=c('red'))
qqnorm(best_model$residuals, main="")
qqline(best_model$residuals)
shapiro.test(best_model$residuals)

plot(a_diversity$weight, a_diversity$Chao_Shannon)
