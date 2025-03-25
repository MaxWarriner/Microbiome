library(phyloseq)
library(glmnet)
library(caret)
library(dplyr)

ps <- readRDS("continuous_data.rds") # read in ps file

#Rarefy Data
set.seed(13)
ps <- rarefy_even_depth(ps, 1e+05)

a_diversity <- estimate_richness(ps) # calculate alpha diversity of each sample

sample_data <- ps@sam_data #extract sample data

sample_data$sample <- rownames(sample_data) #create column with sample

a_diversity$sample <- rownames(sample_data) #create same column in diversity object to be able to merge by

sample_data <- data.frame(sample_data) #convert sample data to data frame for merging

a_diversity <- merge(sample_data, a_diversity, by = c("sample")) |> #merge the data together
  dplyr::select(-sample_data)

a_diversity <- a_diversity %>%       # remove any columns where there is an NA 
  dplyr::select(where(~!any(is.na(.))))

a_diversity <- a_diversity |> # remove sample column now that it's not needed
  dplyr::select(-sample)

a_diversity <- mutate_if(a_diversity, is.character, as.factor) #change all character columns to factor

a_diversity <- a_diversity %>%
  dplyr::select_if(~ !(is.factor(.) & nlevels(.) == 1)) # remove all columns where factor only has one level

a_diversity <- a_diversity %>%
  dplyr::select_if(~ length(unique(.)) > 1) # remove all columns where numeric has only one level

Shannon =lm(Shannon ~ . , data=a_diversity) # generate linear models of data with all the variables
Chao1 =lm(Chao1 ~ . , data=a_diversity)


removeAliased <- function(mod) { # removes variables with colinearity
  ## Retrieve all terms in the model
  X <- attr(mod$terms, "term.label")
  ## Get the aliased coefficients  
  rn <- rownames(alias(mod)$Complete)
  ## remove factor levels from coefficient names to retrieve the terms
  regex.base <- unique(unlist(lapply(mod$model[, sapply(mod$model, is.factor)], levels)))
  aliased <- gsub(paste(regex.base, "$", sep = "", collapse = "|"),  "", gsub(paste(regex.base, ":", sep = "", collapse = "|"), ":", rn))
  uF <- formula(paste(". ~ .", paste(aliased, collapse = "-"), sep = "-"))
  update(mod, uF) 
}

Shannon <- removeAliased(Shannon) # remove colinear variables from data
Chao1 <- removeAliased(Chao1)

#Perform Multivariable Linear Regression
summary(Shannon)
summary(Chao1)

#perform stepwise regression
library(MASS)
Shannon = stepAIC(Shannon, direction = "both", trace = FALSE) # Stepwise regression model with shannon measure
Shannon_summary <- summary(Shannon)
shannon_coeffs <- data.frame(Shannon_summary$coefficients)

Chao1 = stepAIC(Chao1, direction = "both", trace = FALSE) #stepwise regression model with chao1 measure
Chao1_summary <- summary(Chao1)
Chao1_coeffs <- data.frame(Chao1_summary$coefficients)

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1)) #function that rounds all the numbers in the dataframe
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

shannon_coeffs <- round_df(shannon_coeffs, 4) #round each data frame to four decimal places
Chao1_coeffs <- round_df(Chao1_coeffs, 4)


#separate data frames based on which factor group they're apart of

demographic_shannon <- shannon_coeffs[c("Age", "weight", "Height", "ModeofdeliveryVaginal", "EverhadvaccinatedYes", "BCGscarYES", "ChildsfingernailtrimmedYes", "ArechildsfingernailsdirtyYes", "HowoftendoyoutrimyourfingernailsOnce per Week", "HowoftendoyoutrimyourfingernailsOnce per two weeks", "HowoftendoyoutrimyourfingernailsDon't Know", "DidyourparentsteachersorhealthprofessionalsgaveyouadewormingpillNo", "DidyourparentsteachersorhealthprofessionalsgaveyouadewormingpillYes", "DidyourparentsorhealthprofessionalsgaveyouotherantibioticsNo", "DidyourparentsorhealthprofessionalsgaveyouotherantibioticsYes"),]
socioeconomic_shannon <- shannon_coeffs[c("Household.Number","Siblings.Younger.than.12", "HeardALnamebefore", "HeardHIVnamebefore", "HeardInWormnamebefore", "HeardMalanamebefore","HeardTBnamebefore", "HeardSChnamebefore", "familytold", "HPtold", "Teachertold", "Mediatold", "DoyouknowhowintestinalwomstransmittedNo", "DoyouknowhowintestinalwomstransmittedYes", "DoyouknowwhywormsarebadforyourhealthNo", "DoyouknowwhywormsarebadforyourhealthYes", "YourlivingaddressHermata mentina", "YourlivingaddressJiren", "YourlivingaddressMentina", "YourlivingaddressOther", "YourlivingaddressSeto semero", "YourlivingaddressWelda", "FamilyoccupationDailly labor", "FamilyoccupationFarmer", "FamilyoccupationGovernment employee", "FamilyoccupationMerchant", "FamilyoccupationOther", "FamilyoccupationPolice man", "FamilyoccupationSecurity", "FamilyoccupationTeacher", "MaternaleducationalstatusPrimary School", "MaternaleducationalstatusHigh School", "MaternaleducationalstatusHigher Education", "MaternaleducationalstatusDon't Know", "DoyouhaveelectricityinyourhaouseNo", "DoyouhaveelectricityinyourhaouseYes", "DoesyourfamillyownradioNo", "DoesyourfamillyownradioYes", "DoesyourfamillymemberownaphoneNo", "DoesyourfamillymemberownaphoneYes"),]
environmental_shannon <- shannon_coeffs[c("IfthereisatoiletdoesthelatrinehavedoorsYes", "FliesobservedinaroundthelatrineYes", "WhatmaterialyourhousefloormadefromDon't Know", "WhatmaterialyourhousefloormadefromDust", "WhatmaterialyourhousefloormadefromPlastic Covered", "IsyourkitcheninyourhouseorsepatatedSeparated", "IsyourkitcheninyourhouseorsepatatedWithin House", "IfseparetedwhatmaterialsyourkitchenmadefromDon't Know", "IfseparetedwhatmaterialsyourkitchenmadefromDust", "RoofYes", "WallYes", "NoneYes", "GasYes", "CoalYes", "KeroseneYes", "ElectricYes", "CattleYes", "SheepGoatYes", "ChickenYes", "PetYes", "NodomanimalYes", "DoyouhavepotablewaterinyourhouseNo", "DoyouhavepotablewaterinyourhouseYes", "DoyoudrinkdirectlyordoyoutreatDon't Know", "DoyoudrinkdirectlyordoyoutreatTreat", "DoyourfamilyownlatrineYes", "YourlatrineconnectedtoNA/Don't Know", "YourlatrineconnectedtoRiver", "YourlatrineconnectedtoSewage", "YourlatrineconnectedtoWell", "HowoftendoyoubathinriverSometimes", "HowoftendoyoubathinriverAlways", "HowoftendoyoubathinriverNA/Don't Know", "HowoftendoyouwashyourclothesinriverSometimes", "HowoftendoyouwashyourclothesinriverAlways", "DoyoudeficateintheopenfieldSometimes", "DoyoudeficateintheopenfieldNA/Don't Know", "DoyouuseschoollatrineSometimes", "DoyouuseschoollatrineAlways", "DoyouuseschoollatrineNA/Don't Know", "DoyouusetoiletpapertowipeyourbumafterdeficationSometimes", "DoyouusetoiletpapertowipeyourbumafterdeficationAlways", "DoyouusetoiletpapertowipeyourbumafterdeficationNA/Don't Know", "DoyouwashyourhandsaftertoiletSometimes", "DoyouwashyourhandsaftertoiletAlways", "HowdoyouwashyourhandsaftertoiletSoap and Water", "HowdoyouwashyourhandsaftertoiletNA/Don't Know", "HowdoyouwashyourhandsbeforeeatingSoap and Water", "HowdoyouwashyourhandsbeforeeatingNA/Don't Know", "Ifwithsoaphowoftenyouuseit_AAlways", "Ifwithsoaphowoftenyouuseit_ANA/Don't Know", "DoyoueatsoilSometimes", "DoyoueatsoilNA/Don't Know", "DoyouwashyourfruitesbeforeeatingSometimes", "DoyouwashyourfruitesbeforeeatingAlways", "DoyouwashyourfruitesbeforeeatingNA/Don't Know", "IfyeshowoftenyouwashyourvegitablesbeforeeatingSometimes", "IfyeshowoftenyouwashyourvegitablesbeforeeatingAlways", "IfyeshowoftenyouwashyourvegitablesbeforeeatingNA/Don't Know", "DoyouwalkbarefootSometimes", "DoyouwalkbarefootAlways", "DoyouwalkbarefootNA/Don't Know", "whenyouareathomedoyoupreferetousesandalsorshoeSandals", "whenyouareathomedoyoupreferetousesandalsorshoeShoes", "whenyouareathomedoyoupreferetousesandalsorshoeNA/Don't Know"),]

demographic_chao1 <- Chao1_coeffs[c("Age", "SexMale", "weight", "Height", "EverhadvaccinatedYes", "BCGscarYES", "ChildsfingernailtrimmedYes", "HowoftendoyoutrimyourfingernailsOnce per Week", "HowoftendoyoutrimyourfingernailsOnce per two weeks", "HowoftendoyoutrimyourfingernailsDon't Know", "DidyourparentsteachersorhealthprofessionalsgaveyouadewormingpillNo", "DidyourparentsteachersorhealthprofessionalsgaveyouadewormingpillYes", "DidyourparentsorhealthprofessionalsgaveyouotherantibioticsNo", "DidyourparentsorhealthprofessionalsgaveyouotherantibioticsYes"), ]
socioeconomic_chao1 <- Chao1_coeffs[c("Household.Number", "Siblings.Younger.than.12", "HeardALnamebefore", "HeardTTnamebefore", "HeardHWnamebefore", "HeardInWormnamebefore", "HeardMalanamebefore", "HeardTBnamebefore", "familytold", "HPtold", "Teachertold", "Mediatold", "DoyouknowwhywormsarebadforyourhealthNo", "DoyouknowwhywormsarebadforyourhealthYes", "YourlivingaddressHermata mentina", "YourlivingaddressJiren", "YourlivingaddressMentina", "YourlivingaddressOther", "YourlivingaddressSeto semero", "YourlivingaddressWelda", "FamilyoccupationDailly labor", "FamilyoccupationFarmer", "FamilyoccupationGovernment employee", "FamilyoccupationMerchant", "FamilyoccupationOther", "FamilyoccupationPolice man", "FamilyoccupationSecurity", "FamilyoccupationTeacher", "MaternaleducationalstatusPrimary School", "MaternaleducationalstatusHigh School", "MaternaleducationalstatusHigher Education", "MaternaleducationalstatusDon't Know", "DoyouhaveelectricityinyourhaouseNo", "DoyouhaveelectricityinyourhaouseYes", "DoesyourfamillyownradioNo", "DoesyourfamillyownradioYes", "DoesyourfamillyowntelevisionNo", "DoesyourfamillyowntelevisionYes", "DoesyourfamillymemberownaphoneNo", "DoesyourfamillymemberownaphoneYes"),]
environmental_shannon <- Chao1_coeffs[c("FliesobservedinaroundthelatrineYes", "WhatmaterialyourhousefloormadefromDon't Know", "WhatmaterialyourhousefloormadefromDust", "WhatmaterialyourhousefloormadefromPlastic Covered", "IsyourkitcheninyourhouseorsepatatedSeparated", "IsyourkitcheninyourhouseorsepatatedWithin House", "IfseparetedwhatmaterialsyourkitchenmadefromDon't Know", "IfseparetedwhatmaterialsyourkitchenmadefromDust", "RoofYes", "WallYes", "NoneYes", "WoodYes", "GasYes", "CoalYes", "KeroseneYes", "ElectricYes", "CattleYes", "SheepGoatYes", "ChickenYes", "PetYes", "NodomanimalYes", "DoyouhavepotablewaterinyourhouseNo", "DoyouhavepotablewaterinyourhouseYes", "DoyoudrinkdirectlyordoyoutreatDon't Know", "DoyoudrinkdirectlyordoyoutreatTreat", "DoyourfamilyownlatrineYes", "YourlatrineconnectedtoNA/Don't Know", "YourlatrineconnectedtoRiver", "YourlatrineconnectedtoSewage", "YourlatrineconnectedtoWell", "HowoftendoyoubathinriverSometimes", "HowoftendoyoubathinriverAlways", "HowoftendoyoubathinriverNA/Don't Know", "HowoftendoyouwashyourclothesinriverSometimes", "HowoftendoyouwashyourclothesinriverAlways", "DoyoudeficateintheopenfieldSometimes", "DoyoudeficateintheopenfieldNA/Don't Know", "DoyouuseschoollatrineSometimes", "DoyouuseschoollatrineAlways", "DoyouuseschoollatrineNA/Don't Know", "DoyouusetoiletpapertowipeyourbumafterdeficationSometimes", "DoyouusetoiletpapertowipeyourbumafterdeficationAlways", "DoyouusetoiletpapertowipeyourbumafterdeficationNA/Don't Know", "DoyouwashyourhandsaftertoiletSometimes", "DoyouwashyourhandsaftertoiletAlways", "HowdoyouwashyourhandsaftertoiletSoap and Water", "HowdoyouwashyourhandsaftertoiletNA/Don't Know", "HowdoyouwashyourhandsbeforeeatingSoap and Water", "HowdoyouwashyourhandsbeforeeatingNA/Don't Know", "Ifwithsoaphowoftenyouuseit_AAlways", "Ifwithsoaphowoftenyouuseit_ANA/Don't Know", "DoyoueatsoilSometimes", "DoyoueatsoilNA/Don't Know", "DoyouwashyourfruitesbeforeeatingSometimes", "DoyouwashyourfruitesbeforeeatingAlways", "DoyouwashyourfruitesbeforeeatingNA/Don't Know", "IfyeshowoftenyouwashyourvegitablesbeforeeatingSometimes", "IfyeshowoftenyouwashyourvegitablesbeforeeatingAlways", "IfyeshowoftenyouwashyourvegitablesbeforeeatingNA/Don't Know", "DoyouwalkbarefootSometimes", "DoyouwalkbarefootAlways", "DoyouwalkbarefootNA/Don't Know", "whenyouareathomedoyoupreferetousesandalsorshoeSandals", "whenyouareathomedoyoupreferetousesandalsorshoeShoes", "whenyouareathomedoyoupreferetousesandalsorshoeNA/Don't Know"),]
