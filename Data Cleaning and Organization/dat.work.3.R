library(phyloseq)
library(dplyr)

ps <- readRDS("categorized_data.RDS")

sam <- sample_data(ps)

sam <- data.frame(sam)

sam <- sam |>
  dplyr::select(-FinalFamOcc, -OpenFieldDefFinal, -SoilConsumpFinal, -SoilConsumpThinning, -AgeGroupUnderweight, -MatEdThinning, -OpenFieldDefThinning,-WashFruitsThinned, -Foods, -SevStunting, -SevThinning, -logWeight, -logHeight, -logHAZ, -logWAZ, -logBAZ, -Log1Weight, -Log1Height, -Log1WAZ, -Log1HAZ, -Log1BAZ, -lnWeight, -lnHeight, -lnWAZ, -lnBAZ, -sqrtWeight, -sqrtHeight, -sqrtWAZ, -sqrtWAZ, -sqrtHAZ, -sqrtBAZ, -OversqrtWeight, -OversqrtHeight, -OversqrtWAZ, -OversqrtHAZ, -OversqrtBAZ, -sample_data, -lnHAZ)


sam <- sam |>
  mutate(SEI_proxy = ifelse((Maternal_Education_Status == "Illiterate" | Maternal_Education_Status == "Primary School") & House_Floor_Material != "Cement" & Kitchen_Separated.In_House == "Separated" & Underweight == "Underweight", "Lower","Upper"), 
         Cleanliness_proxy = ifelse(How_often_do_you_trim_your_fingernails. != "Once per Week" & Use_of_Toilet_Paper_after_Defecation != "always" & Frequency_of_Using_Soap_After_Using_Toilet != "always" & Frequency_of_Using_Soap_Before_Eating != "always" & Freqency_of_Washing_Raw.Undercooked_Vegetables != "always", "Unclean", "Clean"))

sam <- sam |>
  rename("Method_of_Washing_Hands_After_Toilet" = Howdoyouwashyourhandsaftertoilet)

sample_data(ps) <- sam

write.csv(sam,"sample_data.csv")

write_rds(ps, "categorized_data.RDS")
