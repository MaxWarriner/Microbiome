
ps <- read_rds('categorized_data.RDS')

sam <- ps@sam_data

sam <- data.frame(sam)

sam <- sam |>
  mutate(SEI_proxy = ifelse((Underweight == "Underweight" & (Maternal_Education_Status != "Higher Education" & !is.na(Maternal_Education_Status) & Potable_Water_in_House == 'no')) | Electricity_in_House == 'no', "Low", "High"))

sam <- sam |>
  mutate(Cleanliness_proxy = ifelse((Child_Has_Dirty_Fingernails != "no" & Frequency_of_Bathing_in_River != "never" & Use_of_Toilet_Paper_after_Defecation != "always" & Frequency_of_Hand_Washing_After_Using_Toilet != "always" & Frequency_of_Using_Soap_After_Using_Toilet != "always" & Frequency_of_Using_Soap_Before_Eating != "always") | Use_of_Toilet_Paper_after_Defecation == "never" | Frequency_of_Bathing_in_River == 'always' | Frequency_of_Hand_Washing_After_Using_Toilet == 'never' | Frequency_of_Clothes_Washing_in_River == 'always', "Unclean", "Clean"))

sample_data(ps) <- sam

write_rds(ps, 'categorized_data.RDS')

