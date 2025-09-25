

sam <- data.frame(sam)

sam <- sam |>
  mutate(SEI_proxy = ifelse((Underweight == "Underweight" & (Maternal_Education_Status != "Higher Education" & Maternal_Education_Status != "High School" & !is.na(Maternal_Education_Status)) & House_Floor_Material != "Cement" & !is.na(House_Floor_Material)), "Low", "High"))

sam <- sam |>
  mutate(Cleanliness_proxy = ifelse((Child_Has_Dirty_Fingernails != "no" & Frequency_of_Bathing_in_River != "never" & Use_of_Toilet_Paper_after_Defecation != "always" & Frequency_of_Hand_Washing_After_Using_Toilet != "always" & Frequency_of_Using_Soap_After_Using_Toilet != "always" & Frequency_of_Using_Soap_Before_Eating != "always") | Use_of_Toilet_Paper_after_Defecation == "never", "Unclean", "Clean"))

sample_data(ps) <- sam

write_rds(ps, 'categorized_data.RDS')

