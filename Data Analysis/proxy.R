
library(tidyverse)
library(phyloseq)

ps <- read_rds('categorized_data.RDS')

sam <- ps@sam_data

sam <- data.frame(sam)

sam <- sam |>
  rename(hygiene_index = Cleanliness_proxy, 
         SES_index = SEI_proxy)

sam$hygiene_index <- rep(0, 138)
sam$SES_index <- rep(0,138)
sam$HIQ_index <- rep(0,138)
sam$HSI_index <- rep(0,138)

for (i in 1:138){
  
  if(sam$Frequency_of_Hand_Washing_After_Using_Toilet[i] == 'always'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 2
  }else if(sam$Frequency_of_Hand_Washing_After_Using_Toilet[i] == 'sometimes'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Method_of_Washing_Hands_After_Toilet[i] == 'Soap and Water'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 2
  }else if(sam$Method_of_Washing_Hands_After_Toilet[i] == 'Water'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Frequency_of_Using_School_Latrine[i] == 'always' | sam$Frequency_of_Using_School_Latrine[i] == 'sometimes'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Drink_Water_Directly_or_Treated.[i] == 'Treat'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Defecating_in_Open_Field[i] == 'never'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Child.s_Fingernails_Trimmed[i] == 'yes' & sam$Child_Has_Dirty_Fingernails[i] == 'no'){
    sam$hygiene_index[i] = sam$hygiene_index[i] + 1
  }
  
  if(sam$Maternal_Education_Status[i] == 'Higher Education'){
    sam$SES_index[i] = sam$SES_index[i] + 3
  }else if(sam$Maternal_Education_Status[i] == 'High School'){
    sam$SES_index[i] = sam$SES_index[i] + 2
  }else if(sam$Maternal_Education_Status[i] == 'Primary School'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$Electricity_in_House[i] == 'yes'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$Family_Owns_Radio[i] == 'yes'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$Family_Owns_Television[i] == 'yes'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$Family_Member_with_Phone[i] == 'yes'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$Urban.Rural[i] == 'Urban'){
    sam$SES_index[i] = sam$SES_index[i] + 1
  }
  
  if(sam$House_Floor_Material[i] == 'Cement' | sam$House_Floor_Material[i] == "Plastic Covered"){
    sam$HIQ_index[i] = sam$HIQ_index[i] + 2
  }
  
  if(sam$Connection_to_Latrine[i] == 'Sewage'){
    sam$HIQ_index[i] = sam$HIQ_index[i] + 2
  }else if(sam$Connection_to_Latrine[i] == 'Well'){
    sam$HIQ_index[i] = sam$HIQ_index[i] + 1
  }
  
  if(sam$Potable_Water_in_House[i] == 'yes'){
    sam$HIQ_index[i] = sam$HIQ_index[i] + 1
  }
  
  if(sam$Kitchen_Has_Roof[i] == 'yes' & sam$Kitchen_Has_Wall[i] == 'yes'){
    sam$HIQ_index[i] = sam$HIQ_index[i] + 1
  }
  
  if(sam$Family_Owns_Latrine[i] == 'yes'){
    sam$HSI_index[i] = sam$HSI_index[i] + 1
  }
  
  if(sam$Latrine_Distance_from_House[i] == '>20' | sam$Latrine_Distance_from_House[i] == '10-20'){
    sam$HSI_index[i] = sam$HSI_index[i] + 2
  }else if(sam$Latrine_Distance_from_House[i] == '5-10'){
    sam$HSI_index[i] = sam$HSI_index[i] + 1
  }
  
  if(sam$Defecating_in_Open_Field[i] == 'never'){
    sam$HSI_index[i] = sam$HSI_index[i] + 1
  }
  
  
  
}

sam <- sam |>
  mutate(hygiene_group = ifelse(hygiene_index <= 5, 'low', 'high'), 
         SES_group = ifelse(SES_index <= 4, 'low', 'high'), 
         HSI_group = ifelse(HSI_index <= 2, 'low', 'high'), 
         HIQ_group = ifelse(HIQ_index <= 2, 'low', 'high'))

sample_data(ps) <- sam

write_rds(ps, 'categorized_data.RDS')

