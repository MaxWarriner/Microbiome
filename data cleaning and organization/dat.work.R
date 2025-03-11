
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(stringr)
library(vegan)
library(readr)

dat <- readRDS("ps_clean.rds")

sample_data <- data.frame(dat@sam_data)


############################################################
# Age
############################################################
age_mean_sd <- c(mean(sample_data$Age), sd(sample_data$Age))

age_categories_table <- sample_data |>
  dplyr::select(Age) |>
  mutate(Age = case_when(Age >= 6 & Age <= 9 ~ "6-9", 
                         Age >= 10 & Age <= 13 ~ "10-13", 
                         Age >= 14 & Age <= 17 ~ "14-17")) |>
  group_by(Age) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  mutate(Age = factor(Age, levels = c("6-9", "10-13", "14-17"))) |>
  arrange(Age)

ggplot(data = age_categories_table) +
  geom_col(aes(x = Age, y = n), color = "black") +
  scale_x_continuous(breaks = seq(6, 17)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  ggtitle("Distribution of Ages") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))


############################################################
# Sex
############################################################

sex_table <- sample_data |>
  group_by(Sex) |>
  summarize(n = n()) |>
  mutate(frequency = round(n / 138, 3))

ggplot(data = sample_data, aes(x = Sex)) +
  geom_bar(color = "black") +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  ggtitle("Count of Sex") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))


#####################################################################
# Weight
#####################################################################

weight_mean_sd <- c(mean(sample_data$weight), sd(sample_data$weight))

weight_table <- sample_data |>
  mutate(weight = case_when(weight < 35 ~ "15-35", 
                            weight >= 35 & weight < 55 ~ "35-55", 
                            weight >= 55 ~ ">55")) |>
  group_by(weight) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  mutate(weight = factor(weight, levels = c("15-35","35-55", ">55"))) |>
  arrange(weight)

ggplot(data = sample_data, aes(x = weight)) +
  geom_histogram(binwidth = 5, color = "black") +
  theme_minimal() +
  ggtitle("Distribution of Weight") +
  xlab("Weight (Kg)") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))

underweight_table <- sample_data |>
  group_by(Underweight) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  mutate(Underweight = ifelse(Underweight == 1, "Underweight", "Not Underweight"))
  


#####################################################################
# Height
#####################################################################

height_mean_sd <- c(mean(sample_data$Height), sd(sample_data$Height))

height_table <- sample_data |>
  mutate(Height = case_when(Height < 1.15 ~ "< 1.15", 
                            Height >= 1.15 & Height < 1.35 ~ "1.15-1.35", 
                            Height >= 1.35 & Height < 1.55 ~ "1.35-1.55", 
                            Height >= 1.55 ~ ">1.55")) |>
  group_by(Height) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  mutate(Height = factor(Height, levels = c("< 1.15", "1.15-1.35", "1.35-1.55", ">1.55"))) |>
  arrange(Height)

ggplot(data = sample_data, aes(x = Height)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  theme_minimal() +
  ggtitle("Distribution of Height") +
  xlab("Height (m)") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))

stunted_table <- sample_data |>
  group_by(Stunted) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  mutate(Stunted = ifelse(Stunted == 1, "Stunted Growth", "Not Stunted Growth"))
  


#####################################################################
# Positivity
#####################################################################


positivity_table <- sample_data |>
  dplyr::select(Ascaris, Tricuris, Hookworm, Schistosoma, OtherIPs, AnySTHs, AnyParasite) |>
  mutate(Ascaris = case_when(Ascaris == "Postive" ~ 1, 
                             Ascaris == "Negative" ~ 0), 
         Tricuris = case_when(Tricuris == "Postive" ~ 1, 
                              Tricuris == "Negative" ~ 0), 
         Hookworm = case_when(Hookworm == "Postive" ~ 1, 
                              Hookworm == "Negative" ~ 0), 
         Schistosoma = case_when(Schistosoma == "Postive" ~ 1, 
                                 Schistosoma == "Negative" ~ 0), 
         OtherIPs = case_when(OtherIPs == "Postive" ~ 1, 
                              OtherIPs == "Negative" ~ 0), 
         AnySTHs = case_when(AnySTHs == "Postive" ~ 1, 
                             AnySTHs == "Negative" ~ 0), 
         AnyParasite = case_when(AnyParasite == "Postive" ~ 1, 
                                 AnyParasite == "Negative" ~ 0)) |>
  summarise(Ascaris = mean(Ascaris), 
            Tricuris = mean(Tricuris), 
            Hookworm = mean(Hookworm), 
            Schistosoma = mean(Hookworm), 
            OtherIPs = mean(OtherIPs), 
            AnySTHs = mean(AnySTHs), 
            AnyParasite = mean(AnyParasite)) |>
  pivot_longer(cols = c(Ascaris, Tricuris, Hookworm, OtherIPs, AnySTHs, Schistosoma, AnyParasite), names_to = "disease", values_to = "frequency") |>
  mutate(n = round(frequency*138), 
         frequency = round(frequency, 3))

ggplot(data = positivity_table, aes(x = disease, y = frequency)) +
  geom_col(show.legend = FALSE, color = "black") +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  ylim(0,1) +
  scale_x_discrete(labels = c("Any Parasite", "Any STHs", "Ascaris", "Tricuris", "Schistosoma", "Other IPs", "Any STHs", "Any Parasite")) +
  xlab("Disease") +
  ylab("Frequency") +
  ggtitle("Disease Frequency") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))

############################################################################################
# Number in Household
############################################################################################


household_mean_sd <- c(mean(sample_data$Householdno), sd(sample_data$Householdno))

household_table <- sample_data |>
  mutate(Householdno = case_when(Householdno >= 2 & Householdno <= 3 ~ "2-3", 
                                 Householdno >= 4 & Householdno <= 6 ~ "4-6", 
                                 Householdno >= 7 & Householdno <= 9 ~ "7-9", 
                                 Householdno > 9 ~ ">9")) |>
  group_by(Householdno) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Household Size" = Householdno)

ggplot(data = sample_data) +
  geom_histogram(aes(x = Householdno), 
                 binwidth = 1, color = "black") +
  theme_minimal() +
  ggtitle("Distribution of Number in Household") +
  xlab("Number of People in Household") +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA))


#############################################################################################
# Knowledge on Soil-Transmitted Helminths
#############################################################################################

knowledge_table <- sample_data |>
  dplyr::select(HeardALnamebefore, HeardTTnamebefore, HeardHWnamebefore, HeardHIVnamebefore, HeardInWormnamebefore, HeardMalanamebefore, HeardTBnamebefore, HeardSChnamebefore) |>
  mutate(HeardALnamebefore = ifelse(HeardALnamebefore == 0, 1, 0), 
         HeardTTnamebefore = ifelse(HeardTTnamebefore == 1, 1, 0), 
         HeardHWnamebefore = ifelse(HeardHWnamebefore == 2, 1, 0), 
         HeardHIVnamebefore = ifelse(HeardHIVnamebefore == 3, 1, 0), 
         HeardInWormnamebefore = ifelse(HeardInWormnamebefore == 4, 1, 0), 
         HeardMalanamebefore = ifelse(HeardMalanamebefore == 5, 1, 0), 
         HeardTBnamebefore = ifelse(HeardTBnamebefore == 6, 1, 0), 
         HeardSChnamebefore = ifelse(HeardSChnamebefore == 7, 1, 0)) |>
  mutate(HeardALnamebefore = ifelse(is.na(HeardALnamebefore), 0, HeardALnamebefore), 
         HeardTTnamebefore = ifelse(is.na(HeardTTnamebefore), 0, HeardTTnamebefore), 
         HeardHWnamebefore = ifelse(is.na(HeardHWnamebefore), 0, HeardHWnamebefore), 
         HeardHIVnamebefore = ifelse(is.na(HeardHIVnamebefore), 0, HeardHIVnamebefore), 
         HeardInWormnamebefore = ifelse(is.na(HeardInWormnamebefore), 0, HeardInWormnamebefore), 
         HeardMalanamebefore = ifelse(is.na(HeardMalanamebefore), 0, HeardMalanamebefore), 
         HeardTBnamebefore = ifelse(is.na(HeardTBnamebefore), 0, HeardTBnamebefore), 
         HeardSChnamebefore = ifelse(is.na(HeardSChnamebefore), 0, HeardSChnamebefore)) |>
  mutate(None = ifelse(HeardALnamebefore == 0 & HeardTTnamebefore == 0 & HeardHWnamebefore == 0 & HeardHIVnamebefore == 0 & HeardInWormnamebefore == 0 & HeardMalanamebefore == 0 & HeardTBnamebefore == 0 & HeardSChnamebefore == 0, 1, 0)) |>
  summarise(type = c("Ascharis", "Trichuris", "Hookworms", "HIV/AIDS", "Intestinal Worms", "Malaria", "Tuberculosis", "Schistosoma/Bilharzia", "None"), 
            "Heard Of" = c(round(mean(HeardALnamebefore), 3), round(mean(HeardTTnamebefore), 3), round(mean(HeardHWnamebefore), 3), round(mean(HeardHIVnamebefore), 3), round(mean(HeardInWormnamebefore), 3), round(mean(HeardMalanamebefore), 3), round(mean(HeardTBnamebefore), 3), round(mean(HeardSChnamebefore), 3), round(mean(None), 3))) |>
  mutate(count = round(`Heard Of` * 138))

##############################################################################################
# How do you know about these disease? 
##############################################################################################

told_table <- sample_data |>
  dplyr::select(familytold, HPtold, Teachertold, Mediatold) |>
  mutate(Dont_know = ifelse(familytold == 99, 1, 0),
         familytold = ifelse(familytold == 0, 1, 0), 
         HPtold = ifelse(HPtold == 1, 1, 0), 
         Teachertold = ifelse(Teachertold == 2, 1, 0), 
         Mediatold = ifelse(Mediatold == 3, 1, 0)) |>
  mutate(familytold = ifelse(is.na(familytold), 0, familytold), 
         HPtold = ifelse(is.na(HPtold), 0, HPtold), 
         Teachertold = ifelse(is.na(Teachertold), 0, Teachertold), 
         Mediatold = ifelse(is.na(Mediatold), 0, Mediatold), 
         Dont_know = ifelse(is.na(Dont_know), 0, Dont_know)) |>
  summarise(type = c("Family Told", "Health Professional Told", "Teacher Told", "Media Told", "Don't Know"), 
            frequency = c(round(mean(familytold), 3), round(mean(HPtold), 3), round(mean(Teachertold), 3), round(mean(Mediatold), 3), round(mean(Dont_know), 3))) |>
  mutate(count = round(frequency * 138)) |>
  rename("Who Told You?" = type)





##############################################################################################
# Where do you live?
##############################################################################################

living_table <- sample_data |>
  group_by(WheredoyouliveGROUP) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3))


###############################################################################################
# Maternal Education Status
###############################################################################################

education_table <- sample_data |> 
  mutate(Maternaleducationalstatus = case_when(Maternaleducationalstatus == 0 ~ "Illiterate", 
                                               Maternaleducationalstatus == 1 ~ "Primary School", 
                                               Maternaleducationalstatus == 2 ~ "High School", 
                                               Maternaleducationalstatus == 3 ~ "Higher Education", 
                                               Maternaleducationalstatus == 99 ~ "Don't Know", 
                                               is.na(Maternaleducationalstatus) ~ "Don't Know")) |>
  group_by(Maternaleducationalstatus) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  mutate(Maternaleducationalstatus = factor(Maternaleducationalstatus, levels = c("Illiterate", "Primary School", "High School", "Higher Education", "Don't Know"))) |>
  arrange(Maternaleducationalstatus)


#################################################################################
# De-worming
#################################################################################

deworming <- sample_data |>
  dplyr::slec(Dewormingin6mo, Dewormingin1yr) |>
  mutate(Dewormingin6mo = ifelse(Dewormingin6mo == "YES", 1, 0)) |>
  summarise(Deworming = c("Last 6 Months", "Last 1 Year"), 
            count = c(sum(Dewormingin6mo), sum(Dewormingin1yr)), 
            frequency = c(round(mean(Dewormingin6mo), 3), round(mean(Dewormingin1yr), 3)))



################################################################################
# Mode of Delivery
################################################################################

delivery <- sample_data |>
  mutate(Modeofdelivery = ifelse(Modeofdelivery == 1, "Vaginal", "C-Section")) |>
  group_by(Modeofdelivery) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3))

################################################################################
# Older Brothers / Sisters
################################################################################

older_siblings <- sample_data |>
  mutate(Noofolderbrothres = case_when(Noofolderbrothres >= 1 & Noofolderbrothres <= 2 ~ "1-2", 
                                      Noofolderbrothres >= 3 & Noofolderbrothres <= 4 ~ "3-4", 
                                      Noofolderbrothres > 4 ~ ">4", 
                                      Noofolderbrothres == 0 ~ "0")) |>
  group_by(Noofolderbrothres) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Older Siblings" = Noofolderbrothres)

################################################################################
# Siblings younger than 12
################################################################################

young_siblings <- sample_data |>
  mutate(youngerthan12years = case_when(youngerthan12years >= 1 & youngerthan12years <= 2 ~ "1-2", 
                                       youngerthan12years >= 3 & youngerthan12years <= 4 ~ "3-4", 
                                       youngerthan12years > 4 ~ ">4", 
                                       youngerthan12years == 0 ~ "0")) |>
  group_by(youngerthan12years) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Siblings Younger than 12 Years Old" = youngerthan12years)


################################################################################
# Ever Had Vaccination? 
################################################################################

vaccination <- sample_data |>
  mutate(Everhadvaccinated = ifelse(Everhadvaccinated == 1, "Yes", "No")) |>
  group_by(Everhadvaccinated) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Has Your Child Ever Had a Vaccination?" = Everhadvaccinated)


################################################################################
# BCG Scar
################################################################################

BCG <- sample_data |>
  group_by(BCGscar) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does Your Child Have a BCG Scar" = BCGscar)


################################################################################
# Fever last two weeks
################################################################################

fever <- sample_data |>
  mutate(Inthelasttwoweekschildhadfever = ifelse(Inthelasttwoweekschildhadfever == 1, "Yes", "No")) |> 
  group_by(Inthelasttwoweekschildhadfever) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("In the Last Two Weeks Has Your Child Had a Fever?" = Inthelasttwoweekschildhadfever)

################################################################################
# Diarrhea last two weeks
################################################################################

diarrhea <- sample_data |>
  mutate(InthelasttwoweekschildhadDiahearria = ifelse(InthelasttwoweekschildhadDiahearria == 1, "Yes", "No")) |> 
  group_by(InthelasttwoweekschildhadDiahearria) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("In the Last Two Weeks Has Your Child Had Diahrrea?" = InthelasttwoweekschildhadDiahearria)

################################################################################
# Cough last two weeks
################################################################################

cough <- sample_data |>
  mutate(Inthelasttwoweekschildhadcough = ifelse(Inthelasttwoweekschildhadcough == 1, "Yes", "No")) |> 
  group_by(Inthelasttwoweekschildhadcough) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("In the Last Two Weeks Has Your Child Had a Cough?" = Inthelasttwoweekschildhadcough)

################################################################################
# Fingernails Trimmed
################################################################################

fingernails_trimmed <- sample_data |>
  mutate(Childsfingernailtrimmed = ifelse(Childsfingernailtrimmed == 1, "Yes", "No")) |> 
  group_by(Childsfingernailtrimmed) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Are Your Child's Fingernails Trimmed?" = Childsfingernailtrimmed)

################################################################################
# Fingernails Dirty
################################################################################

fingernails_dirty <- sample_data |>
  mutate(Arechildsfingernailsdirty = ifelse(Arechildsfingernailsdirty == 1, "Yes", "No")) |> 
  group_by(Arechildsfingernailsdirty) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Are Your Child's Fingernails Dirty?" = Arechildsfingernailsdirty)


################################################################################
# Fingernail Trimming Frequency
################################################################################

fingernails_trimmed_frequency <- sample_data |>
  mutate(Howoftendoyoutrimyourfingernails = case_when(Howoftendoyoutrimyourfingernails == 0 ~ "Less Than Once/Month", 
                                                      Howoftendoyoutrimyourfingernails == 1 ~ "Once per Week", 
                                                      Howoftendoyoutrimyourfingernails == 2 ~ "Once per two weeks", 
                                                      is.na(Howoftendoyoutrimyourfingernails) ~ "Don't Know")) |> 
  group_by(Howoftendoyoutrimyourfingernails) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("How Often do you Trim your Child's Fingernails?" = Howoftendoyoutrimyourfingernails)


################################################################################
# Presence of Latrine in School
################################################################################

latrine_in_school <- sample_data |>
  mutate(Istheretoiletintheschool = case_when(Istheretoiletintheschool == 0 ~ "No", 
                                                      Istheretoiletintheschool == 1 ~ "Yes", 
                                                      is.na(Istheretoiletintheschool) ~ "Don't Know")) |> 
  group_by(Istheretoiletintheschool) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Is There a Latrine in Your Child's School?" = Istheretoiletintheschool)


################################################################################
# Doors Around Latrine in School
################################################################################

latrine_doors <- sample_data |>
  mutate(Ifthereisatoiletdoesthelatrinehavedoors = case_when(Ifthereisatoiletdoesthelatrinehavedoors == 0 ~ "No", 
                                                      Ifthereisatoiletdoesthelatrinehavedoors == 1 ~ "Yes", 
                                                      is.na(Ifthereisatoiletdoesthelatrinehavedoors) ~ "Don't Know")) |> 
  group_by(Ifthereisatoiletdoesthelatrinehavedoors) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does the Latrine have Doors?" = Ifthereisatoiletdoesthelatrinehavedoors)


################################################################################
# Flies Around Latrine
################################################################################

latrine_flies <- sample_data |>
  mutate(Fliesobservedinaroundthelatrine = case_when(Fliesobservedinaroundthelatrine == 0 ~ "No", 
                                                             Fliesobservedinaroundthelatrine == 1 ~ "Yes", 
                                                             is.na(Fliesobservedinaroundthelatrine) ~ "Don't Know")) |> 
  group_by(Fliesobservedinaroundthelatrine) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does the Latrine have Flies in or Around?" = Fliesobservedinaroundthelatrine)


################################################################################
# Stool Around Latrine
################################################################################

latrine_stool <- sample_data |>
  mutate(Vissiblestoolobservedonlatrinefloor = case_when(Vissiblestoolobservedonlatrinefloor == 0 ~ "No", 
                                                         Vissiblestoolobservedonlatrinefloor == 1 ~ "Yes", 
                                                         is.na(Vissiblestoolobservedonlatrinefloor) ~ "Don't Know")) |> 
  group_by(Vissiblestoolobservedonlatrinefloor) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does the Latrine have Visible Stool in or Around?" = Vissiblestoolobservedonlatrinefloor)

################################################################################
# Do you Know how Intestinal Worms are Transmitted?
################################################################################

worms_transmitted <- sample_data |>
  mutate(Doyouknowhowintestinalwomstransmitted = case_when(Doyouknowhowintestinalwomstransmitted == 0 ~ "No", 
                                                         Doyouknowhowintestinalwomstransmitted == 1 ~ "Yes", 
                                                         Doyouknowhowintestinalwomstransmitted == 99 ~ "Don't Know")) |> 
  group_by(Doyouknowhowintestinalwomstransmitted) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do you Know How Intestinal Worms are Transmitted?" = Doyouknowhowintestinalwomstransmitted)


################################################################################
# Do you Know why Worms are bad for you health?
################################################################################

worms_health <- sample_data |>
  mutate(Doyouknowwhywormsarebadforyourhealth = case_when(Doyouknowwhywormsarebadforyourhealth == 0 ~ "No", 
                                                           Doyouknowwhywormsarebadforyourhealth == 1 ~ "Yes", 
                                                           Doyouknowwhywormsarebadforyourhealth == 99 ~ "Don't Know")) |> 
  group_by(Doyouknowwhywormsarebadforyourhealth) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do You Know Why Worms are Bad for your Health?" = Doyouknowwhywormsarebadforyourhealth)


################################################################################
# Do you know how you can avoid getting these worms?
################################################################################

worms_avoid <- sample_data |>
  mutate(Doyouknowhowyoucanavoidgettingtheseworms = case_when(Doyouknowhowyoucanavoidgettingtheseworms == 0 ~ "No", 
                                                          Doyouknowhowyoucanavoidgettingtheseworms == 1 ~ "Yes", 
                                                          Doyouknowhowyoucanavoidgettingtheseworms == 99 ~ "Don't Know")) |> 
  group_by(Doyouknowhowyoucanavoidgettingtheseworms) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do You Know How you can Avoid Getting These Worms?" = Doyouknowhowyoucanavoidgettingtheseworms)


################################################################################
# Living Address
################################################################################

living_address <- sample_data |>
  mutate(Yourlivingaddress = ifelse(Yourlivingaddress != "Amenu" & Yourlivingaddress != "Bosa Kito" & Yourlivingaddress != "Hermata mentina" & Yourlivingaddress != "Jiren" & Yourlivingaddress != "Mentina" & Yourlivingaddress != "Seto semero" & Yourlivingaddress != "Welda", "Other", Yourlivingaddress)) |>
  group_by(Yourlivingaddress) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Your Living Address" = Yourlivingaddress) |>
  arrange(n)


################################################################################
# Family Occupation
################################################################################

occupation <- sample_data |>
  group_by(Familyoccupation) |>
  mutate(Familyoccupation = ifelse(Familyoccupation != "Police man" & Familyoccupation != "Teacher" & Familyoccupation != "merchant" & Familyoccupation != "Carprnter" & Familyoccupation != "Dailly labor" & Familyoccupation != "Security" & Familyoccupation != "Farmer" & Familyoccupation != "Government employee" & Familyoccupation != "Merchant", "Other", Familyoccupation)) |>
  mutate(Familyoccupation = ifelse(Familyoccupation == "merchant", "Merchant", Familyoccupation), 
         Familyoccupation = ifelse(is.na(Familyoccupation), "Other", Familyoccupation)) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Family Occupation" = Familyoccupation) |>
  arrange(n)

################################################################################
# House Floor Material
################################################################################

floor <- sample_data |>
  group_by(Whatmaterialyourhousefloormadefrom) |>
  mutate(Whatmaterialyourhousefloormadefrom = case_when(Whatmaterialyourhousefloormadefrom == 0 ~ "Dust", 
                                                        Whatmaterialyourhousefloormadefrom == 1 ~ "Cement", 
                                                        Whatmaterialyourhousefloormadefrom == 2 ~ "Plastic Covered", 
                                                        is.na(Whatmaterialyourhousefloormadefrom) ~ "Don't Know")) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Floor Material" = Whatmaterialyourhousefloormadefrom) |>
  arrange(n)

################################################################################
# Is your Kitchen Separated?
################################################################################

kitchen_separation <- sample_data |>
  group_by(Isyourkitcheninyourhouseorsepatated) |>
  mutate(Isyourkitcheninyourhouseorsepatated = case_when(Isyourkitcheninyourhouseorsepatated == 0 ~ "Within House", 
                                                        Isyourkitcheninyourhouseorsepatated == 2 ~ "Separated", 
                                                        Isyourkitcheninyourhouseorsepatated == 99 | is.na(Isyourkitcheninyourhouseorsepatated) ~ "Don't Know")) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Is Your Kitchen Separated or in Your House?" = Isyourkitcheninyourhouseorsepatated) |>
  arrange(n)

################################################################################
# Kitchen Material if Separated
################################################################################

kitchen_material <- sample_data |>
  group_by(Ifseparetedwhatmaterialsyourkitchenmadefrom) |>
  mutate(Ifseparetedwhatmaterialsyourkitchenmadefrom = case_when(Ifseparetedwhatmaterialsyourkitchenmadefrom == 0 ~ "Dust", 
                                                        Ifseparetedwhatmaterialsyourkitchenmadefrom == 1 ~ "Cement", 
                                                        Ifseparetedwhatmaterialsyourkitchenmadefrom == 2 ~ "Plastic Covered", 
                                                        is.na(Ifseparetedwhatmaterialsyourkitchenmadefrom) | Ifseparetedwhatmaterialsyourkitchenmadefrom == 9 ~ "Don't Know")) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Kitchen Floor Material if Separated" = Ifseparetedwhatmaterialsyourkitchenmadefrom) |>
  arrange(n)

################################################################################
# Kitchen Structures
################################################################################

kitchen_structures <- sample_data |>
  mutate(Roof = ifelse(Roof == 1, 1, 0), 
         Wall = ifelse(Wall == 2, 1, 0), 
         None = ifelse(None == 0, 1, 0)) |>
  mutate(Roof = ifelse(is.na(Roof), 0, Roof), 
         Wall = ifelse(is.na(Wall), 0, Wall), 
         None = ifelse(is.na(None), 0, None)) |>
  summarise(structure = c("Roof", "Wall", "None"), 
            frequency = c(round(mean(Roof), 3), round(mean(Wall), 3), round(mean(None), 3))) |>
  mutate(count = round(frequency*138))


################################################################################
# Cooking Fuel
################################################################################

cooking_fuel <- sample_data |>
  mutate(Wood = ifelse(Wood == 0, 1, 0), 
         Gas = ifelse(Gas == 1, 1, 0), 
         Coal = ifelse(Coal == 2, 1, 0), 
         Kerosisn = ifelse(Kerosisn == 3, 1, 0), 
         Electric = ifelse(Electric == 4, 1, 0)) |>
  mutate(Wood = ifelse(is.na(Wood), 0, Wood), 
         Gas = ifelse(is.na(Gas), 0, Gas), 
         Kerosisn = ifelse(is.na(Kerosisn), 0, Kerosisn), 
         Coal = ifelse(is.na(Coal), 0, Coal), 
         Electric = ifelse(is.na(Electric), 0, Electric)) |>
  summarise(Fuel = c("Wood", "Gas", "Coal", "Kerosene", "Electric"), 
            frequency = c(round(mean(Wood), 3), round(mean(Gas),3),round(mean(Coal), 3), round(mean(Kerosisn), 3), round(mean(Electric), 3))) |>
  mutate(count = round(frequency*138)) |>
  arrange(frequency)


################################################################################
# Electricity in House
################################################################################

electricity <- sample_data |>
  mutate(Doyouhaveelectricityinyourhaouse = case_when(Doyouhaveelectricityinyourhaouse == 0 ~ "No", 
                                                     Doyouhaveelectricityinyourhaouse == 1 ~ "Yes", 
                                                     is.na(Doyouhaveelectricityinyourhaouse) ~ "Don't Know")) |> 
  group_by(Doyouhaveelectricityinyourhaouse) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do you Have Electricity in your House?" = Doyouhaveelectricityinyourhaouse)

################################################################################
# Radio in House
################################################################################

radio <- sample_data |>
  mutate(Doesyourfamillyownradio = case_when(Doesyourfamillyownradio == 0 ~ "No", 
                                                      Doesyourfamillyownradio == 1 ~ "Yes", 
                                                      is.na(Doesyourfamillyownradio) ~ "Don't Know")) |> 
  group_by(Doesyourfamillyownradio) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do you Have a Radio in your House?" = Doesyourfamillyownradio)

################################################################################
# Television in House
################################################################################

tv <- sample_data |>
  mutate(Doesyourfamillyowntelevision = case_when(Doesyourfamillyowntelevision == 0 ~ "No", 
                                             Doesyourfamillyowntelevision == 1 ~ "Yes", 
                                             is.na(Doesyourfamillyowntelevision) ~ "Don't Know")) |> 
  group_by(Doesyourfamillyowntelevision) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do you Have a TV in your House?" = Doesyourfamillyowntelevision)

################################################################################
# Family Member own a phone
################################################################################

phone <- sample_data |>
  mutate(Doesyourfamillymemberownaphone = case_when(Doesyourfamillymemberownaphone == 0 ~ "No", 
                                                  Doesyourfamillymemberownaphone == 1 ~ "Yes", 
                                                  is.na(Doesyourfamillymemberownaphone) ~ "Don't Know")) |> 
  group_by(Doesyourfamillymemberownaphone) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does Your Family Member Have a Phone?" = Doesyourfamillymemberownaphone)

################################################################################
# Use of Phone
################################################################################

phone_use <- sample_data |>
  mutate(Forwhatpurposeyourfamilyuseaphone = case_when(Forwhatpurposeyourfamilyuseaphone == 0 ~ "Call Only", 
                                                    Forwhatpurposeyourfamilyuseaphone == 1 ~ "Call/Radio/Internet", 
                                                    is.na(Forwhatpurposeyourfamilyuseaphone) ~ "Don't Know")) |> 
  group_by(Forwhatpurposeyourfamilyuseaphone) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("For What Purpose Does your Family Use the Phone?" = Forwhatpurposeyourfamilyuseaphone) |>
  arrange(frequency)

################################################################################
# Animal in House/Compound
################################################################################

animals <- sample_data |>
  mutate(Cattle = ifelse(Cattle == 0, 1, 0), 
         SheepGoat = ifelse(SheepGoat == 1, 1, 0), 
         Chicken = ifelse(Chicken == 2, 1, 0), 
         Pet = ifelse(Pet == 3, 1, 0), 
         Nodomanimal = ifelse(Nodomanimal == 4, 1, 0)) |>
  mutate(Cattle = ifelse(is.na(Cattle), 0, Cattle), 
         SheepGoat = ifelse(is.na(SheepGoat), 0, SheepGoat), 
         Chicken = ifelse(is.na(Chicken), 0, Chicken), 
         Pet = ifelse(is.na(Pet), 0, Pet), 
         Nodomanimal = ifelse(is.na(Nodomanimal), 0, Nodomanimal)) |>
  summarise(Animal = c("Cattle", "Sheep/Goats", "Chickens", "Pets", "None"), 
            frequency = c(round(mean(Cattle), 3), round(mean(SheepGoat), 3), round(mean(Chicken), 3), round(mean(Pet), 3), round(mean(Nodomanimal), 3))) |>
  mutate(count = round(frequency*138)) |>
  arrange(count)

################################################################################
# Portable Water
################################################################################

water <- sample_data |>
  mutate(Doyouhavepotablewaterinyourhouse = case_when(Doyouhavepotablewaterinyourhouse == 0 ~ "No", 
                                                       Doyouhavepotablewaterinyourhouse == 1 ~ "Yes", 
                                                       is.na(Doyouhavepotablewaterinyourhouse) ~ "Don't Know")) |> 
  group_by(Doyouhavepotablewaterinyourhouse) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do You Have Portable Water in Your House?" = Doyouhavepotablewaterinyourhouse)


################################################################################
# If not how do you get your water
################################################################################

other_water <- sample_data |>
  mutate(Ifnotwheredoyougetyordrinkingwaterfrom = case_when(Ifnotwheredoyougetyordrinkingwaterfrom == 0 ~ "Neighbor", 
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 1 ~ "River", 
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 2 ~ "Well", 
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 3 ~ "Truck",
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 4 ~ "Tank", 
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 5 ~ "Public Fountain",
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == "Tap water" ~ "Tap Water", 
                                                      Ifnotwheredoyougetyordrinkingwaterfrom == 10 ~ NA)) |> 
  group_by(Ifnotwheredoyougetyordrinkingwaterfrom) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/70,3)) |>
  drop_na() |>
  rename("If not, where do you get your water?" = Ifnotwheredoyougetyordrinkingwaterfrom) |>
  arrange(frequency)


################################################################################
# Do you treat your water?
################################################################################

treat_water <- sample_data |>
  mutate(Doyoudrinkdirectlyordoyoutreat = case_when(Doyoudrinkdirectlyordoyoutreat == 1 ~ "Directly", 
                                                    Doyoudrinkdirectlyordoyoutreat == 0 ~ "Treat", 
                                                    is.na(Doyoudrinkdirectlyordoyoutreat) ~ "Don't Know")) |>
  group_by(Doyoudrinkdirectlyordoyoutreat) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Do you treat your water or drink directly?" = Doyoudrinkdirectlyordoyoutreat) |>
  arrange(frequency)

################################################################################
# How do you treat your water?
################################################################################

treat_water_how <- sample_data |>
  mutate(Ifyoutreathow = case_when(Ifyoutreathow == 0 ~ "Boil", 
                                   Ifyoutreathow == 1 ~ "Chemical", 
                                   Ifyoutreathow == 2 ~ "Filter", 
                                   Ifyoutreathow == 99 ~ "Don't Know")) |>
  group_by(Ifyoutreathow) |>
  summarise(n = n()) |>
  drop_na() |>
  mutate(frequency = round(n/22,3)) |>
  rename("How do you treat your water?" = Ifyoutreathow) |>
  arrange(frequency)

################################################################################
# Does your family own a latrine?
################################################################################

own_latrine <- sample_data |>
  mutate(Doyourfamilyownlatrine = ifelse(Doyourfamilyownlatrine == 1, "Yes", "No")) |>
  group_by(Doyourfamilyownlatrine) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Does your family own a latrine?" = Doyourfamilyownlatrine)

################################################################################
# Is your latrine inside or outside?
################################################################################

latrine_location <- sample_data |>
  mutate(Isyourlatrineinsideoroutside = ifelse(Isyourlatrineinsideoroutside == 1, "Inside", "Outside")) |>
  group_by(Isyourlatrineinsideoroutside) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Is your latrine inside or outside?" = Isyourlatrineinsideoroutside)

################################################################################
# How far from your house?
################################################################################

library(stringr)

latrine_distance_house <- sample_data |>
  mutate(Ifoutsidedistancefromyourhouse = as.numeric(str_sub(Ifoutsidedistancefromyourhouse, 1, -2)), 
         Ifoutsidedistancefromyourhouse = case_when(Ifoutsidedistancefromyourhouse >= 0 & Ifoutsidedistancefromyourhouse < 5 ~ "0-5", 
                                                    Ifoutsidedistancefromyourhouse >= 5 & Ifoutsidedistancefromyourhouse < 10 ~ "5-10", 
                                                    Ifoutsidedistancefromyourhouse >= 10 & Ifoutsidedistancefromyourhouse <= 20 ~ "10-20", 
                                                    Ifoutsidedistancefromyourhouse > 20 ~ ">20"), 
         Ifoutsidedistancefromyourhouse = factor(Ifoutsidedistancefromyourhouse, levels = c("0-5", "5-10", "10-20", ">20"))) |>
  group_by(Ifoutsidedistancefromyourhouse) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Distance from House to Latrine(m)?" = Ifoutsidedistancefromyourhouse)

################################################################################
# How far from your kitchen?
################################################################################

latrine_distance_kitchen <- sample_data |>
  mutate(Distancebetweenlatrineandyourkitchen = as.numeric(str_sub(Distancebetweenlatrineandyourkitchen, 1, -2)), 
         Distancebetweenlatrineandyourkitchen = case_when(Distancebetweenlatrineandyourkitchen >= 0 & Distancebetweenlatrineandyourkitchen < 5 ~ "0-5", 
                                                    Distancebetweenlatrineandyourkitchen >= 5 & Distancebetweenlatrineandyourkitchen < 10 ~ "5-10", 
                                                    Distancebetweenlatrineandyourkitchen >= 10 & Distancebetweenlatrineandyourkitchen <= 20 ~ "10-20", 
                                                    Distancebetweenlatrineandyourkitchen > 20 ~ ">20"), 
         Distancebetweenlatrineandyourkitchen = factor(Distancebetweenlatrineandyourkitchen, levels = c("0-5", "5-10", "10-20", ">20"))) |>
  group_by(Distancebetweenlatrineandyourkitchen) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138,3)) |>
  rename("Distance from Kitchen to Latrine(m)?" = Distancebetweenlatrineandyourkitchen)


################################################################################
# What is your latrine connected to?
################################################################################

latrine_connection <- sample_data |>
  mutate(Yourlatrineconnectedto = case_when(Yourlatrineconnectedto == 0 ~ "Sewage", 
                                            Yourlatrineconnectedto == 1 ~ "Ditch", 
                                            Yourlatrineconnectedto == 2 ~ "River", 
                                            Yourlatrineconnectedto == 3 ~ "Well", 
                                            Yourlatrineconnectedto == 99 | is.na(Yourlatrineconnectedto) ~ "NA/Don't Know")) |>
  group_by(Yourlatrineconnectedto) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("What is your latrine connected to?" = Yourlatrineconnectedto) |>
  arrange(frequency)

################################################################################
# How often do you bathe in river?
################################################################################

bathe_river <- sample_data |>
  mutate(Howoftendoyoubathinriver = case_when(Howoftendoyoubathinriver == 0 ~ "Always", 
                                            Howoftendoyoubathinriver == 1 ~ "Sometimes", 
                                            Howoftendoyoubathinriver == 2 ~ "Never",
                                            Howoftendoyoubathinriver == 99 | is.na(Howoftendoyoubathinriver) ~ "NA/Don't Know"), 
         Howoftendoyoubathinriver = factor(Howoftendoyoubathinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Howoftendoyoubathinriver) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How often do you bathe in the river?" = Howoftendoyoubathinriver)


################################################################################
# How often do you wash your clothes in river?
################################################################################

clothes_river <- sample_data |>
  mutate(Howoftendoyouwashyourclothesinriver = case_when(Howoftendoyouwashyourclothesinriver == 0 ~ "Always", 
                                              Howoftendoyouwashyourclothesinriver == 1 ~ "Sometimes", 
                                              Howoftendoyouwashyourclothesinriver == 2 ~ "Never",
                                              Howoftendoyouwashyourclothesinriver == 99 | is.na(Howoftendoyouwashyourclothesinriver) ~ "NA/Don't Know"), 
         Howoftendoyouwashyourclothesinriver = factor(Howoftendoyouwashyourclothesinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Howoftendoyouwashyourclothesinriver) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How often do you wash clothes in the river?" = Howoftendoyouwashyourclothesinriver)

################################################################################
# How often do you defecate in a open field?
################################################################################

defecate_field <- sample_data |>
  mutate(Doyoudeficateintheopenfield = case_when(Doyoudeficateintheopenfield == 0 ~ "Always", 
                                                         Doyoudeficateintheopenfield == 1 ~ "Sometimes", 
                                                         Doyoudeficateintheopenfield == 2 ~ "Never",
                                                         Doyoudeficateintheopenfield == 99 | is.na(Doyoudeficateintheopenfield) ~ "NA/Don't Know"), 
         Doyoudeficateintheopenfield = factor(Doyoudeficateintheopenfield, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyoudeficateintheopenfield) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How often do you defecate in an open field?" = Doyoudeficateintheopenfield)


################################################################################
# do you use the school latrine?
################################################################################

latrine_school <- sample_data |>
  mutate(Doyouuseschoollatrine = case_when(Doyouuseschoollatrine == 0 ~ "Always", 
                                                 Doyouuseschoollatrine == 1 ~ "Sometimes", 
                                                 Doyouuseschoollatrine == 2 ~ "Never",
                                                 Doyouuseschoollatrine == 99 | is.na(Doyouuseschoollatrine) ~ "NA/Don't Know"), 
         Doyouuseschoollatrine = factor(Doyouuseschoollatrine, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouuseschoollatrine) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you use the school latrine?" = Doyouuseschoollatrine)

################################################################################
# do you use toilet paper after defecating?
################################################################################

toilet_paper <- sample_data |>
  mutate(Doyouusetoiletpapertowipeyourbumafterdefication = case_when(Doyouusetoiletpapertowipeyourbumafterdefication == 0 ~ "Always", 
                                           Doyouusetoiletpapertowipeyourbumafterdefication == 1 ~ "Sometimes", 
                                           Doyouusetoiletpapertowipeyourbumafterdefication == 2 ~ "Never",
                                           Doyouusetoiletpapertowipeyourbumafterdefication == 99 | is.na(Doyouusetoiletpapertowipeyourbumafterdefication) ~ "NA/Don't Know"), 
         Doyouusetoiletpapertowipeyourbumafterdefication = factor(Doyouusetoiletpapertowipeyourbumafterdefication, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouusetoiletpapertowipeyourbumafterdefication) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you use toilet paper after defecation?" = Doyouusetoiletpapertowipeyourbumafterdefication)

################################################################################
# Washing hands - toilet
################################################################################

wash_hands_after_toilet <- sample_data |>
  mutate(Doyouwashyourhandsaftertoilet = case_when(Doyouwashyourhandsaftertoilet == 0 ~ "Always", 
                                                                     Doyouwashyourhandsaftertoilet == 1 ~ "Sometimes", 
                                                                     Doyouwashyourhandsaftertoilet == 2 ~ "Never",
                                                                     Doyouwashyourhandsaftertoilet == 99 | is.na(Doyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
         Doyouwashyourhandsaftertoilet = factor(Doyouwashyourhandsaftertoilet, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouwashyourhandsaftertoilet) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you wash your hands after using the toilet?" = Doyouwashyourhandsaftertoilet)

wash_hands_how_after_toilet <- sample_data |>
  mutate(Howdoyouwashyourhandsaftertoilet = case_when(Howdoyouwashyourhandsaftertoilet == 0 ~ "Water", 
                                                   Howdoyouwashyourhandsaftertoilet == 1 ~ "Soap and Water",
                                                   Howdoyouwashyourhandsaftertoilet == 99 | is.na(Howdoyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
         Howdoyouwashyourhandsaftertoilet = factor(Howdoyouwashyourhandsaftertoilet, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  group_by(Howdoyouwashyourhandsaftertoilet) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How do you wash your hands after using the toilet?" = Howdoyouwashyourhandsaftertoilet)


if_with_soap_after_toilet <- sample_data |>
  mutate(Ifwithsoaphowoftenyouuseit = case_when(Ifwithsoaphowoftenyouuseit == 0 ~ "Always", 
                                                Ifwithsoaphowoftenyouuseit == 1 ~ "Sometimes",
                                                Ifwithsoaphowoftenyouuseit == 2 ~ "Never",
                                                Ifwithsoaphowoftenyouuseit == 99 | is.na(Ifwithsoaphowoftenyouuseit) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit = factor(Ifwithsoaphowoftenyouuseit, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  group_by(Ifwithsoaphowoftenyouuseit) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("If with soap after toilet, how often do you use it?" = Ifwithsoaphowoftenyouuseit)

################################################################################
# Washing hands - Eating
################################################################################

wash_hands_after_eating <- sample_data |>
  mutate(Doyouwashyourhandsbeforeeating = case_when(Doyouwashyourhandsbeforeeating == 0 ~ "Always", 
                                                   Doyouwashyourhandsbeforeeating == 1 ~ "Sometimes", 
                                                   Doyouwashyourhandsbeforeeating == 2 ~ "Never",
                                                   Doyouwashyourhandsbeforeeating == 99 | is.na(Doyouwashyourhandsbeforeeating) ~ "NA/Don't Know"), 
         Doyouwashyourhandsbeforeeating = factor(Doyouwashyourhandsbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouwashyourhandsbeforeeating) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you wash your hands before eating?" = Doyouwashyourhandsbeforeeating)

wash_hands_how_after_eating <- sample_data |>
  mutate(Howdoyouwashyourhandsbeforeeating = case_when(Howdoyouwashyourhandsbeforeeating == 0 ~ "Water", 
                                                      Howdoyouwashyourhandsbeforeeating == 1 ~ "Soap and Water",
                                                      Howdoyouwashyourhandsbeforeeating == 99 | is.na(Howdoyouwashyourhandsbeforeeating) | Howdoyouwashyourhandsbeforeeating == 2 ~ "NA/Don't Know"), 
         Howdoyouwashyourhandsbeforeeating = factor(Howdoyouwashyourhandsbeforeeating, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  group_by(Howdoyouwashyourhandsbeforeeating) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How do you wash your hands before eating?" = Howdoyouwashyourhandsbeforeeating)


if_with_soap_after_eating <- sample_data |>
  mutate(Ifwithsoaphowoftenyouuseit_A = case_when(Ifwithsoaphowoftenyouuseit_A == 0 ~ "Always", 
                                                Ifwithsoaphowoftenyouuseit_A == 1 ~ "Sometimes",
                                                Ifwithsoaphowoftenyouuseit_A == 2 ~ "Never",
                                                Ifwithsoaphowoftenyouuseit_A == 99 | is.na(Ifwithsoaphowoftenyouuseit_A) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit_A = factor(Ifwithsoaphowoftenyouuseit_A, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  group_by(Ifwithsoaphowoftenyouuseit_A) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("If with soap before eating, how often do you use it?" = Ifwithsoaphowoftenyouuseit_A)


################################################################################
# Do you eat soil ?
################################################################################

eat_soil <- sample_data |>
  mutate(Doyoueatsoil = case_when(Doyoueatsoil == 0 ~ "Always", 
                                                    Doyoueatsoil == 1 ~ "Sometimes", 
                                                    Doyoueatsoil == 2 ~ "Never",
                                                    Doyoueatsoil == 99 | is.na(Doyoueatsoil) ~ "NA/Don't Know"), 
         Doyoueatsoil = factor(Doyoueatsoil, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyoueatsoil) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you eat soil?" = Doyoueatsoil)

################################################################################
# Fruit
################################################################################

favorite_fruit = sample_data |>
  mutate(Whatisyourfavoritefruitethatyoueat = ifelse(Whatisyourfavoritefruitethatyoueat != "Mango" & Whatisyourfavoritefruitethatyoueat != "Banana" & Whatisyourfavoritefruitethatyoueat != "Avocado" & Whatisyourfavoritefruitethatyoueat != "Apple" & Whatisyourfavoritefruitethatyoueat != "Pinapple", "Other", Whatisyourfavoritefruitethatyoueat), 
         Whatisyourfavoritefruitethatyoueat = ifelse(is.na(whatfooddoyoueatmost), "Other", Whatisyourfavoritefruitethatyoueat)) |>
  group_by(Whatisyourfavoritefruitethatyoueat) |>
  summarise(n = n()) |>
  rename("What is your favorite fruit?" = Whatisyourfavoritefruitethatyoueat) |>
  mutate(frequency = round(n/138, 3)) |>
  arrange(-n)
  
wash_fruit <- sample_data |>
  mutate(Doyouwashyourfruitesbeforeeating = case_when(Doyouwashyourfruitesbeforeeating == 0 ~ "Always", 
                                                    Doyouwashyourfruitesbeforeeating == 1 ~ "Sometimes", 
                                                    Doyouwashyourfruitesbeforeeating == 2 ~ "Never",
                                                    Doyouwashyourfruitesbeforeeating == 99 | is.na(Doyouwashyourfruitesbeforeeating) ~ "NA/Don't Know"), 
         Doyouwashyourfruitesbeforeeating = factor(Doyouwashyourfruitesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouwashyourfruitesbeforeeating) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you wash your fruit before eating?" = Doyouwashyourfruitesbeforeeating)


################################################################################
# Vegetables
################################################################################

eat_raw_vegetables <- sample_data |>
  mutate(Doyoueatraworundercookedvegitables = case_when(Doyoueatraworundercookedvegitables == 0 ~ "Always", 
                                                      Doyoueatraworundercookedvegitables == 1 ~ "Sometimes", 
                                                      Doyoueatraworundercookedvegitables == 2 ~ "Never",
                                                      Doyoueatraworundercookedvegitables == 99 | is.na(Doyoueatraworundercookedvegitables) ~ "NA/Don't Know"), 
         Doyoueatraworundercookedvegitables = factor(Doyoueatraworundercookedvegitables, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyoueatraworundercookedvegitables) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you eat raw or undercooked vegetables?" = Doyoueatraworundercookedvegitables)

wash_vegetables <- sample_data |>
  mutate(Ifyeshowoftenyouwashyourvegitablesbeforeeating = case_when(Ifyeshowoftenyouwashyourvegitablesbeforeeating == 0 ~ "Always", 
                                                      Ifyeshowoftenyouwashyourvegitablesbeforeeating == 1 ~ "Sometimes", 
                                                      Ifyeshowoftenyouwashyourvegitablesbeforeeating == 2 ~ "Never",
                                                      Ifyeshowoftenyouwashyourvegitablesbeforeeating == 99 | is.na(Ifyeshowoftenyouwashyourvegitablesbeforeeating) ~ "NA/Don't Know"), 
         Ifyeshowoftenyouwashyourvegitablesbeforeeating = factor(Ifyeshowoftenyouwashyourvegitablesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Ifyeshowoftenyouwashyourvegitablesbeforeeating) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How often do you wash your vegetables?" = Ifyeshowoftenyouwashyourvegitablesbeforeeating)


################################################################################
# Footwear
################################################################################

barefoot <- sample_data |>
  mutate(Doyouwalkbarefoot = case_when(Doyouwalkbarefoot == 0 ~ "Always", 
                                                        Doyouwalkbarefoot == 1 ~ "Sometimes", 
                                                        Doyouwalkbarefoot == 2 ~ "Never",
                                                        Doyouwalkbarefoot == 99 | is.na(Doyouwalkbarefoot) ~ "NA/Don't Know"), 
         Doyouwalkbarefoot = factor(Doyouwalkbarefoot, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  group_by(Doyouwalkbarefoot) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Do you walk barefoot?" = Doyouwalkbarefoot)

home_footwear <- sample_data |>
  mutate(whenyouareathomedoyoupreferetousesandalsorshoe = case_when(whenyouareathomedoyoupreferetousesandalsorshoe == 0 ~ "Barefoot", 
                                       whenyouareathomedoyoupreferetousesandalsorshoe == 1 ~ "Sandals", 
                                       whenyouareathomedoyoupreferetousesandalsorshoe == 2 ~ "Shoes",
                                       whenyouareathomedoyoupreferetousesandalsorshoe == 99 | is.na(whenyouareathomedoyoupreferetousesandalsorshoe) ~ "NA/Don't Know"), 
         whenyouareathomedoyoupreferetousesandalsorshoe = factor(whenyouareathomedoyoupreferetousesandalsorshoe, levels = c("Barefoot", "Sandals", "Shoes", "NA/Don't Know"))) |>
  group_by(whenyouareathomedoyoupreferetousesandalsorshoe) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("When you are at home what do you prefer to wear?" = whenyouareathomedoyoupreferetousesandalsorshoe)


#################################################################################
# Deworming pill / antiobiotics / Drugs / Malaria
#################################################################################

deworming_pill <- sample_data |>
  mutate(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill = case_when(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 0 ~ "No", 
                                       Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 1 ~ "Yes",
                                       Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 99 | is.na(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill) ~ "NA/Don't Know")) |>
  group_by(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Did your parents or health professional give you a deworming pill?" = Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill)

deworming_pill_last_time <- sample_data |>
  mutate(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill = case_when(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "1Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "2Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "3month" ~ "1-3 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "5Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "6Month" ~ "4-6 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "7Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "9Month" ~ "7-9 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "11Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "12month" ~ ">9 Months", 
                                                                         is.na(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill) ~ NA)) |>
  group_by(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill) |>
  summarise(n = n()) |>
  drop_na()|>
  mutate(frequency = n/89)

antiobiotics <- sample_data |>
  mutate(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics = case_when(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 0 ~ "No", 
                                                                                      Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 1 ~ "Yes",
                                                                                      Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 99 | is.na(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics) ~ "NA/Don't Know")) |>
  group_by(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Did your parents or health professional give you antibiotics?" = Didyourparentsorhealthprofessionalsgaveyouotherantibiotics)

drugs_illness <- sample_data |>
  mutate(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently = case_when(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 0 ~ "No", 
                                                                                Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 1 ~ "Yes",
                                                                                Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 99 | is.na(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently) ~ "NA/Don't Know")) |>
  group_by(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Are you prescribed a drug for any illness currently?" = Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently)

drug_type <- sample_data |>
  group_by(Ifyesnameandtypeofthedrugs) |>
  summarise(n = n()) |>
  rename("What type of drug?" = Ifyesnameandtypeofthedrugs)

malaria <- sample_data |>
  mutate(Inthepastthreemonthshaveyoutakenanyantimalarialdrug = case_when(Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 0 ~ "No", 
                                                                                Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 1 ~ "Yes",
                                                                                Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 99 | is.na(Inthepastthreemonthshaveyoutakenanyantimalarialdrug) ~ "NA/Don't Know")) |>
  group_by(Inthepastthreemonthshaveyoutakenanyantimalarialdrug) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("In the past three months have you taken any anti-malaria drug?" = Inthepastthreemonthshaveyoutakenanyantimalarialdrug)

################################################################################
# Child Cough / Asthma / Skin
################################################################################

wheezing <- sample_data |>
  mutate(Haveachildwheezingorwhistlinginchest = case_when(Haveachildwheezingorwhistlinginchest == 2 ~ "No", 
                                                                         Haveachildwheezingorwhistlinginchest == 1 ~ "Yes",
                                                                         Haveachildwheezingorwhistlinginchest == 99 | is.na(Haveachildwheezingorwhistlinginchest) ~ "NA/Don't Know")) |>
  group_by(Haveachildwheezingorwhistlinginchest) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Has your child had wheezing or whistling in the chest?" = Haveachildwheezingorwhistlinginchest)

wheezing_attack <- sample_data |>
  mutate(Howmanytimesinthelastyearthechildhadwheeziling = case_when(Howmanytimesinthelastyearthechildhadwheeziling == 2 ~ "1-3", 
                                                          Howmanytimesinthelastyearthechildhadwheeziling == 1 ~ "0",
                                                          Howmanytimesinthelastyearthechildhadwheeziling == 99 | is.na(Howmanytimesinthelastyearthechildhadwheeziling) ~ "NA/Don't Know")) |>
  group_by(Howmanytimesinthelastyearthechildhadwheeziling) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("How many times in the last year has your child had a wheezing attack?" = Howmanytimesinthelastyearthechildhadwheeziling)

asthma <- sample_data |>
  mutate(Hasachildeverhadasthma = case_when(Hasachildeverhadasthma == 2 ~ "No", 
                                                          Hasachildeverhadasthma == 1 ~ "Yes",
                                                          Hasachildeverhadasthma == 99 | is.na(Hasachildeverhadasthma) ~ "NA/Don't Know")) |>
  group_by(Hasachildeverhadasthma) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Has your child ever had asthma?" = Hasachildeverhadasthma)

itchy_skin <- sample_data |>
  mutate(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases = case_when(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 2 ~ "No", 
                                                          Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 1 ~ "Yes",
                                                          Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 99 | is.na(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases) ~ "NA/Don't Know")) |>
  group_by(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Has your child ever had itchy skin or rash which affected the skin creases?" = Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases)

itchy_skin_area <- sample_data |>
  mutate(Ifyeshasthisrashaffectedtheelbowfolds = ifelse(is.na(Ifyeshasthisrashaffectedtheelbowfolds), 0, Ifyeshasthisrashaffectedtheelbowfolds), 
         Ifyeshasthisrashaffectedbehindtheknees = ifelse(is.na(Ifyeshasthisrashaffectedbehindtheknees), 0, Ifyeshasthisrashaffectedbehindtheknees), 
         Ifyeshasthisrashaffectedinfrontoftheankls = ifelse(is.na(Ifyeshasthisrashaffectedinfrontoftheankls), 0, Ifyeshasthisrashaffectedinfrontoftheankls), 
         Ifyeshasthisrashaffectedunderthebuttucks = ifelse(is.na(Ifyeshasthisrashaffectedunderthebuttucks), 0, Ifyeshasthisrashaffectedunderthebuttucks), 
         Ifyeshasthisrashaffectedaroundtheneck = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheneck), 0, Ifyeshasthisrashaffectedaroundtheneck), 
         Ifyeshasthisrashaffectedaroundtheeyesears = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheeyesears), 0, Ifyeshasthisrashaffectedaroundtheeyesears)) |>
  summarise(`Body Area Affected` = c("Elbow folds", "Behind the knees", "Front of the ankles", "Under the buttocks", "Around the neck", "Around the eyes/ears"), 
            frequency = c(round(mean(Ifyeshasthisrashaffectedtheelbowfolds), 3), round(mean(Ifyeshasthisrashaffectedbehindtheknees), 3), round(mean(Ifyeshasthisrashaffectedinfrontoftheankls), 3), round(mean(Ifyeshasthisrashaffectedunderthebuttucks), 3), round(mean(Ifyeshasthisrashaffectedaroundtheneck), 3), round(mean(Ifyeshasthisrashaffectedaroundtheeyesears), 3))) |>
  mutate(n = round(frequency * 138))

hay_fever_sneezing <- sample_data |>
  mutate(Doesthechildhadhayfeverorpersistantsneezing = case_when(Doesthechildhadhayfeverorpersistantsneezing == 2 ~ "No", 
                                                                                    Doesthechildhadhayfeverorpersistantsneezing == 1 ~ "Yes",
                                                                                    Doesthechildhadhayfeverorpersistantsneezing == 99 | is.na(Doesthechildhadhayfeverorpersistantsneezing) ~ "NA/Don't Know")) |>
  group_by(Doesthechildhadhayfeverorpersistantsneezing) |>
  summarise(n = n()) |>
  mutate(frequency = round(n/138, 3)) |>
  rename("Has your child had hay fever or persistent sneezing?" = Doesthechildhadhayfeverorpersistantsneezing)


##########################################################################################
# Big data table
##########################################################################################

cleaned <- sample_data |>
  mutate(Age = case_when(Age >= 6 & Age <= 12 ~ "6-12", 
                         Age >= 13 & Age <= 17 ~ "13-17")) |>
  mutate(Age = factor(Age, levels = c("6-12", "13-17"))) |>
  mutate(weight = case_when(weight < 35 ~ "15-35", 
                            weight >= 35 & weight < 55 ~ "35-55", 
                            weight >= 55 ~ ">55")) |>
  mutate(weight = factor(weight, levels = c("15-35","35-55", ">55"))) |>
  mutate(Underweight = ifelse(Underweight == 1, "Underweight", "Not Underweight")) |>
  mutate(Householdno = case_when(Householdno >= 2 & Householdno <= 3 ~ "2-3", 
                                 Householdno >= 4 & Householdno <= 6 ~ "4-6", 
                                 Householdno >= 7 & Householdno <= 9 ~ "7-9", 
                                 Householdno > 9 ~ ">9"), 
         Householdno = factor(Householdno, levels = c("2-3", "4-6", "7-9", ">9"))) |>
  rename("Household Number" = Householdno) |>
  mutate(HeardALnamebefore = ifelse(HeardALnamebefore == 0, 1, 0), 
         HeardTTnamebefore = ifelse(HeardTTnamebefore == 1, 1, 0), 
         HeardHWnamebefore = ifelse(HeardHWnamebefore == 2, 1, 0), 
         HeardHIVnamebefore = ifelse(HeardHIVnamebefore == 3, 1, 0), 
         HeardInWormnamebefore = ifelse(HeardInWormnamebefore == 4, 1, 0), 
         HeardMalanamebefore = ifelse(HeardMalanamebefore == 5, 1, 0), 
         HeardTBnamebefore = ifelse(HeardTBnamebefore == 6, 1, 0), 
         HeardSChnamebefore = ifelse(HeardSChnamebefore == 7, 1, 0)) |>
  mutate(HeardALnamebefore = ifelse(is.na(HeardALnamebefore), 0, HeardALnamebefore), 
         HeardTTnamebefore = ifelse(is.na(HeardTTnamebefore), 0, HeardTTnamebefore), 
         HeardHWnamebefore = ifelse(is.na(HeardHWnamebefore), 0, HeardHWnamebefore), 
         HeardHIVnamebefore = ifelse(is.na(HeardHIVnamebefore), 0, HeardHIVnamebefore), 
         HeardInWormnamebefore = ifelse(is.na(HeardInWormnamebefore), 0, HeardInWormnamebefore), 
         HeardMalanamebefore = ifelse(is.na(HeardMalanamebefore), 0, HeardMalanamebefore), 
         HeardTBnamebefore = ifelse(is.na(HeardTBnamebefore), 0, HeardTBnamebefore), 
         HeardSChnamebefore = ifelse(is.na(HeardSChnamebefore), 0, HeardSChnamebefore)) |>
  mutate(familytold = ifelse(familytold == 0, 1, 0), 
         HPtold = ifelse(HPtold == 1, 1, 0), 
         Teachertold = ifelse(Teachertold == 2, 1, 0), 
         Mediatold = ifelse(Mediatold == 3, 1, 0)) |>
  mutate(familytold = ifelse(is.na(familytold), 0, familytold), 
         HPtold = ifelse(is.na(HPtold), 0, HPtold), 
         Teachertold = ifelse(is.na(Teachertold), 0, Teachertold), 
         Mediatold = ifelse(is.na(Mediatold), 0, Mediatold)) |>
  mutate(Maternaleducationalstatus = case_when(Maternaleducationalstatus == 0 ~ "Illiterate", 
                                               Maternaleducationalstatus == 1 ~ "Primary School", 
                                               Maternaleducationalstatus == 2 ~ "High School", 
                                               Maternaleducationalstatus == 3 ~ "Higher Education", 
                                               Maternaleducationalstatus == 99 ~ "Don't Know", 
                                               is.na(Maternaleducationalstatus) ~ "Don't Know")) |>
  mutate(Maternaleducationalstatus = factor(Maternaleducationalstatus, levels = c("Illiterate", "Primary School", "High School", "Higher Education", "Don't Know"))) |>
  mutate(Dewormingin1yr = ifelse(Dewormingin1yr == 1, "Yes", "No")) |>
  mutate(Modeofdelivery = ifelse(Modeofdelivery == 1, "Vaginal", "C-Section")) |>
  mutate(Noofolderbrothres = case_when(Noofolderbrothres >= 1 & Noofolderbrothres <= 2 ~ "1-2", 
                                       Noofolderbrothres >= 3 & Noofolderbrothres <= 4 ~ "3-4", 
                                       Noofolderbrothres > 4 ~ ">4", 
                                       Noofolderbrothres == 0 ~ "0"), 
         Noofolderbrothres = factor(Noofolderbrothres, levels = c("0", "1-2", "3-4", ">4"))) |>
  rename(`Older Siblings` = Noofolderbrothres) |>
  mutate(youngerthan12years = case_when(youngerthan12years >= 1 & youngerthan12years <= 2 ~ "1-2", 
                                        youngerthan12years >= 3 & youngerthan12years <= 4 ~ "3-4", 
                                        youngerthan12years > 4 ~ ">4", 
                                        youngerthan12years == 0 ~ "0"), 
         youngerthan12years = factor(youngerthan12years, levels = c("0", "1-2", "3-4", ">4"))) |>
  rename("Siblings Younger than 12" = youngerthan12years) |>
  mutate(Everhadvaccinated = ifelse(Everhadvaccinated == 1, "Yes", "No")) |>
  mutate(Inthelasttwoweekschildhadfever = ifelse(Inthelasttwoweekschildhadfever == 1, "Yes", "No")) |>
  mutate(InthelasttwoweekschildhadDiahearria = ifelse(InthelasttwoweekschildhadDiahearria == 1, "Yes", "No")) |>
  mutate(Inthelasttwoweekschildhadcough = ifelse(Inthelasttwoweekschildhadcough == 1, "Yes", "No")) |>
  mutate(Childsfingernailtrimmed = ifelse(Childsfingernailtrimmed == 1, "Yes", "No")) |>
  mutate(Howoftendoyoutrimyourfingernails = case_when(Howoftendoyoutrimyourfingernails == 0 ~ "Less Than Once/Month", 
                                                      Howoftendoyoutrimyourfingernails == 1 ~ "Once per Week", 
                                                      Howoftendoyoutrimyourfingernails == 2 ~ "Once per two weeks", 
                                                      is.na(Howoftendoyoutrimyourfingernails) ~ "Don't Know"), 
         Howoftendoyoutrimyourfingernails = factor(Howoftendoyoutrimyourfingernails, levels = c("Less Than Once/Month", "Once per Week", "Once per two weeks", "Don't Know"))) |>
  mutate(Istheretoiletintheschool = case_when(Istheretoiletintheschool == 0 ~ "No", 
                                              Istheretoiletintheschool == 1 ~ "Yes", 
                                              is.na(Istheretoiletintheschool) ~ "Don't Know")) |>
  mutate(Ifthereisatoiletdoesthelatrinehavedoors = case_when(Ifthereisatoiletdoesthelatrinehavedoors == 0 ~ "No", 
                                                             Ifthereisatoiletdoesthelatrinehavedoors == 1 ~ "Yes", 
                                                             is.na(Ifthereisatoiletdoesthelatrinehavedoors) ~ "Don't Know")) |>
  mutate(Fliesobservedinaroundthelatrine = case_when(Fliesobservedinaroundthelatrine == 0 ~ "No", 
                                                     Fliesobservedinaroundthelatrine == 1 ~ "Yes", 
                                                     is.na(Fliesobservedinaroundthelatrine) ~ "Don't Know")) |>
  mutate(Vissiblestoolobservedonlatrinefloor = case_when(Vissiblestoolobservedonlatrinefloor == 0 ~ "No", 
                                                         Vissiblestoolobservedonlatrinefloor == 1 ~ "Yes", 
                                                         is.na(Vissiblestoolobservedonlatrinefloor) ~ "Don't Know")) |>
  mutate(Doyouknowhowintestinalwomstransmitted = case_when(Doyouknowhowintestinalwomstransmitted == 0 ~ "No", 
                                                           Doyouknowhowintestinalwomstransmitted == 1 ~ "Yes", 
                                                           Doyouknowhowintestinalwomstransmitted == 99 ~ "Don't Know")) |>
  mutate(Doyouknowwhywormsarebadforyourhealth = case_when(Doyouknowwhywormsarebadforyourhealth == 0 ~ "No", 
                                                          Doyouknowwhywormsarebadforyourhealth == 1 ~ "Yes", 
                                                          Doyouknowwhywormsarebadforyourhealth == 99 ~ "Don't Know")) |>
  mutate(Doyouknowhowyoucanavoidgettingtheseworms = case_when(Doyouknowhowyoucanavoidgettingtheseworms == 0 ~ "No", 
                                                              Doyouknowhowyoucanavoidgettingtheseworms == 1 ~ "Yes", 
                                                              Doyouknowhowyoucanavoidgettingtheseworms == 99 ~ "Don't Know")) |>
  mutate(Yourlivingaddress = ifelse(Yourlivingaddress != "Amenu" & Yourlivingaddress != "Bosa Kito" & Yourlivingaddress != "Hermata mentina" & Yourlivingaddress != "Jiren" & Yourlivingaddress != "Mentina" & Yourlivingaddress != "Seto semero" & Yourlivingaddress != "Welda", "Other", Yourlivingaddress)) |>
  mutate(Familyoccupation = ifelse(Familyoccupation != "Police man" & Familyoccupation != "Teacher" & Familyoccupation != "merchant" & Familyoccupation != "Carprnter" & Familyoccupation != "Dailly labor" & Familyoccupation != "Security" & Familyoccupation != "Farmer" & Familyoccupation != "Government employee" & Familyoccupation != "Merchant", "Other", Familyoccupation)) |>
  mutate(Familyoccupation = ifelse(Familyoccupation == "merchant", "Merchant", Familyoccupation), 
         Familyoccupation = ifelse(is.na(Familyoccupation), "Other", Familyoccupation)) |>
  mutate(Whatmaterialyourhousefloormadefrom = case_when(Whatmaterialyourhousefloormadefrom == 0 ~ "Dust", 
                                                        Whatmaterialyourhousefloormadefrom == 1 ~ "Cement", 
                                                        Whatmaterialyourhousefloormadefrom == 2 ~ "Plastic Covered", 
                                                        is.na(Whatmaterialyourhousefloormadefrom) ~ "Don't Know")) |>
  mutate(Isyourkitcheninyourhouseorsepatated = case_when(Isyourkitcheninyourhouseorsepatated == 0 ~ "Within House", 
                                                         Isyourkitcheninyourhouseorsepatated == 2 ~ "Separated", 
                                                         Isyourkitcheninyourhouseorsepatated == 99 | is.na(Isyourkitcheninyourhouseorsepatated) ~ "Don't Know")) |>
  mutate(Ifseparetedwhatmaterialsyourkitchenmadefrom = case_when(Ifseparetedwhatmaterialsyourkitchenmadefrom == 0 ~ "Dust", 
                                                                 Ifseparetedwhatmaterialsyourkitchenmadefrom == 1 ~ "Cement", 
                                                                 Ifseparetedwhatmaterialsyourkitchenmadefrom == 2 ~ "Plastic Covered", 
                                                                 is.na(Ifseparetedwhatmaterialsyourkitchenmadefrom) | Ifseparetedwhatmaterialsyourkitchenmadefrom == 9 ~ "Don't Know")) |>
  mutate(Roof = ifelse(Roof == 1, "Yes", "No"), 
         Wall = ifelse(Wall == 2, "Yes", "No"), 
         None = ifelse(None == 0, "Yes", "No")) |>
  mutate(Roof = ifelse(is.na(Roof), "No", Roof), 
         Wall = ifelse(is.na(Wall), "No", Wall), 
         None = ifelse(is.na(None), "No", None)) |>
  mutate(Wood = ifelse(Wood == 0, "Yes", "No"), 
         Gas = ifelse(Gas == 1, "Yes", "No"), 
         Coal = ifelse(Coal == 2, "Yes", "No"), 
         Kerosisn = ifelse(Kerosisn == 3, "Yes", "No"), 
         Electric = ifelse(Electric == 4, "Yes", "No")) |>
  mutate(Wood = ifelse(is.na(Wood), "No", Wood), 
         Gas = ifelse(is.na(Gas), "No", Gas), 
         Kerosisn = ifelse(is.na(Kerosisn), "No", Kerosisn), 
         Coal = ifelse(is.na(Coal), "No", Coal), 
         Electric = ifelse(is.na(Electric), "No", Electric)) |>
  mutate(Doyouhaveelectricityinyourhaouse = case_when(Doyouhaveelectricityinyourhaouse == 0 ~ "No", 
                                                      Doyouhaveelectricityinyourhaouse == 1 ~ "Yes", 
                                                      is.na(Doyouhaveelectricityinyourhaouse) ~ "Don't Know")) |>
  mutate(Doesyourfamillyownradio = case_when(Doesyourfamillyownradio == 0 ~ "No", 
                                             Doesyourfamillyownradio == 1 ~ "Yes", 
                                             is.na(Doesyourfamillyownradio) ~ "Don't Know")) |>
  mutate(Doesyourfamillyowntelevision = case_when(Doesyourfamillyowntelevision == 0 ~ "No", 
                                                  Doesyourfamillyowntelevision == 1 ~ "Yes", 
                                                  is.na(Doesyourfamillyowntelevision) ~ "Don't Know")) |>
  mutate(Doesyourfamillymemberownaphone = case_when(Doesyourfamillymemberownaphone == 0 ~ "No", 
                                                    Doesyourfamillymemberownaphone == 1 ~ "Yes", 
                                                    is.na(Doesyourfamillymemberownaphone) ~ "Don't Know")) |>
  mutate(Forwhatpurposeyourfamilyuseaphone = case_when(Forwhatpurposeyourfamilyuseaphone == 0 ~ "Call Only", 
                                                       Forwhatpurposeyourfamilyuseaphone == 1 ~ "Call/Radio/Internet", 
                                                       is.na(Forwhatpurposeyourfamilyuseaphone) ~ "Don't Know")) |>
  mutate(Cattle = ifelse(Cattle == 0, "Yes", "No"), 
         SheepGoat = ifelse(SheepGoat == 1, "Yes", "No"), 
         Chicken = ifelse(Chicken == 2, "Yes", "No"), 
         Pet = ifelse(Pet == 3, "Yes", "No"), 
         Nodomanimal = ifelse(Nodomanimal == 4, "Yes", "No")) |>
  mutate(Cattle = ifelse(is.na(Cattle), "No", Cattle), 
         SheepGoat = ifelse(is.na(SheepGoat), "No", SheepGoat), 
         Chicken = ifelse(is.na(Chicken), "No", Chicken), 
         Pet = ifelse(is.na(Pet), "No", Pet), 
         Nodomanimal = ifelse(is.na(Nodomanimal), "No", Nodomanimal)) |>
  mutate(Doyouhavepotablewaterinyourhouse = case_when(Doyouhavepotablewaterinyourhouse == 0 ~ "No", 
                                                      Doyouhavepotablewaterinyourhouse == 1 ~ "Yes", 
                                                      is.na(Doyouhavepotablewaterinyourhouse) ~ "Don't Know")) |>
  mutate(Ifnotwheredoyougetyordrinkingwaterfrom = case_when(Ifnotwheredoyougetyordrinkingwaterfrom == 0 ~ "Neighbor", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 1 ~ "River", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 2 ~ "Well", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 3 ~ "Truck",
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 4 ~ "Tank", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 5 ~ "Public Fountain",
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == "Tap water" ~ "Tap Water", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 10 ~ NA)) |>
  mutate(Doyoudrinkdirectlyordoyoutreat = case_when(Doyoudrinkdirectlyordoyoutreat == 1 ~ "Directly", 
                                                    Doyoudrinkdirectlyordoyoutreat == 0 ~ "Treat", 
                                                    is.na(Doyoudrinkdirectlyordoyoutreat) ~ "Don't Know")) |>
  mutate(Ifyoutreathow = case_when(Ifyoutreathow == 0 ~ "Boil", 
                                   Ifyoutreathow == 1 ~ "Chemical", 
                                   Ifyoutreathow == 2 ~ "Filter", 
                                   Ifyoutreathow == 99 ~ "Don't Know")) |>
  mutate(Doyourfamilyownlatrine = ifelse(Doyourfamilyownlatrine == 1, "Yes", "No")) |>
  mutate(Isyourlatrineinsideoroutside = ifelse(Isyourlatrineinsideoroutside == 1, "Inside", "Outside")) |>
  mutate(Ifoutsidedistancefromyourhouse = as.numeric(str_sub(Ifoutsidedistancefromyourhouse, 1, -2)), 
         Ifoutsidedistancefromyourhouse = case_when(Ifoutsidedistancefromyourhouse >= 0 & Ifoutsidedistancefromyourhouse < 5 ~ "0-5", 
                                                    Ifoutsidedistancefromyourhouse >= 5 & Ifoutsidedistancefromyourhouse < 10 ~ "5-10", 
                                                    Ifoutsidedistancefromyourhouse >= 10 & Ifoutsidedistancefromyourhouse <= 20 ~ "10-20", 
                                                    Ifoutsidedistancefromyourhouse > 20 ~ ">20"), 
         Ifoutsidedistancefromyourhouse = factor(Ifoutsidedistancefromyourhouse, levels = c("0-5", "5-10", "10-20", ">20"))) |>
  mutate(Distancebetweenlatrineandyourkitchen = as.numeric(str_sub(Distancebetweenlatrineandyourkitchen, 1, -2)), 
         Distancebetweenlatrineandyourkitchen = case_when(Distancebetweenlatrineandyourkitchen >= 0 & Distancebetweenlatrineandyourkitchen < 5 ~ "0-5", 
                                                          Distancebetweenlatrineandyourkitchen >= 5 & Distancebetweenlatrineandyourkitchen < 10 ~ "5-10", 
                                                          Distancebetweenlatrineandyourkitchen >= 10 & Distancebetweenlatrineandyourkitchen <= 20 ~ "10-20", 
                                                          Distancebetweenlatrineandyourkitchen > 20 ~ ">20"), 
         Distancebetweenlatrineandyourkitchen = factor(Distancebetweenlatrineandyourkitchen, levels = c("0-5", "5-10", "10-20", ">20"))) |>
  mutate(Yourlatrineconnectedto = case_when(Yourlatrineconnectedto == 0 ~ "Sewage", 
                                            Yourlatrineconnectedto == 1 ~ "Ditch", 
                                            Yourlatrineconnectedto == 2 ~ "River", 
                                            Yourlatrineconnectedto == 3 ~ "Well", 
                                            Yourlatrineconnectedto == 99 | is.na(Yourlatrineconnectedto) ~ "NA/Don't Know")) |>
  mutate(Howoftendoyoubathinriver = case_when(Howoftendoyoubathinriver == 0 ~ "Always", 
                                              Howoftendoyoubathinriver == 1 ~ "Sometimes", 
                                              Howoftendoyoubathinriver == 2 ~ "Never",
                                              Howoftendoyoubathinriver == 99 | is.na(Howoftendoyoubathinriver) ~ "NA/Don't Know"), 
         Howoftendoyoubathinriver = factor(Howoftendoyoubathinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howoftendoyouwashyourclothesinriver = case_when(Howoftendoyouwashyourclothesinriver == 0 ~ "Always", 
                                                         Howoftendoyouwashyourclothesinriver == 1 ~ "Sometimes", 
                                                         Howoftendoyouwashyourclothesinriver == 2 ~ "Never",
                                                         Howoftendoyouwashyourclothesinriver == 99 | is.na(Howoftendoyouwashyourclothesinriver) ~ "NA/Don't Know"), 
         Howoftendoyouwashyourclothesinriver = factor(Howoftendoyouwashyourclothesinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyoudeficateintheopenfield = case_when(Doyoudeficateintheopenfield == 0 ~ "Always", 
                                                 Doyoudeficateintheopenfield == 1 ~ "Sometimes", 
                                                 Doyoudeficateintheopenfield == 2 ~ "Never",
                                                 Doyoudeficateintheopenfield == 99 | is.na(Doyoudeficateintheopenfield) ~ "NA/Don't Know"), 
         Doyoudeficateintheopenfield = factor(Doyoudeficateintheopenfield, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouuseschoollatrine = case_when(Doyouuseschoollatrine == 0 ~ "Always", 
                                           Doyouuseschoollatrine == 1 ~ "Sometimes", 
                                           Doyouuseschoollatrine == 2 ~ "Never",
                                           Doyouuseschoollatrine == 99 | is.na(Doyouuseschoollatrine) ~ "NA/Don't Know"), 
         Doyouuseschoollatrine = factor(Doyouuseschoollatrine, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouusetoiletpapertowipeyourbumafterdefication = case_when(Doyouusetoiletpapertowipeyourbumafterdefication == 0 ~ "Always", 
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 1 ~ "Sometimes", 
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 2 ~ "Never",
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 99 | is.na(Doyouusetoiletpapertowipeyourbumafterdefication) ~ "NA/Don't Know"), 
         Doyouusetoiletpapertowipeyourbumafterdefication = factor(Doyouusetoiletpapertowipeyourbumafterdefication, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouwashyourhandsaftertoilet = case_when(Doyouwashyourhandsaftertoilet == 0 ~ "Always", 
                                                   Doyouwashyourhandsaftertoilet == 1 ~ "Sometimes", 
                                                   Doyouwashyourhandsaftertoilet == 2 ~ "Never",
                                                   Doyouwashyourhandsaftertoilet == 99 | is.na(Doyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
         Doyouwashyourhandsaftertoilet = factor(Doyouwashyourhandsaftertoilet, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howdoyouwashyourhandsaftertoilet = case_when(Howdoyouwashyourhandsaftertoilet == 0 ~ "Water", 
                                               Howdoyouwashyourhandsaftertoilet == 1 ~ "Soap and Water",
                                               Howdoyouwashyourhandsaftertoilet == 99 | is.na(Howdoyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
        Howdoyouwashyourhandsaftertoilet = factor(Howdoyouwashyourhandsaftertoilet, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  mutate(Ifwithsoaphowoftenyouuseit = case_when(Ifwithsoaphowoftenyouuseit == 0 ~ "Always", 
                                                Ifwithsoaphowoftenyouuseit == 1 ~ "Sometimes",
                                                Ifwithsoaphowoftenyouuseit == 2 ~ "Never",
                                                Ifwithsoaphowoftenyouuseit == 99 | is.na(Ifwithsoaphowoftenyouuseit) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit = factor(Ifwithsoaphowoftenyouuseit, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  mutate(Doyouwashyourhandsbeforeeating = case_when(Doyouwashyourhandsbeforeeating == 0 ~ "Always", 
                                                  Doyouwashyourhandsbeforeeating == 1 ~ "Sometimes", 
                                                  Doyouwashyourhandsbeforeeating == 2 ~ "Never",
                                                  Doyouwashyourhandsbeforeeating == 99 | is.na(Doyouwashyourhandsbeforeeating) ~ "NA/Don't Know"), 
       Doyouwashyourhandsbeforeeating = factor(Doyouwashyourhandsbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howdoyouwashyourhandsbeforeeating = case_when(Howdoyouwashyourhandsbeforeeating == 0 ~ "Water", 
                                                       Howdoyouwashyourhandsbeforeeating == 1 ~ "Soap and Water",
                                                       Howdoyouwashyourhandsbeforeeating == 99 | is.na(Howdoyouwashyourhandsbeforeeating) | Howdoyouwashyourhandsbeforeeating == 2 ~ "NA/Don't Know"), 
         Howdoyouwashyourhandsbeforeeating = factor(Howdoyouwashyourhandsbeforeeating, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  mutate(Ifwithsoaphowoftenyouuseit_A = case_when(Ifwithsoaphowoftenyouuseit_A == 0 ~ "Always", 
                                                  Ifwithsoaphowoftenyouuseit_A == 1 ~ "Sometimes",
                                                  Ifwithsoaphowoftenyouuseit_A == 2 ~ "Never",
                                                  Ifwithsoaphowoftenyouuseit_A == 99 | is.na(Ifwithsoaphowoftenyouuseit_A) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit_A = factor(Ifwithsoaphowoftenyouuseit_A, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  mutate(Doyoueatsoil = case_when(Doyoueatsoil == 0 ~ "Always", 
                                  Doyoueatsoil == 1 ~ "Sometimes", 
                                  Doyoueatsoil == 2 ~ "Never",
                                  Doyoueatsoil == 99 | is.na(Doyoueatsoil) ~ "NA/Don't Know"), 
         Doyoueatsoil = factor(Doyoueatsoil, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Whatisyourfavoritefruitethatyoueat = ifelse(Whatisyourfavoritefruitethatyoueat != "Mango" & Whatisyourfavoritefruitethatyoueat != "Banana" & Whatisyourfavoritefruitethatyoueat != "Avocado" & Whatisyourfavoritefruitethatyoueat != "Apple" & Whatisyourfavoritefruitethatyoueat != "Pinapple", "Other", Whatisyourfavoritefruitethatyoueat), 
         Whatisyourfavoritefruitethatyoueat = ifelse(is.na(whatfooddoyoueatmost), "Other", Whatisyourfavoritefruitethatyoueat)) |>
  mutate(Doyouwashyourfruitesbeforeeating = case_when(Doyouwashyourfruitesbeforeeating == 0 ~ "Always", 
                                                      Doyouwashyourfruitesbeforeeating == 1 ~ "Sometimes", 
                                                      Doyouwashyourfruitesbeforeeating == 2 ~ "Never",
                                                      Doyouwashyourfruitesbeforeeating == 99 | is.na(Doyouwashyourfruitesbeforeeating) ~ "NA/Don't Know"), 
         Doyouwashyourfruitesbeforeeating = factor(Doyouwashyourfruitesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyoueatraworundercookedvegitables = case_when(Doyoueatraworundercookedvegitables == 0 ~ "Always", 
                                                        Doyoueatraworundercookedvegitables == 1 ~ "Sometimes", 
                                                        Doyoueatraworundercookedvegitables == 2 ~ "Never",
                                                        Doyoueatraworundercookedvegitables == 99 | is.na(Doyoueatraworundercookedvegitables) ~ "NA/Don't Know"), 
         Doyoueatraworundercookedvegitables = factor(Doyoueatraworundercookedvegitables, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Ifyeshowoftenyouwashyourvegitablesbeforeeating = case_when(Ifyeshowoftenyouwashyourvegitablesbeforeeating == 0 ~ "Always", 
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 1 ~ "Sometimes", 
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 2 ~ "Never",
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 99 | is.na(Ifyeshowoftenyouwashyourvegitablesbeforeeating) ~ "NA/Don't Know"), 
         Ifyeshowoftenyouwashyourvegitablesbeforeeating = factor(Ifyeshowoftenyouwashyourvegitablesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouwalkbarefoot = case_when(Doyouwalkbarefoot == 0 ~ "Always", 
                                       Doyouwalkbarefoot == 1 ~ "Sometimes", 
                                       Doyouwalkbarefoot == 2 ~ "Never",
                                       Doyouwalkbarefoot == 99 | is.na(Doyouwalkbarefoot) ~ "NA/Don't Know"), 
         Doyouwalkbarefoot = factor(Doyouwalkbarefoot, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(whenyouareathomedoyoupreferetousesandalsorshoe = case_when(whenyouareathomedoyoupreferetousesandalsorshoe == 0 ~ "Barefoot", 
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 1 ~ "Sandals", 
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 2 ~ "Shoes",
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 99 | is.na(whenyouareathomedoyoupreferetousesandalsorshoe) ~ "NA/Don't Know"), 
         whenyouareathomedoyoupreferetousesandalsorshoe = factor(whenyouareathomedoyoupreferetousesandalsorshoe, levels = c("Barefoot", "Sandals", "Shoes", "NA/Don't Know"))) |>
  mutate(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill = case_when(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 0 ~ "No", 
                                                                                      Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 1 ~ "Yes",
                                                                                      Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 99 | is.na(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill) ~ "NA/Don't Know"))|>
  mutate(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill = case_when(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "1Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "2Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "3month" ~ "1-3 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "5Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "6Month" ~ "4-6 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "7Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "9Month" ~ "7-9 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "11Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "12month" ~ ">9 Months", 
                                                                         is.na(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill) ~ NA)) |>
  mutate(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics = case_when(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 0 ~ "No", 
                                                                                Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 1 ~ "Yes",
                                                                                Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 99 | is.na(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics) ~ "NA/Don't Know")) |>
  mutate(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently = case_when(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 0 ~ "No", 
                                                                                  Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 1 ~ "Yes",
                                                                                  Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 99 | is.na(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently) ~ "NA/Don't Know")) |>
  mutate(Inthepastthreemonthshaveyoutakenanyantimalarialdrug = case_when(Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 0 ~ "No", 
                                                                         Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 1 ~ "Yes",
                                                                         Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 99 | is.na(Inthepastthreemonthshaveyoutakenanyantimalarialdrug) ~ "NA/Don't Know")) |>
  mutate(Haveachildwheezingorwhistlinginchest = case_when(Haveachildwheezingorwhistlinginchest == 2 ~ "No", 
                                                          Haveachildwheezingorwhistlinginchest == 1 ~ "Yes",
                                                          Haveachildwheezingorwhistlinginchest == 99 | is.na(Haveachildwheezingorwhistlinginchest) ~ "NA/Don't Know")) |>
  mutate(Howmanytimesinthelastyearthechildhadwheeziling = case_when(Howmanytimesinthelastyearthechildhadwheeziling == 2 ~ "1-3", 
                                                                    Howmanytimesinthelastyearthechildhadwheeziling == 1 ~ "0",
                                                                    Howmanytimesinthelastyearthechildhadwheeziling == 99 | is.na(Howmanytimesinthelastyearthechildhadwheeziling) ~ "NA/Don't Know")) |>
  mutate(Hasachildeverhadasthma = case_when(Hasachildeverhadasthma == 2 ~ "No", 
                                            Hasachildeverhadasthma == 1 ~ "Yes",
                                            Hasachildeverhadasthma == 99 | is.na(Hasachildeverhadasthma) ~ "NA/Don't Know")) |>
  mutate(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases = case_when(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 2 ~ "No", 
                                                                                    Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 1 ~ "Yes",
                                                                                    Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 99 | is.na(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases) ~ "NA/Don't Know")) |>
  mutate(Ifyeshasthisrashaffectedtheelbowfolds = ifelse(is.na(Ifyeshasthisrashaffectedtheelbowfolds), 0, Ifyeshasthisrashaffectedtheelbowfolds), 
         Ifyeshasthisrashaffectedbehindtheknees = ifelse(is.na(Ifyeshasthisrashaffectedbehindtheknees), 0, Ifyeshasthisrashaffectedbehindtheknees), 
         Ifyeshasthisrashaffectedinfrontoftheankls = ifelse(is.na(Ifyeshasthisrashaffectedinfrontoftheankls), 0, Ifyeshasthisrashaffectedinfrontoftheankls), 
         Ifyeshasthisrashaffectedunderthebuttucks = ifelse(is.na(Ifyeshasthisrashaffectedunderthebuttucks), 0, Ifyeshasthisrashaffectedunderthebuttucks), 
         Ifyeshasthisrashaffectedaroundtheneck = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheneck), 0, Ifyeshasthisrashaffectedaroundtheneck), 
         Ifyeshasthisrashaffectedaroundtheeyesears = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheeyesears), 0, Ifyeshasthisrashaffectedaroundtheeyesears)) |>
  dplyr::select(-AL, -TT, -HW, -SM, -Others, -Date, -parentconsent, -childassent, -Grade, -GradeLet, -WAZ, -HAZ, -BAZ, -Wheredoyoulive, -Doesachildasthmainthelastyear, -Doesachildhadasthmainthelast2years, -Hasthisbeenconfirmedbyadoctor, -Doesthechildhadhayfeverorsneezingwithrunningnoseinthelast2years, -Doesthechildhadhayfeverorsneezingwithrunningnoseinthelastyear, -V137, -AgeGroup, -HouseholdGroup) |>
  mutate(Doesthechildhadhayfeverorpersistantsneezing = case_when(Doesthechildhadhayfeverorpersistantsneezing == 2 ~ "No", 
                                                                 Doesthechildhadhayfeverorpersistantsneezing == 1 ~ "Yes",
                                                                 Doesthechildhadhayfeverorpersistantsneezing == 99 | is.na(Doesthechildhadhayfeverorpersistantsneezing) ~ "NA/Don't Know")) |>
  mutate(Arechildsfingernailsdirty = ifelse(Arechildsfingernailsdirty == 1, "Yes", "No")) |>
  rename(Kerosene = Kerosisn) |>
  dplyr::select(-ID, -Seq_ID) |>
  mutate(Height = case_when(Height < 1.15 ~ "< 1.15", 
                            Height >= 1.15 & Height < 1.35 ~ "1.15-1.35", 
                            Height >= 1.35 & Height < 1.55 ~ "1.35-1.55", 
                            Height >= 1.55 ~ ">1.55")) |>
  mutate(Height = factor(Height, levels = c("< 1.15", "1.15-1.35", "1.35-1.55", ">1.55")))
  






  

analysis_dat <- sample_data |>
  mutate(Underweight = ifelse(Underweight == 1, "Underweight", "Not Underweight")) |>
  mutate(HeardALnamebefore = ifelse(HeardALnamebefore == 0, 1, 0), 
         HeardTTnamebefore = ifelse(HeardTTnamebefore == 1, 1, 0), 
         HeardHWnamebefore = ifelse(HeardHWnamebefore == 2, 1, 0), 
         HeardHIVnamebefore = ifelse(HeardHIVnamebefore == 3, 1, 0), 
         HeardInWormnamebefore = ifelse(HeardInWormnamebefore == 4, 1, 0), 
         HeardMalanamebefore = ifelse(HeardMalanamebefore == 5, 1, 0), 
         HeardTBnamebefore = ifelse(HeardTBnamebefore == 6, 1, 0), 
         HeardSChnamebefore = ifelse(HeardSChnamebefore == 7, 1, 0)) |>
  mutate(HeardALnamebefore = ifelse(is.na(HeardALnamebefore), 0, HeardALnamebefore), 
         HeardTTnamebefore = ifelse(is.na(HeardTTnamebefore), 0, HeardTTnamebefore), 
         HeardHWnamebefore = ifelse(is.na(HeardHWnamebefore), 0, HeardHWnamebefore), 
         HeardHIVnamebefore = ifelse(is.na(HeardHIVnamebefore), 0, HeardHIVnamebefore), 
         HeardInWormnamebefore = ifelse(is.na(HeardInWormnamebefore), 0, HeardInWormnamebefore), 
         HeardMalanamebefore = ifelse(is.na(HeardMalanamebefore), 0, HeardMalanamebefore), 
         HeardTBnamebefore = ifelse(is.na(HeardTBnamebefore), 0, HeardTBnamebefore), 
         HeardSChnamebefore = ifelse(is.na(HeardSChnamebefore), 0, HeardSChnamebefore)) |>
  mutate(familytold = ifelse(familytold == 0, 1, 0), 
         HPtold = ifelse(HPtold == 1, 1, 0), 
         Teachertold = ifelse(Teachertold == 2, 1, 0), 
         Mediatold = ifelse(Mediatold == 3, 1, 0)) |>
  mutate(familytold = ifelse(is.na(familytold), 0, familytold), 
         HPtold = ifelse(is.na(HPtold), 0, HPtold), 
         Teachertold = ifelse(is.na(Teachertold), 0, Teachertold), 
         Mediatold = ifelse(is.na(Mediatold), 0, Mediatold)) |>
  mutate(Maternaleducationalstatus = case_when(Maternaleducationalstatus == 0 ~ "Illiterate", 
                                               Maternaleducationalstatus == 1 ~ "Primary School", 
                                               Maternaleducationalstatus == 2 ~ "High School", 
                                               Maternaleducationalstatus == 3 ~ "Higher Education", 
                                               Maternaleducationalstatus == 99 ~ "Don't Know", 
                                               is.na(Maternaleducationalstatus) ~ "Don't Know")) |>
  mutate(Maternaleducationalstatus = factor(Maternaleducationalstatus, levels = c("Illiterate", "Primary School", "High School", "Higher Education", "Don't Know"))) |>
  mutate(Dewormingin1yr = ifelse(Dewormingin1yr == 1, "Yes", "No")) |>
  mutate(Modeofdelivery = ifelse(Modeofdelivery == 1, "Vaginal", "C-Section")) |>
  rename(`Older Siblings` = Noofolderbrothres) |>
  rename("Household Number" = Householdno) |>
  rename("Siblings Younger than 12" = youngerthan12years) |>
  mutate(Everhadvaccinated = ifelse(Everhadvaccinated == 1, "Yes", "No")) |>
  mutate(Inthelasttwoweekschildhadfever = ifelse(Inthelasttwoweekschildhadfever == 1, "Yes", "No")) |>
  mutate(InthelasttwoweekschildhadDiahearria = ifelse(InthelasttwoweekschildhadDiahearria == 1, "Yes", "No")) |>
  mutate(Inthelasttwoweekschildhadcough = ifelse(Inthelasttwoweekschildhadcough == 1, "Yes", "No")) |>
  mutate(Childsfingernailtrimmed = ifelse(Childsfingernailtrimmed == 1, "Yes", "No")) |>
  mutate(Howoftendoyoutrimyourfingernails = case_when(Howoftendoyoutrimyourfingernails == 0 ~ "Less Than Once/Month", 
                                                      Howoftendoyoutrimyourfingernails == 1 ~ "Once per Week", 
                                                      Howoftendoyoutrimyourfingernails == 2 ~ "Once per two weeks", 
                                                      is.na(Howoftendoyoutrimyourfingernails) ~ "Don't Know"), 
         Howoftendoyoutrimyourfingernails = factor(Howoftendoyoutrimyourfingernails, levels = c("Less Than Once/Month", "Once per Week", "Once per two weeks", "Don't Know"))) |>
  mutate(Istheretoiletintheschool = case_when(Istheretoiletintheschool == 0 ~ "No", 
                                              Istheretoiletintheschool == 1 ~ "Yes", 
                                              is.na(Istheretoiletintheschool) ~ "Don't Know")) |>
  mutate(Ifthereisatoiletdoesthelatrinehavedoors = case_when(Ifthereisatoiletdoesthelatrinehavedoors == 0 ~ "No", 
                                                             Ifthereisatoiletdoesthelatrinehavedoors == 1 ~ "Yes", 
                                                             is.na(Ifthereisatoiletdoesthelatrinehavedoors) ~ "Don't Know")) |>
  mutate(Fliesobservedinaroundthelatrine = case_when(Fliesobservedinaroundthelatrine == 0 ~ "No", 
                                                     Fliesobservedinaroundthelatrine == 1 ~ "Yes", 
                                                     is.na(Fliesobservedinaroundthelatrine) ~ "Don't Know")) |>
  mutate(Vissiblestoolobservedonlatrinefloor = case_when(Vissiblestoolobservedonlatrinefloor == 0 ~ "No", 
                                                         Vissiblestoolobservedonlatrinefloor == 1 ~ "Yes", 
                                                         is.na(Vissiblestoolobservedonlatrinefloor) ~ "Don't Know")) |>
  mutate(Doyouknowhowintestinalwomstransmitted = case_when(Doyouknowhowintestinalwomstransmitted == 0 ~ "No", 
                                                           Doyouknowhowintestinalwomstransmitted == 1 ~ "Yes", 
                                                           Doyouknowhowintestinalwomstransmitted == 99 ~ "Don't Know")) |>
  mutate(Doyouknowwhywormsarebadforyourhealth = case_when(Doyouknowwhywormsarebadforyourhealth == 0 ~ "No", 
                                                          Doyouknowwhywormsarebadforyourhealth == 1 ~ "Yes", 
                                                          Doyouknowwhywormsarebadforyourhealth == 99 ~ "Don't Know")) |>
  mutate(Doyouknowhowyoucanavoidgettingtheseworms = case_when(Doyouknowhowyoucanavoidgettingtheseworms == 0 ~ "No", 
                                                              Doyouknowhowyoucanavoidgettingtheseworms == 1 ~ "Yes", 
                                                              Doyouknowhowyoucanavoidgettingtheseworms == 99 ~ "Don't Know")) |>
  mutate(Yourlivingaddress = ifelse(Yourlivingaddress != "Amenu" & Yourlivingaddress != "Bosa Kito" & Yourlivingaddress != "Hermata mentina" & Yourlivingaddress != "Jiren" & Yourlivingaddress != "Mentina" & Yourlivingaddress != "Seto semero" & Yourlivingaddress != "Welda", "Other", Yourlivingaddress)) |>
  mutate(Familyoccupation = ifelse(Familyoccupation != "Police man" & Familyoccupation != "Teacher" & Familyoccupation != "merchant" & Familyoccupation != "Carprnter" & Familyoccupation != "Dailly labor" & Familyoccupation != "Security" & Familyoccupation != "Farmer" & Familyoccupation != "Government employee" & Familyoccupation != "Merchant", "Other", Familyoccupation)) |>
  mutate(Familyoccupation = ifelse(Familyoccupation == "merchant", "Merchant", Familyoccupation), 
         Familyoccupation = ifelse(is.na(Familyoccupation), "Other", Familyoccupation)) |>
  mutate(Whatmaterialyourhousefloormadefrom = case_when(Whatmaterialyourhousefloormadefrom == 0 ~ "Dust", 
                                                        Whatmaterialyourhousefloormadefrom == 1 ~ "Cement", 
                                                        Whatmaterialyourhousefloormadefrom == 2 ~ "Plastic Covered", 
                                                        is.na(Whatmaterialyourhousefloormadefrom) ~ "Don't Know")) |>
  mutate(Isyourkitcheninyourhouseorsepatated = case_when(Isyourkitcheninyourhouseorsepatated == 0 ~ "Within House", 
                                                         Isyourkitcheninyourhouseorsepatated == 2 ~ "Separated", 
                                                         Isyourkitcheninyourhouseorsepatated == 99 | is.na(Isyourkitcheninyourhouseorsepatated) ~ "Don't Know")) |>
  mutate(Ifseparetedwhatmaterialsyourkitchenmadefrom = case_when(Ifseparetedwhatmaterialsyourkitchenmadefrom == 0 ~ "Dust", 
                                                                 Ifseparetedwhatmaterialsyourkitchenmadefrom == 1 ~ "Cement", 
                                                                 Ifseparetedwhatmaterialsyourkitchenmadefrom == 2 ~ "Plastic Covered", 
                                                                 is.na(Ifseparetedwhatmaterialsyourkitchenmadefrom) | Ifseparetedwhatmaterialsyourkitchenmadefrom == 9 ~ "Don't Know")) |>
  mutate(Roof = ifelse(Roof == 1, "Yes", "No"), 
         Wall = ifelse(Wall == 2, "Yes", "No"), 
         None = ifelse(None == 0, "Yes", "No")) |>
  mutate(Roof = ifelse(is.na(Roof), "No", Roof), 
         Wall = ifelse(is.na(Wall), "No", Wall), 
         None = ifelse(is.na(None), "No", None)) |>
  mutate(Wood = ifelse(Wood == 0, "Yes", "No"), 
         Gas = ifelse(Gas == 1, "Yes", "No"), 
         Coal = ifelse(Coal == 2, "Yes", "No"), 
         Kerosisn = ifelse(Kerosisn == 3, "Yes", "No"), 
         Electric = ifelse(Electric == 4, "Yes", "No")) |>
  mutate(Wood = ifelse(is.na(Wood), "No", Wood), 
         Gas = ifelse(is.na(Gas), "No", Gas), 
         Kerosisn = ifelse(is.na(Kerosisn), "No", Kerosisn), 
         Coal = ifelse(is.na(Coal), "No", Coal), 
         Electric = ifelse(is.na(Electric), "No", Electric)) |>
  mutate(Doyouhaveelectricityinyourhaouse = case_when(Doyouhaveelectricityinyourhaouse == 0 ~ "No", 
                                                      Doyouhaveelectricityinyourhaouse == 1 ~ "Yes", 
                                                      is.na(Doyouhaveelectricityinyourhaouse) ~ "Don't Know")) |>
  mutate(Doesyourfamillyownradio = case_when(Doesyourfamillyownradio == 0 ~ "No", 
                                             Doesyourfamillyownradio == 1 ~ "Yes", 
                                             is.na(Doesyourfamillyownradio) ~ "Don't Know")) |>
  mutate(Doesyourfamillyowntelevision = case_when(Doesyourfamillyowntelevision == 0 ~ "No", 
                                                  Doesyourfamillyowntelevision == 1 ~ "Yes", 
                                                  is.na(Doesyourfamillyowntelevision) ~ "Don't Know")) |>
  mutate(Doesyourfamillymemberownaphone = case_when(Doesyourfamillymemberownaphone == 0 ~ "No", 
                                                    Doesyourfamillymemberownaphone == 1 ~ "Yes", 
                                                    is.na(Doesyourfamillymemberownaphone) ~ "Don't Know")) |>
  mutate(Forwhatpurposeyourfamilyuseaphone = case_when(Forwhatpurposeyourfamilyuseaphone == 0 ~ "Call Only", 
                                                       Forwhatpurposeyourfamilyuseaphone == 1 ~ "Call/Radio/Internet", 
                                                       is.na(Forwhatpurposeyourfamilyuseaphone) ~ "Don't Know")) |>
  mutate(Cattle = ifelse(Cattle == 0, "Yes", "No"), 
         SheepGoat = ifelse(SheepGoat == 1, "Yes", "No"), 
         Chicken = ifelse(Chicken == 2, "Yes", "No"), 
         Pet = ifelse(Pet == 3, "Yes", "No"), 
         Nodomanimal = ifelse(Nodomanimal == 4, "Yes", "No")) |>
  mutate(Cattle = ifelse(is.na(Cattle), "No", Cattle), 
         SheepGoat = ifelse(is.na(SheepGoat), "No", SheepGoat), 
         Chicken = ifelse(is.na(Chicken), "No", Chicken), 
         Pet = ifelse(is.na(Pet), "No", Pet), 
         Nodomanimal = ifelse(is.na(Nodomanimal), "No", Nodomanimal)) |>
  mutate(Doyouhavepotablewaterinyourhouse = case_when(Doyouhavepotablewaterinyourhouse == 0 ~ "No", 
                                                      Doyouhavepotablewaterinyourhouse == 1 ~ "Yes", 
                                                      is.na(Doyouhavepotablewaterinyourhouse) ~ "Don't Know")) |>
  mutate(Ifnotwheredoyougetyordrinkingwaterfrom = case_when(Ifnotwheredoyougetyordrinkingwaterfrom == 0 ~ "Neighbor", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 1 ~ "River", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 2 ~ "Well", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 3 ~ "Truck",
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 4 ~ "Tank", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 5 ~ "Public Fountain",
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == "Tap water" ~ "Tap Water", 
                                                            Ifnotwheredoyougetyordrinkingwaterfrom == 10 ~ NA)) |>
  mutate(Doyoudrinkdirectlyordoyoutreat = case_when(Doyoudrinkdirectlyordoyoutreat == 1 ~ "Directly", 
                                                    Doyoudrinkdirectlyordoyoutreat == 0 ~ "Treat", 
                                                    is.na(Doyoudrinkdirectlyordoyoutreat) ~ "Don't Know")) |>
  mutate(Ifyoutreathow = case_when(Ifyoutreathow == 0 ~ "Boil", 
                                   Ifyoutreathow == 1 ~ "Chemical", 
                                   Ifyoutreathow == 2 ~ "Filter", 
                                   Ifyoutreathow == 99 ~ "Don't Know")) |>
  mutate(Doyourfamilyownlatrine = ifelse(Doyourfamilyownlatrine == 1, "Yes", "No")) |>
  mutate(Isyourlatrineinsideoroutside = ifelse(Isyourlatrineinsideoroutside == 1, "Inside", "Outside")) |>
  mutate(Ifoutsidedistancefromyourhouse = as.numeric(str_sub(Ifoutsidedistancefromyourhouse, 1, -2))) |>
  mutate(Distancebetweenlatrineandyourkitchen = as.numeric(str_sub(Distancebetweenlatrineandyourkitchen, 1, -2))) |>
  mutate(Yourlatrineconnectedto = case_when(Yourlatrineconnectedto == 0 ~ "Sewage", 
                                            Yourlatrineconnectedto == 1 ~ "Ditch", 
                                            Yourlatrineconnectedto == 2 ~ "River", 
                                            Yourlatrineconnectedto == 3 ~ "Well", 
                                            Yourlatrineconnectedto == 99 | is.na(Yourlatrineconnectedto) ~ "NA/Don't Know")) |>
  mutate(Howoftendoyoubathinriver = case_when(Howoftendoyoubathinriver == 0 ~ "Always", 
                                              Howoftendoyoubathinriver == 1 ~ "Sometimes", 
                                              Howoftendoyoubathinriver == 2 ~ "Never",
                                              Howoftendoyoubathinriver == 99 | is.na(Howoftendoyoubathinriver) ~ "NA/Don't Know"), 
         Howoftendoyoubathinriver = factor(Howoftendoyoubathinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howoftendoyouwashyourclothesinriver = case_when(Howoftendoyouwashyourclothesinriver == 0 ~ "Always", 
                                                         Howoftendoyouwashyourclothesinriver == 1 ~ "Sometimes", 
                                                         Howoftendoyouwashyourclothesinriver == 2 ~ "Never",
                                                         Howoftendoyouwashyourclothesinriver == 99 | is.na(Howoftendoyouwashyourclothesinriver) ~ "NA/Don't Know"), 
         Howoftendoyouwashyourclothesinriver = factor(Howoftendoyouwashyourclothesinriver, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyoudeficateintheopenfield = case_when(Doyoudeficateintheopenfield == 0 ~ "Always", 
                                                 Doyoudeficateintheopenfield == 1 ~ "Sometimes", 
                                                 Doyoudeficateintheopenfield == 2 ~ "Never",
                                                 Doyoudeficateintheopenfield == 99 | is.na(Doyoudeficateintheopenfield) ~ "NA/Don't Know"), 
         Doyoudeficateintheopenfield = factor(Doyoudeficateintheopenfield, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouuseschoollatrine = case_when(Doyouuseschoollatrine == 0 ~ "Always", 
                                           Doyouuseschoollatrine == 1 ~ "Sometimes", 
                                           Doyouuseschoollatrine == 2 ~ "Never",
                                           Doyouuseschoollatrine == 99 | is.na(Doyouuseschoollatrine) ~ "NA/Don't Know"), 
         Doyouuseschoollatrine = factor(Doyouuseschoollatrine, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouusetoiletpapertowipeyourbumafterdefication = case_when(Doyouusetoiletpapertowipeyourbumafterdefication == 0 ~ "Always", 
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 1 ~ "Sometimes", 
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 2 ~ "Never",
                                                                     Doyouusetoiletpapertowipeyourbumafterdefication == 99 | is.na(Doyouusetoiletpapertowipeyourbumafterdefication) ~ "NA/Don't Know"), 
         Doyouusetoiletpapertowipeyourbumafterdefication = factor(Doyouusetoiletpapertowipeyourbumafterdefication, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouwashyourhandsaftertoilet = case_when(Doyouwashyourhandsaftertoilet == 0 ~ "Always", 
                                                   Doyouwashyourhandsaftertoilet == 1 ~ "Sometimes", 
                                                   Doyouwashyourhandsaftertoilet == 2 ~ "Never",
                                                   Doyouwashyourhandsaftertoilet == 99 | is.na(Doyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
         Doyouwashyourhandsaftertoilet = factor(Doyouwashyourhandsaftertoilet, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howdoyouwashyourhandsaftertoilet = case_when(Howdoyouwashyourhandsaftertoilet == 0 ~ "Water", 
                                                      Howdoyouwashyourhandsaftertoilet == 1 ~ "Soap and Water",
                                                      Howdoyouwashyourhandsaftertoilet == 99 | is.na(Howdoyouwashyourhandsaftertoilet) ~ "NA/Don't Know"), 
         Howdoyouwashyourhandsaftertoilet = factor(Howdoyouwashyourhandsaftertoilet, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  mutate(Ifwithsoaphowoftenyouuseit = case_when(Ifwithsoaphowoftenyouuseit == 0 ~ "Always", 
                                                Ifwithsoaphowoftenyouuseit == 1 ~ "Sometimes",
                                                Ifwithsoaphowoftenyouuseit == 2 ~ "Never",
                                                Ifwithsoaphowoftenyouuseit == 99 | is.na(Ifwithsoaphowoftenyouuseit) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit = factor(Ifwithsoaphowoftenyouuseit, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  mutate(Doyouwashyourhandsbeforeeating = case_when(Doyouwashyourhandsbeforeeating == 0 ~ "Always", 
                                                    Doyouwashyourhandsbeforeeating == 1 ~ "Sometimes", 
                                                    Doyouwashyourhandsbeforeeating == 2 ~ "Never",
                                                    Doyouwashyourhandsbeforeeating == 99 | is.na(Doyouwashyourhandsbeforeeating) ~ "NA/Don't Know"), 
         Doyouwashyourhandsbeforeeating = factor(Doyouwashyourhandsbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Howdoyouwashyourhandsbeforeeating = case_when(Howdoyouwashyourhandsbeforeeating == 0 ~ "Water", 
                                                       Howdoyouwashyourhandsbeforeeating == 1 ~ "Soap and Water",
                                                       Howdoyouwashyourhandsbeforeeating == 99 | is.na(Howdoyouwashyourhandsbeforeeating) | Howdoyouwashyourhandsbeforeeating == 2 ~ "NA/Don't Know"), 
         Howdoyouwashyourhandsbeforeeating = factor(Howdoyouwashyourhandsbeforeeating, levels = c("Water", "Soap and Water", "NA/Don't Know"))) |>
  mutate(Ifwithsoaphowoftenyouuseit_A = case_when(Ifwithsoaphowoftenyouuseit_A == 0 ~ "Always", 
                                                  Ifwithsoaphowoftenyouuseit_A == 1 ~ "Sometimes",
                                                  Ifwithsoaphowoftenyouuseit_A == 2 ~ "Never",
                                                  Ifwithsoaphowoftenyouuseit_A == 99 | is.na(Ifwithsoaphowoftenyouuseit_A) ~ "NA/Don't Know"), 
         Ifwithsoaphowoftenyouuseit_A = factor(Ifwithsoaphowoftenyouuseit_A, levels = c("Never", "Sometimes","Always", "NA/Don't Know"))) |>
  mutate(Doyoueatsoil = case_when(Doyoueatsoil == 0 ~ "Always", 
                                  Doyoueatsoil == 1 ~ "Sometimes", 
                                  Doyoueatsoil == 2 ~ "Never",
                                  Doyoueatsoil == 99 | is.na(Doyoueatsoil) ~ "NA/Don't Know"), 
         Doyoueatsoil = factor(Doyoueatsoil, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Whatisyourfavoritefruitethatyoueat = ifelse(Whatisyourfavoritefruitethatyoueat != "Mango" & Whatisyourfavoritefruitethatyoueat != "Banana" & Whatisyourfavoritefruitethatyoueat != "Avocado" & Whatisyourfavoritefruitethatyoueat != "Apple" & Whatisyourfavoritefruitethatyoueat != "Pinapple", "Other", Whatisyourfavoritefruitethatyoueat), 
         Whatisyourfavoritefruitethatyoueat = ifelse(is.na(whatfooddoyoueatmost), "Other", Whatisyourfavoritefruitethatyoueat)) |>
  mutate(Doyouwashyourfruitesbeforeeating = case_when(Doyouwashyourfruitesbeforeeating == 0 ~ "Always", 
                                                      Doyouwashyourfruitesbeforeeating == 1 ~ "Sometimes", 
                                                      Doyouwashyourfruitesbeforeeating == 2 ~ "Never",
                                                      Doyouwashyourfruitesbeforeeating == 99 | is.na(Doyouwashyourfruitesbeforeeating) ~ "NA/Don't Know"), 
         Doyouwashyourfruitesbeforeeating = factor(Doyouwashyourfruitesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyoueatraworundercookedvegitables = case_when(Doyoueatraworundercookedvegitables == 0 ~ "Always", 
                                                        Doyoueatraworundercookedvegitables == 1 ~ "Sometimes", 
                                                        Doyoueatraworundercookedvegitables == 2 ~ "Never",
                                                        Doyoueatraworundercookedvegitables == 99 | is.na(Doyoueatraworundercookedvegitables) ~ "NA/Don't Know"), 
         Doyoueatraworundercookedvegitables = factor(Doyoueatraworundercookedvegitables, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Ifyeshowoftenyouwashyourvegitablesbeforeeating = case_when(Ifyeshowoftenyouwashyourvegitablesbeforeeating == 0 ~ "Always", 
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 1 ~ "Sometimes", 
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 2 ~ "Never",
                                                                    Ifyeshowoftenyouwashyourvegitablesbeforeeating == 99 | is.na(Ifyeshowoftenyouwashyourvegitablesbeforeeating) ~ "NA/Don't Know"), 
         Ifyeshowoftenyouwashyourvegitablesbeforeeating = factor(Ifyeshowoftenyouwashyourvegitablesbeforeeating, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(Doyouwalkbarefoot = case_when(Doyouwalkbarefoot == 0 ~ "Always", 
                                       Doyouwalkbarefoot == 1 ~ "Sometimes", 
                                       Doyouwalkbarefoot == 2 ~ "Never",
                                       Doyouwalkbarefoot == 99 | is.na(Doyouwalkbarefoot) ~ "NA/Don't Know"), 
         Doyouwalkbarefoot = factor(Doyouwalkbarefoot, levels = c("Never", "Sometimes", "Always", "NA/Don't Know"))) |>
  mutate(whenyouareathomedoyoupreferetousesandalsorshoe = case_when(whenyouareathomedoyoupreferetousesandalsorshoe == 0 ~ "Barefoot", 
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 1 ~ "Sandals", 
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 2 ~ "Shoes",
                                                                    whenyouareathomedoyoupreferetousesandalsorshoe == 99 | is.na(whenyouareathomedoyoupreferetousesandalsorshoe) ~ "NA/Don't Know"), 
         whenyouareathomedoyoupreferetousesandalsorshoe = factor(whenyouareathomedoyoupreferetousesandalsorshoe, levels = c("Barefoot", "Sandals", "Shoes", "NA/Don't Know"))) |>
  mutate(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill = case_when(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 0 ~ "No", 
                                                                                      Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 1 ~ "Yes",
                                                                                      Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill == 99 | is.na(Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill) ~ "NA/Don't Know"))|>
  mutate(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill = case_when(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "1Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "2Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "3month" ~ "1-3 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "5Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "6Month" ~ "4-6 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "7Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "9Month" ~ "7-9 Months", 
                                                                         Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "11Month" | Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill == "12month" ~ ">9 Months", 
                                                                         is.na(Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill) ~ NA)) |>
  mutate(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics = case_when(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 0 ~ "No", 
                                                                                Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 1 ~ "Yes",
                                                                                Didyourparentsorhealthprofessionalsgaveyouotherantibiotics == 99 | is.na(Didyourparentsorhealthprofessionalsgaveyouotherantibiotics) ~ "NA/Don't Know")) |>
  mutate(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently = case_when(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 0 ~ "No", 
                                                                                  Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 1 ~ "Yes",
                                                                                  Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently == 99 | is.na(Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently) ~ "NA/Don't Know")) |>
  mutate(Inthepastthreemonthshaveyoutakenanyantimalarialdrug = case_when(Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 0 ~ "No", 
                                                                         Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 1 ~ "Yes",
                                                                         Inthepastthreemonthshaveyoutakenanyantimalarialdrug == 99 | is.na(Inthepastthreemonthshaveyoutakenanyantimalarialdrug) ~ "NA/Don't Know")) |>
  mutate(Haveachildwheezingorwhistlinginchest = case_when(Haveachildwheezingorwhistlinginchest == 2 ~ "No", 
                                                          Haveachildwheezingorwhistlinginchest == 1 ~ "Yes",
                                                          Haveachildwheezingorwhistlinginchest == 99 | is.na(Haveachildwheezingorwhistlinginchest) ~ "NA/Don't Know")) |>
  mutate(Howmanytimesinthelastyearthechildhadwheeziling = case_when(Howmanytimesinthelastyearthechildhadwheeziling == 2 ~ "1-3", 
                                                                    Howmanytimesinthelastyearthechildhadwheeziling == 1 ~ "0",
                                                                    Howmanytimesinthelastyearthechildhadwheeziling == 99 | is.na(Howmanytimesinthelastyearthechildhadwheeziling) ~ "NA/Don't Know")) |>
  mutate(Hasachildeverhadasthma = case_when(Hasachildeverhadasthma == 2 ~ "No", 
                                            Hasachildeverhadasthma == 1 ~ "Yes",
                                            Hasachildeverhadasthma == 99 | is.na(Hasachildeverhadasthma) ~ "NA/Don't Know")) |>
  mutate(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases = case_when(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 2 ~ "No", 
                                                                                    Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 1 ~ "Yes",
                                                                                    Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases == 99 | is.na(Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases) ~ "NA/Don't Know")) |>
  mutate(Ifyeshasthisrashaffectedtheelbowfolds = ifelse(is.na(Ifyeshasthisrashaffectedtheelbowfolds), 0, Ifyeshasthisrashaffectedtheelbowfolds), 
         Ifyeshasthisrashaffectedbehindtheknees = ifelse(is.na(Ifyeshasthisrashaffectedbehindtheknees), 0, Ifyeshasthisrashaffectedbehindtheknees), 
         Ifyeshasthisrashaffectedinfrontoftheankls = ifelse(is.na(Ifyeshasthisrashaffectedinfrontoftheankls), 0, Ifyeshasthisrashaffectedinfrontoftheankls), 
         Ifyeshasthisrashaffectedunderthebuttucks = ifelse(is.na(Ifyeshasthisrashaffectedunderthebuttucks), 0, Ifyeshasthisrashaffectedunderthebuttucks), 
         Ifyeshasthisrashaffectedaroundtheneck = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheneck), 0, Ifyeshasthisrashaffectedaroundtheneck), 
         Ifyeshasthisrashaffectedaroundtheeyesears = ifelse(is.na(Ifyeshasthisrashaffectedaroundtheeyesears), 0, Ifyeshasthisrashaffectedaroundtheeyesears)) |>
  dplyr::select(-AL, -TT, -HW, -SM, -Others, -Date, -parentconsent, -childassent, -Grade, -GradeLet, -WAZ, -HAZ, -BAZ, -Wheredoyoulive, -Doesachildasthmainthelastyear, -Doesachildhadasthmainthelast2years, -Hasthisbeenconfirmedbyadoctor, -Doesthechildhadhayfeverorsneezingwithrunningnoseinthelast2years, -Doesthechildhadhayfeverorsneezingwithrunningnoseinthelastyear, -V137, -AgeGroup, -HouseholdGroup) |>
  mutate(Doesthechildhadhayfeverorpersistantsneezing = case_when(Doesthechildhadhayfeverorpersistantsneezing == 2 ~ "No", 
                                                                 Doesthechildhadhayfeverorpersistantsneezing == 1 ~ "Yes",
                                                                 Doesthechildhadhayfeverorpersistantsneezing == 99 | is.na(Doesthechildhadhayfeverorpersistantsneezing) ~ "NA/Don't Know")) |>
  mutate(Arechildsfingernailsdirty = ifelse(Arechildsfingernailsdirty == 1, "Yes", "No")) |>
  rename(Kerosene = Kerosisn) |>
  dplyr::select(-ID, -Seq_ID)




sample_data(dat)  <- cleaned

write_rds(dat, "categorized_data.rds")  



  
  
  
