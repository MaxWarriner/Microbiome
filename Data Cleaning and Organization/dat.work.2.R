library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(stringr)
library(vegan)
library(readr)


dat <- readRDS("categorized_data.RDS")

sam <- sample_data(dat)

sam <- data.frame(sam)

sam <- sam |>
  rename(Weight = weight, 
         "Mode of Delivery" = Modeofdelivery,
         "Vaccinated" = Everhadvaccinated, 
         "BCG Scar" = BCGscar, 
         "Fever in Last Two Weeks" = Inthelasttwoweekschildhadfever, 
         "Number of Fevers in Last Two Weeks" = howmanytimesperweek, 
         "Diarrhea in Last Two Weeks" = InthelasttwoweekschildhadDiahearria, 
         "Number of Diarrhea in Last Two Weeks" = Howmanytimesperweek_A, 
         "Cough in Last Two Weeks" = Inthelasttwoweekschildhadcough, 
         "Number of Cough in Last Two Weeks" = Howmanytimesperweek_B, 
         "Child's Fingernails Trimmed" = Childsfingernailtrimmed, 
         "Child Has Dirty Fingernails" = Arechildsfingernailsdirty, 
         "How often do you trim your fingernails?" = Howoftendoyoutrimyourfingernails, 
         "Toilet in School" = Istheretoiletintheschool, 
         "Does Toilet in School Have Doors?" = Ifthereisatoiletdoesthelatrinehavedoors, 
         "Flies Observed Around Latrine" = Fliesobservedinaroundthelatrine, 
         "Visible Stool on Latrine Floor" = Vissiblestoolobservedonlatrinefloor, 
         "Heard of Ascharis" = HeardALnamebefore, 
         "Heard of Trichuris" = HeardTTnamebefore, 
         "Heard of Hookworm" = HeardHWnamebefore, 
         "Heard of HIV" = HeardHIVnamebefore, 
         "Heard of Intestinal Worms" = HeardInWormnamebefore, 
         "Heard of Malaria" = HeardMalanamebefore, 
         "Heard of Tuberculosis" = HeardTBnamebefore, 
         "Heard of Schistosoma" = HeardSChnamebefore, 
         "Told by Family" = familytold, 
         "Told by Health Professional" = HPtold, 
         "Told by Teacher" = Teachertold, 
         "Told by Media" = Mediatold, 
         "Knowledge of Intestinal Worm Transmission" = Doyouknowhowintestinalwomstransmitted, 
         "How Are Intestinal Worms Transmitted?" = IfyesHow, 
         "Knowledge of Health Effects of Worms" = Doyouknowwhywormsarebadforyourhealth, 
         "How are Worms Bad for Your Health?" = IfyesHow_A, 
         "Knowledge of How To Avoid Intestinal Worms" = Doyouknowhowyoucanavoidgettingtheseworms, 
         "How Can You Avoid Getting Intestinal Worms?" = Ifyeshow_B, 
         "Living Address" = Yourlivingaddress, 
         "Family Occupation" = Familyoccupation, 
         "Maternal Education Status" = Maternaleducationalstatus, 
         "House Floor Material" = Whatmaterialyourhousefloormadefrom, 
         "Kitchen Separated/In House" = Isyourkitcheninyourhouseorsepatated, 
         "Kitchen Material" = Ifseparetedwhatmaterialsyourkitchenmadefrom, 
         "Kitchen Has Roof" = Roof, 
         "Kitchen Has Wall" = Wall, 
         "Kitchen Has No Roof or Wall" = None, 
         "Wood as Cooking Fuel" = Wood, 
         "Gas as Cooking Fuel" = Gas, 
         "Coal as Cooking Fuel" = Coal, 
         "Kerosene as Cooking Fuel" = Kerosene, 
         "Electricity for Cooking" = Electric, 
         "Electricity in House" = Doyouhaveelectricityinyourhaouse, 
         "Family Owns Radio" = Doesyourfamillyownradio, 
         "Family Owns Television" = Doesyourfamillyowntelevision, 
         "Family Member with Phone" = Doesyourfamillymemberownaphone, 
         "Purpose of Phone Use" = Forwhatpurposeyourfamilyuseaphone, 
         "Cattle in House/Compound" = Cattle, 
         "Sheep or Goat in House/Compound" = SheepGoat, 
         "Chicken in House/Compound" = Chicken, 
         "Pet in House/Compound" = Pet, 
         "No Animals in House/Compound" = Nodomanimal, 
         "Potable Water in House" = Doyouhavepotablewaterinyourhouse, 
         "Source of Water if not House" = Ifnotwheredoyougetyordrinkingwaterfrom, 
         "Drink Water Directly or Treated?" = Doyoudrinkdirectlyordoyoutreat, 
         "Method of Treating Water" = Ifyoutreathow, 
         "Family Owns Latrine" = Doyourfamilyownlatrine, 
         "Latrine Inside or Outside" = Isyourlatrineinsideoroutside, 
         "Latrine Distance from House" = Ifoutsidedistancefromyourhouse, 
         "Distance Between Latrine and Kitchen" = Distancebetweenlatrineandyourkitchen, 
         "Connection to Latrine" = Yourlatrineconnectedto, 
         "Frequency of Bathing in River" = Howoftendoyoubathinriver, 
         "Frequency of Clothes Washing in River" = Howoftendoyouwashyourclothesinriver, 
         "Defecating in Open Field" = Doyoudeficateintheopenfield, 
         "Frequency of Using School Latrine" = Doyouuseschoollatrine, 
         "Use of Toilet Paper after Defecation" = Doyouusetoiletpapertowipeyourbumafterdefication, 
         "Frequency of Hand Washing After Using Toilet" = Doyouwashyourhandsaftertoilet, 
         "Frequency of Using Soap After Using Toilet" = Ifwithsoaphowoftenyouuseit, 
         "Frequency of Hand Washing Before Eating" = Doyouwashyourfruitesbeforeeating, 
         "Method of Washing Hands Before Eating" = Howdoyouwashyourhandsbeforeeating,
         "Frequency of Using Soap Before Eating" = Ifwithsoaphowoftenyouuseit_A, 
         "Frequency of Eating Soil" = Doyoueatsoil, 
         "Favorite Fruit" = Whatisyourfavoritefruitethatyoueat, 
         "Frequency of Washing Fruit Before Eating" = Doyouwashyourfruitesbeforeeating, 
         "Frequency of Eating Raw/Undercooked Vegetables" = Doyoueatraworundercookedvegitables, 
         "Freqency of Washing Raw/Undercooked Vegetables" = Ifyeshowoftenyouwashyourvegitablesbeforeeating, 
         "Frequency of Walking Barefoot" = Doyouwalkbarefoot, 
         "Prefer Sandals or Barefoot in House" = whenyouareathomedoyoupreferetousesandalsorshoe, 
         "Activities Where Barefoot" = Inwhichactivitiesofthedayareyoubarefoot, 
         "Deworming Pill" = Didyourparentsteachersorhealthprofessionalsgaveyouadewormingpill, 
         "Last Deworming Pill" = Ifyeswhenwasthelasttimetheygaveyousuchdewormingpill, 
         "Antibiotics" = Didyourparentsorhealthprofessionalsgaveyouotherantibiotics, 
         "Most Eaten Food" = whatfooddoyoueatmost, 
         "Drug Prescribed For Illness Currently" = Isanydrugsprescribedbyhealthinstituionforanyillnesscurrently, 
         "Name and Type of Drug Currently Prescribed" = Ifyesnameandtypeofthedrugs, 
         "Anti-Malarial Drug" = Inthepastthreemonthshaveyoutakenanyantimalarialdrug, 
         "Wheezing or Whistling in Chest" = Haveachildwheezingorwhistlinginchest, 
         "Wheezing or Whistling in Last 2 Years" = wheezilingorwhistlinginthelast2years, 
         "Wheezing or Whistling in Last 1 Year" = wheezilingorwhistlinginthelast1year, 
         "How Many Times Wheezing or Whistling" = Howmanytimesinthelastyearthechildhadwheeziling, 
         "Ever Had Asthma" = Hasachildeverhadasthma, 
         "Rash Affecting Skin Creases" = Doesachildeverhadanitchyskinrashwhichhasaffectedtheskincreases, 
         "Rash Affecting Elbow Folds" = Ifyeshasthisrashaffectedtheelbowfolds, 
         "Rash Behind the Knees" = Ifyeshasthisrashaffectedbehindtheknees, 
         "Rash In Front of Ankles" = Ifyeshasthisrashaffectedinfrontoftheankls, 
         "Rash Under the Buttocks" = Ifyeshasthisrashaffectedunderthebuttucks, 
         "Rash Affecting the Neck" = Ifyeshasthisrashaffectedaroundtheneck, 
         "Rash Affecting Eyes or Ears" = Ifyeshasthisrashaffectedaroundtheeyesears, 
         "Hay Fever or Persistent Sneezing" = Doesthechildhadhayfeverorpersistantsneezing, 
         "Urban/Rural" = WheredoyouliveGROUP, 
         "Stunted Growth" = Stunted, 
         "Household Number" = Household.Number, 
         "Number of Older Siblings" = Older.Siblings, 
         "Number of Siblings Younger than 12" = Siblings.Younger.than.12)

clean_yes_no_columns <- function(df, yes_values = c("yes", "Yes", "YES"), 
                                 no_values = c("no", "No", "NO"), 
                                 ignore_case = TRUE) {
  # Iterate through each column of the data frame
  for (col_name in names(df)) {
    column <- df[[col_name]]
    
    # Skip if column is not character or factor
    if (!is.character(column)) {
      if (is.factor(column)) {
        column <- as.character(column)
      } else {
        next
      }
    }
    
    # Standardize case if ignore_case is TRUE
    if (ignore_case) {
      column_lower <- tolower(column)
      yes_values_lower <- tolower(yes_values)
      no_values_lower <- tolower(no_values)
      
      # Replace values
      clean_column <- ifelse(column_lower %in% yes_values_lower, "yes",
                             ifelse(column_lower %in% no_values_lower, "no", NA))
    } else {
      # Case-sensitive comparison
      clean_column <- ifelse(column %in% yes_values, "yes",
                             ifelse(column %in% no_values, "no", NA))
    }
    
    # Only update the column if we found at least some yes/no values
    if (any(clean_column %in% c("yes", "no"), na.rm = TRUE)) {
      df[[col_name]] <- clean_column
    }
  }
  
  return(df)
}

clean_frequency_columns <- function(df, 
                                    always_values = c("always", "Always", "ALWAYS"),
                                    sometimes_values = c("sometimes", "Sometimes", "SOMETIMES"),
                                    never_values = c("never", "Never", "NEVER"),
                                    ignore_case = TRUE) {
  
  # Iterate through each column of the data frame
  for (col_name in names(df)) {
    column <- df[[col_name]]
    
    # Skip if column is not character or factor
    if (!is.character(column)) {
      if (is.factor(column)) {
        column <- as.character(column)
      } else {
        next
      }
    }
    
    # Standardize case if ignore_case is TRUE
    if (ignore_case) {
      column_lower <- tolower(column)
      always_lower <- tolower(always_values)
      sometimes_lower <- tolower(sometimes_values)
      never_lower <- tolower(never_values)
      
      # Replace values
      clean_column <- ifelse(column_lower %in% always_lower, "always",
                             ifelse(column_lower %in% sometimes_lower, "sometimes",
                                    ifelse(column_lower %in% never_lower, "never", NA)))
    } else {
      # Case-sensitive comparison
      clean_column <- ifelse(column %in% always_values, "always",
                             ifelse(column %in% sometimes_values, "sometimes",
                                    ifelse(column %in% never_values, "never", NA)))
    }
    
    # Only update the column if we found at least some frequency values
    if (any(clean_column %in% c("always", "sometimes", "never"), na.rm = TRUE)) {
      df[[col_name]] <- clean_column
    }
  }
  
  return(df)
}

clean_dont_know_responses <- function(df, 
                                      dont_know_values = c("don't know", "do not know", "dont know", 
                                                           "Don't Know", "Do Not Know", "Dont Know",
                                                           "DON'T KNOW", "DO NOT KNOW", "DONT KNOW",
                                                           "unknown", "Unknown", "UNKNOWN",
                                                           "not sure", "Not Sure", "NOT SURE",
                                                           "unsure", "Unsure", "UNSURE",
                                                           "no answer", "No Answer", "NO ANSWER",
                                                           "refused", "Refused", "REFUSED"),
                                      ignore_case = TRUE) {
  
  # Iterate through each column of the data frame
  for (col_name in names(df)) {
    column <- df[[col_name]]
    
    # Skip if column is not character or factor
    if (!is.character(column)) {
      if (is.factor(column)) {
        column <- as.character(column)
      } else {
        next
      }
    }
    
    # Standardize case if ignore_case is TRUE
    if (ignore_case) {
      column_lower <- tolower(column)
      dont_know_lower <- tolower(dont_know_values)
      
      # Replace values
      clean_column <- ifelse(column_lower %in% dont_know_lower, NA, column)
    } else {
      # Case-sensitive comparison
      clean_column <- ifelse(column %in% dont_know_values, NA, column)
    }
    
    # Update the column
    df[[col_name]] <- clean_column
  }
  
  return(df)
}

sam <- clean_yes_no_columns(sam)

sam <- clean_frequency_columns(sam)

sam <- clean_dont_know_responses(sam)

clean_names_simple <- function(df) {
  colnames(df) <- gsub(" ", "_", colnames(df))
  return(df)
}

sam <- clean_names_simple(sam)


sample_data(dat) <- sam

write_rds(dat, "categorized_data.RDS")

