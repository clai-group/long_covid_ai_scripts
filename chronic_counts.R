library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(data.table)

choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

# Input data ------
site <- "" ## UCLA, UTH
cohort <- "longhaulers"
longhaulers_final_file <- file.choose()  ## longhaulers_final_site.RData under output folder
outputDirectory <- choose_directory(caption = "select output data directory")

# Load data ------
load(longhaulers_file)

chronic_sizes <- longhaulers_final %>%
  group_by(Chronic_Status) %>%
  summarise(count_pt = n_distinct(patient_num), 
            count_condition = n(),
            .groups = "drop"  )%>%  
  bind_rows(tibble(Chronic_Status = "Total_from_lh", 
                     count_pt = n_distinct(longhaulers_final$patient_num),
                     count_condition = nrow(longhaulers_final)))

write.csv(chronic_sizes, file= paste0(outputDirectory,"/chronic_lh_ptcount_",site, ".csv"), row.names = FALSE)
