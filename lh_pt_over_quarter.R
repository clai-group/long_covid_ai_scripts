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
cov_pats_file <- file.choose() ## cov_pats.RData under output folder
outputDirectory <- choose_directory(caption = "select output data directory")

# Load data ------
load(longhaulers_file)
load(cov_pats_file)

all_quarters <- tibble(
  year = c(rep(2020:2024, each = 4), rep(2025, 3)),
  quarter = c(rep(1:4, 5), 1:3)
) %>%
  mutate(
    quarter_label = paste0(year, "-Q", quarter),
    quarter_start = ymd(paste0(year, "-", (quarter - 1) * 3 + 1, "-01")),
    quarter_end = quarter_start + months(3) - days(1)
  )

cov_pt_counts <- all_quarters %>%
  left_join(
    cov_pats %>%
      crossing(all_quarters) %>%
      filter(as.Date(start_date) <= quarter_end) %>%
      group_by(quarter_label) %>%
      summarise(
        cov_pt_count = n_distinct(patient_num),
        .groups = "drop"
      ),
    by = "quarter_label"
  ) %>%
  mutate(cov_pt_count = replace_na(cov_pt_count, 0)) %>%
  select(quarter_label, cov_pt_count)

lh_pt_counts <- all_quarters %>%
  left_join(
    longhaulers_final %>%
      crossing(all_quarters) %>%
      filter(as.Date(cov_date) <= quarter_end) %>%
      group_by(quarter_label) %>%
      summarise(
        lh_pt_count = n_distinct(patient_num),
        .groups = "drop"
      ),
    by = "quarter_label"
  ) %>%
  mutate(lh_pt_count = replace_na(lh_pt_count, 0)) %>%
  select(quarter_label, lh_pt_count)

pt_count <- cov_pt_counts %>%
  left_join(lh_pt_counts, by = "quarter_label") %>%
  mutate(percentage = round((lh_pt_count / cov_pt_count) * 100, 3))

write.csv(pt_count, file= paste0(outputDirectory,"/pt_counts_overquarter_",site, ".csv"), row.names = FALSE)

df <- longhaulers_final %>% 
  ungroup()%>%
  select("organ", "combo","duration", "subcombo", "patient_num") %>%
  distinct()

subcombo <- df %>% group_by(organ, combo, subcombo) %>%
  summarise(count = n())

write.csv(subcombo, file= paste0(outputDirectory,"/subcombo_count_",site, ".csv"), row.names = FALSE)
