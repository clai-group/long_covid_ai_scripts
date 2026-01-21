library(readr)
library(dplyr)
library(tidyr)
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
longhaulers_file <- file.choose()  ## longhaulers_site.RData under output folder
cov_pats_file <- file.choose() ## cov_pats.RData under output folder
dems_file_path <- file.choose() ## dems_cases.csv under data/cases folder
icds_map_file <- file.choose() ## icds_map_likelihood.csv
ref_combo_organ_file <- file.choose() ## combo_update_202511.csv
outputDirectory <- choose_directory(caption = "select output data directory")

# Load data ------
load(longhaulers_file)
load(cov_pats_file)
dems_file <-  data.table::fread(dems_file_path)
icds_map <- data.table::fread(icds_map_file)
ref_combo_organ <-  data.table::fread(ref_combo_organ_file)


if(site == "UTH"){
  #### For UTH
  map_lh <- longhaulers %>%
    rename(concept_cd_old= concept_cd) %>%
    mutate(concept_cd = str_replace(concept_cd_old, "^ICD10CM:", "ICD10:")) %>%
    left_join(icds_map, by = c("concept_cd"))
}else if (site == "UCLA"){
  #### For UCLA
  map_lh <- longhaulers %>%
    rename(concept_cd_old= concept_cd) %>%
    mutate(concept_cd = str_extract(concept_cd_old, "(?<=\\|).*")) %>%
    left_join(icds_map, by = c("concept_cd"))
}else{
  ## For MGB
  map_lh <- longhaulers %>%
    left_join(icds_map, by = c("concept_cd"))
}

map_lh_filter <- map_lh[map_lh$Long_COVID_Likelihood_Speculative!='Unlikely',]

organ_ln <- map_lh_filter %>%
  select(-c(organ, combo, subcombo, `Unnamed: 0`)) %>%
  left_join(ref_combo_organ, by = c("phenx", "icd10_desc"))

longhaulers_final <- organ_ln[!is.na(organ_ln$organ), ]
save(longhaulers_final, file= paste0(outputDirectory,"/longhaulers_final_",site, ".RData"))

summary <- function(dems_all, df) {
  dems_all <- dems_all[dems_all$patient_num %in% df$patient_num, ]
  dems <- dems_all[, c("patient_num", "age", "sex_cd", "CHARLSON_INDEX", "race_cd", "ethnicity_cd")]
  rm(dems_all)
  
  # If the patient has more than 1 ages, select the oldest.
  dems <- unique(dems) %>% group_by(patient_num) %>% 
    filter(age == max(age)) %>% 
    ungroup()
  
  # If the patient has more than 1 CHARLSON INDEX, select the largest.
  dems <- unique(dems) %>% group_by(patient_num) %>% 
    filter(CHARLSON_INDEX == max(CHARLSON_INDEX)) %>% 
    ungroup()
  
  # Prepare df (lowercase column names)
  names(df) <- tolower(names(df))
  dat.sum <- subset(df, df$patient_num %in% dems$patient_num)
  colnames(dat.sum) <- tolower(colnames(dat.sum))
  
  # Ensure dems only contains patients in df
  dems <- subset(dems, dems$patient_num %in% df$patient_num)
  
  # Create summary statistics
  dems_temp1 <- data.frame(rbind(c("patients", nrow(dems))))
  dems_temp2 <- data.frame(rbind(c("mean age", mean(dems$age))))
  dems_sd <- data.frame(rbind(c("sd age", sd(dems$age))))
  dems_temp3 <- data.frame(rbind(c("percent female", (nrow(subset(dems, dems$sex_cd == "F"))/nrow(dems))*100)))
  dems_temp4 <- data.frame(rbind(c("mean charlson", mean(dems$CHARLSON_INDEX))))
  
  tb1 <- rbind(dems_temp1, dems_temp2, dems_sd, dems_temp3, dems_temp4) 
  tb2 <- data.frame(rbind(c("unique CCSR concepts", length(unique(dat.sum$phenx)))))
  rm(dems_temp1, dems_temp2, dems_temp3, dems_temp4)
  
  
  table1 <- rbind(tb1, tb2)
  table1$site <- site
  table1$cohort <- cohort
  rm(tb1, tb2)
  
  # Race and ethnicity statistics
  race <- data.frame(table(dems$race_cd))
  race <- subset(race, race$Freq > 30)
  race$site <- site
  race$cohort <- cohort 
  
  eth <- data.frame(table(dems$ethnicity_cd))
  eth$site <- site
  eth$cohort <- cohort
  
  return(list(dems_stat = table1, race_stat = race, eth_stat = eth)) 
}

# Output Longhaulers summary statistics ---------------
stat <- summary(dems_all = dems_file, 
                df = longhaulers_final)

lapply(1:length(stat), function(i) write.csv(stat[[i]], 
                                             file = paste0(outputDirectory, "/", cohort, "_", names(stat[i]), "_", site, ".csv"),
                                             row.names = FALSE))

# Output Longhaulers organ combo count ---------------
df <- longhaulers_final %>% 
  ungroup()%>%
  select("organ", "combo","duration", "subcombo", "patient_num")

duration_count = df %>% group_by(duration, organ, combo) %>%
  summarise(count = n())

duration_organ_count = df %>% group_by(organ) %>%
  summarise(count_organ = n())

duration_combo_count = df %>% group_by(organ, combo) %>%
  summarise(count_combo = n())

subcombo <- df %>% group_by(organ, combo, subcombo) %>%
  summarise(count = n())

# save(duration_count, file= paste0(outputDirectory,"/longhaulers_duration_organ_combo_count_",site, ".RData"))
# save(duration_organ_count, file= paste0(outputDirectory,"/longhaulers_organ_count_",site, ".RData"))
# save(duration_combo_count, file= paste0(outputDirectory,"/longhaulers_organ_combo_count_",site, ".RData"))
write.csv(duration_count, file= paste0(outputDirectory,"/longhaulers_duration_organ_combo_count_",site, ".csv"), row.names = FALSE)
write.csv(duration_organ_count, file= paste0(outputDirectory,"/longhaulers_organ_count_",site, ".csv"), row.names = FALSE)
write.csv(duration_combo_count, file= paste0(outputDirectory,"/longhaulers_organ_combo_count_",site, ".csv"), row.names = FALSE)
write.csv(subcombo, file= paste0(outputDirectory,"/longhaulers_subcombo_count_",site, ".csv"), row.names = FALSE)


# Output Prevalence over quarters ---------------
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
      filter(start_date <= quarter_end) %>%
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
      filter(as.Date(cov_date)) <= quarter_end) %>%
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

# save(pt_count, file= paste0(outputDirectory,"/pt_counts_overquarter_",site, ".RData"))
write.csv(pt_count, file= paste0(outputDirectory,"/pt_counts_overquarter_",site, ".csv"), row.names = FALSE)

# Output chronic counts ---------------
chronic_sizes <- longhaulers_final %>%
  group_by(Chronic_Status) %>%
  summarise(count_pt = n_distinct(patient_num), 
            count_condition = n(),
            .groups = "drop"  )%>%  
  bind_rows(tibble(Chronic_Status = "Total_from_lh", 
                     count_pt = n_distinct(longhaulers_final$patient_num),
                     count_condition = nrow(longhaulers_final)))

write.csv(chronic_sizes, file= paste0(outputDirectory,"/chronic_lh_ptcount_",site, ".csv"), row.names = FALSE)
