###this script maps the phenx back to ICD10 description
###adds organ, modified CCSR category, and clinical problem for results analysis.

library(readr)
library(dplyr)
library(lubridate)
library(data.table)

choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

site <- "MGB"
longhaulers_file <- file.choose()  ##the longCOVID_patients_site_ref_thresholds0.05.csv file
cov_pat_incident_FileName <- file.choose() ##the cov_pats.RData file
dbmartCases_FileName <-   file.choose() ##CCSR-mapped cases
ref_phenxlookup <- file.choose() ##the pre-computed delivered with the docker container -- ref_phenxlookup.RData
ref_combo_organ <- file.choose() ##the pre-computed delivered with the docker container -- combo_updated.csv
ccsr_icd_mapback <- file.choose() ##the pre-computed delivered with the docker container -- ccsr_icddesc_mapback.csv
outputDirectory <- choose_directory(caption = "select output data directory")


##load data 
longhaulers <- data.table::fread(longhaulers_file)
ref_combo_organ <-  data.table::fread(ref_combo_organ)
cases <- data.table::fread(dbmartCases_FileName)
ccsr <- data.table::fread(ccsr_icd_mapback)
load(cov_pat_incident_FileName)
load(ref_phenxlookup)



##map to find the phenx
colnames(phenxlookup) <- c("infection", "startPhen")
df_infect <- merge(longhaulers, phenxlookup, by ="startPhen", all.x = TRUE)
df_infect$infection_seq <- as.numeric(gsub("[^0-9.]+", "", df_infect$infection))

##map to find the cov_date
df_infect_covpts <- merge(df_infect,cov_pats, by = c("patient_num", "infection_seq"), all.x=TRUE)
rm(df_infect)

##find the ICD10 description
colnames(cases) <- tolower(colnames(cases))
# cases <- select(cases, -c_fullname)
# names(cases)[names(cases) == 'concept_path'] <- 'c_fullname'
cases$patient_num <- as.integer(cases$patient_num)
colnames(cases)[colnames(cases) == "start_date"] <- "phenx_date"
cases$phenx_date <- as.Date(cases$phenx_date, format="%Y-%m-%d")
cases_s <- merge(cases, ccsr, by = c("c_fullname"), all.x = TRUE)
rm(cases)
# 
# ccsr <-  read_csv("P:/PASC/scripts_on_git/CCSR_PASC_ACT_Mapping_022024.csv")
# colnames(ccsr) <- tolower(colnames(ccsr))
# ccsr <- ccsr[,c("concept_path", "ccsr_key", "phenx", "icd10_desc")]
# names(ccsr)[names(ccsr) == 'concept_path'] <- 'c_fullname'
# cases_s <- merge(cases, ccsr, by = c("c_fullname", "ccsr_key", "phenx"), all.x = TRUE)


colnames(df_infect_covpts)[colnames(df_infect_covpts) == "start_date"] <- "cov_date"
map <- merge(df_infect_covpts, cases_s, by = c("patient_num", "phenx"), all.x=TRUE)

# 0 = [0,14], 1= (14, 30], 2=(30,60],3=(60,90], 4=(90,120], 5=(120,150]
map_date <- map %>% filter(
  (duration == 0 & phenx_date >= cov_date & phenx_date <= cov_date + days(14)) |
    (duration == 1 & phenx_date > cov_date + days(14) &  phenx_date <= cov_date + days(30)) |
    (duration >= 2 & phenx_date > cov_date + days(30 * (duration - 1)) & phenx_date <= cov_date + days(30 * duration))
) %>% distinct()

rm(map)

colnames(map_date) <- tolower(colnames(map_date))
map_date <- map_date[, c("patient_num", "phenx","phenx_date","duration","durationbucket", "infection", 
                         "concept_cd", "cov_date", "icd10_desc")]

##add organ, modified CCSR category, and clinical problem
colnames(ref_combo_organ) <- tolower(colnames(ref_combo_organ))
ref_combo_organ <- ref_combo_organ[, -c("count")]
map_co <- merge(map_date, ref_combo_organ, by = c("phenx", "icd10_desc"), all.x=TRUE) %>% distinct()
rm(map_date)

##remove duplicates by filtering the earliest dates
map_back <- map_co %>% 
  group_by(patient_num, duration,phenx, icd10_desc) %>% 
  filter(phenx_date == min(phenx_date)) 

##remove cancers
longhaulers <- subset(map_back, !(icd10_desc %in% c("Mammographic calcification found on diagnostic imaging of breast",
                                                   "Other abnormal and inconclusive findings on diagnostic imaging of breast",
                                                   "Inconclusive mammogram",
                                                   "Mammographic microcalcification found on diagnostic imaging of breast",
                                                   "Acquired absence of left breast and nipple",
                                                   "Acquired absence of unspecified breast and nipple",
                                                   "Acquired absence of right breast and nipple",
                                                   "Acquired absence of bilateral breasts and nipples",
                                                   "Benign neoplasm of right breast",
                                                   "Benign neoplasm of left breast")))

##output the file
save(longhaulers, file= paste0(outputDirectory,"/longhaulers_",site, ".RData"))


##summary statistics of organ,combo, subcombo
df <- longhaulers %>% 
  select("organ", "combo","duration", "subcombo", "patient_num") %>%
  # mutate(combo = ifelse(!is.na(subcombo),subcombo, organ)) %>%
  distinct()

na <- df[is.na(df$organ),]

duration_count = df %>% group_by(duration, organ, combo) %>%
  summarise(count = n())

new_duration_all_count = df %>% group_by(organ) %>%
  summarise(count_organ = n())

# df %>% group_by(Organ) %>%
#   summarise(count_organ = n_distinct(patient_num))
# the number in that duration / the number of cases in all durations within just that combo
duration_combo_count = df %>% group_by(organ, combo) %>%
  summarise(count_combo = n())

save(duration_count, file= paste0(outputDirectory,"/longhaulers_duration_organ_combo_count_",site, ".RData"))
save(duration_all_count, file= paste0(outputDirectory,"/longhaulers_organ_count_",site, ".RData"))
save(duration_combo_count, file= paste0(outputDirectory,"/longhaulers_organ_combo_count_",site, ".RData"))
