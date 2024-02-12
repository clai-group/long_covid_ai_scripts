### this script extracts the positive incidences data 
### and generates the summary statistics for the study population.
 
# Load packages ------------------------------------
library(dplyr)
library(stringr)
library(data.table)

choose_directory = function(caption = 'Select directory') {
  if(commandArgs()[[1L]] == "RStudio"){
    rstudioapi::selectDirectory(caption = caption)
  }else if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else{
    tcltk::tk_choose.dir(caption = caption)
  }
}

# Modify the following lines ------------------------------------ 
site <- "MGB"
cohort <- "cases" #cases, controls, controls_pre
cohortDirectory <- choose_directory(caption = "select cohort data directory")  ## where the inputs come from
outputDirectory <- choose_directory(caption = "select output data directory") ## where the outputs are saved


# Extract cov_pats.RData  ------------------------------------ 
df <- data.table::fread(paste0(cohortDirectory, "/",cohort, ".csv")) 
names(df) <- toupper(names(df))
# Patient data should have at least these five columns
target_cols <- c("PATIENT_NUM", "START_DATE", "CONCEPT_CD","C_FULLNAME", "PHENX")
if (all(target_cols %in% colnames(df))) {
  df <- subset(df, select = target_cols)
} else {
  print("The input data doesn't include 'PATIENT_NUM', 'START_DATE', 'CONCEPT_CD', 'C_FULLNAME' or 'PHENX' columns.")
}

df_rmdup <- unique(df)
df_rmdup$START_DATE <- as.POSIXct(df_rmdup$START_DATE, "%Y-%m-%d")
df_rmdup$PATIENT_NUM <- as.integer(df_rmdup$PATIENT_NUM)

rm(df)

if(cohort == 'cases'){
  rdx <-  df_rmdup[(df_rmdup$CONCEPT_CD=="ICD10:U07.1"), ]
  rp1 <-  df_rmdup[startsWith(df_rmdup$C_FULLNAME, "\\ACT\\UMLS_C0031437\\SNOMED_3947185011\\UMLS_C0037088\\SNOMED_3947183016\\"), ]
  rp2 <-  df_rmdup[startsWith(df_rmdup$C_FULLNAME, "\\ACT\\UMLS_C0031437\\SNOMED_3947185011\\UMLS_C0022885\\UMLS_C1335447\\"), ]
  rp3 <- df_rmdup[startsWith(df_rmdup$C_FULLNAME, "\\ACT\\UMLS_C0031437\\SNOMED_3947185011\\UMLS_C0022885\\ACT_LOCAL_LAB_ANY_POSITIVE\\"), ] 
  cases_encs <- rbind(rdx, rp1, rp2, rp3) %>% distinct()
  
  df_rmdup <- anti_join(df_rmdup, cases_encs) 
  
  last_dates <- df_rmdup %>%
    dplyr::group_by(PATIENT_NUM) %>%
    dplyr::summarise(max_date=max(START_DATE))
  last_dates$max_date_possib <- as.Date(last_dates$max_date)-365
  
  cases_encs <- merge(cases_encs, last_dates, by ="PATIENT_NUM")
  cases_encs <- cases_encs[as.Date(cases_encs$START_DATE) <= cases_encs$max_date_possib, ]
  cases_encs <- dplyr::select(cases_encs, -c(max_date, max_date_possib))
  
  save(cases_encs, file= paste0(cohortDirectory, "/cov_pats.RData")) 
  
  df_rmdup <- df_rmdup[df_rmdup$PATIENT_NUM %in% cases_encs$PATIENT_NUM, ]
  
  rm(rdx, rp1, rp2, rp3, cases_encs, last_dates)
  gc()
  
}
# else{
#   print("Continue to the next step.")
# }

# Summary Statistics  ------------------------------------
summary <- function(cohort){
  dems_all <-  data.table::fread(paste0(cohortDirectory, "/dems_",cohort, ".csv"))
  dems <- dems_all[, c("patient_num", "age", "sex_cd", "CHARLSON_INDEX", "race_cd", "ethnicity_cd")]
  rm(dems_all)
  
  # If the patient has more than 1 ages, select the oldest.
  dems <- unique(dems) %>% group_by(patient_num) %>% 
    filter(age == max(age)) %>% 
    ungroup()
  # If the patient has more than CHARLSON INDEX, select the largest.
  dems <- unique(dems) %>% group_by(patient_num) %>% 
    filter(CHARLSON_INDEX == max(CHARLSON_INDEX)) %>% 
    ungroup()
  
  names(df_rmdup) <- tolower(names(df_rmdup))
  dat.sum <- subset(df_rmdup,df_rmdup$patient_num %in% dems$patient_num)
  colnames(dat.sum) <- tolower(colnames(dat.sum))
  
  dems_temp1 <- data.frame(rbind(c("patients",nrow(dems))))
  dems_temp2 <- data.frame(rbind(c("mean age",mean(dems$age))))
  dems_sd <- data.frame(rbind(c("sd age",sd(dems$age))))
  dems_temp3 <- data.frame(rbind(c("percent female",(nrow(subset(dems,dems$sex_cd == "F"))/nrow(dems))*100)))
  dems_temp4 <- data.frame(rbind(c("mean charlson",mean(dems$CHARLSON_INDEX))))
  
  tb1 <- rbind(dems_temp1,dems_temp2,dems_sd,dems_temp3,dems_temp4) 
  tb2 <- data.frame(rbind(c("unique CCSR concepts",length(unique(dat.sum$phenx)))))
  rm(dems_temp1, dems_temp2, dems_temp3, dems_temp4)
  
  length.all <- dat.sum %>% group_by(patient_num) %>% summarise(max=max(start_date),min=min(start_date))
  length.all$durate=as.numeric(difftime(length.all$max, length.all$min,unit="weeks"))/52.25
  rm(dat.sum)
  
  tb3 <- data.frame(rbind(c("mean data depth",mean(length.all$durate))))
  length.all$max <- as.Date(length.all$max)
  tb4 <- data.frame(X1 = c("data to"), X2 = max(length.all$max))
  tb4$X2 <- as.character(tb4$X2)
  length.all$min <- as.Date(length.all$min)
  tb5 <- data.frame(X1 = c("data from"), X2 = min(length.all$min))
  tb5$X2 <- as.character(tb5$X2)
  
  table1 <- rbind(tb1,tb2,tb3,tb4,tb5)
  table1$site <- site
  table1$cohort <- cohort
  rm(tb1,tb2,tb3,tb4,tb5)
  
  race <- data.frame(table(dems$race_cd))
  race <- subset(race,race$Freq > 30)
  race$site <- site
  race$cohort <- cohort 
  
  eth <- data.frame(table(dems$ethnicity_cd))
  eth$site <- site
  eth$cohort <- cohort

  return(list(dems_stat = table1, race_stat =race, eth_stat=eth))
}

stat <- summary(cohort)

# Output  ------------------------------------
lapply(1:length(stat), function(i) write.csv(stat[[i]], 
                                                file = paste0(outputDirectory,"/",cohort, "_", names(stat[i]),"_", site,".csv"),
                                                row.names = FALSE))
fwrite(df_rmdup, paste0(cohortDirectory,"/", cohort, "_map_CCSR_",site, ".csv"))
