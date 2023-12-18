###this script maps the ICD codes to CCSR to obtain the phenx for further analysis
###uses ICD10 codes and codes domain (either ICD10CM or ICD10PCS) to map
###extracts the positive incidences data
###and generates the summary statistics for the study population.

# Load packages ------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

# Modify the following lines ------------------------------------ 
site <- "MGB" 
labtest_code <- "" # positive lab test code 
groupDirectory <- choose_directory(caption = "select group data directory")  ## where the inputs come from 
outputDirectory <- choose_directory(caption = "select output data directory") ## where the outputs are saved
group <- "cases" # cases, controls, controls_pre

# Load Data ------------------------------------ 
input_file <- paste0(groupDirectory, "/",group, ".csv")
df <-  data.table::fread(input_file) 
ccsr <- data.table::fread("P:/CCSR/CCSR_PASC_ICD.csv")

# Extract cov_pats.RData ------------------------------------ 
names(df) <- toupper(names(df))
df_rmdup <- unique(df)
gc()

if(group == 'cases'){
  
  rdx <-  df_rmdup[(df_rmdup$ICD10=="U07.1"), ]
  rlab <- df_rmdup[(df_rmdup$ICD10==labtest_code), ]
  cases_encs <- rbind(rdx, rlab)
  
  df_rmdup <- anti_join(df_rmdup, cases_encs)
  
  last_dates <- df_rmdup %>%
    dplyr::group_by(PATIENT_NUM) %>%
    dplyr::summarise(max_date=max(START_DATE))
  last_dates$max_date_possib <- as.Date(last_dates$max_date)-365
  
  cases_encs <- merge(cases_encs, last_dates, by ="PATIENT_NUM")
  cases_encs <- cases_encs[as.Date(cases_encs$START_DATE) <= cases_encs$max_date_possib, ]
  cases_encs <- cases_encs[,-c("max_date", "max_date_possib")]
  
  save(cases_encs, file= paste0(groupDirectory, "/cov_pats.RData")) 
  
  df_rmdup <- df_rmdup[df_rmdup$PATIENT_NUM %in% cases_encs$PATIENT_NUM, ]
  
  rm(rdx,rlab, cases_encs, last_dates)
  gc()
  
}

# Mapping based on ICD10 and DOMAIN------------------------------------ 
map_code<- merge(df_rmdup, ccsr, by = c("ICD10", "DOMAIN"), all.x = TRUE)

ms_rate <- n_distinct(map_code$ICD10[!is.na(map_code$CCSR_Key)])/n_distinct(map_code$ICD10)

if ((ms_rate) < 0.9){
  print("Need to manually check the cases data.")
} else {
  df_mapped <- map_code[!is.na(map_code$CCSR_Key),]
}

# Mapping results ------------------------------------ 
smr_df_ICD10 <- df_rmdup %>% group_by(DOMAIN) %>% summarise(nunique=n_distinct(ICD10), total = n())
smr_mapped <- df_mapped %>% group_by(DOMAIN) %>% summarise(nunique=n_distinct(ICD10), total=n())

result <-  cbind(smr_mapped, round(smr_mapped[-1]/smr_df_ICD10[-1]*100, digits = 3))
names(result)[-1] <- c("mapped_unique_codes", "mapped_total_records", "mapped_unique_codes_percentage", "mapped_total_records_percentage")

# Dems results ------------------------------------ 
dems_all <-  data.table::fread(paste0(groupDirectory, "/dems_",group, ".csv"))
dems <- dems_all[dems_all$patient_num %in% df_mapped$PATIENT_NUM, 
                 c("patient_num", "age", "sex_cd", "CHARLSON_INDEX", "race_cd", "ethnicity_cd")]
rm(dems_all)

# When the patient has more than 1 ages, select the oldest.
dems <- unique(dems) %>% group_by(patient_num) %>% 
  filter(age == max(age)) %>% 
  ungroup()
# When the patient has more than CHARLSON INDEX, select the largest.
dems <- unique(dems) %>% group_by(patient_num) %>% 
  filter(CHARLSON_INDEX == max(CHARLSON_INDEX)) %>% 
  ungroup()

dat.sum <- subset(df_mapped,df_mapped$PATIENT_NUM %in% dems$patient_num)
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
tb4 <- data.frame(X1 = c("data to"), X2 = max(length.all$max))
tb4$X2 <- as.character(tb4$X2)
tb5 <- data.frame(X1 = c("data from"), X2 = min(length.all$min))
tb5$X2 <- as.character(tb5$X2)

table1 <- rbind(tb1,tb2,tb3,tb4,tb5)
table1$site <- site
table1$group <- group
rm(tb1,tb2,tb3,tb4,tb5)

race <- data.frame(table(dems$race_cd))
race <- subset(race,race$Freq > 30)
race$site <- site
race$group <- group

eth <- data.frame(table(dems$ethnicity_cd))
eth$site <- site
eth$group <- group

# Output ------------------------------------ 
fwrite(table1, paste0(outputDirectory,"/", group,"_dems_", site, ".csv" ))
fwrite(race, paste0(outputDirectory, "/",group,"_race_", site, ".csv" ))
fwrite(eth, paste0(outputDirectory, "/",group,"_eth_", site, ".csv" ))
fwrite(result, paste0(outputDirectory, "/",group,"_map_stats_", site, ".csv" ))
fwrite(df_mapped, paste0(groupDirectory,"/", group, "_map_CCSR_",site, ".csv"))

