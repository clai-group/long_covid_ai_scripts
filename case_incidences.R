###this script creates incident-level data from patient encounters for COVID infections
### rule is to cluster infections dates and recognize an infection if a cluster is 90 days or longer apart from another
if(!require(pacman)) install.packages("pacman")

pacman::p_load(data.table, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,Rmisc,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,
               ggridges, forcats, stats)



#### load the encounters data
### we will need 3 columns patient_num, start_date, ENCOUNTER_NUM
## better not have duplicated rows!

load("~/clai_share/shared_workspace/ACT/outputs/cases/encounters_cases.RData")


##remove potential duplicates within a month
#break down dates by quarter
cases_encs$date_ym <- format(cases_encs$start_date,format="%y%m")
# ##if there is multiple records within a month, let's only take the first one

##this is to merge source later
source_temp <- dplyr::select(cases_encs,patient_num,start_date)

cases_encs <- cases_encs %>%
  dplyr::group_by(patient_num,date_ym) %>%
  dplyr::summarise(start_date=min(start_date)) %>%
  dplyr::select(patient_num,start_date)
#
cases_encs <- dplyr::distinct(cases_encs, .keep_all = TRUE)


length(unique(cases_encs$patient_num))

##let's clean up the data so it only captures the first PCR in 90 days
reps <- data.frame(table(cases_encs$patient_num))
reps <- subset(reps,reps$Freq > 1)


uniqpats <- c(as.character(unique((reps$Var1))))


cases_encs <- setDT(cases_encs)

###boost loop speed by 100 times
dat.mi <- subset(cases_encs,cases_encs$patient_num %in% uniqpats)
dat.mi$start_date <- as.POSIXct(dat.mi$start_date, "%Y-%m-%d")

# 
# cores<-detectCores()
# cl <- makeCluster(cores[1]-92,"FORK") ### 2 is minimum suggested cores to be removed from parallelization. it also has memory implications!
# registerDoParallel(cl)
# 
# 

# 
# cases_encs_mi <- foreach(i = 1:length(uniqpats),#
#                           .combine='rbind', 
#                           .packages = c("data.table","dplyr", "tidyr","plyr","scales","DT", "lubridate", "tidyverse","reshape2")
#                          ) %dopar% {
#                             tryCatch({

cases_encs_mi <- list()
for (i in 1:length(uniqpats)) {
  tryCatch({
#     
    print(paste0(i," of ",length(uniqpats)))
    
    pat.dat.i <- dat.mi %>%
      filter(patient_num == uniqpats[i]) %>%
      arrange(start_date) %>%
      mutate(since_piorinf = c(NA,as.numeric(diff(start_date), units="days")))  %>%
      filter(since_piorinf > 90 | is.na(since_piorinf))%>%
      mutate(infection_seq = dense_rank((start_date)))
    
    # pat.dat.i
cases_encs_mi[[i]] <- pat.dat.i
  },
  error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
}
cases_encs_mi <- data.table::rbindlist(cases_encs_mi)


# },
# error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
#                           }


##those who do not have more than 1 positive signal
cases_encs <- subset(cases_encs,!(cases_encs$patient_num %in% cases_encs_mi$patient_num))
cases_encs <- dplyr::select(cases_encs,patient_num,start_date)
cases_encs$infection_seq <- 1
cases_encs$since_piorinf <- NA

reps <- data.frame(table(cases_encs_mi$patient_num))
reps <- subset(reps,reps$Freq > 1)

##those who have multiple labs but within 90 days
cases_encs_mi_no <- subset(cases_encs_mi,!(cases_encs_mi$patient_num %in% reps$Var1))
##those who truly have more than infections
cases_encs_mi <- subset(cases_encs_mi,(cases_encs_mi$patient_num %in% reps$Var1))


cov_pats <- rbind(cases_encs,cases_encs_mi,cases_encs_mi_no)

rm(cases_encs_mi,cases_encs_mi_no,cases_encs)

stopCluster(cl)
rm(list=setdiff(ls(), "cov_pats"))
save.image("~/clai_share/shared_workspace/ACT/outputs/cases/cov_pats.RData")