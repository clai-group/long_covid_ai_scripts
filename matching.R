# matching script for pulling data from the control group


####  Install and load the required packages
if(!require(pacman)) install.packages("pacman")



require(dplyr)
require(tidyr)
require(DT)
pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,RcppParallel,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools)

##utils::choose.dir is a windows functionality use tk on other systems
choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

###select demoraphic file for cases and controls
dems_cases <- file.choose() 
dems_controls <- file.choose() 

###select the output directory:  data/controls/ 
outputDirectory <- choose_directory(caption = "select output data directory") ## where outputs are saved


####now you can match by demographics and comorbidity
dems = rbind(dems_cases, dems_controls)
dems$case = factor(c(rep(1,dim(dems_cases)[1]), rep(0, dim(dems_controls)[1])))

library("MatchIt")
match = matchit(case ~ sex_cd + race_cd + age + CHARLSON_INDEX, data = dems,
                method = "nearest", distance = "glm", ratio = 1, replace = F)
mdat = match.data(match)
write.csv(mdat$patient_num[mdat$case == 0],file=paste0(outputDirectory,"/patient_nums.csv"))

patients <- c(mdat$patient_num[mdat$case == 0],"x")

dems_controls <- subset(dems_controls,dems_controls$patient_num %in% patients)
write.csv(dems_controls,file=paste0(outputDirectory,"/dems_controls.csv"))

