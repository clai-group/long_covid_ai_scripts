###this script implements WHO definition of long COVID
###this is step 3 for diagnosis of exclusion 
### we will use cases

## initial settings
Sys.setenv(R_MAX_NUM_DLLS = 999)
options("scipen"=100, "digits"=4)
seed <- 10000
set.seed(seed)

options(java.parameters = "-Xmx8048m")

####  Install and load the required packages
if(!require(pacman)) install.packages("pacman")

require(dplyr)
require(tidyr)
#devtools::install_github("hestiri/mlho")
require(mlho)
require(DT)
pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,RcppParallel,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools,tcltk)

##request parameters:
mem_buffer <- 5 #in GB. just a buffer to make sure the computer wont crash
cores_buffer <- 5 # choose the number of cores to free up make sure not to overload your computer!


##utils::choose.dir is a windows functionality use tk on other systems
choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}



# cov_pat incident level data
cov_pat_incident_FileName <- file.choose() ##the cov_pats.RData file
dbmartCases_FileName <-   file.choose() ##CCSR-mapped cases
corrsFileName <= file.choose() ## EITHER the file computed in the previous step or the precomputed delivered with the docker container
outputDirectory <- choose_directory(caption = "select output data directory") ## where outputs are saved
outputDirectory

#   ###### NON-INTERACTIVE MODE ### CHANGE THIS VARIABLES TO THE CORRECT PATH
# cov_pat_incident_FileName <- "set Path"
# dbmartCases_FileName <-   "setPath"
# dbmartControlls_FileName <-   "setPath"
# outputDirectory <- "select output data directory"
# corrsFileName <-  paste0(outputDirectory, "/corrs.RData")

numOfChunksFileName <- paste0(outputDirectory,"/num_of_case_chunks.RData")
phenxlookup_FileName <- paste0(outputDirectory, "/phenxlookup.RData")
apdativeDbFilenName <- paste0(outputDirectory,"/adpativeDBMart.RData")
#base file names will be completed in the loop
jBaseFileName <- paste0(outputDirectory,"/J_chunk_")
dbBaseFileName <- paste0(outputDirectory, "/db_longhauler_chunk_")
resultsFileName <- paste0(outputDirectory, "/point5_ccsr_mod_longCOVID.csv")




sparsity = 0.001
numOfThreads = detectCores()-cores_buffer

##load Js and correlations
load(phenxlookup_FileName)
load(numOfChunksFileName)
load(apdativeDbFilenName)
load(corrsFileName)

for(i in seq(1:numOfChunks)){
  
  jFileName = paste0(jBaseFileName, i, ".RData")
  load(file=jFileName)

  dbmart_num <- dbmart_adapt$chunks[[i]]
  
  gc()
  ###covid codes
  cov_cods <- c(subset(phenxlookup$num_Phenx,phenxlookup$phenx %like% "COVID")) 
  virals_cod <- c(subset(phenxlookup$num_Phenx,phenxlookup$phenx %in% c("Viral infection")))
  
  ################ identify correlated Js
  ###we want Js thathave significant correlations and also correlated to covid
  corrs_cov_J_sig <- subset(corrs,corrs$sequence >0 & corrs$p.adjust <=0.05 & corrs$rho >=0.3 & 
                              corrs$startPhen %in% cov_cods & !(corrs$endPhenx %in% cov_cods))
  unique(corrs_cov_J_sig$endPhenx)
  J <- subset(J,J$endPhenx %in% corrs_cov_J_sig$endPhenx)
  

  
  endPhenx = c(J$endPhenx,cov_cods) #TODO look up id for covid phenx in db$phenxLookUp
  temporalBucket =  c(0,1,3)
  minDuration = 0 #techical parameter, ignore for now
  bitShift = 0 #techical parameter, ignore for now
  lengthOfPhenx = 7 #techical parameter
  storeSequencesDuringCreation = FALSE #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity
  
  ###get all the sequences that end with
  corseq <- tSPMPlus::getSequencesWithEndPhenx(dbmart_num,
                                               bitShift,
                                               lengthOfPhenx,
                                               temporalBucket,
                                               endPhenx,
                                               includeCorBuckets=TRUE,
                                               minDuration,
                                               storeSequencesDuringCreation,
                                               numOfThreads = numOfThreads,
                                               sparsityValue = sparsity)
  corseq <- dplyr::distinct(corseq, .keep_all = TRUE)
  gc()
  
  ##########################################################################
  ###########################################################################
  #######################  let the EXCLUSIONS begin!
  ##### ########## ########## ########## ########## ########## ########## #####
  ##### ########## ########## ########## ########## ########## ########## #####
  ##### ########## ########## ########## ########## ########## ########## ########## #####
  ##### ########## ########## ########## ########## ########## ########## ########## #####
  
  ##construct covid --> J sequences
  Js <- list()
  for(j in 1:length(cov_cods)) {
    Js.j <- J
    Js.j$sequence <- ifelse(Js.j$endPhenx<10,paste0(cov_cods[j],"000000",Js.j$endPhenx),NA)
    Js.j$sequence <- ifelse(Js.j$endPhenx<100 & Js.j$endPhenx>9,paste0(cov_cods[j],"00000",Js.j$endPhenx),Js.j$sequence)
    Js.j$sequence <- ifelse(Js.j$endPhenx<1000 & Js.j$endPhenx>99,paste0(cov_cods[j],"0000",Js.j$endPhenx),Js.j$sequence)
    
    Js[[j]] <- Js.j
  }
  Js <- data.table::rbindlist(Js)
  Js$sequence <- as.numeric(Js$sequence)
  
  Js <- merge(Js[,"sequence"],corrs,by="sequence")
  ##sort by sequence number to make sure we start from infection 1
  Js <- Js[order(startPhen),]
  rm(Js.j);gc()
  
  corrs$key_corr <- paste0(corrs$endPhenx,"--",corrs$startPhen_dur)
  corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))
  corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)
  corseq$key_corr <- paste0(corseq$endPhenx,"--",corseq$startPhen_dur)
  
  ########################################################################
  ####### now incorporating correlations ########################################################
  ########################################################################
  
  ######## HYPERPARAMETER ALERT!
  exlusion_corrs <- subset(corrs,corrs$p.adjust <= 0.005 & corrs$rho >= 0.5) 
  
  ###J->J sequences
  # we will only remove Js when J->J with longer than 3 months is significant and also the patient has a J-->COVID sequence
  JJ3 <- exlusion_corrs %>%
    filter(endPhenx == startPhen & startPhen_dur %like% "-3" & !(endPhenx %in% cov_cods)) 
  # we will only remove j when X->J is significant and X !=J
  XJ <- exlusion_corrs %>%
    filter(endPhenx != startPhen & !(endPhenx %in% cov_cods)) %>%
    filter(!(startPhen_dur %like% "-0")) #### removed co-0ccurences for now!!!!!!!

  ##### ########## ########## ########## ########## ########## #####
  ### create a database to cut from for identifying long COVID patients
  # db_longhaulers <- subset(corseq,corseq$sequence %in% c(Js$sequence))
  ##### ########## ########## ########## ########## ########## #####

  #setup parallel backend to use many processors
  cores<-detectCores()
  cl <- parallel::makeCluster(cores_buffer) 
  doParallel::registerDoParallel(cl)
  
  # PARALELIZING THE longhauler db creation
  db_longhaulers <- foreach(i = 1:nrow(Js),#
                            .combine='rbind', 
                            .packages = c("data.table","dplyr", "tidyr","plyr","scales","DT", "lubridate", "tidyverse","reshape2"),
                            .multicombine=TRUE) %dopar% {
                              tryCatch({
                                sequence.i <- as.numeric(Js[i,"sequence"])
                                J.i <- as.numeric(Js[i,"endPhenx"])
                                cov.i <- as.numeric(Js[i,"startPhen"])
                                
                                ### data from all patients have the sequence
                                corseq.i <- subset(corseq,corseq$patient_num %in% 
                                                     unique(subset(corseq$patient_num,corseq$sequence == sequence.i)))
                                
                                ##################  ################## ##################
                                ###exclusion by explanation
                                ### see if COVID-->J can be explained by an event prior
                                
                                ###J->J sequences
                                # we will only remove Js when J->J with longer than 3 months is significant and also the patient has a J-->COVID sequence
                                JJ3.i <- subset(JJ3,JJ3$startPhen %in% J.i)
                                
                                if(nrow(JJ3.i) > 0) {
                                  JJ3.i$JCOV <- ifelse(cov.i<10,paste0(JJ3.i$startPhen,"000000",cov.i),NA)
                                  JJ3.i$JCOV <- ifelse(cov.i<100 & cov.i>9,paste0(JJ3.i$startPhen,"00000",cov.i),JJ3.i$JCOV)
                                  JJ3.i$JCOV <- ifelse(cov.i<1000 & cov.i>99,paste0(JJ3.i$startPhen,"0000",cov.i),JJ3.i$JCOV)
                                  
                                  JJ3.i$JCOV <- as.numeric(JJ3.i$JCOV)
                                  ###patients who have the J-J sequences before COVI->J, cannot have the given COV->J
                                  #### RTI?!?!
                                  exclude.i.jj3 <- subset(corseq.i,corseq.i$patient_num %in% 
                                                            unique(subset(corseq.i$patient_num,corseq.i$key_corr == JJ3.i$key_corr)))
                                  
                                  if (nrow(exclude.i.jj3) > 0) {
                                    exclude.i.jj3 <-subset(exclude.i.jj3,exclude.i.jj3$patient_num %in% 
                                                             unique(subset(exclude.i.jj3$patient_num,exclude.i.jj3$sequence %in% JJ3.i$JCOV)))
                                  }
                                  
                                  if (nrow(exclude.i.jj3) == 0) {exclude.i.jj3 <- c(999999999999)}
                                  if (nrow(exclude.i.jj3) > 0) { exclude.i.jj3 <- c(unique(exclude.i.jj3$patient_num))}
                                }
                                if (nrow(JJ3.i) == 0) {exclude.i.jj3 <- c(999999999999)
                                }
                                
                                # exclude.i.jj3s cannot have sequence.i
                                
                                ### X->J sequences
                                # we will only remove j when X->J is significant and X !=J
                                XJ.i <- subset(XJ,XJ$endPhenx %in% J.i)
                                
                                if (nrow(XJ.i)) {
                                  ##construct J --> covid sequences
                                  startphen.i <- c(unique(XJ.i$startPhen))
                                  XJ.is <- list()
                                  for(o in 1:length(startphen.i)) {
                                    XJ.i.o <- subset(XJ.i,XJ.i$startPhen == startphen.i[o])
                                    
                                    XJ.i.o$JCOV <- ifelse(cov.i<10,paste0(startphen.i[o],"000000",cov.i),NA)
                                    XJ.i.o$JCOV <- ifelse(cov.i<100 & cov.i>9,paste0(startphen.i[o],"00000",cov.i),XJ.i.o$JCOV)
                                    XJ.i.o$JCOV <- ifelse(cov.i<1000 & cov.i>99,paste0(startphen.i[o],"0000",cov.i),XJ.i.o$JCOV)
                                    XJ.i.o$JCOV <- as.numeric(XJ.i.o$JCOV)
                                    
                                    XJ.is[[o]] <- XJ.i.o
                                    rm(XJ.i.o)
                                  }
                                  XJ.i <- data.table::rbindlist(XJ.is)
                                  rm(XJ.is,startphen.i)
                                  
                                  ###patients who have the X-J sequences before COVI->J, cannot have the given COV->J
                                  #### RTI?!?!
                                  exclude.i.xj <- subset(corseq.i,corseq.i$patient_num %in% 
                                                           unique(subset(corseq.i$patient_num,corseq.i$key_corr %in% c(XJ.i$key_corr))))
                                  
                                  if (nrow(exclude.i.xj) > 0) {
                                    exclude.i.xj <-subset(exclude.i.xj,exclude.i.xj$patient_num %in% 
                                                            unique(subset(exclude.i.xj$patient_num,exclude.i.xj$sequence %in% c(XJ.i$JCOV))))
                                  }
                                  
                                  if (nrow(exclude.i.xj) == 0) {exclude.i.xj <- c(999999999999)}
                                  if (nrow(exclude.i.xj) > 0) { exclude.i.xj <- c(unique(exclude.i.xj$patient_num))}
                                }
                                if (nrow(XJ.i) == 0) {exclude.i.xj <- c(999999999999)}

                                #### remove by correlation
                                exclude.by.corr <- c(unique(c(exclude.i.jj3,exclude.i.xj)))
                                
                                db_longhaulers.i <- subset(corseq,corseq$sequence == sequence.i & !(corseq$patient_num %in% exclude.by.corr))
                                
                                ##### now working on recurrence with the db_longhaulers.i
                                
                                ########################################################################################
                                ######## start recurrence screening
                                ########################################################################################################
                                
                                ###remove sequences with duration less than 1 months or longer than 12 months
                                db_longhaulers.i <- db_longhaulers.i %>%
                                  filter(duration >= 1) %>%
                                  filter(duration<=12)
                                
                                ##if a patient does not have a second sequence of COVID->J, sequence has to be dropped from that patient because the 2+ months does not apply
                                
                                exlcude.i_count1 <- db_longhaulers.i %>%
                                  dplyr::group_by(patient_num) %>%
                                  dplyr::summarise(count=length((patient_num))) %>%
                                  filter(count == 1)
                                
                                db_longhaulers.i <- subset(db_longhaulers.i,!(db_longhaulers.i$patient_num %in% exlcude.i_count1$patient_num ))
                                
                                # recalculate and merge min date to corseq
                                db_longhaulers.i_min <- db_longhaulers.i %>%
                                  dplyr::group_by(patient_num) %>%
                                  dplyr::summarise(min_duration=min((duration)))
                                
                                db_longhaulers.i <- merge(db_longhaulers.i,db_longhaulers.i_min,by=c("patient_num"))
                                db_longhaulers.i$delta_dur <- db_longhaulers.i$duration-db_longhaulers.i$min_duration
                                
                                exclude.i.recurrence <- db_longhaulers.i %>%
                                  dplyr::group_by(patient_num) %>%
                                  dplyr::summarise(max_delta_dur=max((delta_dur))) %>%
                                  filter(max_delta_dur < 2) 
                                
                                db_longhaulers.i <- subset(db_longhaulers.i,!(db_longhaulers.i$patient_num %in% exclude.i.recurrence$patient_num ))
                                
                                rm(sequence.i,J.i,cov.i,JJ3.i,exclude.i.jj3.pats,startphen.i,exclude.i.xj,XJ.i,
                                   db_longhaulers.i_min,exclude.i.recurrence,exlcude.i_count1);gc()
                                
                                db_longhaulers.i <- subset(db_longhaulers.i,db_longhaulers.i$delta_dur == 0)
                                
                                db_longhaulers.i
                                
                              },
                              error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
                            }

  db_longhaulers <- merge(db_longhaulers,J,by="endPhenx")
  ### refresh the Js remaining
  J_updated <- c(unique(db_longhaulers$endPhenx))
  J$update <- ifelse(J$endPhenx %in% J_updated,1,0)
  
  #TODO: save Js from chunks
  dbFileName = paste0(dbBaseFileName, i, ".RData")
  save(db_longhaulers,file = dbFileName)
  #jFileName = paste0("P:/PASC/data/J_notExcluded_chunk_", i, ".RData")
  #save(J, file=jFileName)
  rm(J,J_updated, db_longhaulers,Js, corseq,exlusion_corrs,XJ,JJ3)
  gc()
}


#load chunks and merge

for(i in seq(1:numOfChunks_J)){
  dbFileName = paste0(dbBaseFileName, i, ".RData")
  #jFileName = paste0("P:/PASC/data/J_notExcluded_chunk_", i, ".RData")
  
  if(i == 1){
    db_longhaulers <- load(dbFileName)
    #J <- load(jFileName)
  }else{
    db_longhaulers <- rbind(db_longhaulers, load(dbFileName))
    #J <- rbind(J, load(jFileName))
  }
}


###get descriptive statistics
length(unique(db_longhaulers$patient_num))
long_COVID <- data.frame(table(db_longhaulers$phenx))
write.csv(long_COVID_poin5,file=resultsFileName)