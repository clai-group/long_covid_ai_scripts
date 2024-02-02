###this script implements WHO definition of long COVID using the reference J thresholds from the MGB study
###this is step 3 implements diagnosis of exclusion
### we will use cases only to perform exclusion by correlation

### use this script for implementing the pre-trained algorithm, skipping local tuning steps

###hyper-parameters
site <- "MGB"
param1 <- "ref_" ###
param2 <- "thresholds" ### 
p <- 0.05 ### p-value for correlations
sparsity <- 0.001 ##sparsity screening
mem_buffer <- 5 #in GB. just a buffer to make sure the computer wont crash
cores_buffer <- 1 # choose the number of cores to free up make sure not to overload your computer!
### this number will be taken from available cores. Set a number that get's you  3-5 cores max, based on the memory available


### load the needed libraries
{
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




##utils::choose.dir is a windows functionality use tk on other systems
choose_directory = function(caption = 'Select directory') {
  if(commandArgs()[[1L]] == "RStudio"){
    rstudioapi::selectDirectory(caption = caption)
  }else if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else{
    tcltk::tk_choose.dir(caption = caption)
  }
}

}

# cov_pat incident level data
cov_pat_incident_FileName <- file.choose() ##the cov_pats.RData file
dbmartCases_FileName <-   file.choose() ##CCSR-mapped cases
corrsFileName <- file.choose() ## the pre-computed delivered with the docker container -- ref_corrs.RData
refJFileName <- file.choose() ## the pre-computed delivered with the docker container -- ref_J_thresholds.RData
JFileName <- file.choose() ## the pre-computed delivered with the docker container -- ref_J.RData

outputDirectory <- choose_directory(caption = "select output data directory") ## where outputs are saved


numOfChunksFileName <- paste0(outputDirectory,"/num_of_case_chunks.RData")
phenxlookup_FileName <- paste0(outputDirectory, "/phenxlookup.RData")
patlookup_FileName <- paste0(outputDirectory, "/patlookup.RData")
apdativeDbFilenName <- paste0(outputDirectory,"/adpativeDBMart.RData")
  #base file names will be completed in the loop
jBaseFileName <- paste0(outputDirectory,"/J_chunk_")
corrsBaseFileName <-  paste0(outputDirectory, "/corrs_chunk_")
dbBaseFileName <- paste0(outputDirectory, "/db_longhauler_chunk_")
resultsFileName_sum <- paste0(outputDirectory, "/longCOVID_summary_",site,param1,param2,p,".csv")
resultsFileName_longCOVID_patients <- paste0(outputDirectory, "/longCOVID_patients_",site,param1,param2,p,".csv")


####load the cases data
dbmart_cases_map_ccsr <- data.table::fread(dbmartCases_FileName)
##load the cov_pat incident level data
load(cov_pat_incident_FileName)
load(corrsFileName)
load(refJFileName)
load(JFileName)

colnames(J_threshold_correlations)[1] <- "phenx"

### IS EVERYONE GOING TO HAVE THE RELEVANCE COLUMN OR ARE WE CUTTING IT FROM THE QUERY?
# dbmart_cases_map_ccsr <- subset(dbmart_cases_map_ccsr,is.na(dbmart_cases_map_ccsr$relevance))


colnames(dbmart_cases_map_ccsr) <- tolower(colnames(dbmart_cases_map_ccsr))
length(unique(dbmart_cases_map_ccsr$phenx))
colnames(cov_pats) <- tolower(colnames(cov_pats))


dbmart <- dplyr::select(dbmart_cases_map_ccsr,patient_num,start_date,phenx)
rm(dbmart_cases_map_ccsr);gc()

cov_pats <- subset(cov_pats,cov_pats$infection_seq <= 5)
cov_pats$phenx <- paste0("COVID",cov_pats$infection_seq)
dbmart <- rbind(dbmart,dplyr::select(cov_pats,patient_num,start_date,phenx))

dbmart <- dplyr::distinct(dbmart, .keep_all = TRUE)

db <- tSPMPlus::transformDbMartToNumeric(dbmart)

phenxlookup <- db$phenxLookUp
save(phenxlookup,file=phenxlookup_FileName)
patlookup <- db$patientLookUp
save(patlookup, file = patlookup_FileName)

  ### define sequencing parameters

numOfThreads = detectCores()-cores_buffer


J_threshold_correlations <- merge(J_threshold_correlations,phenxlookup,by="phenx")
cov_cods <- c(subset(phenxlookup$num_Phenx,phenxlookup$phenx %like% "COVID\\d+"))


endPhenx = c(J_threshold_correlations$num_Phenx,cov_cods) #TODO look up id for covid phenx in db$phenxLookUp
temporalBucket =  c(0,1,3)
minDuration = 0 #technical parameter, ignore for now
bitShift = 0 #technical parameter, ignore for now
lengthOfPhenx = 7 #technical parameter
storeSequencesDuringCreation = FALSE #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity


### adaptive chunk sizes
buffer<- 100000000*mem_buffer #extra buffer in bytes find good value *a higher value here decrease the chunk size, but left more memory to work with after the sequencing, the chunk size is calculated on the memory consumption of the non-sparse sequences, reduce buffer
dbmart_adapt <- tSPMPlus::splitdbMartInChunks(db$dbMart, includeCorSeq = TRUE, buffer = buffer)
numOfChunks = length(dbmart_adapt$chunks)

# save(numOfChunks,file=numOfChunksFileName)
# save(dbmart_adapt, file=apdativeDbFilenName )
names(J_threshold_correlations)[3] <- "endPhenx"

J <- subset(J,J$phenx %in% J_threshold_correlations$phenx)
J_base <- J


for(i in seq(1:numOfChunks)){

  # jFileName = paste0(jBaseFileName, i, ".RData")
  # load(file=jFileName)

dbmart_num <- dbmart_adapt$chunks[[i]]

gc()


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
  #######################  personalized EXCLUSION by temporal association
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

 
  exlusion_corrs <- merge(J_threshold_correlations,corrs,by="endPhenx")
  exlusion_corrs <- subset(exlusion_corrs,exlusion_corrs$p.adjust <= p & exlusion_corrs$rho >= exlusion_corrs$`cut-off threshold`)
  # exlusion_corrs <- subset(exlusion_corrs,exlusion_corrs$endPhenx %in% J$endPhenx)
  
  ###J->J sequences
  # we will only remove Js when the patient has a J-->COVID sequence
  JJ3 <- exlusion_corrs %>%
    filter(endPhenx == startPhen & !(endPhenx %in% cov_cods))
  # we will only remove j when X->J is significant and X !=J
  XJ <- exlusion_corrs %>%
    filter(endPhenx != startPhen & !(endPhenx %in% cov_cods)) %>%
    filter(!(startPhen_dur %like% "-0")) #### removed co-0ccurences for now!!!!!!!
  

  #setup parallel backend to use many processors
  cl <- parallel::makeCluster(numOfThreads)
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

                                if (nrow(XJ.i) > 0) {
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
  rm(J_updated, db_longhaulers,Js, corseq,exlusion_corrs,XJ,JJ3)
  stopCluster(cl)
  gc()
}


#load chunks and merge
load(apdativeDbFilenName)

for(i in seq(1:numOfChunks)){
  dbFileName = paste0(dbBaseFileName, i, ".RData")
  #jFileName = paste0("P:/PASC/data/J_notExcluded_chunk_", i, ".RData")
  load(dbFileName)
  nextChunk <- db_longhaulers;rm(db_longhaulers)
  colnames(nextChunk)[colnames(nextChunk)== "patient_num"] = "chunk_pat_num"
  nextChunk <-  dplyr::left_join(nextChunk, dbmart_adapt$lookUps[[i]], by="chunk_pat_num") %>%
    dplyr::left_join(patlookup, by="num_pat_num") %>%
    dplyr::select(-"num_pat_num", -"chunk_pat_num")

  if(i == 1){
    longhaulers <- nextChunk
    #J <- load(jFileName)
  }else{
    nextChunk <-  nextChunk[,c(colnames(longhaulers))]
    longhaulers <- rbind(longhaulers, nextChunk)
    #J <- rbind(J, load(jFileName))
  }
}


###get descriptive statistics
temp <- data.frame(length(unique(longhaulers$patient_num)))
temp$Var1 <- "unique long haulers"
colnames(temp)[1] <- "Freq"
temp <- temp[c("Var1", "Freq")]
long_COVID_list <- longhaulers %>%
  dplyr::group_by(phenx) %>%
  dplyr::summarise(unique_patients=length(unique(patient_num)))
colnames(long_COVID_list) <- c("Var1", "Freq")
long_COVID_list <- rbind(long_COVID_list,temp)

long_COVID_list <- cbind(long_COVID_list, param1 = c(param1), param2 = c(param2))
longhaulers <- cbind(longhaulers, param1 = c(param1), param2 = c(param2))

write.csv(long_COVID_list,file=resultsFileName_sum)
write.csv(longhaulers,file=resultsFileName_longCOVID_patients)


