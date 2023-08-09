###this script implements WHO definition of long COVID
###this is step 1 to identify potential candidates or Js

##request parameters:
mem_buffer <- 1 #in GB. just a buffer to make sure the computer wont crash
cores_buffer <- 1 # choose the number of cores to free up - make sure not to overload your computer!


###libraries and stuff!
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
require(mlho)
require(DT)
pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,RcppParallel,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools, tcltk)





  ##utils::choose.dir is a windows functionality use tk on other systems
  choose_directory = function(caption = 'Select directory') {
    if (exists('utils::choose.dir')) {
      choose.dir(caption = caption) 
    } else {
      tcltk::tk_choose.dir(caption = caption)
    }
  }
  
}
  
    # cov_pat incident level data
    cov_pat_incident_FileName <- file.choose() ##the cov_pats.RData file
    dbmartCases_FileName <-   file.choose() ##CCSR-mapped cases
    outputDirectory <- choose_directory(caption = "select output data directory") ## where outputs are saved
    # dbmartControlls_FileName <-   file.choose()

    
    
    
### run the following chunk
{
    #   ###### NON-INTERACTIVE MODE ### CHANGE THIS VARIABLES TO THE CORRECT PATH
    # cov_pat_incident_FileName <- "set Path"
    # dbmartCases_FileName <-   "setPath"
    # dbmartControlls_FileName <-   "setPath"
    # outputDirectory <- "select output data directory"

  numOfChunksFileName <- paste0(outputDirectory,"/num_of_case_chunks.RData")
  phenxlookup_FileName <- paste0(outputDirectory, "/phenxlookup.RData")
  patlookup_FileName <- paste0(outputDirectory, "/patlookup.RData")
  apdativeDbFilenName <- paste0(outputDirectory,"/adpativeDBMart.RData")
  #base file names will be completed in the loop
  jBaseFileName <- paste0(outputDirectory,"/J_chunk_")
  corrsBaseFileName <-  paste0(outputDirectory, "/corrs_chunk_")
  dbBaseFileName <- paste0(outputDirectory, "/db_longhauler_chunk_")
  resultsFileName <- paste0(outputDirectory, "/point5_ccsr_mod_longCOVID.csv")



####load the cases data
dbmart_cases_map_ccsr <- data.table::fread(dbmartCases_FileName)
##load the cov_pat incident level data
load(cov_pat_incident_FileName)
dbmart_cases_map_ccsr <- subset(dbmart_cases_map_ccsr,is.na(dbmart_cases_map_ccsr$relevance))

length(unique(dbmart_cases_map_ccsr$phenx))

colnames(dbmart_cases_map_ccsr) <- tolower(colnames(dbmart_cases_map_ccsr))
colnames(cov_pats) <- tolower(colnames(cov_pats))


dbmart <- dplyr::select(dbmart_cases_map_ccsr,patient_num,start_date,phenx)
rm(dbmart_cases_map_ccsr);gc()

cov_pats <- subset(cov_pats,cov_pats$infection_seq <= 7)
cov_pats$phenx <- paste0("COVID",cov_pats$infection_seq)
dbmart <- rbind(dbmart,dplyr::select(cov_pats,patient_num,start_date,phenx))

dbmart <- dplyr::distinct(dbmart, .keep_all = TRUE)

db <- tSPMPlus::transformDbMartToNumeric(dbmart)

phenxlookup <- db$phenxLookUp
save(phenxlookup,file=phenxlookup_FileName)
patlookup <- db$patientLookUp
save(patlookup, file = patlookup_FileName)

### define sequencing parameters
sparsity = 0.001
numOfThreads = detectCores()-cores_buffer

phenxOfInterest = c(as.numeric(db$phenxLookUp[(db$phenxLookUp$phenx %like% "COVID"),"num_Phenx"]$num_Phenx)) #TODO look up id for covid phenx in db$phenxLookUp
temporalBucket =  c(0,1,3)
minDuration = 0 #techical parameter, ignore for now ##TODO if not working with 0 use 1
bitShift = 0 #techical parameter, ignore for now
lengthOfPhenx = 7 #techical parameter
storeSequencesDuringCreation = FALSE #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity


### adaptive chunk sizes
buffer<- 100000000*mem_buffer #extra buffer in bytes find good value *a higher value here decrease the chunk size, but left more memory to work with after the sequencing, the chunk size is calculated on the memory consumption of the non-sparse sequences, reduce buffer 
dbmart_adapt <- tSPMPlus::splitdbMartInChunks(db$dbMart, includeCorSeq = TRUE, buffer = buffer)
numOfChunks = length(dbmart_adapt$chunks)

save(numOfChunks,file=numOfChunksFileName)
save(dbmart_adapt, file=apdativeDbFilenName )
J_loop <- list()
for (i in seq(1:numOfChunks)) {
  tryCatch({
  dbmart_num <- dbmart_adapt$chunks[[i]]
  corseq <- tSPMPlus::getCandidateSequencesForPOI(dbmart_num,
                                                  minDuration,
                                                  bitShift,
                                                  lengthOfPhenx,
                                                  temporalBucket,
                                                  phenxOfInterest,
                                                  storeSequencesDuringCreation,
                                                  numOfThreads = numOfThreads,
                                                  sparsityValue = sparsity)
  
  
  corseq <- dplyr::distinct(corseq, .keep_all = TRUE)
  gc()
  
  ##construct covid --> J sequences
  J <- data.frame(unique(corseq$endPhenx))
  colnames(J) <- "endPhenx"
  Js <- list()
  for(j in 1:length(phenxOfInterest)) {
    Js.j <- J
    
    Js.j$candidate_sequences <- ifelse(Js.j$endPhenx<10,paste0(phenxOfInterest[j],"000000",Js.j$endPhenx),NA)
    Js.j$candidate_sequences <- ifelse(Js.j$endPhenx<100 & Js.j$endPhenx>9,paste0(phenxOfInterest[j],"00000",Js.j$endPhenx),Js.j$candidate_sequences)
    Js.j$candidate_sequences <- ifelse(Js.j$endPhenx<1000 & Js.j$endPhenx>99,paste0(phenxOfInterest[j],"0000",Js.j$endPhenx),Js.j$candidate_sequences)
    
    Js[[j]] <- Js.j
  }
  J <- data.table::rbindlist(Js)
  J$candidate_sequences <- as.numeric(J$candidate_sequences)
  rm(Js,Js.j)
  ###
  corseq_sub <- subset(corseq,corseq$sequence %in% c(J$candidate_sequences))
  
  ##now start looping through with the last infection
  infections <- data.frame(db$phenxLookUp[(db$phenxLookUp$phenx %like% "COVID"),c("phenx","num_Phenx")])
  infections$num <- as.numeric(sub('.*COVID', '', infections$phenx))
  
  J_low <- list()
  
  for (inf in max(infections$num):1){
    tryCatch({
      gc()
      phen_num_inf <- infections[infections$num <- inf,"num_Phenx"]
      corseq_sub_inf <- corseq_sub %>%
        dplyr::filter(sequence %like% paste0("^",phen_num_inf))
      
      pat_num_inf <- length(unique(corseq_sub_inf$patient_num))
      
      if (nrow(corseq_sub_inf) > 0) {
        ### now we want to make sure that patients have a data point 12 months after this last infection
        ### if not, they are out for this round of sparsity calculations
        
        pat_sub_inf <- corseq_sub_inf %>%
          dplyr::group_by(patient_num) %>%
          dplyr::summarise(duration_max=max(duration)) %>%
          dplyr::filter(duration_max >= 12)
        
        
        if (nrow(pat_sub_inf) > 0){
          corseq_sub_inf <- corseq_sub_inf %>%
            ### sequences longer than 12 months are not of interest!
            dplyr::filter(duration <= 12) %>%
            dplyr::filter(patient_num %in% pat_sub_inf$patient_num)
          
          ##### these patients have follow up data to consider
          
          ###### now make sure the patient has more than 1 of the sequences post that infection
          corseq_sub_count_inf <- corseq_sub_inf %>%
            dplyr::group_by(patient_num,sequence) %>%
            dplyr::summarise(count=length((patient_num))) #%>%
          # filter(count > sparsity*length(unique(cov_pats$PATIENT_NUM)))
          
          ##remove patients who only had 1 sequence to redo the sparsity screen
          corseq_sub_count_inft_single <- corseq_sub_count_inf %>%
            filter(count==1)
          
          corseq_sub_inf <- subset(corseq_sub_inf,!(paste0(corseq_sub_inf$patient_num,corseq_sub_inf$sequence) %in%
                                                      paste0(corseq_sub_count_inft_single$patient_num,corseq_sub_count_inft_single$sequence)))
 
          ###now sparsity screening
          ##sparsity screen
          corseq_sub_inf_count <- corseq_sub_inf %>%
            dplyr::group_by(sequence) %>%
            dplyr::summarise(count=length(unique(patient_num)))%>%
            filter(count > sparsity*pat_num_inf)
          
          
          corseq_sub_inf <- subset(corseq_sub_inf,corseq_sub_inf$sequence %in% c(unique(corseq_sub_inf_count$sequence)))
          
          
          #### the 2-month requirement
          
          
          ###calculate the min duration for each patient and sequence
          corseq_sub_inf_min <- corseq_sub_inf %>%
            dplyr::group_by(patient_num,sequence) %>%
            dplyr::summarise(min_duration=min((duration)))
          
          
          corseq_sub_inf <- merge(corseq_sub_inf,corseq_sub_inf_min,by=c("patient_num","sequence"))
          
          corseq_sub_inf <- subset(corseq_sub_inf,corseq_sub_inf$duration == corseq_sub_inf$min_duration | 
                                     corseq_sub_inf$duration >= corseq_sub_inf$min_duration+2)
          
          corseq_sub_inf_count <- corseq_sub_inf %>%
            dplyr::group_by(patient_num,sequence) %>%
            dplyr::summarise(count=length((patient_num)))
          
          ##remove patients who only had 1 sequences to redo the sparsity screen
          corseq_sub_inf_count_single <- corseq_sub_inf_count %>%
            filter(count==1)
          
          
          corseq_sub_inf <- subset(corseq_sub_inf,!(paste0(corseq_sub_inf$patient_num,corseq_sub_inf$sequence) %in%
                                                      paste0(corseq_sub_inf_count_single$patient_num,corseq_sub_inf_count_single$sequence)))
          
          ###see if we pass the sparsity screen
          
          corseq_sub_inf_count <- corseq_sub_inf %>%
            dplyr::group_by(sequence) %>%
            dplyr::summarise(count=length(unique(patient_num)))%>%
            filter(count > sparsity*pat_num_inf)
          
          
          corseq_sub_inf <- subset(corseq_sub_inf,corseq_sub_inf$sequence %in% c(unique(corseq_sub_inf_count$sequence)))
          
          
          ###the remaining sequences are the ones of interest. Now extract J from them
          J_inf <- subset(J,J$candidate_sequences  %in% c(unique(corseq_sub_inf_count$sequence)))
          J_inf <- merge(J_inf,db$phenxLookUp,by.x = "endPhenx",by.y ="num_Phenx" )
          
          J_low[[inf]] <- J_inf
          
          rm(J_inf,pat_sub_inf,corseq_sub_count_inf,corseq_sub_count_inft_single,
             corseq_sub_inf_count_single,corseq_sub_inf_count,corseq_sub_inf_min,pat_num_inf)
          gc()
        }
      }
      
      
      
    }, 
    error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
  }
  
  
  J <- data.table::rbindlist(J_low)
  J <- dplyr::distinct(dplyr::select(J,endPhenx,phenx), .keep_all = TRUE)
  J_loop[[i]] <- J

  jFileName <- paste0(jBaseFileName, i, ".RData")
  save(J,file=jFileName)
  rm(j,jFileName, infections, corseq_sub)
  gc();
  
  },
  error = function(foll) {cat("ERROR in chunk ",i, ": ", conditionMessage(foll), "\n")})
}
  
J <- data.table::rbindlist(J_loop)

###we will limit candidate Js to those that came up in half of the chunks
J <- J %>%
  dplyr::group_by(endPhenx,phenx) %>%
  dplyr::summarise(count_perc=length((endPhenx))/numOfChunks)%>%
  filter(count_perc > 0.5) %>%
  dplyr::select(endPhenx,phenx)

save(J,file=paste0(outputDirectory,"/J.RData"))
}