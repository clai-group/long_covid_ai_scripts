
###this script implements WHO definition of long COVID
###this is step 2 to calculate temporal correlations for candidates Js
### we will use both  cases and controls



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
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools)

#setup parallel backend to use many processors
cores<-detectCores()
cl <- parallel::makeCluster(3) #not to overload your computer
doParallel::registerDoParallel(cl)

####load the cases data
dbmart_cases_map_ccsr <- read_csv("P:/PASC/data/dbmart_cases_map_ccsr_modified.csv")
dbmart_control_map_ccsr <- read_csv("P:/PASC/data/dbmart_control_map_ccsr_modified.csv")
dbmart_cases_map_ccsr <- rbind(dbmart_cases_map_ccsr,dbmart_control_map_ccsr)
rm(dbmart_control_map_ccsr)
dbmart_cases_map_ccsr <- subset(dbmart_cases_map_ccsr,!(dbmart_cases_map_ccsr$type %in% c("unrelated") | is.na(dbmart_cases_map_ccsr$type)))

##load the cov_pat incident level data
load("P:/PASC/data/cov_pats.RData")

colnames(dbmart_cases_map_ccsr)[13] <- "phenx"

length(unique(dbmart_cases_map_ccsr$phenx))

dbmart <- dplyr::select(dbmart_cases_map_ccsr,PATIENT_NUM,START_DATE,phenx)
rm(dbmart_cases_map_ccsr);gc()

# cov_pats$phenx <- "COVID1"
cov_pats$phenx <- paste0("COVID",cov_pats$infection_seq)
dbmart <- rbind(dbmart,dplyr::select(cov_pats,PATIENT_NUM,START_DATE,phenx))
colnames(dbmart) <- c("patient_num","start_date","phenx")

load("P:/PASC/data/phenxlookup.RData")

### we need to make sure that the data elements

dbmart <- subset(dbmart,dbmart$phenx %in% phenxlookup$phenx)
db <- tSPMPlus::transformDbMartToNumeric(dbmart)

##moved loading Js
#load number of chunks 
load(file="P:/PASC/data/numOfChunks_J.RData")
load(file="P:/PASC/data/numOfChunks_1.RData")




# phenxlookup <- db$phenxLookUp
# save(phenxlookup,file="P:/PASC/data/phenxlookup.RData")


######################################
######################################
## create adaptive chunks
######################################

buffer<- 100000000 #extra buffer in bytes find good value *a higher value here decrease the chunk size, but left more memory to work with after the sequencing, the chunk size is calculated on the memory consumption of the non-sparse sequences, reduce buffer 
dbmart_adapt <- tSPMPlus::splitdbMartInChunks(db$dbMart, includeCorSeq = TRUE, buffer = buffer)
dbmart_num <- dbmart_adapt$chunks[[1]]

#dbmart_num <- db$dbMart
sparsity = 0.005
numOfThreads = detectCores()


endPhenx = c(J$endPhenx) #TODO look up id for covid phenx in db$phenxLookUp
temporalBucket =  c(0,1,3)
minDuration = 1 #techical parameter, ignore for now
bitShift = 0 #techical parameter, ignore for now
lengthOfPhenx = 7 #techical parameter
storeSequencesDuringCreation = FALSE #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity


##old version!
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

#####calculate correlations
corseq <- corseq %>%
  dplyr::group_by(patient_num,sequence,endPhenx,durationBucket) %>%
  dplyr::summarise(count=length((patient_num)))

##currently buggy, use old version!
# corseq <- tSPMPlus::getSequencesWithEndPhenx(dbmart_num,
#                                              bitShift,
#                                              lengthOfPhenx,
#                                              temporalBucket,
#                                              endPhenx,
#                                              includeCorBuckets = TRUE,
#                                              minDuration,
#                                              storeSequencesDuringCreation,
#                                              numOfThreads = numOfThreads,
#                                              sparsityValue = sparsity,
#                                              returnSummary = TRUE)

corseq$value.var <- 1

corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))

corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)

dat <- data.table(corseq[c("patient_num","endPhenx","value.var","startPhen_dur")])
wide.dat.start <- dcast.data.table(dat, patient_num ~ startPhen_dur, value.var="value.var", fun=length)



for (i in seq(1:numOfChunks_J)) {
  gc()
  tryCatch({
    ##load Js
    jFileName = paste0("P:/PASC/data/J_chunk_", i, ".RData")
    load(file=jFileName)
    end <- c(unique(J$endPhenx))

    #setup parallel backend to use many processors
    # PARALELIZING THE CORRELATION CALCULATIONS
    corrs <- foreach(j = 1:length(end),#
                     .combine='rbind', 
                     .multicombine=TRUE,
                     .packages = c("data.table")) %dopar% {
                       tryCatch({
                         dat.j <- data.table(corseq[corseq$endPhenx %in% end[j],c("patient_num","endPhenx","value.var")])
                         lab.j <- dcast.data.table(dat.j, patient_num ~ endPhenx, value.var="value.var", fun=length)
                         colnames(lab.j)[2] <- "label"
                         wide.j <- merge(wide.dat.start,lab.j,by="patient_num",all.x = T)
                         wide.j[is.na(wide.j)] <- 0
                         wide.j$patient_num <- NULL
                         
                         
                         out <- apply(wide.j[, -("label")], 2, cor.test, wide.j$label, method="spearman")
                         p <- data.frame(sapply(out, "[[", "p.value"))
                         p$startPhen_dur <- rownames(p)
                         rownames(p) <- NULL
                         colnames(p)[1] <- "p.value"
                         p$p.adjust <- p.adjust(p$p.value, method = "bonferroni", n = nrow(p))# consider "holm" or "hochberg"
                         
                         rho <- data.frame(sapply(out, "[[", "estimate"))
                         rho$startPhen_dur <- rownames(rho)
                         rownames(rho) <- NULL
                         colnames(rho)[1] <- "rho"
                         rho$startPhen_dur <- substr(rho$startPhen_dur,1,nchar(rho$startPhen_dur)-4)
                         cor.j <- merge(rho,p,by="startPhen_dur")
                         cor.j$rho.abs <- abs(cor.j$rho)
                         cor.j$endPhenx <- end[j]
                         
                         
                         rm(rho,p);gc()
                         
                         # write.csv(cor.j,file = paste0("P:/PASC/data/cors/",j,".csv"))
                         
                         
                         cor.j
                         
                       },
                       error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
                     }
    
    
    corrs <- merge(corrs,db$phenxLookUp,by.x = "endPhenx",by.y ="num_Phenx" )
    
    corrs$startPhen <- sub("\\-.*", "", corrs$startPhen_dur)
    corrs$sequence <- ifelse (corrs$endPhenx <10,paste0(corrs$startPhen,"000000",corrs$endPhenx),NA)
    corrs$sequence <- ifelse(corrs$endPhenx <100 & corrs$endPhenx >9 ,paste0(corrs$startPhen,"00000",corrs$endPhenx),corrs$sequence)
    corrs$sequence <- ifelse(corrs$endPhenx <1000 & corrs$endPhenx >99 ,paste0(corrs$startPhen,"0000",corrs$endPhenx),corrs$sequence)
    
    corrs$sequence <- as.numeric(corrs$sequence)
    
    corrsFileName <- paste0("P:/PASC/data/corrs_chunk_", i, ".RData")
    save(corrs,file=corrsFileName)
    rm(corrs, end)
    gc()
  },
  error = function(foll) {cat("ERROR in chunk ",i, ": ", conditionMessage(foll), "\n")})
  
}



