
library(data.table)
library(bnlearn)
library(arules)
library(missForest)
library(parallel)
library(caret)
library(ggplot2)
library(gtools)
library(plyr)
library(stringr)
library(doMC)
library(readxl)
library(dplyr)
library(missForest)

##extract all the convertors and the de-novo subjects##
aggDiagnosis <- ddply(adnimerge, .(PTID), summarize, DX = toString(DX))
patientid <- c()
for(i in 1:nrow(aggDiagnosis)){
  if(grepl("Dementia", aggDiagnosis[i,2])){
    patientid <- c(patientid, aggDiagnosis[i,1])
  }
}

##dataFrame with these convertors and denovo ids##
convDenovoDf <- subset(adnimerge, adnimerge$PTID %in% patientid)

##generalized function to extract demographic columns from a data frame##
func.demog.id <- function(input_data){
  ids <- grep(".*ID$", names(input_data), value=TRUE)
  if(length(ids)>1){
    id_idx1 <- which(colnames(input_data) %in% ids[1])
    id_idx2 <- which(colnames(input_data) %in% ids[2])
  }
  time_period <- c("M", "Month", "month", "Time", "time", "Years")
  time <-  colnames(input_data)[which(names(input_data) %in% time_period)]
  if(length(time)>1){
    time_idx <- which(colnames(input_data) %in% time[1])
  }
  visit_id <- c("VISCODE", "visit", "VISIT")
  visit <-  colnames(input_data)[which(names(input_data) %in% visit_id)]
  visit_idx <- which(colnames(input_data) %in% visit)
  prm_key <- c("primary_key", "primary_key1", "primary_key2")
  pk <- colnames(input_data)[which(names(input_data) %in% prm_key)]
  pk_id <- which(colnames(input_data) %in% pk)
  age <- c("Age", "age", "AGE")
  age_col <- colnames(input_data)[which(names(input_data) %in% age)]
  age_idx <- which(colnames(input_data) %in% age_col)
  diagnosis_bl <- c("DX.bl", "Diagnosis.bl", "Stage.bl")
  diag_bl <-  colnames(input_data)[which(names(input_data) %in% diagnosis_bl)]
  diag_bl_idx <- which(colnames(input_data) %in% diag_bl)
  diagnosis <- c("DX", "Diagnosis", "Stage")
  diag <-  colnames(input_data)[which(names(input_data) %in% diagnosis)]
  diag_idx <- which(colnames(input_data) %in% diag)
  education <- c("Education", "EDUCATION", "PTEDUCAT")
  edu <-  colnames(input_data)[which(names(input_data) %in% education)]
  edu_idx <- which(colnames(input_data) %in% edu)
  gender <- c("GENDER", "PTGENDER", "Gender") 
  gen <-  colnames(input_data)[which(names(input_data) %in% gender)]
  gen_idx <- which(colnames(input_data) %in% gen)
  ethnicity <- c("PTETHCAT", "ethinicity", "ethics")
  eth <-  colnames(input_data)[which(names(input_data) %in% ethnicity)]
  eth_idx <- which(colnames(input_data) %in% eth)
  race <- c("PTRACCAT", "RACE", "race", "Race")
  race_col <-  colnames(input_data)[which(names(input_data) %in% race)]
  race_idx <- which(colnames(input_data) %in% race_col)
  marital_status <- c("PTMARRY", "MARRIED", "Married", "married", "MARRY")
  marital <-  colnames(input_data)[which(names(input_data) %in% marital_status)]
  marital_idx <- which(colnames(input_data) %in% marital)
  demog_cols <- c(ids, time,visit, pk, age_col, diag_bl, diag, edu, gen, eth, race_col, marital)
  if(length(pk_id == 2)){
    pk_upd <- c(pk_id[1], pk_id[2])
  }
  else{
    pk_upd <- pk_id
  }
  demog_col_id <- c(id_idx1, id_idx2, time_idx,visit_idx,pk_upd, age_idx, 
                    diag_bl_idx, diag_idx, edu_idx, gen_idx, eth_idx, race_idx, marital_idx)
  return(list(demog_col_id, demog_cols))
}

demog_list <- func.demog.id(convDenovoDf)  
demog_cols <- unlist(demog_list[2])

##Vector for each type of biomarkers##
mri <- c("FDG", "AV45","Ventricles", "Hippocampus", "WholeBrain", "Entorhinal", "Fusiform", "MidTemp", "ICV")
ravlt <- grep("^RAVLT.*[^bl]$", colnames(adnimerge), value = TRUE)
ecog <- grep("^Ecog.*[^bl]$",  colnames(adnimerge), value = TRUE)
cognitionFunctional <- c("CDRSB", "MMSE", "MOCA", "ADAS11", "ADAS13", "FAQ")
csf <- c("ABETA", "TAU", "PTAU")
apoe <- c("APOE4")
fdgav <- c("FDG","AV45")
mostRelevantCols <- c(ptid,mri, ravlt, ecog, cognitionFunctional, csf, apoe,fdgav)

##subset the data frame according to the most relevant columns you need##
selectedFeatures <- convDenovoDf[,c(demog_cols, mostRelevantCols)]

# CSF biomarkers values are represented as characters due to presence of "<" or ">" symbols
#convert them to numeric
#print("converting to Numeric")
selectedFeatures$ABETA <- gsub('(<|>)','',selectedFeatures$ABETA)
selectedFeatures$PTAU <- gsub('(<|>)','',selectedFeatures$PTAU)
selectedFeatures$TAU <- gsub('(<|>)','',selectedFeatures$TAU)

selectedFeatures$ABETA <- as.numeric(selectedFeatures$ABETA)
selectedFeatures$PTAU <- as.numeric(selectedFeatures$PTAU)
selectedFeatures$TAU <- as.numeric(selectedFeatures$TAU)

##Removing rows with no diagnostic information
rem.rid <- c()
for(i in 1:nrow(selectedFeatures)){
  if(is.na(selectedFeatures[i,"DX"]) == TRUE){
    if(sum(is.na(selectedFeatures[i,])) != (ncol(selectedFeatures)-10)){
      rem.rid <- c(rem.rid, i)
    }
  }
}
selectedFeatures <- selectedFeatures[-rem.rid,] 

##removing baseline columns as they are repititive##
blCols <- grep("bl", colnames(selectedFeatures), value = TRUE)
selectedFeatures[,blCols] <- NULL

##taking out demographic features##
demogs <- grep("(PTID$|AGE|^PT[^ID|AU].*)", colnames(selectedFeatures), value = TRUE)
demoCols <- selectedFeatures[,demogs]
demoCols <- unique(demoCols)

##converting the data frame to wider format##
##extracting the columns that need to be converted to wider format##
selectedCols <- setdiff(colnames(selectedFeatures), demogsAndOthers)
convDenovoSelectedDf <- selectedFeatures[,selectedCols]
reshapeDf = reshape(convDenovoSelectedDf,idvar='PTID.1',timevar="VISCODE",dir='w')
names(reshapeDf)[names(reshapeDf) == 'PTID.1'] <- 'PTID'
copyreshape <- reshapeDf

# baseline features with less than 50% missing value
baselineFeatures <- reshapeDf[,grep("bl", colnames(reshapeDf), value = TRUE)]
baselineFeaturesUpd = baselineFeatures[ , -which(colMeans(is.na(baselineFeatures)) > 0.5)]

##removing the columns having more than 50% missing data from the main data frame##
remove50PercMissingCols <- setdiff(colnames(baselineFeatures),colnames(baselineFeaturesUpd))
toMatch <- (unlist(strsplit(remove50PercMissingCols, "\\.bl")))
matches <- unique(grep(paste(toMatch,collapse="|"), colnames(reshapeDf), value=TRUE))
reshapeDf[,matches] <- NULL

#####create a variable, store similar features in it and add auxiliary column

#working with Brain regions
brainReg <- c("Entorhinal","Fusiform","MidTemp","Hippocampus","ICV","Ventricles","WholeBrain")

#get names of the 7 brain region columns from reshapeDf
brainVolCols <-  grep("PTID|Entorhinal|Fusiform|MidTemp|Hippocampus|ICV|Ventricles|WholeBrain",colnames(reshapeDf), value = TRUE)


##subset the columns for brain regions
brainVolDf <- reshapeDf[,brainVolCols]

#Get only those months that have less than 50% missing data, same applies for other longitudinal data. 
#This missingness pattern had already been studied before

brVolmons <- grep("(PTID|bl|m06|m12$|m24)",colnames(brainVolDf),value=TRUE)
brainVolDf <- brainVolDf[,brVolmons]

#split the column name. i.e. Ventricles.m24 will be separated as "Ventricles" and "m24"
monthsForBrainVol <- as.vector(strsplit((grep("Entorhinal|Fusiform|MidTemp|Hippocampus|ICV|Ventricles|WholeBrain",colnames(brainVolDf), value = TRUE)), "\\."))

bvolmons <- c()
for(i in 1:length(monthsForBrainVol)){
  bvolmons <- c(bvolmons, monthsForBrainVol[[i]][2])
}

bvolmons <- unique(bvolmons)

#adding auxiliary variable, the code goes through all the 7 brain regions of a given visit and checks where all values are missing or not. 
#if all values are missing flag is set to 1 else 0
for (i in 1:length(bvolmons)){
  
  ent <- paste(brainReg[1],".",bvolmons[i],sep="")  #Entorhinal.bl
  fus <- paste(brainReg[2],".",bvolmons[i],sep="")
  mid <- paste(brainReg[3],".",bvolmons[i],sep="")
  hip <- paste(brainReg[4],".",bvolmons[i],sep="")
  icv <- paste(brainReg[5],".",bvolmons[i],sep="")
  ven <- paste(brainReg[6],".",bvolmons[i],sep="")
  wbr <- paste(brainReg[7],".",bvolmons[i],sep="")
  
  allvol <- c(ent,fus,mid,hip,icv,ven,wbr)
  
  bvolaux <- paste("allbrainvol",".",bvolmons[i],".aux",sep="")
  
  brainVolDf[,bvolaux] <- 0
  
  count = 0
  for (j in 1:nrow(brainVolDf)){
    
    for (k in 1:length(allvol)){
      
      if(is.na(brainVolDf[,allvol[k]][j])){
        count = count + 1
      }
      
    }  
    if (count == 7){
      brainVolDf[,bvolaux][j] <- 1
      #for(i in 1:length(allvol)){
      #dfDXs[,allvol[i]][j] <- 0
      #}
    }  
    count = 0  
    
  }
}

######################

#Working with Cognition tests and adding auxiliary variables

#acquire all columns names of Cog. test
allCogTests <- grep("(CDRSB|ADAS11|ADAS13|MMSE|FAQ|RAVLT)",colnames(reshapeDf),value = TRUE)

##subset the columns for brain regions
CogScore <- reshapeDf[,allCogTests]

MonCog <- grep("(PTID|bl|m06|m12$|m18|m24|m36|m48|m60)",colnames(CogScore),value=TRUE)
CogScore <- CogScore[,MonCog]

#get cog Names

#remove all the months i.e avoid everything after the last "."
CogNames <- sub(".[^.]*$","",allCogTests)
CogNames <- unique(CogNames)

#get months
CogMon <- c("bl","m06","m12","m18","m24","m36","m48","m60")



for (i in 1:length(CogMon)){
  
  allCogTests <- grep("(CDRSB|ADAS11|ADAS13|MMSE|FAQ|RAVLT)",colnames(CogScore),value = TRUE)
  CogNames <- sub(".[^.]*$","",allCogTests)
  CogNames <- unique(CogNames)
  
  adas11 <- paste(CogNames[1],".",CogMon[i],sep="")  #"ADAS11.bl" 
  adas13 <- paste(CogNames[2],".",CogMon[i],sep="") 
  cdrsb <- paste(CogNames[3],".",CogMon[i],sep="") 
  faq <- paste(CogNames[4],".",CogMon[i],sep="") 
  mmse <- paste(CogNames[5],".",CogMon[i],sep="") 
  ravltFor <- paste(CogNames[6],".",CogMon[i],sep="") 
  ravltImm <- paste(CogNames[7],".",CogMon[i],sep="") 
  ravltLearn <- paste(CogNames[8],".",CogMon[i],sep="") 
  ravltPerFor <- paste(CogNames[9],".",CogMon[i],sep="") 
  
  CogNames <- c(adas11,adas13,cdrsb,faq,mmse,ravltFor,ravltImm,ravltLearn,ravltPerFor)
  
  CogAux <- paste("CogScore",".",CogMon[i],".aux",sep="")
  
  CogScore[,CogAux] <- 0
  
  count = 0
  for (j in 1:nrow(CogScore)){
    
    for (k in 1:length(CogNames)){
      
      if(is.na(CogScore[,CogNames[k]][j])){
        count = count + 1
        #print(count)
      }
      
    }  
    if (count == 9){
      
      CogScore[,CogAux][j] <- 1
      for(i in 1:length(CogNames)){
        CogScore[,CogNames[i]][j] <- 0
      }
    }
    count = 0  
    
  }
}

CogScore <- cbind(reshapeDf$PTID,CogScore)
names(CogScore)[names(CogScore) == 'reshapeDf$PTID'] <- 'PTID'
####

#csf

csfBioM <- c("ABETA","PTAU","TAU")
getCSF <- grep("PTID|ABETA.bl|PTAU.bl|TAU.bl",colnames(reshapeDf),value = TRUE)
csfmonths <- c("bl")
getCSF <- reshapeDf[,getCSF] 

for (i in 1:length(csfmonths)){
  abeta <- paste(csfBioM[1],".",csfmonths[i],sep="")  #ABETA.bl
  #print(abeta)
  ptau <- paste(csfBioM[2],".",csfmonths[i],sep="")
  #print(ptau)
  tau <- paste(csfBioM[3],".",csfmonths[i],sep="")
  #print(tau)
  
  CSFaux <- paste("csf",".",csfmonths[i],".aux",sep="")
  #print(CSFaux)
  
  getCSF[,CSFaux] <- 0
  
  for (j in 1:nrow(getCSF)){
    if ((is.na(getCSF[,abeta][j]))&(is.na(getCSF[,ptau][j]))&(is.na(getCSF[,tau][j]))){
      getCSF[,CSFaux][j] <- 1
      
    }  
  }
}



################
#add FDG and AV45 from baseline
FDGAV45 <- grep("(PTID|FDG.bl|AV45.bl)",colnames(reshapeDf),value = TRUE)
FDGAV45 <- reshapeDf[,FDGAV45]
colnames(FDGAV45) <- c("PTID","FDG.bl","AV45.bl")

######
###################


####other 68 brain regions 
#import 68 brain region values
#EMC_ADNI_FS60_Phenotypes_Desikan_20180219 <- read_excel("~/Documents/PhDWork/BayesianNetworkAD/EMC_ADNI_FS60_Phenotypes_Desikan_20180219.xlsx")
volumeRegion <- as.data.frame(EMC_ADNI_FS60_Phenotypes_Desikan_20180219)
names(volumeRegion)[1] <- "PTID"
brainRegion <- c(colnames(volumeRegion)[grep("Left|Right", colnames(volumeRegion))])
volumeRegion <- volumeRegion[,c("PTID",brainRegion)]
#volumeRegion <- setDT(volumeRegion, keep.rownames = TRUE)[]


snps.dat.full = c(colnames(dat.full)[grep("(PTID|rs[0-9]+|APOE4.0)", colnames(dat.full))])
pathways = c(colnames(dat.full)[grep("PTID|^.*Homo.sapiens|\\.[a-z]+$|\\.[A-Z][a-z]+$|^[A-Z][a-z]+$",colnames(dat.full))])
snp.pathway.data <- dat.full[,c(snps.dat.full, pathways)]
snp.pathway.data <- setDT(snp.pathway.data, keep.rownames = TRUE)[]
names(snp.pathway.data)[1] <- "PTID"


#merge all dataframes
BeforeImp <- merge(brainVolDf,CogScore,all.x = TRUE)
BeforeImp <- merge(BeforeImp,getCSF,all.x = TRUE)
BeforeImp <- merge(BeforeImp,snp.pathway.data,all.x = TRUE)
BeforeImp <- merge(BeforeImp,volumeRegion,all.x = TRUE)
BeforeImp <- merge(BeforeImp,FDGAV45,all.x = TRUE)
BeforeImp <- merge(BeforeImp,demoCols,all.x = TRUE)

#convert all aux columns to factors
AllAux <- grep("aux",colnames(BeforeImp),value=TRUE)
BeforeImp[,AllAux] <- lapply(BeforeImp[,AllAux],factor)

##Imputation step
Before.Imp <- missForest(BeforeImp[2:ncol(BeforeImp)], ntree=100, parallelize = "no",verbose=TRUE,maxiter=3) 

#extract imputed data from the missForest object
ExtImpADNI <- Before.Imp$ximp
ExtImpADNI <- cbind(BeforeImp$PTID,ExtImpADNI)
names(ExtImpADNI)[names(ExtImpADNI)=="PrepImp$PTID"] <- "PTID"

save.image(file = "ADNI_prepare.RData")
