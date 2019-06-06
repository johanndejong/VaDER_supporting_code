library(data.table)
# library(bnlearn)
# library(arules)
library(missForest)
library(parallel)
library(caret)
library(ggplot2)
library(gtools)
library(plyr)
library(stringr)
library(ADNIMERGE)
library(dummy)
library(varhandle)
# library(doMC)

missingness_allowed <- 0.8
include_cols <- c("m06", "m12", "m24", "m36", "m48", "m60")

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
mostRelevantCols <- c(mri, ravlt, ecog, cognitionFunctional, csf, apoe)

##subset the data frame according to the most relevant columns you need##
selectedFeatures <- convDenovoDf[,c(demog_cols, mostRelevantCols)]

# CSF biomarkers values are represented as characters due to presence of "<" or ">" symbols
#convert them to numeric
print("converting to Numeric")
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
demogs <- grep("^PT[^ID|AU].*", colnames(selectedFeatures), value = TRUE)
demogsAndOthers <- c(demogs, "AGE", "M", "PTID", "Month", "APOE4")
##converting the data frame to wider format##
##extracting the columns that need to be converted to wider format##
selectedCols <- setdiff(colnames(selectedFeatures), demogsAndOthers)
convDenovoSelectedDf <- selectedFeatures[,selectedCols]
reshapeDf = reshape(convDenovoSelectedDf,idvar='RID',timevar="VISCODE",dir='w')

# baseline features with less than 50% missing value
baselineFeatures <- reshapeDf[,grep("bl", colnames(reshapeDf), value = TRUE)]
baselineFeaturesUpd = baselineFeatures[ , -which(
  colMeans(is.na(baselineFeatures)) > 1 - missingness_allowed | 
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(baselineFeatures))
)]

##removing the columns having more than 50% missing data from the main data frame##
remove50PercMissingCols <- setdiff(colnames(baselineFeatures),colnames(baselineFeaturesUpd))
toMatch <- (unlist(strsplit(remove50PercMissingCols, "\\.bl")))
matches <- unique(grep(paste(toMatch,collapse="|"), colnames(reshapeDf), value=TRUE))
reshapeDf[,matches] <- NULL
##1.csf data frame##
csfDf <- reshapeDf[,grep("ABETA|TAU|PTAU", colnames(reshapeDf), value = TRUE)]

##Volumetric data frame##
volumeDf <- reshapeDf[,grep("Entorhinal|Fusiform|MidTemp|Hippocampus|ICV|Ventricles|WholeBrain", colnames(reshapeDf),
                            value = TRUE)]

##Cognitive test data frame##
cogTestDf <- reshapeDf[,grep("CDRSB|ADAS11|ADAS13|MMSE|FAQ|RAVLT", colnames(reshapeDf), value = TRUE)]

##fdg##
fdgDf <- reshapeDf[,grep("FDG", colnames(reshapeDf), value = TRUE)]

##Diagnostic data frame##
dxDf <- reshapeDf[,grep("PTID|DX", colnames(reshapeDf), value = TRUE)]

# Remove columns with less than 50% missing value
csfDf = csfDf[ , which(
  colMeans(is.na(csfDf)) < missingness_allowed |
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(csfDf))
)]
volumeDf = volumeDf[ , which(
  colMeans(is.na(volumeDf)) < missingness_allowed |
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(volumeDf))
)]
cogTestDf = cogTestDf[ , which(
  colMeans(is.na(cogTestDf)) < missingness_allowed | 
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(cogTestDf))
)]
fdgDf = fdgDf[ , which(
  colMeans(is.na(fdgDf)) < missingness_allowed |
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(fdgDf))
)]
dxDf = dxDf[ , which(
  colMeans(is.na(dxDf)) < missingness_allowed |
  grepl(paste(sprintf("\\.%s$", include_cols), collapse = "|"), colnames(dxDf))
)]

##extract data related to above features##
demogsAndOthersDf <- convDenovoDf[rownames(dxDf),demogsAndOthers]

##Add bl extension to the demogsAndOthers features##
colnames(demogsAndOthersDf) <- paste(colnames(demogsAndOthersDf) , ".bl", sep = "")

##Combine all the data frames##
allFeatures <- cbind.data.frame(demogsAndOthersDf, csfDf, volumeDf, cogTestDf, fdgDf, dxDf)

# # Add group name to columns
# colnames(csfDf) = paste0("csf_",colnames(csfDf))
# colnames(volumeDf) = paste0("brain_",colnames(volumeDf))
# colnames(cogTestDf) = paste0("Cog_",colnames(cogTestDf))

dt <- data.table(allFeatures)

# m <- colSums(is.na(dt)) / nrow(dt)
# m[m > 0.8]

# remove all-NA columns, and keep only bl if the all-NA column is part of 
# a series
remove <- colnames(dt)[colSums(is.na(dt)) == nrow(dt)]
remove <- unlist(lapply(gsub("\\.m[0-9]+$", "", remove), function(pattern) {
  grep(sprintf("^%s\\.m[0-9]+$", pattern), colnames(dt))
}))
if (length(remove) > 0) {
  dt <- dt[, -remove, with = FALSE]
}

colnames(dt) <- gsub("\\.bl", ".m00", colnames(dt))
# colnames(dt)[colnames(dt) == "fdgDf"] <- "fdgDf.m00"

dt$Month.m00 <- dt$M.m00 <- NULL


# visit-wise imputation per group
visits <- gsub("^\\.", "", str_extract(colnames(dt), "\\.bl$|\\.m[0-9]+$"))
groups <- rep("other", ncol(dt))
ii <- match(colnames(csfDf), colnames(dt))
groups[ii] <- "csf"
ii <- match(colnames(volumeDf), colnames(dt))
groups[ii] <- "volume"
ii <- match(colnames(cogTestDf), colnames(dt))
groups[ii] <- "cog"

# the auxiliary variables that reflect MNAR
dt_aux <- do.call("cbind", lapply(split(1:ncol(dt), list(groups, visits)), function(ii) {
  if (length(ii) > 0) {
    rowSums(is.na(dt[, ii, with = FALSE])) == length(ii)
  } else {
    rep(FALSE, nrow(dt))
  }
}))
dt_aux <- dt_aux[,colSums(dt_aux) != 0 & colSums(dt_aux) != nrow(dt_aux)]
# group "other" has too few time vars: treat them as MAR (not MNAR)
dt_aux <- dt_aux[, !grepl("^other\\.", colnames(dt_aux))]
dt_aux <- data.table(dt_aux)
colnames(dt_aux) <- sprintf("aux_%s", colnames(dt_aux))

# impute
dt_raw <- dt
ptid <- dt$PTID.m00
dt$PTID.m00 <- NULL

ii <- which(lapply(lapply(dt, class), tail, 1) %in% c("factor", "character"))
dt_nondummy <- dt[, ii, with = FALSE]
dt_dummy <- dummy::dummy(dt_nondummy, int = TRUE)
colnames(dt_dummy) <- sprintf("dummy_%s", colnames(dt_dummy))
m <- as.matrix(cbind(dt[, -ii, with = FALSE], dt_dummy, dt_aux))

visits <- gsub("^\\.", "", str_extract(colnames(m), "\\.m[0-9]+"))
I <- split(1:ncol(m), visits)
dt_imp <- do.call("cbind", lapply(1:length(I), function(i) {
  ii <- I[[i]]
  cat(sprintf("########\n# i = %i\n########\n", i))
  missForest(m[, ii], ntree = 1000)$ximp
}))
dt_imp <- data.table(dt_imp[, !grepl("dummy_", colnames(dt_imp))])
dt_imp <- cbind(dt_nondummy, dt_imp)
dt_imp$PTID.m00 <- ptid
dt_imp <- dt_imp[, !grepl("^aux_", colnames(dt_imp)), with = FALSE]

dt_imp <- dt_imp[, !grepl("^DX\\.", colnames(dt_imp)), with = FALSE]
dt_raw <- dt_raw[, !grepl("^DX\\.", colnames(dt_raw)), with = FALSE]

# organize and sort the column names
s <- do.call("rbind", strsplit(colnames(dt_imp), "\\.m"))
I <- split(1:ncol(dt_imp), s[,1])
bl_only <- unlist(I[lengths(I) == 1])
dt1 <- dt_imp[, bl_only, with = FALSE]
colnames(dt1) <- gsub("\\.m00", "", colnames(dt1))
dt2 <- dt_imp[, -bl_only, with = FALSE]
dt2 <- dt2[, order(s[-bl_only,1], as.numeric(s[-bl_only,2])), with = FALSE]
dt_imp <- cbind(dt1, dt2)
s <- do.call("rbind", strsplit(colnames(dt_raw), "\\.m"))
I <- split(1:ncol(dt_raw), s[,1])
bl_only <- unlist(I[lengths(I) == 1])
dt1 <- dt_raw[, bl_only, with = FALSE]
colnames(dt1) <- gsub("\\.m00", "", colnames(dt1))
dt2 <- dt_raw[, -bl_only, with = FALSE]
dt2 <- dt2[, order(s[-bl_only,1], as.numeric(s[-bl_only,2])), with = FALSE]
dt_raw <- cbind(dt1, dt2)

# normalize w.r.t. baseline
dt_norm <- dt_imp
s <- do.call("rbind", strsplit(colnames(dt_norm), "\\.m"))
I <- split(which(s[,1] != s[,2]), s[s[,1] != s[,2], 1])
for (ii in I) {
  bl <- grep("\\.m00", colnames(dt_norm)[ii])
  d <- dt_norm[[ii[bl]]]
  for (i in ii) {
    dt_norm[[i]] <- (dt_norm[[i]] - d) / sd(d, na.rm = TRUE)
  }
}

save(dt_raw, dt_imp, dt_aux, file = file.path(dirout, "..", "ADNI.RData"))

