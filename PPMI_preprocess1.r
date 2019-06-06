####################################################################
########## Setting up the PPMI data set ############################
########## Getting the important variables and imputing them #######
rm(list = ls())
setwd("C:/Users/ashar/Dropbox/PPMI")
library(PPMI)

################### Choosing patients with enrolled status #########
pd <- getCohort(cohort="PD", enrolledOnly=TRUE)


######## Extracting Variables ####################
##### UPDRS and Motor Scores ####################
updrs <- findVarByGroup("UPDRS")
var.updrs <- c('TD','PIGD','UPDRS1','UPDRS2','UPDRS3','UPDRS')
pd.updrs <- extractVariables(
  patients= pd,
  variables= var.updrs,
  events = c('BL','V01','V03','V04','V05','V06','V07','V08','V09','V10')
)
table.updrs <- do.call(cbind, pd.updrs) 


##### Non Motor Scores ####################
nm <- findVarByGroup("Non-Motor")
var.nm <- c('BJLOT','ESS','GDS','HVLTIR','HVLTDR','QUIP','RBD','SCOPA','SFT','STA','STAI.State','STAI.Trait')
pd.nm <- extractVariables(
  patients= pd,
  variables= var.nm,
  events = c('BL','V04','V06','V08','V10')
)
table.nm <- do.call(cbind, pd.nm) 

######### Imaging features ##############
##### Imaging features that we are interested in
img <- findVarByGroup("Imaging")
img.rel <- rownames(img)[18:28]
pd.img <- extractVariables(
  patients= pd,
  variables= img.rel,
  events = c('V04','V06','V10')
)
table.img <- do.call(cbind, pd.img) 

####### Extract CSF biomarkers ########
var.csf <- c('ABeta 1-42','CSF Alpha-synuclein','pTau','tTau')
pd.csf <- extractVariables(
  patients= pd,
  variables= var.csf,
  events = c('BL','V02','V04')
)
table.csf <- do.call(cbind, pd.csf)

##### Extract RNA data #############
var.rna <- c('DHPR','DJ-1','FBXO7-001','FBXO7-005','FBXO7-007','FBXO7-008','FBXO7-010','GLT25D1','GUSB','MON1B','RPL13','SNCA-007','SNCA-3UTR-1','SNCA-3UTR-2','SNCA-E3E4','SNCA-E4E6','SOD2','SRCAP','UBC','ZNF160','ZNF746')
pd.rna <- extractVariables(
  patients= pd,
  variables= var.rna,
  events = c('BL')
)
table.rna <- do.call(cbind, pd.rna)

######### Cognitive scores #####################

####### MOCA scores ############################
mc.rel <- c('MCATOT')
pd.mc <-  extractVariables(
  patients= pd,
  variables= mc.rel,
  events = c('V04','V06','V08','V10')
)

table.moca <- do.call(cbind, pd.mc)
  
  
############################################################
####### Hoehn Yahr Scale, Modified Schwab England Scale ####
###########################################################
hm.rel <- c('NHY','MSEADLG')
pd.hm <-  extractVariables(
  patients= pd,
  variables= hm.rel,
  events = c('BL','V01','V02','V03','V04','V05','V06','V07','V08','V09','V10')
)

table.hm <- do.call(cbind, pd.hm)

######## Storing Data ################################
##### Combine the different variables ###########
###### Join the tables ############################
table.overall <- cbind(table.updrs, table.nm, table.img, table.hm, table.moca, table.csf, table.rna)

####### Missing data imputation #####

###### Use the R-package missForest ##############
###### Set seed ##################################
set.seed(81)
library('missForest')


pd.mis <- as.matrix(table.overall)

###### Carry out the imputation############
pd.imp <- missForest(pd.mis, maxiter = 100)
######## Assign the matrix to a new matrix ###
pd.com <- pd.imp$ximp

################### Extract relavent tables ######
################ Main Variables #################
table.main <- cbind(pd, table.updrs, table.nm)
table.main.com <- cbind(pd, pd.com[,1:120])
####### Save the results ##################
write.csv(table.main, 'PPMI_data.csv')
write.csv(table.main.com, 'PPMI_complete.csv')

################ Other Variables #################
table.other <- cbind(pd, table.img, table.csf, table.moca, table.nm, table.rna, table.hm)
table.other.com <- cbind(pd, pd.com[,121:153], pd.com[,180:191], pd.com[,176:179],pd.com[,61:120], pd.com[,192:212], pd.com[,154:175])
####### Save the results ##################
write.csv(table.other.com, 'PPMI_othervariables.csv')
write.csv(table.other, 'PPMI_OtherVariables_NAsincl.csv')






