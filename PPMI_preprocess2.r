library(abind)
library(data.table)
library(matrixStats)

f_out <- file.path(dir_ppmi, "PPMI.RData")

# the locations of the files used in this script
dir_ppmi <- file.path("..", "data", "PPMI")
f_ppmi_imp <- file.path(dir_ppmi, "PPMI_complete.csv")
f_ppmi_ori <- file.path(dir_ppmi, "PPMI_data.csv")
f_other_imp <- file.path(dir_ppmi, "PPMI_othervariables.csv")
f_other_ori <- file.path(dir_ppmi, "PPMI_OtherVariables_NAsincl.csv")
f_genotype <- file.path(dir_ppmi, "from_Ping", "PPMI_combPD_genotype.csv")
f_cadd <- file.path(dir_ppmi, "from_Ping", "PPMI_combPD_pass_CADD.csv")
f_disgenet <- file.path(dir_ppmi, "from_holger", "CuratedSNPs_DisGeNET_PD.txt")
f_desc <- file.path(dir_ppmi, "OtherVariableDescription.csv")
f_patchar <- file.path(dir_ppmi, "PatientCharacteristics_PPMI.csv")

func <- function(f_in, impute) {
  pd.complete <- read.csv(file = f_in)
  n = dim(pd.complete)[1]
  colnames(pd.complete)
  
  time <- c(0,1,3,4,6,8,10,12,14,16)
  m = length(time)
  
  time2 <- c(0,4,8,12,16)
  m2 <- length(time2)
  
  ############### Get all the scores ####################
  pid <- pd.complete[[1]]
  td <- pd.complete[1:n,2:11]
  pigd <- pd.complete[1:n,12:21]
  updrs2 <- pd.complete[1:n,32:41]
  updrs3 <- pd.complete[1:n,42:51]
  updrs <-  pd.complete[1:n,52:61]
  updrs1 <- pd.complete[1:n,22:31]
  rbd <- pd.complete[1:n,92:96]
  ess <- pd.complete[1:n,67:71]
  scopa <- pd.complete[1:n,97:101]
  
  ######## Calculate the normalized version of the scores ########
  ############ Motor scores ###############
  z.td <- matrix(NA, nrow = n, ncol = m)
  z.pigd <- matrix(NA, nrow = n, ncol = m)
  z.updrs1 <- matrix(NA, nrow = n, ncol = m)
  z.updrs2 <- matrix(NA, nrow = n, ncol = m)
  z.updrs3 <- matrix(NA, nrow = n, ncol = m)
  z.updrs <- matrix(NA, nrow = n, ncol = m)
  ############ Non-motor scores #############
  z.rbd <-  matrix(NA, nrow = n, ncol = m2)
  z.ess <-  matrix(NA, nrow = n, ncol = m2)
  z.scopa <-  matrix(NA, nrow = n, ncol = m2)
  ########### Two normalized scores ###################################
  ############ First for the UPDRS scores #############################
  for ( i in 1:m){
    
    z.td[,i] <-  (td[,i] - td[,1])/sd(td[,1], na.rm = TRUE)
    z.pigd[,i] <-  (pigd[,i] - pigd[,1])/sd(pigd[,1], na.rm = TRUE)
    z.updrs1[,i] <-  (updrs1[,i] - updrs1[,1])/sd(updrs1[,1], na.rm = TRUE)
    z.updrs2[,i] <-  (updrs2[,i] - updrs2[,1])/sd(updrs2[,1], na.rm = TRUE)
    z.updrs3[,i] <-  (updrs3[,i] - updrs3[,1])/sd(updrs3[,1], na.rm = TRUE)
    z.updrs[,i] <-  (updrs[,i] - updrs[,1])/sd(updrs[,1], na.rm = TRUE)
    
  }
  ############ Second the Non-motor scores##############################
  for ( i in 1:m2){
    z.rbd[,i] <-  (rbd[,i] - rbd[,1])/sd(rbd[,1], na.rm = TRUE)
    z.ess[,i] <-  (ess[,i] - ess[,1])/sd(ess[,1], na.rm = TRUE)
    z.scopa[,i] <-  (scopa[,i] - scopa[,1])/sd(scopa[,1], na.rm = TRUE)
  }
  
  ###########################################################
  # Combine the motor and non-motor datasets, and interpolate
  motor <- list(z.td, z.pigd, z.updrs1, z.updrs2, z.updrs3, z.updrs)
  nonmotor <- list(z.rbd, z.ess, z.scopa)
  # Linearly interpolate nonmotor to cover same time points as motor
  # NOTE: the exact values do not matter, because they will zero-weighted
  # during model fitting.
  if (impute) {
    nonmotor <- lapply(nonmotor, function(Y) {
      t(apply(Y, 1, function(y) {
        approx(time2, y, xout = time, method = "linear")$y
      }))
    })    
  } else {
    nonmotor <- lapply(nonmotor, function(Y) {
      y <- matrix(NA, nrow = nrow(Y), ncol = length(time))
      y[,match(time2, time)] <- Y
      y
    })
  }
  # construct data array
  X <- do.call(abind, list(c(motor, nonmotor), along = 3))
  dimnames(X) <- list(
    pid,
    time,
    c("td", "pigd", "updrs1", "updrs2", "updrs3", "updrs", "rbd", "ess", "scopa")
  )
  X
}

X <- func(f_ppmi_imp, impute = TRUE)
MNAR <- func(f_ppmi_ori, impute = FALSE)
time2 <- c(0,4,8,12,16)
for (i in 1:dim(MNAR)[3]) {
  # MNAR[, colSums(!is.na(MNAR[,,i])) == 0, i] <- Inf # dummy
  mask <- is.na(MNAR[,,i])
  MNAR[,,i][mask] <- 1
  MNAR[,,i][!mask] <- 0
}


############################################3
############################################3
############################################3
############################################3
############################################3
############################################3

# # construct weight array
# mw <- array(1, c(dim(X)[1], dim(X)[2], length(motor)))
# nmw <- array(1, c(dim(X)[1], dim(X)[2], length(nonmotor)))
# nmw[,is.na(match(time, time2)),] <- 0
# W <- abind(mw, nmw, along = 3)
# dimnames(W) <- dimnames(X)

# combine with genetics data
g <- fread(f_genotype)
g$V1 <- gsub("_[A-Z].*$", "", g$V1)
cadd <- fread(f_cadd)
# g <- g[V1 %in% cadd$ID[cadd$CADD_PHRED > 10], ]
pid <- colnames(g)[-1]
rs <- g$V1
g$V1 <- NULL
g <- t(as.matrix(g))
colnames(g) <- rs
rownames(g) <- gsub("PPMISI", "", rownames(g))
G <- g
colnames(G) <- gsub("_[ACGT]+$", "", colnames(G))
rs <- scan(
  f_disgenet,
  what = "",
  sep = "\n"
)
G <- G[, colnames(G) %in% rs]

# other variables
V <- fread(f_other_imp)
N <- fread(f_other_ori)
V[is.na(N)] <- NA
pid <- V$V1
V$V1 <- NULL
V <- as.matrix(V)
rownames(V) <- pid
D <- fread(f_desc, header = TRUE)
postfix <- stringr::str_extract(colnames(V), "\\.V[0-9]+$|\\.BL$")
postfix[is.na(postfix)] <- ""
postfix <- gsub("\\.", "_", postfix)
nms <- gsub("\\.V[0-9]+$|\\.BL$", "", colnames(V))
ii <- match(nms, D$`Abbreviation USED`)
colnames(V) <- D$`Variable Description`[ii]
colnames(V) <- gsub(" ", "_", colnames(V))
colnames(V) <- gsub("-", "", colnames(V))
colnames(V) <- gsub("\\(|\\)", "", colnames(V))
colnames(V) <- gsub("\\/", "", colnames(V))
colnames(V) <- sprintf("%s%s", colnames(V), postfix)

C <- fread(f_patchar)
names(C)[names(C) == "V1"] <- "PID"
C$Medfree <- NULL

save(X, MNAR, G, C, V, file = f_out)
