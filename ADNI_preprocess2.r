library(missForest)
library(data.table)
library(dummy)
library(matrixStats)
library(doParallel)
library(gplots)
library(ADNIMERGE)
library(stringr)

dir_out <- file.path("..", "data", "ADNI", "plots")

# the locations of the files used in this script
dir_karki <- file.path(dir_out, "..", "from_karki")
dir_asif <- file.path(dir_out, "..", "from_asif")
dir_holger <- file.path(dir_out, "..", "from_holger")
f_ADNI_imp <- file.path(dir_karki, "ImputedADNI.csv")
f_ADNI_ori <- file.path(dir_karki, "BeforeImpADNI.csv")
f_genetics_adni1 <- file.path(dir_asif, "adadni1rec.raw")
f_genetics_adni2 <- file.path(dir_asif, "adadni2rec.raw")
f_disgenet <- file.path(dir_holger, "CuratedSNPs_DisGeNET_AD.txt")
f_subcortical <- file.path(dir_from_holger, "EMC_ADNI_FS60_Phenotypes_Subcortical_plus_Destrieux_20180219.csv")


ad <- data.table(adnimerge)
t2AD <- ad[, .SD$M[min(which(.SD$DX == "Dementia"))], by = PTID]
colnames(t2AD) <- c("PTID", "M")
t2AD$M[is.na(t2AD$M)] <- Inf
# dx <- dcast(ad, PTID ~ M, value.var = "DX")
# ab <- data.table(plasmaabeta)

dir.create(dir_out)

dt <- fread(f_ADNI_imp)
dt_raw <- fread(f_ADNI_ori)
dt_raw <- dt_raw[, match(colnames(dt), colnames(dt_raw)), with = FALSE]
dt$TIME2AD.bl <- dt_raw$TIME2AD.bl <- t2AD$M[match(dt$PTID,t2AD$PTID)]

##########################
# cap the RAVLT outliers #
##########################
cols <- grep("RAVLT.perc.forgetting", colnames(dt))
for (col in cols) {
  dt[[col]] <- pmax(-100, dt[[col]])
}

pdf(file.path(dir_out, "distributions_before_normalization.pdf"))
for (i in 1:ncol(dt)) {
  if (is.numeric(dt[[i]])) {
    hist(
      dt[[i]], 
      col = "red",
      xlab = colnames(dt)[i],
      main = sprintf("%s: dt_raw", colnames(dt)[i])
    )
    hist(
      dt_raw[[i]], 
      col = "red",
      xlab = colnames(dt)[i],
      main = sprintf("%s: imputed", colnames(dt)[i])
    )
  }
}
dev.off()

aux <- dt[, grepl("\\.aux$", colnames(dt)), with = FALSE]
dt <- dt[, !grepl("\\.aux$", colnames(dt)), with = FALSE]
dt_raw <- dt_raw[, !grepl("\\.aux$", colnames(dt_raw)), with = FALSE]
sort(unique(gsub("^.*\\.", "", colnames(dt)[grep("\\.m[0-9]{2}", colnames(dt))])))

colnames(dt) <- gsub("\\.bl", ".m00", colnames(dt))
colnames(dt_raw) <- gsub("\\.bl", ".m00", colnames(dt_raw))
selected <- grepl("\\.m[0-9]{2}", colnames(dt))

ptid <- dt$PTID
annot <- dt[, !selected, with = FALSE]
dt <- dt[, selected, with = FALSE]
annot_raw <- dt_raw[, !selected, with = FALSE]
dt_raw <- dt_raw[, selected, with = FALSE]
I <- split(1:ncol(dt), gsub("\\.m[0-9]{2}", "", colnames(dt)))
annot <- dt[, unlist(I[lengths(I) == 1]), with = FALSE]
dt <- dt[, unlist(I[lengths(I) > 1]), with = FALSE]
annot_raw <- dt_raw[, unlist(I[lengths(I) == 1]), with = FALSE]
dt_raw <- dt_raw[, unlist(I[lengths(I) > 1]), with = FALSE]
I <- split(1:ncol(dt), gsub("\\.m[0-9]{2}", "", colnames(dt)))

pdf(file.path(dir_out, "trajectories_before_normalization.pdf"))
for (i in 1:length(I)) {
  y <- t(as.matrix(dt_raw[, I[[i]], with = FALSE]))
  x <- as.integer(gsub("^.*\\.m", "", rownames(y)))
  x <- matrix(rep(x, ncol(y)), ncol = ncol(y))
  matplot(x, y,
    type = "l",
    ylab = "score",
    xlab = "Month",
    main = sprintf("%s: dt_raw", names(I)[i])
  )

  y <- t(as.matrix(dt[, I[[i]], with = FALSE]))
  x <- as.integer(gsub("^.*\\.m", "", rownames(y)))
  x <- matrix(rep(x, ncol(y)), ncol = ncol(y))
  matplot(x, y,
    type = "l",
    ylab = "score",
    xlab = "Month",
    main = sprintf("%s: imputed", names(I)[i])
  )

}
dev.off()

pdf(file.path(dir_out, "trajectories_after_normalization.pdf"))
for (i in 1:length(I)) {
  y <- t(as.matrix(dt[, I[[i]], with = FALSE]))
  y <- t((t(y) - y[1,]) / sd(y[1,]))
  x <- as.integer(gsub("^.*\\.m", "", rownames(y)))
  x <- matrix(rep(x, ncol(y)), ncol = ncol(y))
  matplot(x, y,
    type = "l",
    ylab = "score",
    xlab = "Month",
    main = names(I)[i]
  )
}
dev.off()

# Extract cognitive assessments:
# Cognitive assessments can be identified by having 8 time points
# annot <- cbind(annot, dt[, unlist(I[lengths(I) == 4]), with = FALSE])
# annot_raw <- cbind(annot_raw, dt_raw[, unlist(I[lengths(I) == 4]), with = FALSE])
dt <- dt[, unlist(I[lengths(I) %in% c(4, 8)]), with = FALSE]
dt_raw <- dt_raw[, unlist(I[lengths(I) %in% c(4, 8)]), with = FALSE]
I <- split(1:ncol(dt), gsub("\\.m[0-9]{2}", "", colnames(dt)))

dt_unnormalized <- dt

# normalize w.r.t. baseline
dt <- do.call("cbind", lapply(1:length(I), function(i) {
  d <- dt[, I[[i]], with = FALSE]
  dr <- dt_raw[, I[[i]], with = FALSE]
  # if ( grepl("RAVLT", names(I)[i]) ) {
  #   d <- (d - d[[1]])[, -1, drop = FALSE]
  # } else {
  #   d <- log2( d / d[[1]] )[, -1, drop = FALSE]
  # }
  
  # d <- (d - d[[1]])[, -1, drop = FALSE]

  d <- (d - d[[1]]) / sd(dr[[1]], na.rm = TRUE)

  d
}))
colnames(dt) <- colnames(dt_unnormalized)


# visit-wise imputation per group
visits <- gsub("^\\.", "", str_extract(colnames(dt), "\\.m[0-9]+$"))
cog <- c("ADAS11", "ADAS13", "CDRSB", "FAQ", "MMSE", "RAVLT")
brain <- c("Entorhinal", "Fusiform", "Hippocampus", "ICV", "MidTemp", "Ventricles", "WholeBrain")
groups <- rep(NA, ncol(dt))
ii <- which(gsub("\\..*$", "", colnames(dt)) %in% cog)
groups[ii] <- "cog"
ii <- which(gsub("\\..*$", "", colnames(dt)) %in% brain)
groups[ii] <- "brain"

# the auxiliary variables that reflect MNAR
mnar <- matrix(FALSE, nrow = nrow(dt), ncol = ncol(dt), dimnames = dimnames(dt))
I <- split(1:ncol(dt), list(groups, visits))
for (ii in I) {
  if (length(ii) > 0) {
    b <- rowSums(is.na(dt_raw[, ii, with = FALSE])) == length(ii)
    for (i in ii) {
      mnar[,i] <- b
    }
  } else {
    rep(FALSE, nrow(dt))
  }
}
mnar <- data.table(mnar)

G1 <- fread(f_genetics_adni1)
G2 <- fread(f_genetics_adni2)
colnames(G1) <- gsub("_[ACGT]+$", "", colnames(G1))
colnames(G2) <- gsub("_[ACGT]+$", "", colnames(G2))
rs <- scan(
  f_disgenet,
  what = "",
  sep = "\n"
)
pid <- c(G1$FID, G2$FID)
G1 <- G1[, colnames(G1) %in% rs, with = FALSE]
G2 <- G2[, colnames(G2) %in% rs, with = FALSE]
G <- rbind(G1, G2, fill = TRUE)
G <- as.matrix(G)
rownames(G) <- pid

remove <- grep("^rs[0-9]+$", colnames(annot))
if (length(remove) > 0) {
  snps <- colnames(annot)[remove]
  annot <- annot[,-remove,with=FALSE]
  annot_raw <- annot_raw[,-remove,with=FALSE]
} else {
  snps <- character(0)
}

annot$PTID <- annot_raw$PTID <- ptid

d <- fread(f_subcortical)
d <- d[d$ergoid.x %in% ptid,]
p <- d[1:nrow(annot),]
p[!is.na(p)] <- NA
colnames(p) <- colnames(d)
ii <- match(d$ergoid.x, ptid)
p[ii,] <- d
p$ergoid.x <- NULL
brain_subsub <- as.matrix(p)
rownames(brain_subsub) <- ptid

# annot$time2AD <- t2AD$M[match(, dt$)]


save(dt, dt_unnormalized, mnar, brain_subsub, dt_raw, aux, ptid, annot, annot_raw, G, file = file.path(dir_out, "..", "ADNI.RData"))



