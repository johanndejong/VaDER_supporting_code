# some input file dependencies
dir_in <- file.path("..", "data", "ADNI")
dir_out <- file.path("..", "results", "ADNI", "analyze_clusters")
files <- file.path(c(
  file.path("..", "results", "ADNI", "vader", "hyperparameter_optimization", "20190417124750")
), "clustering.RData")
f_ADNI <- file.path(dir_in, "ADNI.RData")
f_adni_ori <- file.path(dir_in, "from_karki", "BeforeImpADNI.csv")

library(matrixStats)
library(data.table)
library(nnet)
library(mlogit)

get_weights <- function(y) {
  weights <- 1 / table(y)
  weights <- weights[match(y, as.numeric(names(weights)))]
  weights <- length(y) * weights / sum(weights)
  weights
}

my_scale <- function(data) {
  for (i in 1:ncol(data)) {
    if (is.numeric(data[[i]])) {
      (data[[i]] - min(data[[i]])) / diff(range(data[[i]]))
    }
  }
  data
}

dir.create(dir_out)

Y <- do.call("cbind", lapply(files, function(file) {
  get(load(file))
}))
Y <- cbind(Y, sample(0:1, nrow(Y), replace = TRUE))
colnames(Y) <- c("VADER", "RVAE")

load(f_ADNI)

ptid <- intersect(rownames(Y), rownames(G))
G <- G[match(ptid, rownames(G)),]
brain_subsub <- brain_subsub[match(ptid, rownames(brain_subsub)),]
colnames(brain_subsub) <- gsub("4th", "Fourth", gsub("3rd", "Third", gsub("[-&]", "_", colnames(brain_subsub))))
Y <- Y[match(ptid, rownames(Y)),]

dt <- fread(f_adni_ori)
dt <- dt[, !grepl("aux$", colnames(dt)), with = FALSE]
dt$TIME2AD.bl <- annot$TIME2AD.m00[match(annot$PTID, dt$PTID)]
ii <- match(ptid, dt$PTID)
patient_characteristic <- dt[ii, c(3:8, ncol(dt)), with = FALSE]
brain_main <- dt[ii, 9:36, with = FALSE]
marker <- dt[ii, c(37:39, 243, 354, 355), with = FALSE]
# G <- dt[ , 244:353, with = FALSE]
cognitive_assessment <- dt[ii, c(80:111), with = FALSE]
brain_sub <- dt[ii, 356:423, with = FALSE]
colnames(Y) <- sprintf("clustering____%s", colnames(Y))
colnames(G) <- sprintf("snp____%s", colnames(G))
colnames(patient_characteristic) <- sprintf("patient_characteristic____%s", colnames(patient_characteristic))
colnames(brain_main) <- sprintf("brain_level1____%s", colnames(brain_main))
colnames(marker) <- sprintf("marker____%s", colnames(marker))
colnames(brain_sub) <- sprintf("brain_level2____%s", colnames(brain_sub))
colnames(brain_subsub) <- sprintf("brain_level3____%s", colnames(brain_subsub))

DT <- data.table(cbind(
  # Y, G, patient_characteristic, brain_main, marker, brain_sub, brain_subsub
  Y, patient_characteristic, brain_main, marker, brain_sub, brain_subsub
))
for (i in 1:ncol(DT)) {
  if (is.numeric(DT[[i]])) {
    ii <- which(!is.finite(DT[[i]]))
    if (length(ii) > 0) {
      DT[[i]][ii] <- NA
    }
  }
}
groups <- c(
  rep("clustering", ncol(Y)),
  # rep("G", ncol(G)),
  rep("patient_characteristic", ncol(patient_characteristic)),
  rep("brain_level1", ncol(brain_main)),
  rep("marker", ncol(marker)),
  rep("brain_level2", ncol(brain_sub)),
  rep("brain_level3", ncol(brain_subsub))
)
keep <- unlist(lapply(1:length(DT), function(i) {
  x <- DT[[i]]
  if (is.character(x)) {
    length(unique(x, na.rm = TRUE)) > 1
  } else {
    var(x, na.rm = TRUE) > 0
  }
}))
DT <- DT[, keep, with = FALSE]
groups <- groups[keep]
# remove <- colnames(DT) %in% c("marker____APOE4.bl")
# DT <- DT[, !remove, with = FALSE]
# groups <- groups[!remove]
colnames(DT) <- gsub("^.*____", "", colnames(DT))
cnames <- colnames(DT)

confounders <- c("AGE.bl", "PTEDUCAT.bl")
clusterings <- c("RVAE", "VADER")
mask <- !colnames(DT) %in% confounders & !colnames(DT) %in% clusterings
CNAMES <- cnames[mask]
GROUPS <- groups[mask]

tab <- table(gsub("\\.bl|\\.m[0-9]{2}$", "", CNAMES))
isseries <- names(tab[tab > 1])
isseries <- grepl(paste(sprintf("^%s\\.", isseries), collapse = "|"), CNAMES)
isbaseline <- !grepl("\\.m[0-9]{2}$", CNAMES)

pdf(file.path(dir_out, "confounders.pdf"))
p <- numeric(0)
for (confounder in confounders) {
  if (length(unique(DT[[confounder]])) < 5 | is.character(DT[[confounder]])) {
    tab <- table(DT[[confounder]],DT$VADER)
    p[confounder] <- chisq.test(tab)$p.value
    tab <- t(t(tab) / colSums(tab))
    b <- barplot(
      tab,
      ylab = sprintf("%s fraction", confounder),
      legend = TRUE
    )
    text(
      mean(b),
      0.5,
      label = sprintf("p = %.2g", p[confounder]),
      col = "red"
    )
  } else {
    # kw <- coin::kruskal_test(DT[[confounder]] ~ factor(DT$VADER))
    a <- summary(aov(DT[[confounder]] ~ factor(DT$VADER)))
    p[confounder] <- a[[1]][["Pr(>F)"]][1]
    boxplot(
      DT[[confounder]] ~ DT$VADER,
      ylab = confounder
    )
    text(
      2.5,
      mean(DT[[confounder]]),
      label = sprintf("p = %.2g", p[confounder]),
      col = "red"
    )
  }
}
dev.off()

# only use those confounders that are significant
confounders <- confounders[p < 0.05]

#####################################
# analyze the baseline measurements #
#####################################
dt <- DT
cnames <- CNAMES[isbaseline]
groups <- GROUPS[isbaseline]
for (clustering in clusterings) {
  P <- data.table(
    variable = cnames,
    group = groups,
    pval = rep(1.0, length(cnames))
  )
  for (i in 1:length(cnames)) {
    cname <- cnames[i]
    cat(sprintf("%i of %i\n", i, length(cnames)))
    cols <- c(clustering, cname, confounders)
    data <- dt[, cols, with = FALSE]
    data <- data[rowSums(is.na(data)) == 0,]
    if (length(unique(data[[cname]])) > 1) {
      fit1 <- multinom(
        formula = formula(sprintf(
          "%s ~ %s", 
          clustering, paste(c(cname, confounders), collapse = " + ")
        )),
        data = my_scale(data),
        # weights = get_weights(data[[clustering]]),
        trace = FALSE
      )
      
      fit0 <- multinom(
        formula = formula(sprintf(
          "%s ~ %s", 
          clustering, if (length(confounders) > 0) paste(confounders, collapse = " + ") else "1"
        )),
        data = my_scale(data),
        # weights = get_weights(data[[clustering]]),
        trace = FALSE
      )
      
      lr <- as.data.frame(anova(fit1, fit0))
      lr_stat <- lr$`LR stat.`[2]
      p_value <- lr$`Pr(Chi)`[2]
      
      P$pval[which(P$variable == cname)] <- p_value
    }
  }

  # P$qval <- p.adjust(P$pval, "fdr")
  I <- split(1:nrow(P), P$group)
  P$qval <- rep(1.0, nrow(P))
  for (ii in I) {
    P$qval[ii] <- p.adjust(P$pval[ii], "fdr")
  }
  P <- P[order(P$qval, P$pval, decreasing = FALSE),]
  fwrite(P, file = file.path(dir_out, sprintf("%s.csv", clustering)))
}

####################
# Do some plotting #
####################

plotfunc <- function(confounder, pval = TRUE, label = "p", description = NULL) {
  dt <- dt[ !is.na(dt[[confounder]]),]
  main <- if (is.null(description)) confounder else description
  main <- paste(strwrap(gsub("_", " ", main), 40), collapse = "\n")
  if (nrow(dt) > 0) {
    if (!is.numeric(dt[[confounder]]) || length(unique(dt[[confounder]])) < 4) {
      tab <- table(dt[[confounder]], dt$VADER)
      b <- barplot(
        t(t(tab) / colSums(tab)),
        ylim = 0:1,
        ylab = "Fraction",
        xlab = "Cluster",
        main = main,
        legend = TRUE
      )
      if (is.logical(pval) && pval) {
        pval <- chisq.test(tab)$p.value
        frac <- head(table(dt[[confounder]]), -1) / nrow(dt)
        abline(
          h = frac,
          col = "red",
          lty = 2
        )
        text(
          mean(b),
          frac,
          labels = "Expected",
          col = "red",
          pos = 3
        )
      }
      if (is.numeric(pval)) {
        text(
          mean(b),
          0,
          labels = sprintf("%s = %.2g", label, pval),
          col = "red",
          pos = 3
        )
      }
    } else {
      boxplot(
        dt[[confounder]] ~ dt$VADER,
        ylab = "Value",
        xlab = "",
        outline = FALSE,
        main = main,
        ylim = range(dt[[confounder]])
      )
      stripchart(
        dt[[confounder]] ~ dt$VADER,
        add = TRUE,
        vertical = TRUE,
        axes = FALSE,
        col = "gray",
        pch = 20,
        cex = 0.5,
        method = "jitter",
        jitter = 0.2
      )
      if (is.logical(pval) && pval) {
        pval <- coin::pvalue(coin::kruskal_test(dt[[confounder]] ~ factor(dt$VADER)))
      }
      if (is.numeric(pval)) {
        text(
          length(unique(dt$VADER)) / 2 + 0.5,
          mean(range(dt[[confounder]])),
          labels = sprintf("%s = %.2g", label, pval),
          col = "red",
          pos = 3
        )
      }
    }
  }
}

P <- fread(file.path(dir_out, "VADER.csv"))

# pdf(file.path(dir_out, "variables.pdf"), width = 50, height = 40)
# par(mfrow = c(14 , 20), mar = c(2, 4, 4, 2) + .1)
# selected_vars <- CNAMES[isbaseline & !CNAMES %in% c(confounders, "RVAE", "VADER", "PD.progression.HDDC")]
# selected_vars <- selected_vars[order(P$qval[match(selected_vars, P$variable)])]
# for (selected_var in selected_vars) {
#   plotfunc(selected_var, P$qval[P$variable == selected_var], label = "q")
# }
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + .1)
# dev.off()

selected_vars <- CNAMES[isbaseline & !CNAMES %in% c(confounders, "RVAE", "VADER", "PD.progression.HDDC")]
selected_vars <- selected_vars[order(P$qval[match(selected_vars, P$variable)])]
g <- groups[match(selected_vars, cnames)]
G <- split(selected_vars, g)
for (i in 1:length(G)) {
  subdir_out <- file.path(dir_out, sprintf("variables_%s", names(G)[i]))
  dir.create(subdir_out)
  sel_vars <- G[[i]]
  n <- ceiling(sqrt(length(sel_vars)))
  m <- ceiling(sqrt(length(sel_vars)))

  # for (sel_var in sel_vars) {
  #   pdf(
  #     file.path(subdir_out, sprintf("%s.pdf", sel_var)),
  #     width = 3,
  #     height = 5
  #   )
  #   plotfunc(sel_var, P$qval[P$variable == sel_var], label = "q")
  #   dev.off()
  # }

  pdf(
    file.path(dir_out, sprintf("variables_%s.pdf", names(G)[i])), 
    width = n * 2.5, 
    height = m * 1.75
  )
  par(mfrow = c(n, m), mar = c(2, 4, 4, 2) + .1)
  for (sel_var in sel_vars) {
    plotfunc(
      sel_var, 
      P$qval[P$variable == sel_var], 
      label = "q"
    )
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + .1)
  dev.off()
}

######################
# analyze the series #
######################
dt <- DT
cnames <- gsub("\\.bl$|\\.m[0-9]{2}", "", CNAMES[isseries])
ii <- sapply(split(1:length(cnames), cnames), head, 1)
cnames <- cnames[ii]
groups <- GROUPS[isseries][ii]
for (clustering in clusterings) {
  P <- data.table(
    variable = cnames,
    group = groups,
    pval_marker_time_markertime = rep(1.0, length(cnames)),
    pval_marker = rep(1.0, length(cnames)),
    pval_time = rep(1.0, length(cnames)),
    pval_markertime = rep(1.0, length(cnames))
  )
  for (i in 1:length(cnames)) {
    cname <- cnames[i]
    test_vars <- colnames(dt)[grepl(sprintf("%s\\.bl|%s\\.m[0-9]{2}", cname, cname), colnames(dt))]
    cat(sprintf("%i of %i\n", i, length(cnames)))
    cols <- c(
      clustering, 
      test_vars, 
      confounders
    )
    data <- dt[, cols, with = FALSE]
    data <- data[rowSums(is.na(data)) == 0,]
    data <- melt(data, measure.vars = test_vars)
    data$time <- as.integer(
      gsub("^.*\\.m", "", gsub("\\.bl", ".m00", data$variable))
    )
    data$variable <- NULL
    colnames(data)[colnames(data) == "value"] <- "marker"

    fit1 <- multinom(
      formula = formula(sprintf(
        "%s ~ %s", 
        clustering, 
        paste(c("marker", "time", "time*marker", confounders), collapse = " + ")
      )),
      data = my_scale(data),
      trace = FALSE
    )
    fit0 <- multinom(
      formula = formula(sprintf(
        "%s ~ %s", 
        clustering,
        if (length(confounders) > 0) paste(confounders, collapse = " + ") else "1"
      )),
      data = my_scale(data),
      trace = FALSE
    )
    lr <- as.data.frame(anova(fit1, fit0))
    P$pval_marker_time_markertime[which(P$variable == cname)] <- lr$`Pr(Chi)`[2]
    
    fit1 <- multinom(
      formula = formula(sprintf(
        "%s ~ %s", 
        clustering, 
        paste(c("marker", confounders), collapse = " + ")
        
      )),
      data = my_scale(data),
      trace = FALSE
    )
    lr <- as.data.frame(anova(fit1, fit0))
    P$pval_marker[which(P$variable == cname)] <- lr$`Pr(Chi)`[2]
    
    fit1 <- multinom(
      formula = formula(sprintf(
        "%s ~ %s", 
        clustering, 
        paste(c("time", confounders), collapse = " + ")
      )),
      data = my_scale(data),
      trace = FALSE
    )
    lr <- as.data.frame(anova(fit1, fit0))
    P$pval_time[which(P$variable == cname)] <- lr$`Pr(Chi)`[2]
    
    fit1 <- multinom(
      formula = formula(sprintf(
        "%s ~ %s", 
        clustering, 
        paste(c("time*marker", confounders), collapse = " + ")
      )),
      data = my_scale(data),
      trace = FALSE
    )
    lr <- as.data.frame(anova(fit1, fit0))
    P$pval_markertime[which(P$variable == cname)] <- lr$`Pr(Chi)`[2]
    
  }
  
  P$qval_marker_time_markertime <- p.adjust(P$pval_marker_time_markertime, "fdr")

  # P$qval_marker_time_markertime <- p.adjust(P$qval_marker_time_markertime, "fdr")
  I <- split(1:nrow(P), P$group)
  P$qval_marker_time_markertime <- rep(1.0, nrow(P))
  for (ii in I) {
    P$qval_marker_time_markertime[ii] <- p.adjust(P$pval_marker_time_markertime[ii], "fdr")
  }
  # P <- P[order(P$qval, P$pval, decreasing = FALSE),]
  fwrite(P, file = file.path(dir_out, sprintf("%s_series.csv", clustering)))
}

plotfunc <- function(cname, pval = TRUE, label = "q", description = NULL) {
  vars <- colnames(DT)[grepl(sprintf("%s\\.bl$|%s\\.m[0-9]{2}", cname, cname), colnames(DT))]
  varname <- unique(gsub("\\.bl$|\\.m[0-9]{2}", "", vars))
  main <- if (is.null(description)) varname else description
  main <- paste(strwrap(gsub("_", " ", main), 30), collapse = "\n")
  X <- as.integer(gsub("^.*\\.m", "", gsub("\\.bl$", ".m00", vars)))
  dt <- DT[, c(vars, "VADER"), with = FALSE]
  Y <- dt[, lapply(.SD, mean, na.rm = TRUE), by = VADER]
  ii <- order(as.numeric(unique(Y$VADER)))
  YCI <- dt[, lapply(.SD, function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(x) * 2) * 1.96
  }), by = VADER]
  clusters <- Y$VADER
  Y$VADER <- YCI$VADER <- NULL
  Y <- Y[order(clusters),]
  YCI <- YCI[order(clusters),]
  clusters <- sort(clusters)
  Y <- cbind(t(Y), t(Y - YCI), t(Y + YCI))
  matplot(
    X,
    Y,
    type = "l",
    lty = rep(c(1, 2, 2), each = length(ii)),
    lwd = rep(c(3, 1, 1), each = length(ii)),
    col = rainbow(length(ii)),
    xlab = "Month",
    ylab = "Value",
    main = main
  )
  text(
    mean(range(X)),
    mean(range(Y)),
    labels = paste(sprintf("%s = %.2g", label, pval), collapse = "\n"),
    col = "black",
    pos = 1
  )
  legend(
    "topright",
    lty = c(rep(1, length(ii)), 2),
    lwd = c(rep(3, length(ii)), 1),
    col = c(rainbow(length(ii)), "black"),
    legend = c(clusters, "95% CI")
  )
}

P <- fread(file.path(dir_out, "VADER_series.csv"))

pdf(file.path(dir_out, "variables_series.pdf"), width = 10, height = 8)
par(mfrow = c(3, 3))
for (selected_var in cnames) {
  q <- P$qval_marker_time_markertime[P$variable == selected_var]
  if (q < 0.05) {
    pvals <- P[P$variable == selected_var, c(
      "qval_marker_time_markertime", "pval_marker", "pval_time", "pval_markertime"
    ), with = FALSE]
    cols <- c("q(marker + time + marker*time)", "p(marker)", "p(time)", "p(marker*time)")
  } else {
    pvals <- P[P$variable == selected_var, c("qval_marker_time_markertime"), with = FALSE]
    cols <- "q(marker + time + marker*time)"
  }
  pdf(file.path(subdir_out, sprintf("%s.pdf", selected_var)), width = 4, height = 7)
  plotfunc(selected_var, pval = pvals, label = cols)
  dev.off()
  plotfunc(selected_var, pval = pvals, label = cols)
}
par(mfrow = c(1, 1))
dev.off()

subdir_out <- file.path(dir_out, "variables_series")
dir.create(subdir_out)
for (selected_var in cnames) {
  q <- P$qval_marker_time_markertime[P$variable == selected_var]
  if (q < 0.05) {
    pvals <- P[P$variable == selected_var, c(
      "pval_marker", "pval_markertime"
    ), with = FALSE]
    cols <- c("p(marker)", "p(marker*time)")
  } else {
    pvals <- P[P$variable == selected_var, c("qval_marker_time_markertime"), with = FALSE]
    cols <- "q(marker + time + marker*time)"
  }
  pdf(file.path(subdir_out, sprintf("%s.pdf", selected_var)), width = 5, height = 5)
  plotfunc(selected_var, pval = pvals, label = cols)
  dev.off()
}
