library(matrixStats)
library(data.table)
library(abind)
library(reticulate)

parse_results <- function(results_dir, metric, decreasing = TRUE) {
  files <- list.files(results_dir, "_seed[0-9]+\\.RData", full.names = TRUE)
  P <- lapply(files, function(file) {
    p <- get(load(file))
    p$rank <- order(order(p[[metric]], sample(nrow(p)), decreasing = decreasing, na.last = TRUE), decreasing = FALSE, na.last = TRUE)
    p$top_ranking <- p$rank == 1
    as.matrix(p)
  })
  
  # check if all are properly sorted
  param_cols <- 1:tail(grep("n_hidden", colnames(P[[1]])), 1)
  C <- lapply(P, function(p) {
    do.call(paste, lapply(param_cols, function(param_col) { p[,param_col] }))
  })
  b <- all(sapply(C, function(p) {
    all(sapply(C, function(q) {
      all(p == q)
    }))
  }))
  if (b) {
    cat("Everything is properly sorted\n")
  } else {
    cat("Something is wrong with the sorting!\n")
  }
  
  perf_cols <- (tail(grep("n_hidden", colnames(P[[1]])), 1) + 1):ncol(P[[1]])
  perf <- do.call("cbind", lapply(perf_cols, function(perf_col) {
    p <- do.call("cbind", lapply(1:length(P), function(i) {
      q <- P[[i]][,perf_col]
    }))
    p <- cbind(
      rowMeans(p, na.rm = TRUE),
      rowSdDiffs(p, na.rm = TRUE) / sqrt(rowSums(!is.na(p)) * 2) * 1.96 # 95% CI
    )
    colnames(p) <- sprintf("%s%s", colnames(P[[1]])[perf_col], c("", "__95CI"))
    p
  }))
  perf <- data.table(cbind(P[[1]][,param_cols], perf))
  
  setorderv(perf, cols = metric, na.last = TRUE, order = if (decreasing) -1 else 1)
  perf
}

###############################################
# parse results from                          #
# non-variational hyperparameter optimization #
###############################################

# First run ADNI_hyperparameter_optimization.r, which 
# writes its output to a directory for which you need
# to provide the name below in the variable dir_out. 
# Given the hyperparameter optimization results, the
# code below is runnable and can be uncommented

time_stamp <- "20190412210934"
metric <- "test_reconstruction_loss"
decreasing <- FALSE
results_dir <- file.path("..", "results", "ADNI", "vader", "hyperparameter_optimization", time_stamp)
perf <- parse_results(results_dir, metric, decreasing = decreasing)
fwrite(perf, file = file.path(results_dir, "grid_search.csv"))


############################################
# parse results from                       #
# variational models with optimal          #
# hyperparameters for each k: choose best k#
############################################

# First run ADNI_hyperparameter_optimization.r, which 
# writes its output to a directory for which you need
# to provide the name below in the variable dir_out. 
# Given the hyperparameter optimization results, the
# code below is runnable and can be uncommented

time_stamp <- "20190417124750"
metric <- "prediction_strength"
decreasing <- TRUE
results_dir <- file.path("..", "results", "ADNI", "vader", "hyperparameter_optimization", time_stamp)
p <- parse_results(results_dir, metric, decreasing = decreasing)
fwrite(p, file = file.path(results_dir, "number_of_clusters.csv"))


pdf(file.path(results_dir, "number_of_clusters.pdf"))

mu <- p[[metric]]
sigma <- p[[sprintf("%s__95CI", metric)]]
mu_null <- p[[sprintf("%s_null", metric)]]
sigma_null <- p[[sprintf("%s_null__95CI", metric)]]
ii <- order(p$k, decreasing = FALSE)
x <- matrix(rep(p$k, 6), ncol = 6)[ii,]
y <- cbind(
  mu - sigma, mu, mu + sigma,
  mu_null - sigma_null, mu_null, mu_null + sigma_null
)[ii,]
matplot(
  x = x,
  y = y,
  xlab = "k",
  ylab = metric,
  type = "l",
  lty = rep(c(2, 1, 2), 2),
  lwd = rep(c(1, 2, 1), 2),
  col = c(rep("blue", 3), rep("red", 3)),
  main = ""
)
legend(
  "topright",
  legend = c("model", "null", "95% CI"),
  col = c("blue", "red", "gray"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 1)
)

mu <- p$effective_k
sigma <- p$effective_k__95CI
ii <- order(p$k, decreasing = FALSE)
x <- matrix(rep(p$k, 3), ncol = 3)[ii,]
y <- cbind(mu - sigma, mu, mu + sigma)[ii,]
matplot(
  x = x,
  y = y,
  xlab = "k",
  ylab = "effective k",
  type = "l",
  lty = c(2, 1, 2),
  lwd = c(1, 2, 1),
  col = "red",
  main = ""
)
abline(0, 1, col = "blue", lty = 3)
legend(
  "topright",
  legend = c("model", "null", "95% CI"),
  col = c("blue", "red", "gray"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 1)
)

dev.off()

############################################
# train final model with optimal           #
# hyperparameters and optimal k            #
############################################

k <- 3
learning_rate <- 1E-3
batch_size <- 16
n_layer <- 2
n_hidden1 <- 64
n_hidden2 <- 32
code_dir <- "."
code_file <- "ADNI_hyperparameter_optimization.r"
f_in <- file.path("..", "data", "ADNI", "ADNI.RData")
# f_in <- file.path("ADNI_artificial_data.RData")
vader_path <- "../VaDER/"
vars <- c("ADAS13", "CDRSB", "MMSE", "FAQ")

s <- parse(file.path(code_dir, code_file))
# makes available the functions load_data(), ...
eval(s[grep("<- function\\(", s)]) 

L <- load_data(f_in, vars)
X_train <- L$X
W_train <- L$W
ptid <- L$ptid

# X <- X_train
# X[W_train == 0] <- NA
# X <- X[,8,]
# C <- cor(X, use = "pairwise.complete.obs")
# heatmap(abs(C))

VADER <- reticulate::import_from_path("vader", path = vader_path)$VADER

# The final model, on all data
save_path <- file.path("vader", "vader.ckpt")
vader <- VADER(
  X_train = X_train,
  save_path = save_path,
  n_hidden = lapply(c(n_hidden1, n_hidden2), as.integer),
  k = as.integer(k),
  learning_rate = learning_rate,
  batch_size = as.integer(batch_size),
  output_activation = NULL,
  alpha = 1.0,
  seed = 123L, # makes it reproducible
  recurrent = TRUE,
  # n_thread = 1L,
  W_train = W_train
  # phi = rep(1/k, k)
)
system.time(vader$pre_fit(n_epoch = as.integer(100), verbose = TRUE))
system.time(vader$fit(n_epoch = as.integer(100), verbose = TRUE))
lat <- vader$map_to_latent(X_train, W_train)
yhat <- vader$cluster(X_train, W_train)
print(table(yhat))

names(yhat) <- ptid
save(yhat, file = file.path(results_dir, "clustering.RData"))


pdf(file.path(results_dir, "cluster_means.pdf"), width = 7, height = 7)

# Only plot values that are non-missing (missing values were assigned an
# arbitrary numerical value of 0 for the purpose of model fitting.)
X_train[W_train == 0] <- NA

n <- ceiling(sqrt(dim(X_train)[3]))
par(mfrow = c(n, n))
YLIM <- list()
for (i in 1:dim(X_train)[3]) {
  dt <- data.table(X_train[,,i])
  means <- as.matrix(dt[ , lapply(.SD, mean, na.rm = TRUE), by = yhat])
  means <- means[order(means[,"yhat"]),]
  y <- t(means[,-1,drop=FALSE])
  clusters <- means[,1,drop=FALSE]
  nc <- ncol(y)
  # yci <- aggregate(as.matrix(dt), list(yhat), function(x) {
  #   colSdDiffs(x, na.rm = TRUE) / sqrt(colSums(!is.na(x)) * 2) * 1.96
  # })
  yci <- as.matrix(dt[ , lapply(.SD, function(x) {
    sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)) * 2) * 1.96
  }), by = yhat])
  yci <- yci[order(yci[,"yhat"]),]
  yci <- t(yci[,-1,drop=FALSE])
  isna <- is.na(yci) | is.na(y)
  if (any(isna)) {
    yci[isna] <- y[isna] <- NA
  }
  y <- cbind(y, y - yci, y + yci)
  x <- as.numeric(rownames(y))
  YLIM[[i]] <- range(y, na.rm = TRUE)
  matplot(x, y,
          type = "l",
          col = rainbow(nc),
          lty = rep(c(1, 2, 2), each = nc),
          lwd = rep(c(3, 1, 1), each = nc),
          xlab = "Month",
          ylim = YLIM[[i]],
          ylab = "z-score relative to baseline",
          main = sprintf("%s", dimnames(X_train)[[3]][i])
  )
  if (i == dim(X_train)[3]) {
    tab <- table(yhat)
    legend(
      "bottomleft",
      legend = c(sprintf("%s (n = %i)", clusters, tab[match(clusters, names(tab))]), "95% CI"),
      col = c(rainbow(nc), "black"),
      lty = c(rep(1, nc), 2),
      lwd = c(rep(3, nc), 1)
    )
  }
}
par(mfrow = c(1, 1))
dev.off()

#############################################################
# Use VaDER to simulate patient data from the optimal model #
#############################################################

Xy <- vader$generate(as.integer(dim(X_train)[1]))
y_sim <- Xy$clusters
X_sim <- Xy$samples
dimnames(X_sim) <- dimnames(X_train)
W_sim <- W_train # generate_missingness_pattern(X_sim, k = 0.8)
mode(W_sim) <- "integer"
dimnames(W_sim) <- dimnames(X_train)

save_path <- file.path("vader_sim", "vader.ckpt")
vader_sim <- VADER(
  X_train = X_sim,
  save_path = save_path,
  n_hidden = lapply(c(n_hidden1, n_hidden2), as.integer),
  k = as.integer(k),
  learning_rate = learning_rate,
  batch_size = as.integer(batch_size),
  output_activation = NULL,
  alpha = 1.0,
  seed = 123L, # makes it reproducible
  y_train = y_sim,
  recurrent = TRUE,
  W_train = W_sim
  # n_thread = 1L,
  # phi = rep(1/k, k)
)
system.time(vader_sim$pre_fit(n_epoch = as.integer(100), verbose = TRUE))
system.time(vader_sim$fit(n_epoch = as.integer(100), verbose = TRUE))

pdf(file.path(results_dir, "cluster_means_simulated.pdf"), width = 7, height = 7)
yhat <- y_sim
n <- ceiling(sqrt(dim(X_sim)[3]))
par(mfrow = c(n, n))
YLIM <- list()
for (i in 1:dim(X_sim)[3]) {
  dt <- data.table(X_sim[,,i])
  means <- as.matrix(dt[ , lapply(.SD, mean, na.rm = TRUE), by = yhat])
  means <- means[order(means[,"yhat"]),]
  y <- t(means[,-1,drop=FALSE])
  clusters <- means[,1,drop=FALSE]
  nc <- ncol(y)
  # yci <- aggregate(as.matrix(dt), list(yhat), function(x) {
  #   colSdDiffs(x, na.rm = TRUE) / sqrt(colSums(!is.na(x)) * 2) * 1.96
  # })
  yci <- as.matrix(dt[ , lapply(.SD, function(x) {
    sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)) * 2) * 1.96
  }), by = yhat])
  yci <- yci[order(yci[,"yhat"]),]
  yci <- t(yci[,-1,drop=FALSE])
  isna <- is.na(yci) | is.na(y)
  if (any(isna)) {
    yci[isna] <- y[isna] <- NA
  }
  y <- cbind(y, y - yci, y + yci)
  x <- as.numeric(rownames(y))
  YLIM[[i]] <- range(y, na.rm = TRUE)
  matplot(x, y,
          type = "l",
          col = rainbow(nc),
          lty = rep(c(1, 2, 2), each = nc),
          lwd = rep(c(3, 1, 1), each = nc),
          xlab = "Month",
          ylim = YLIM[[i]],
          ylab = "z-score relative to baseline",
          main = sprintf("%s", dimnames(X_sim)[[3]][i])
  )
  if (i == dim(X_sim)[3]) {
    tab <- table(yhat)
    legend(
      "bottomleft",
      legend = c(sprintf("%s (n = %i)", clusters, tab[match(clusters, names(tab))]), "95% CI"),
      col = c(rainbow(nc), "black"),
      lty = c(rep(1, nc), 2),
      lwd = c(rep(3, nc), 1)
    )
  }
}
par(mfrow = c(1, 1))
dev.off()

X <- X_sim
W <- W_sim
X[W == 0] <- NA
y <- y_sim
save(X, W, y, file = file.path(results_dir, "ADNI_artificial_data.RData"))
X <- do.call(cbind, lapply(1:dim(X)[3], function(i) {
  x <- X[,,i]
  colnames(x) <- sprintf("%s_%sm", dimnames(X)[[3]][i], dimnames(X)[[2]])
  x
}))
rownames(X) <- 1:nrow(X)
dt <- data.table(cbind(cluster = y, X))
fwrite(dt, file = file.path(results_dir, "ADNI_artificial_data.csv"))
