library(matrixStats)
library(data.table)
library(abind)
library(reticulate)

parse_results <- function(dir_out, metric, decreasing = TRUE) {
  files <- list.files(dir_out, "_seed[0-9]+\\.RData", full.names = TRUE)
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

time_stamp <- "20190412210958"
metric <- "test_reconstruction_loss"
decreasing <- FALSE
dir_out <- file.path("..", "results", "PPMI", "vader", "hyperparameter_optimization", time_stamp)
perf <- parse_results(dir_out, metric, decreasing = decreasing)
fwrite(perf, file = file.path(dir_out, "grid_search.csv"))


############################################
# parse results from                       #
# variational models with optimal          #
# hyperparameters for each k: choose best k#
############################################

time_stamp <- "20190417134522"
metric <- "prediction_strength"
decreasing <- TRUE
dir_out <- file.path("..", "results", "PPMI", "vader", "hyperparameter_optimization", time_stamp)
p <- parse_results(dir_out, metric, decreasing = decreasing)
fwrite(p, file = file.path(dir_out, "number_of_clusters.csv"))

# x <- (p$prediction_strength - p$prediction_strength__95CI) - 
#   (p$prediction_strength_null - p$prediction_strength_null__95CI)
# y <- q$prediction_strength / q$prediction_strength_null
# plot(x, y, log = "y", col = "white")
# text(x, y, label = p$k)

# plot(p$test_latent_loss, p$adj_rand_index, col = "white")
# text(p$test_latent_loss, p$adj_rand_index, label = p$k)

# plot(q$prediction_strength, q$rand_index, col = "white")
# text(q$prediction_strength, q$rand_index, label = p$k)



pdf(file.path(dir_out, "number_of_clusters.pdf"))

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
dev.off()

############################################
# final model with optimal hyperparameters #
# and optimal k                            #
############################################

k <- 3
learning_rate <- 1E-3
batch_size <- 16
n_layer <- 2
n_hidden1 <- 128
n_hidden2 <- 32
# f_in <- file.path("..", "data", "PPMI", "PPMI.RData")
f_in <- file.path("PPMI_artificial_data.RData")
vader_path <- "../VaDER/"

load(f_in)
W <- 1 - MNAR
dimnames(X)[[3]] <- dimnames(W)[[3]] <- toupper(dimnames(X)[[3]])
X[W == 0] <- 0 # set missing values arbitrarily to 0 (can be any value)

# The final model, on all data
VADER <- import_from_path("vader", path = vader_path)$VADER
save_path <- file.path(
  "vader", "vader.ckpt"
)
vader <- VADER(
  X_train = X,
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
  weights = W # represents the missing values
)
vader$pre_fit(n_epoch = as.integer(100), verbose = TRUE)
vader$fit(n_epoch = as.integer(100), verbose = TRUE)
yhat <- vader$cluster(X, W)
print(table(yhat))

save(yhat, file = file.path(dir_out, "clusters.RData"))

pdf(file.path(dir_out, "cluster_means.pdf"), width = 7, height = 7)

par(mfrow = c(3, 3))
YLIM <- list()
for (i in 1:dim(X)[3]) {
  dt <- data.table(X[,,i])
  means <- as.matrix(dt[ , lapply(.SD, mean), by = yhat])
  means <- means[order(means[,"yhat"]),]
  y <- t(means[,-1,drop=FALSE])
  clusters <- means[,1,drop=FALSE]
  nc <- ncol(y)
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
  YLIM[[i]] <- range(y)
  matplot(x, y,
    type = "l",
    col = rainbow(nc),
    lty = rep(c(1, 2, 2), each = nc),
    lwd = rep(c(3, 1, 1), each = nc),
    xlab = "Month",
    ylim = YLIM[[i]],
    ylab = "z-score relative to baseline",
    main = sprintf("%s", dimnames(X)[[3]][i])
  )
  if (i == dim(X)[3]) {
    tab <- table(yhat)
    legend(
      "topleft",
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

X_train <- X
W_train <- W
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
  weights = W_sim
  # n_thread = 1L,
  # phi = rep(1/k, k)
)
system.time(vader_sim$pre_fit(n_epoch = as.integer(100), verbose = TRUE))
system.time(vader_sim$fit(n_epoch = as.integer(100), verbose = TRUE))

pdf(file.path(dir_out, "cluster_means_simulated.pdf"), width = 7, height = 7)
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
save(X, W, y, file = file.path(dir_out, "PPMI_artificial_data.RData"))
X <- do.call(cbind, lapply(1:dim(X)[3], function(i) {
  x <- X[,,i]
  colnames(x) <- sprintf("%s_%sm", dimnames(X)[[3]][i], dimnames(X)[[2]])
  x
}))
rownames(X) <- 1:nrow(X)
dt <- data.table(cbind(cluster = y, X))
fwrite(dt, file = file.path(dir_out, "PPMI_artificial_data.csv"))

