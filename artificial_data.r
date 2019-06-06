vader_path <- file.path("..", "VaDER")
n_proc <- 120
n_repeat <- 100

n_samp <- 1e4
n_t <- 8
n_seq <- 4
n_clust <- 3

# variance of the VAR coefficients:
# the higher the variance, the more different the clusters
varscales <- seq(0.001, 0.25, length.out = 10)

k <- n_clust
learning_rate <- 1E-3
batch_size <- 64
n_layer <- 2
n_hidden1 <- 36
n_hidden2 <- 4
alpha <- 0.2
n_epoch <- 1e2
n_epoch_pretrain <- n_epoch

dir_out <- file.path("..", "results", "vader", "artificial_data")


library(MASS)
library(reticulate)
library(tsDyn)
library(abind)
library(caret)
library(parallel)
library(matrixStats)
library(latex2exp)
library(data.table)
library(fossil)
use_python(vader_path, required = TRUE)
VADER <- reticulate::import_from_path("vader", path = vader_path)$VADER


dir.create(dir_out, recursive = TRUE)


# lag_max <- n_t
# b <- .2
# B <- matrix(runif(n_seq ^ 2, -b, b), nrow = n_seq)
# if (lag_max > 1) for (lag in 2:lag_max) {
#   B <- cbind(B, matrix(runif(n_seq ^ 2, -b, b), nrow = n_seq))
# }
# B <- cbind(runif(n_seq, -1, 1), B) * 1
# A <- matrix(runif(n_seq ^ 2, -1, 1), ncol = n_seq)
# Sigma <- A %*% t(A)
# X <- array(dim = c(n_samp, n_t, n_seq))
# for (i in 1:n_samp) {
#   X[i,,] <- VAR.sim(B = B, n = n_t, lag = lag_max, include = "const", varcov = Sigma)
# }
# for (i in 1:n_seq) {
#   X[,,i] <- (X[,,i] - X[,1,i]) / sd(X[,1,i])
# }
# matplot(
#   t(X[,,1]),
#   type = "l",
#   lty = 1,
#   lwd = 1,
#   xlab = "time",
#   ylab = "value"
# )

is_stable <- function(B) {
  B1 <- tail(B, -1)
  B2 <- tail(B, 1)[[1]]
  n <- sum(sapply(B1, ncol))
  A <- cbind(
    rbind(
      do.call("cbind", B1),
      diag(n)
    ),
    rbind(
      B2,
      matrix(0, nrow = n, ncol = ncol(B2))
    )
  )
  all(Mod(eigen(A)$values) < 1)
}

random_matrix_with_eigenvalues <- function(lambda) {
  Q <- randortho(length(lambda))
  t(Q) %*% diag(lambda) %*% Q
}

generate_data <- function(n_samp, n_t, n_seq, n_clust, varscale = 1) {
  # theta <- rnorm(n_clust)
  # theta <- exp(theta) / sum(exp(theta))
  theta <- rep(1 / n_clust, n_clust)
  
  y <- sample(n_clust, n_samp, prob = theta, replace = TRUE)
  X <- array(dim = c(n_samp, n_t, n_seq))
  for (k in 1:n_clust) {
    ii <- which(y == k)
    if (length(ii) > 0) {
      B <- lapply(1:n_t, function(i) {
        matrix(runif(n_seq ^ 2, -1, 1) * varscale, nrow = n_seq)
      })
      while(!is_stable(B)) {
        B <- lapply(1:n_t, function(i) {
          matrix(runif(n_seq ^ 2, -1, 1) * varscale, nrow = n_seq)
        })
      }
      B <- do.call("cbind", B)
      B <- cbind(runif(n_seq, -1, 1), B) # const
      A <- matrix(runif(n_seq ^ 2, -.1, .1), ncol = n_seq)
      Sigma <- A %*% t(A)
      for (i in ii) {
        X[i,,] <- VAR.sim(B = B, n = n_t, lag = n_t, include = "const", varcov = Sigma)
      }
      for (i in 1:n_seq) {
        X[ii,,i] <- (X[ii,,i] - mean(X[ii,,i])) / sd(X[ii,,i])
      }
      # X[ii,,] <- replicate(n_seq, {
      #   # one separate GP for each sequence
      #   mu <- runif(n_t, -1, 1)
      #   A <- matrix(runif(n_t ^ 2, -1, 1), ncol = n_t)
      #   Sigma <- A %*% t(A)
      #   # normalize by the max-norm
      #   Sigma <- varscale * Sigma / max(abs(Sigma))
      #   # Sigma <- cov(matrix(runif(n_t ^ 2, -1, 1) * varscale, ncol = n_t))
      #   mvrnorm(sum(y == i), mu = mu, Sigma = Sigma)
      # })
    }
  }
  # # z-normalize
  # for (i in 1:n_seq) {
  #   X[,,i] <- (X[,,i] - X[,1,i]) / sd(X[,1,i])
  # }
  y <- y - 1
  mode(y) <- "integer"
  list(X = X, y = y)
}

##################################
#
# Complete data (VaDER + hierarchical)
#
##################################

set.seed(12345)

D <- lapply(varscales, function(varscale) {
  generate_data(n_samp, n_t, n_seq, n_clust, varscale = varscale)
})
pdf(file.path(dir_out, "example_data.pdf"), width = 4 * 2.5, height = n_seq * 2.5)
par(mfrow = c(n_seq, 4))
for (i in 1:n_seq) {
  for (j in seq(1, length(varscales), length.out = 4)) {
    varscale <- varscales[j]
    ii <- sample(n_samp, 5e2)
    X <- D[[j]]$X[ii,,]
    y <- D[[j]]$y[ii]
    matplot(
      t(X[,,i]),
      type = "l",
      lty = 1,
      lwd = 1,
      col = rainbow(n_clust)[y + 1],
      xlab = "time point",
      ylab = sprintf("variable %i", i),
      main = TeX(sprintf("$\\lambda = %.2g$", varscale))
    )
  }
}
par(mfrow = c(1, 1))
dev.off()
# I <- split(1:n_samp, y)
# for (i in 1:length(I)) {
#   ii <- I[[i]]
#   C <- unlist(lapply(1:n_seq, function(i) {
#     lapply(1:n_seq, function(j) {
#       d <- cor(X[ii,,i], X[ii,,j])
#       d[upper.tri(d)]
#     })
#   }), recursive = FALSE)
#   names(C) <- apply(expand.grid(1:n_seq, 1:n_seq), 1, paste, collapse = ",")
#   boxplot(
#     C, 
#     col = "red", 
#     las = 2, 
#     main = sprintf("Correlations between features\nwithin cluster %i", i),
#     ylab = "Pearson correlation coefficient"
#   )
# }
# C <- unlist(lapply(1:n_seq, function(i) {
#   lapply(1:n_seq, function(j) {
#     d <- cor(X[,,i], X[,,j])
#     d[upper.tri(d)]
#   })
# }), recursive = FALSE)
# names(C) <- apply(expand.grid(1:n_seq, 1:n_seq), 1, paste, collapse = ",")
# boxplot(
#   C, 
#   col = "red", 
#   las = 2, 
#   main = "Correlations between features\nacross all clusters",
#   ylab = "Pearson correlation coefficient"
# )

# tx: vector of length ncol(x) - the time points
STS <- function(x, tx, eps = 1e-8) {
  d <- matrix(nrow = nrow(x), ncol = nrow(x))
  rownames(d) <- colnames(d) <- rownames(x)
  dx <- rowDiffs(x) + eps
  dtx <- diff(tx) + eps
  for (i in 1:nrow(dx)) {
    d[i,] <- sqrt(colSums( (dx[i,] / dtx - t(dx) / dtx )^2 ))
  }
  d
}

hierarchical <- function(X, y, method = "euc") {
  if (method == "euc") {
    dim(X) <- c(dim(X)[1], dim(X)[2] * dim(X)[3])
    D <- dist(X)
  } else if (method == "cor") {
    dim(X) <- c(dim(X)[1], dim(X)[2] * dim(X)[3])
    D <- as.dist(1 - cor(t(X)))
  } else if (method == "sts") {
    D <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[1])
    for (i in 1:dim(X)[3]) {
      D <- D + STS(X[,,1], 1:dim(X)[2])
    }
    D <- as.dist(D)
  }
  hc <- hclust(D)
  cutree(hc, k = n_clust)
}

vader <- function(X, y) {

  save_path <- file.path(
    "vader", sprintf("%s.ckpt", gsub("\\.", "", sprintf("%.21f", as.numeric(Sys.time()))))
  )
  vader <- VADER(
    X_train = X,
    y_train = y,
    save_path = save_path,
    n_hidden = lapply(c(n_hidden1, n_hidden2), as.integer),
    k = as.integer(k),
    learning_rate = learning_rate,
    batch_size = as.integer(batch_size),
    output_activation = NULL,
    # alpha = alpha,
    # seed = 123L,
    recurrent = TRUE,
    n_thread = 1L
  )
  vader$pre_fit(n_epoch = as.integer(n_epoch_pretrain), verbose = FALSE)
  vader$fit(n_epoch = as.integer(n_epoch), verbose = FALSE)
  vader$cluster(X)
}


perf <- do.call("rbind", mclapply(rep(varscales, each = n_repeat), function(varscale) {
  # print(varscale)
  data <- generate_data(n_samp, n_t, n_seq, n_clust, varscale = varscale)
  X <- data$X
  y <- data$y
  mode(y) <- "integer"
  res <- c(
    varscale = varscale,
    vader = adj.rand.index(vader(X, y), y),
    hierarchical_euclidean = adj.rand.index(hierarchical(X, y, "euc"), y),
    hierarchical_correlation = adj.rand.index(hierarchical(X, y, "cor"), y),
    hierarchical_sts = adj.rand.index(hierarchical(X, y, "sts"), y)
  )
  print(res)
  res
}, mc.cores = n_proc))
save(perf, file = file.path(dir_out, "vader_vs_hierarchical.RData"))


# load(file.path(dir_out, "vader_vs_hierarchical_n1E4_ep1E2.RData"))

mu <- aggregate(perf, by = list(perf[,"varscale"]), mean, na.rm = TRUE)
x <- mu$varscale
mu <- as.matrix(mu[,-(1:2)])
ci <- as.matrix(aggregate(perf, by = list(perf[,"varscale"]), sd, na.rm = TRUE)[,-c(1:2)]) / sqrt(n_repeat) * 1.96



pdf(file.path(dir_out, "vader_vs_hierarchical.pdf"), width = 4.5, height = 4.5)
par(mar = c(5, 4, 1, 1) + .1)
matplot(
  x = x,
  y = cbind(mu - ci, mu, mu + ci),
  xlab = expression(lambda),
  ylab = "adjusted Rand index",
  type = "l",
  lty = rep(c(2, 1, 2), each = ncol(mu)),
  lwd = rep(c(1, 3, 1), each = ncol(mu)),
  col = rep(rainbow(ncol(mu)), 3)
)
legend(
  "bottomright",
  lty = c(rep(1, ncol(mu)), 2),
  lwd = c(rep(3, ncol(mu)), 1),
  col = c(rainbow(ncol(mu)), "black"),
  cex = 0.8,
  legend = c("VaDER", "Euclidean", "Correlation", "STS", "95% CI")
)
par(mar = c(5, 4, 4, 2) + .1)
dev.off()


##################################
#
# MAR (only VaDER)
#
##################################

set.seed(12345)

thetas <- seq(0, 0.8, 0.2)
ks <- thetas * 2
varscales <- c(0.001, 0.025, 0.05, 0.25)

# theta: "MAR" - parameter to Bernoulli distribution
# k: "MNAR" - logistic function growth rate
generate_missingness_pattern <- function(X, theta = NULL, k = NULL) {
  x0 <- (1 + dim(X)[2]) / 2
  fun <- function(x, k) {
    k / (1 + exp(x0 - x))
  }
  W <- rep(1, prod(dim(X)))
  dim(W) <- dim(X)
  if (!is.null(theta)) { # "MAR"
    W[runif(length(W)) < theta] <- 0
  }
  if (!is.null(k)) { # "MNAR"
    for (i in 1:dim(X)[2]) {
      W[,i,][which(runif(length(W[,i,])) < fun(i, k))] <- 0
    }
  }
  mode(W) <- "integer"
  W
}

# 460 seconds for one run
vader <- function(X, y, W) {
  save_path <- file.path(
    "vader", sprintf("%s.ckpt", gsub("\\.", "", sprintf("%.21f", as.numeric(Sys.time()))))
  )
  mode(W) <- "integer"
  vader <- VADER(
    X_train = X,
    y_train = y,
    weights = W,
    save_path = save_path,
    n_hidden = lapply(c(n_hidden1, n_hidden2), as.integer),
    k = as.integer(k),
    learning_rate = learning_rate,
    batch_size = as.integer(batch_size),
    output_activation = NULL,
    # alpha = alpha,
    seed = 123L,
    recurrent = TRUE,
    n_thread = 1L
  )
  vader$pre_fit(n_epoch = as.integer(n_epoch_pretrain), verbose = FALSE)
  vader$fit(n_epoch = as.integer(n_epoch), verbose = FALSE)
  vader$cluster(X)
}

params <- expand.grid(1:n_repeat, varscales, thetas)
names(params) <- c("repeat", "lambda", "theta")
params$ari <- do.call("rbind", mclapply(1:nrow(params), function(i) {
  cat(sprintf("%i of %i\n", i, nrow(params)))
  varscale <- params$lambda[i]
  theta <- params$theta[i]
  data <- generate_data(n_samp, n_t, n_seq, n_clust, varscale = varscale)
  X <- data$X
  y <- data$y
  W <- generate_missingness_pattern(X, theta = theta)
  mode(y) <- "integer"
  ari <- adj.rand.index(vader(X, y, W), y)
}, mc.cores = n_proc))
MAR <- params

params <- expand.grid(1:n_repeat, varscales, ks)
names(params) <- c("repeat", "lambda", "k")
params$ari <- do.call("rbind", mclapply(1:nrow(params), function(i) {
  cat(sprintf("%i of %i\n", i, nrow(params)))
  varscale <- params$lambda[i]
  k <- params$k[i]
  data <- generate_data(n_samp, n_t, n_seq, n_clust, varscale = varscale)
  X <- data$X
  y <- data$y
  W <- generate_missingness_pattern(X, k = k)
  mode(y) <- "integer"
  ari <- adj.rand.index(vader(X, y, W), y)
}, mc.cores = n_proc))
MNAR <- params

save(MAR, MNAR, file = file.path(dir_out, "missingness.RData"))

load(file.path(dir_out, "missingness.RData"))
L <- list(MAR = MAR, MNAR = MNAR)
for (i in 1:length(L)) {
  dt <- data.table(L[[i]])
  if (names(L)[i] == "MAR") {
    mu <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(theta, lambda)]
    sigma <- dt[, lapply(.SD, sd, na.rm = TRUE), by = .(theta, lambda)]
    n <- dt[, lapply(.SD, function(col) {
      sum(!is.nan(col) & !is.na(col))
    }), by = .(theta, lambda)]
    mu$`repeat` <- sigma$`repeat` <- n$`repeat` <- NULL
    rnms <- unique(mu$lambda)
    cnms <- unique(mu$theta)
    perc <- cnms
  } else {
    mu <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(k, lambda)]
    sigma <- dt[, lapply(.SD, sd, na.rm = TRUE), by = .(k, lambda)]
    n <- dt[, lapply(.SD, function(col) {
      sum(!is.nan(col) & !is.na(col))
    }), by = .(k, lambda)]
    mu$`repeat` <- sigma$`repeat` <- n$`repeat` <- NULL
    rnms <- unique(mu$lambda)
    cnms <- unique(mu$k)
    perc <- cnms / 2
  }
  mu <- t(matrix(mu$ari, nrow = length(rnms)))
  sigma <- t(matrix(sigma$ari, nrow = length(rnms)))
  n <- t(matrix(n$ari, nrow = length(rnms)))
  ci <- sigma / sqrt(n * 2) * 1.96
  pdf(file.path(dir_out, sprintf("%s.pdf", names(L)[i])), width = 4.5, height = 4.5)
  par(mar = c(5, 4, 1, 1) + .1)
  matplot(
    perc,
    mu,
    type = "l",
    lwd = 3,
    lty = 1,
    xlab = TeX("$\\theta$"),
    ylab = "Adjusted Rand index",
    ylim = range(as.numeric(mu), as.numeric(mu+ci), as.numeric(mu - ci)),
    col = rainbow(ncol(mu))
  )
  matlines(
    perc,
    cbind(mu - ci, mu + ci),
    lwd = 1,
    lty = 2,
    col = rep(rainbow(ncol(mu)), 2)
  )
  legend(
    "bottomleft",
    lwd = c(rep(3, ncol(mu)), 1),
    lty = c(rep(1, ncol(mu)), 2),
    col = c(rainbow(ncol(mu)), "black"),
    legend = c(TeX(sprintf("$\\lambda = %.2g$", rnms)), "95% CI")
  )
  par(mar = c(5, 4, 4, 2) + .1)
  dev.off()
} 
