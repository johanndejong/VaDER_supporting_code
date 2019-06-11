args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  N_PROC <- 2
  SEED <- 12345
  TIME_STAMP <- system("date +'%Y%m%d%H%M%S'", intern = TRUE)
} else {
  N_PROC <- as.integer(args[1])
  SEED <- as.integer(args[2])
  TIME_STAMP <- as.character(args[3])
}

USE_PYTHON <- NULL # "/nmitapps/lib/python/anaconda/3.6/bin/python3.6"
VADER_PATH <- file.path("..", "VaDER")
PRINT_OUT <- "ADNI_hyperparameter_optimization.out"
DIR_OUT <- file.path("..", "results", "ADNI", "vader", "hyperparameter_optimization", TIME_STAMP)
F_OUT <- file.path(DIR_OUT, sprintf("grid_search_seed%i.RData", SEED))
DIR_IN <- file.path("..", "data", "ADNI")
# F_IN <- file.path(DIR_IN, "ADNI.RData")
F_IN <- file.path("ADNI_artificial_data.RData")
VARS <- c("ADAS13", "CDRSB", "MMSE", "FAQ")
N_SAMPLE <- 1e4 # Inf; 100 takes 90 minutes using 3x100 cores
N_PERM <- 1e3
PARAM_SEED <- 12345
COMPUTE_PREDICTION_STRENGTH <- TRUE
if (COMPUTE_PREDICTION_STRENGTH) {
  N_FOLD <- 2
  # test one model for different k
  PARAM_GRID <- list( 
    k = 2:15,
    n_layer = 2,
    alpha = 1.0,
    n_node = c(32, 64),
    learning_rate = 1E-3,
    batch_size = 16
  )
} else {
  N_FOLD <- 10
  # hyperparameter optimization using non-variational AEs
  PARAM_GRID <- list(
    k = 1, # dummy if alpha = 0.0
    n_layer = 1:2,
    alpha = 0.0, # hyperparameter optimization using non-variational AEs
    n_node = c(1, 2, 4, 8, 16, 32, 64),
    learning_rate = 10^(-4:-1),
    batch_size = 2^(4:7)
  )
}

library(caret)
library(data.table)
library(abind)
library(parallel)
library(caret)
library(gplots)
library(matrixStats)
library(fossil)
library(reticulate)
# call directly after loading reticulate
if (!is.null(USE_PYTHON)) {
  use_python(USE_PYTHON, required = TRUE)
}
file.remove(PRINT_OUT)
dir.create(DIR_OUT, recursive = TRUE)

load_data <- function(f_in, vars) {
	s <- load(f_in)
	if ("y" %in% s) { # cluster label is present: artificial data
	  # set missing values arbitrarily to 0 (can be any value)
	  X[W == 0] <- 0 
	  list(X = X, W = W, ptid = 1:dim(X)[1], y = y)
	} else {
	  dt <- dt[, grepl(paste(vars, collapse = "|"), colnames(dt)), with = FALSE]
	  # dt <- dt[, grepl("m00$|m06$|m12$|m24$", colnames(dt)), with = FALSE]
	  dt_raw <- !is.na(dt_raw[, match(colnames(dt), colnames(dt_raw)), with = FALSE])
	  mnar <- 1 - mnar[, match(colnames(dt), colnames(mnar)), with = FALSE]
	  I <- split(1:ncol(dt), gsub("\\.m[0-9]{2}", "", colnames(dt)))
	  colnames(dt) <- as.character(as.integer(gsub("^.*\\.m", "", colnames(dt))))
	  X <- lapply(I, function(ii) { as.matrix(dt[, ii, with = FALSE]) })
	  names(X) <- names(I)
	  X_train <- do.call(abind, c(X, along = 3))
	  # X_train <- X_train
	  W <- mnar # mnar, dt_raw
	  I <- split(1:ncol(W), gsub("\\.m[0-9]{2}", "", colnames(W)))
	  colnames(W) <- as.character(as.integer(gsub("^.*\\.m", "", colnames(W))))
	  W <- lapply(I, function(ii) {
	    as.matrix(W[, ii, with = FALSE])
	  })
	  names(W) <- names(I)
	  W_train <- do.call(abind, c(W, along = 3))
	  mode(W_train) <- "numeric"
	  # set missing values arbitrarily to 0 (can be any value)
	  X_train[W_train == 0] <- 0 
	  list(X = X_train, W = W_train, ptid = annot$PTID)
	}
}

cross_validate <- function(params, data, weights, n_fold, n_perm, seed = NULL) {
  # cat("begin cross_validate\n", file = "vader_hyperparameter_optimization.out", append = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  adj_rand_index <- mclust::adjustedRandIndex
  
  # unadjusted rand index
  rand_index <- function(
    p, # clusters predicted on the test data using the model that was trained on the training data
    q # clusters inferred directly from the test data
  ) {
    n <- length(p)
    
    # y: cluster assignment vector
    # returns: length(y) * length(y) binary matrix indicating 
    # which points fall into the same cluster. 
    f <- function(y) {
      m <- matrix(rep(y, length(y)), ncol = length(y))
      m == t(m)
    }
    mp <- f(p)
    mq <- f(q)
    a <- (sum(mp & mq) - n) / 2
    b <- sum(!(mp | mq)) / 2
    (a + b) / choose(n, 2)
  }
  prediction_strength <- function(
    p, # clusters predicted on the test data using the model that was trained on the training data
    q # clusters inferred directly from the test data
  ) {
    n <- length(p)
    
    # y: cluster assignment vector
    # returns: length(y) * length(y) binary matrix indicating 
    # which points fall into the same cluster. 
    f <- function(y) {
      m <- matrix(rep(y, length(y)), ncol = length(y))
      m == t(m)
    }
    mp <- f(p)
    mq <- f(q)
    # mpq <- !mq | mp
    mpq <- mq & mp
    min(tapply(1:n, q, function(ii) {
      n_ii <- length(ii)
      (sum(mpq[,ii]) - n_ii) / n_ii / (n - 1)
    }))
  }
  
  VADER <- reticulate::import_from_path("vader", path = VADER_PATH)$VADER
  
  # train the model
  folds <- caret::createFolds(1:nrow(data), n_fold)
  perf <- do.call("rbind", lapply(1:length(folds), function(i) {
    fold <- folds[[i]]
    save_dir <- file.path(
      "temp", 
      paste(sample(letters, 64, replace = TRUE), collapse = ""),
      gsub("-", "minus", paste(unlist(sapply(params, as.character)), collapse = "_"))
    )
    dir.create(save_dir, recursive = TRUE)
    # save path for VADER
    save_path <- file.path(save_dir, "vader.ckpt")
    
    n_hidden <- lapply(params[grep("n_hidden", names(params))], as.integer)
    names(n_hidden) <- NULL
    
    vader <- VADER(
      X_train = data[-fold,,, drop = FALSE],
      save_path = save_path,
      n_hidden = n_hidden,
      k = as.integer(params$k),
      learning_rate = params$learning_rate,
      batch_size = as.integer(params$batch_size), 
      alpha = params$alpha,
      output_activation = NULL,
      seed = if (is.null(seed)) NULL else as.integer(seed),
      recurrent = TRUE,
      W_train = weights[-fold,,, drop = FALSE],
      n_thread = 1L
    )
    if (COMPUTE_PREDICTION_STRENGTH) {
      vader$pre_fit(n_epoch = as.integer(100), verbose = FALSE)
    }
    test_loss <- vader$get_loss(data[fold,,, drop = FALSE], weights[fold,,, drop = FALSE])
    loss <- c(
      vader$reconstruction_loss, 
      vader$latent_loss, 
      test_loss$reconstruction_loss, 
      test_loss$latent_loss
    )
    names(loss) <- c(
      "train_reconstruction_loss", 
      "train_latent_loss", 
      "test_reconstruction_loss", 
      "test_latent_loss"
    )
    if (COMPUTE_PREDICTION_STRENGTH) {
      # how many mixture components are effectively used?
      effective_k <- length(table(vader$cluster(data[-fold,,, drop = FALSE])))
      y_pred <- vader$cluster(
        data[fold,,, drop = FALSE],
        weights[fold,,, drop = FALSE]
      )
      
      vader <- VADER(
        X_train = data[fold,,, drop = FALSE],
        save_path = save_path,
        n_hidden = n_hidden,
        k = as.integer(params$k),
        learning_rate = params$learning_rate,
        batch_size = as.integer(params$batch_size), 
        alpha = params$alpha,
        output_activation = NULL,
        seed = if (is.null(seed)) NULL else as.integer(seed),
        recurrent = TRUE,
        W_train = weights[fold,,, drop = FALSE],
        n_thread = 1L
      )
      vader$pre_fit(n_epoch = as.integer(100), verbose = FALSE)
      vader$fit(n_epoch = as.integer(100), verbose = FALSE)
      y_true <- vader$cluster(
        data[fold,,, drop = FALSE], 
        weights[fold,,, drop = FALSE]
      )
      
      arindex <- adj_rand_index(y_pred, y_true)
      rindex <- rand_index(y_pred, y_true)
      pstrength <- prediction_strength(y_pred, y_true)
      null <- t(replicate(n_perm, {
        sample_y_pred <- sample(y_pred)
        c(
          rindex = rand_index(sample_y_pred, y_true),
          arindex = adj_rand_index(sample_y_pred, y_true),
          pstrength = prediction_strength(sample_y_pred, y_true)
        )
      }))
      res <- c(
        loss, 
        effective_k = effective_k, 
        rand_index = rindex, 
        rand_index_null = mean(null[,"rindex"]),
        adj_rand_index = arindex, 
        adj_rand_index_null = mean(null[,"arindex"]),
        prediction_strength = pstrength, 
        prediction_strength_null = mean(null[,"pstrength"])
      )
    } else {
      res <- loss
    }
    # delete model-related files
    unlink(save_dir, recursive=TRUE)
    
    res
  }))
  colMeans(perf)
}

explore_grid <- function(
  data,
  weights = NULL,
  param_grid,
  n_sample, # how many random samples to take from the grid?
  n_fold,
  n_proc,
  n_perm,
  seed = NULL,
  param_seed = NULL
) {
  
  nms <- names(param_grid)[names(param_grid) != "n_node"]
  paramspace <- data.table(expand.grid(param_grid[nms]))
  layer_configs <- lapply(param_grid$n_layer, function(n_layer) {
    p <- data.table(do.call(expand.grid, lapply(1:n_layer, function(i) {
      param_grid$n_node
    })))
    names(p) <- sprintf("n_hidden%i", 1:n_layer)
    p
  })
  names(layer_configs) <- param_grid$n_layer

  paramspace <- unlist(lapply(1:nrow(paramspace), function(i) {
    p <- as.matrix(cbind(
      paramspace[i,],
      layer_configs[[as.character(paramspace$n_layer[i])]]
    ))
    nms <- colnames(p)
    p <- split(p, row(p))
    p <- lapply(p, function(pi) {
      names(pi) <- nms
      pi
    })
  }), recursive = FALSE)
  keep <- unlist(lapply(paramspace, function(p) {
    b <- TRUE
    if ("n_hidden2" %in% names(p)) {
      b <- b && p["n_hidden2"] < p["n_hidden1"]
    }
    b
  }))
  paramspace <- paramspace[keep]
  if (n_sample < length(paramspace)) {
    if (!is.null(param_seed)) {
      set.seed(param_seed)
    }
    paramspace <- paramspace[sample(length(paramspace), n_sample)]
  }
  for (i in 1:length(paramspace)) {
    paramspace[[i]] <- as.list(paramspace[[i]])
  }

  cl <- makePSOCKcluster(min(nrow(paramspace), n_proc))
  is_first_iteration <- TRUE
  clusterExport(cl, envir = environment(), varlist = c(
    "data", "paramspace", "cross_validate", "n_fold", 
    "is_first_iteration", "n_perm", "seed", "weights", "n_proc", "PRINT_OUT",
    "COMPUTE_PREDICTION_STRENGTH", "VADER_PATH", "USE_PYTHON"
  ))
  
  clusterEvalQ(cl, {
    library(reticulate)
    if (!is.null(USE_PYTHON)) {
      use_python(USE_PYTHON, required = TRUE)
    }
  })
  cv_func <- function(i) {
    if (is_first_iteration) {
      Sys.sleep(10 * (i %% n_proc)) # to avoid high CPU loads simultaneously from all processes
      is_first_iteration <- FALSE
    }
    params <- paramspace[[i]]
    # # make sure each worker is seeded differently, but deterministically depends on the master thread
    # if (!is.null(seed)) { 
    #   seed <- sample(.Machine$integer.max, i)[i]
    # }
    # cat("begin call cross_validate\n", file = "vader_hyperparameter_optimization.out", append = TRUE)
    res <- cross_validate(params, data, weights, n_fold, n_perm, seed = seed)
    cat(
      sprintf(
        "%i of %i: \teffective_k=%.2g \t[%s]\n",
        i, length(paramspace), res["effective_k"],
        paste(sprintf("%s=%.4f", names(params), as.numeric(params)), collapse = "; ")
      ),
      file = PRINT_OUT,
      append = TRUE
    )
    res
  }
  environment(cv_func) <- .GlobalEnv 
  perf <- do.call("rbind", clusterApplyLB(cl, 1:length(paramspace), cv_func))
  stopCluster(cl)
  paramspace <- do.call("rbind", c(lapply(paramspace, function(params) {
    data.table(t(unlist(params)))
  }), fill = TRUE))
  dt <- data.table(cbind(paramspace, perf))
  setorderv(dt, cols = colnames(dt)[1:tail(grep("n_hidden", colnames(dt)), 1)])
}

L <- load_data(f_in = F_IN, vars = VARS)
perf <- explore_grid(
  data = L$X,
  weights = L$W,
  param_grid = PARAM_GRID,
  n_sample = N_SAMPLE,
  n_fold = N_FOLD,
  n_proc = N_PROC,
  n_perm = N_PERM,
  seed = SEED,
  param_seed = PARAM_SEED
)
save(perf, file = F_OUT)
