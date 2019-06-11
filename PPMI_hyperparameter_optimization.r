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

# some input file dependencies and parameters
USE_PYTHON <- NULL # "/nmitapps/lib/python/anaconda/3.6/bin/python3.6"
VADER_PATH <- file.path("..", "VaDER")
PRINT_OUT <- "PPMI_hyperparameter_optimization.out"
CODE_DIR <- "."
CODE_FILE <- "ADNI_hyperparameter_optimization.r"
PARAM_SEED <- 12345
DIR_IN <- file.path("..", "data", "PPMI")
DIR_OUT <- file.path("..", "results", "PPMI", "rgmvae", "hyperparameter_optimization", TIME_STAMP)
F_OUT <- file.path(DIR_OUT, sprintf("grid_search_seed%i.RData", SEED))
# F_IN <- file.path(DIR_IN, "PPMI.RData")
F_IN <- file.path("ADNI_artificial_data.RData")
N_SAMPLE <- 4e3 # Inf; 100 takes 90 minutes using 3x100 cores
N_PERM  <- 1e3
COMPUTE_PREDICTION_STRENGTH <- TRUE
if (COMPUTE_PREDICTION_STRENGTH) {
  N_FOLD <- 2
  # test one model for different k
  PARAM_GRID <- list(
    k = 2:15,
    n_layer = 2,
    alpha = 1.0,
    n_node = c(32, 128),
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
    n_node = c(1, 2, 4, 8, 16, 32, 64, 128),
    learning_rate = 10^(-4:-1),
    batch_size = 2^(4:7)
  )
}

library(abind)
library(caret)
library(data.table)
library(parallel)
library(caret)
library(gplots)
library(matrixStats)
library(reticulate)
# call directly after loading reticulate
if (!is.null(USE_PYTHON)) {
  use_python(USE_PYTHON, required = TRUE)
}
file.remove(PRINT_OUT)
dir.create(DIR_OUT, recursive = TRUE)

s <- parse(file.path(CODE_DIR, CODE_FILE))
# makes available the functions explore_grid and cross_validate 
# that we also used for ADNI
eval(s[grep("<- function\\(", s)])

load_data <- function(f_in) {
  s <- load(f_in)
  if ("y" %in% s) { # cluster label is present: artificial data
    X[W == 0] <- 0 # set missing values arbitrarily to 0 (can be any value)
    list(X = X, W = W, y = y)
  } else {
    W <- 1 - MNAR
    dimnames(X)[[3]] <- dimnames(W)[[3]] <- toupper(dimnames(X)[[3]])
    X[W == 0] <- 0 # set missing values arbitrarily to 0 (can be any value)
    list(X = X, W = W)
  }
}

L <- load_data(F_IN)

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

