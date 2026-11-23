## -----------------------------------------------------------------------------
## Parallel back end
##
## All parallel loops of the package go through `enetApply()`. It chooses the
## back end once, so that the rest of the code does not need to know about
## operating systems or clusters:
##
##   * `ncores == 1`            plain `lapply()` (no overhead)
##   * Unix-alikes (default)    forked workers via `parallel::mclapply()`
##   * Windows (default)        a PSOCK cluster, started once per top-level call
##                              by `startWorkers()` and stopped by `stopWorkers()`
##
## The back end can be forced with `options(enetLTS.backend = "fork" | "psock")`.
##
## None of the parallelized tasks draws random numbers inside a worker, so the
## results do not depend on the number of cores (see `cvFoldIds()` and
## `localSeed()`).
## -----------------------------------------------------------------------------

.enetEnv <- new.env(parent = emptyenv())


## Validate the `ncores` argument; `NA` means "use all available cores".
checkNcores <- function(ncores) {
  ncores <- rep(ncores, length.out = 1)
  if (is.na(ncores)) {
    ncores <- parallel::detectCores()
    if (is.na(ncores)) ncores <- 1L
  }
  if (!is.numeric(ncores) || !is.finite(ncores) || ncores < 1) {
    warning("invalid value of 'ncores'; using default value")
    return(1L)
  }
  as.integer(ncores)
}


## Name of the back end used when more than one core is requested.
parallelBackend <- function() {
  default <- if (.Platform$OS.type == "windows") "psock" else "fork"
  match.arg(getOption("enetLTS.backend", default), c("fork", "psock"))
}


## Start the PSOCK workers if (and only if) they are needed. The cluster is
## kept in the package environment so that nested internal functions can reach
## it without passing it around. Call `stopWorkers()` in `on.exit()`.
startWorkers <- function(ncores) {
  stopWorkers()
  if (ncores > 1L && parallelBackend() == "psock") {
    cl <- parallel::makePSOCKcluster(ncores)
    ok <- tryCatch({
      # workers must find the package wherever the main process found it
      # (`.libPaths` keeps its state in a closure, so it must be *called* on the
      # worker rather than shipped there)
      parallel::clusterCall(cl, eval, call(".libPaths", .libPaths()), baseenv())
      parallel::clusterCall(cl, loadNamespace, "enetLTS")
      TRUE
    }, error = function(e) {
      parallel::stopCluster(cl)
      stop("could not load package 'enetLTS' on the workers: ",
           conditionMessage(e), call. = FALSE)
    })
    if (ok) .enetEnv$cl <- cl
  }
  invisible(NULL)
}

stopWorkers <- function() {
  cl <- .enetEnv$cl
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    .enetEnv$cl <- NULL
  }
  invisible(NULL)
}


## `lapply()` that runs in parallel when `ncores > 1`. Errors raised inside a
## worker are re-thrown in the main process. Extra arguments in `...` must not
## be called `X`, `x`, `FUN`, `fun` or `cl` (they would clash with the
## arguments of the cluster functions); capture such objects in a closure.
enetApply <- function(X, FUN, ..., ncores = 1L) {
  if (ncores <= 1L || length(X) <= 1L) return(lapply(X, FUN, ...))

  cl <- .enetEnv$cl
  res <- if (!is.null(cl)) {
    parallel::parLapply(cl, X, FUN, ...)
  } else if (parallelBackend() == "fork" && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, ..., mc.cores = min(ncores, length(X)),
                       mc.allow.recursive = FALSE)
  } else {
    lapply(X, FUN, ...)
  }

  failed <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("a parallel worker failed: ",
         conditionMessage(attr(res[[which(failed)[1]]], "condition")),
         call. = FALSE)
  }
  res
}


## Evaluate `expr` with a fixed random seed and restore the previous state of
## the random number generator afterwards. Used wherever an algorithm needs
## random numbers inside a parallel task (e.g. the MCD in the multinomial
## C-step): the draw then depends only on the task, not on the worker.
localSeed <- function(seed, expr) {
  env <- globalenv()
  had.seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had.seed) old <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had.seed) assign(".Random.seed", old, envir = env)
    else if (exists(".Random.seed", envir = env, inherits = FALSE))
      rm(".Random.seed", envir = env)
  })
  set.seed(seed)
  expr
}

## A reproducible integer seed derived from a subset of observations.
subsetSeed <- function(index) {
  as.integer(sum(as.numeric(index) * seq_along(index)) %% 2147483647)
}
