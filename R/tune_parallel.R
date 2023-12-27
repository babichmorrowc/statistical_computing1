tune_parallel <- function(METHOD, train.x, train.y = NULL, data = list(),
                          validation.x = NULL, validation.y = NULL,
                          ranges = NULL, predict.func = predict,
                          tunecontrol = tune.control(),
                          ...
) {
  call <- match.call()
  
  n_cores <- detectCores()
  registerDoParallel(n_cores * 0.8)
  
  ## internal helper functions
  resp <- function(formula, data) {
    
    model.response(model.frame(formula, data))
  }
  
  classAgreement <- function (tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev])) / n
    } else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m])) / n
    }
    p0
  }
  
  # Function for extracting the data for validation
  getValidationData <- function(train.x, validation.x, useFormula, data, train_ind_sample) {
    if (!is.null(validation.x)) {
      return(validation.x)
    } else if (useFormula) {
      return(data[-train_ind_sample,, drop = FALSE])
    } else if (inherits(train.x, "matrix.csr")) {
      return(train.x[-train_ind_sample, ])
    } else {
      return(train.x[-train_ind_sample,, drop = FALSE])
    }
  }
  
  # Function for extracting true y values
  getTrueY <- function(validation.y, useFormula, train.x, data, train_ind_sample) {
    if (!is.null(validation.y)) {
      return(validation.y)
    } else if (useFormula) {
      if (!is.null(validation.x)) {
        return(resp(train.x, validation.x))
      } else {
        return(resp(train.x, data[-train_ind_sample,]))
      }
    } else {
      return(train.y[-train_ind_sample])
    }
  }
  
  # Function for computing performance metric
  computeError <- function(true.y, pred, error_fun) {
    if (!is.null(error_fun)) {
      return(error_fun(true.y, pred))
    } else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) {
      # Classification error
      return(1 - classAgreement(table(pred, true.y)))
    } else if (is.numeric(true.y) && is.numeric(pred)) {
      # Mean squared error
      return(crossprod(pred - true.y) / length(pred))
    } else {
      stop("Dependent variable has the wrong type!")
    }
  }
  
  ## parameter handling
  if (tunecontrol$sampling == "cross")
    validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0))
    data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x))
  if (is.data.frame(train.y))
    train.y <- as.matrix(train.y)
  
  ## prepare training indices
  if (!is.null(validation.x)) tunecontrol$fix <- 1
  n <- nrow(if (useFormula) data else train.x)
  perm.ind <- sample(n)
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n)
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1)
      stop(sQuote("cross"), " must be greater than 1!")
  }
  train.ind <- if (tunecontrol$sampling == "cross")
    tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
  else if (tunecontrol$sampling == "fix")
    list(perm.ind[1:trunc(n * tunecontrol$fix)])
  else ## bootstrap
    lapply(1:tunecontrol$nboot,
           function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))
  
  ## find best model
  parameters <- if (is.null(ranges))
    data.frame(dummyparameter = 0)
  else
    expand.grid(ranges)
  p <- nrow(parameters)
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1)
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random),]
    p <- nrow(parameters)
  }
  model.variances <- model.errors <- c()
  
  sampling.errors <- foreach(para.set = 1:p, .combine = "cbind") %:%
    foreach(sample = 1:length(train.ind), .combine = "c") %dopar% {
      repeat.errors <- numeric(tunecontrol$nrepeat)
      
      pars <- if (is.null(ranges)) {
        NULL
      } else {
        lapply(parameters[para.set,,drop = FALSE], unlist)
      }
      
      train_data <- if (useFormula) {
        list(train.x, data = data, subset = train.ind[[sample]], ...)
      } else {
        list(train.x[train.ind[[sample]],], y = train.y[train.ind[[sample]]], ...)
      }
      
      # Train models and predict validation set
      models <- replicate(tunecontrol$nrepeat, {
        model <- do.call(METHOD, c(train_data, pars))
        predict.func(model, getValidationData(train.x, validation.x, useFormula, data, train.ind[[sample]]))
      }, simplify = FALSE)
      
      ## compute performance measure
      true.y <- getTrueY(validation.y, useFormula, train.x, data, train.ind[[sample]])
      
      if (is.null(true.y)) true.y <- rep(TRUE, length(models[[1]]))
      
      repeat.errors <- sapply(models, function(pred) {
        computeError(true.y, pred, tunecontrol$error.fun)
      })
      
      tunecontrol$repeat.aggregate(repeat.errors)
    }
  
  # clean up the cluster
  stopImplicitCluster()
  
  sampling.errors <- matrix(unlist(sampling.errors), nrow = length(train.ind), ncol = p)
  model.errors <- apply(sampling.errors, 2, tunecontrol$sampling.aggregate)
  model.variances <- apply(sampling.errors, 2, tunecontrol$sampling.dispersion)
  
  
  ## return results
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))
    NULL
  else
    lapply(parameters[best,,drop = FALSE], unlist)
  structure(list(best.parameters  = parameters[best,,drop = FALSE],
                 best.performance = model.errors[best],
                 method           = if (!is.character(METHOD))
                   deparse(substitute(METHOD)) else METHOD,
                 nparcomb         = nrow(parameters),
                 train.ind        = train.ind,
                 sampling         = switch(tunecontrol$sampling,
                                           fix = "fixed training/validation set",
                                           bootstrap = "bootstrapping",
                                           cross = if (tunecontrol$cross == n) "leave-one-out" else
                                             paste(tunecontrol$cross,"-fold cross validation", sep="")
                 ),
                 performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                 best.model       = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula)
                     do.call(METHOD, c(list(train.x, data = data),
                                       pars, list(...)))
                   else
                     do.call(METHOD, c(list(x = train.x,
                                            y = train.y),
                                       pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }
  ),
  class = "tune"
  )
}
