## hsize = "optimal" (experimental): see ?enetLTS

hsizeCriterion <- function (fit, trim = 0.25, nfold = 5, repl = 2) 
{
    inp <- fit$inputs
    family <- inp$family
    xx <- inp$xx
    xx <- xx[, colnames(xx) != "(Intercept)", drop = FALSE]
    x <- standardize.x(xx, centerFun = median, scaleFun = mad)
    y <- if (family == "gaussian") 
        center.y(inp$yy, centerFun = median)
    else inp$yy
    spec <- cvSpec(family)
    groups <- spec$groups(y)
    H <- fit$best
    n <- nrow(x)
    loss <- localSeed(subsetSeed(H), {
        sapply(seq_len(repl), function(q) {
            fold <- cvFoldIds(groups, nfold)
            l <- rep(NA_real_, n)
            for (f in seq_len(nfold)) {
                test <- which(fold == f)
                train <- setdiff(H, test)
                g <- tryCatch(glmnet(x[train, , drop = FALSE], 
                  rowsOf(y, train), family = spec$glm.family, 
                  alpha = fit$alpha, lambda = spec$lambda(fit$lambda, 
                    length(train)), standardize = FALSE, intercept = FALSE), 
                  error = function(e) NULL)
                if (is.null(g)) 
                  return(rep(NA_real_, n))
                l[test] <- spec$loss(g, x[test, , drop = FALSE], 
                  rowsOf(y, test))
            }
            if (family == "gaussian") 
                l^2
            else l
        })
    })
    l <- rowMeans(matrix(loss, nrow = n))
    if (anyNA(l)) 
        return(Inf)
    mean(sort(l)[seq_len(floor(n * (1 - trim)))])
}

enetLTSoptimal <- function (call, env, hgrid) 
{
    if (!is.numeric(hgrid) || length(hgrid) < 2 || any(!is.finite(hgrid)) || 
        any(hgrid < 0.5 | hgrid > 1)) {
        stop("'hgrid' must contain at least two values between 0.5 and 1")
    }
    hgrid <- sort(unique(hgrid))
    fits <- lapply(hgrid, function(h) {
        cl <- call
        cl[[1]] <- enetLTS
        cl$hsize <- h
        cl$hgrid <- NULL
        eval(cl, env)
    })
    crit <- vapply(fits, hsizeCriterion, numeric(1))
    if (all(!is.finite(crit))) 
        stop("the criterion for 'hsize' could not be computed")
    crit.ok <- is.finite(crit)
    best <- max(which(crit.ok & crit <= min(crit[crit.ok]) * 
        1.1))
    if (hgrid[best] < 0.7) {
        warning("the optimal 'hsize' is ", hgrid[best], ": more than ", 
            round(100 * (1 - hgrid[best])), "% of the observations are trimmed. ", 
            "The data may contain a large fraction of outliers; interpret the ", 
            "results with care.", call. = FALSE)
    }
    out <- fits[[best]]
    out$hsize.search <- data.frame(hsize = hgrid, criterion = crit, 
        selected = seq_along(hgrid) == best)
    out$call <- call
    out
}

