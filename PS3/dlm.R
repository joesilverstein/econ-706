function (...) 
{
    if (nargs() == 1) 
        x <- as.list(...)
    else x <- list(...)
    nm <- c("m0", "C0", "FF", "V", "GG", "W")
    nmInd <- match(nm, names(x))
    if (any(is.na(nmInd))) 
        stop(paste("Component(s)", paste(nm[is.na(nmInd)], collapse = ", "), 
            "is (are) missing"))
    x[nmInd[-1]] <- lapply(x[nmInd[-1]], as.matrix)
    if (!is.numeric(x$FF)) 
        stop("Component FF must be numeric")
    m <- nrow(x$FF)
    p <- ncol(x$FF)
    if (!is.numeric(x$V)) 
        stop("Component V must be numeric")
    if (!(nrow(x$V) == m && ncol(x$V) == m)) 
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$GG)) 
        stop("Component GG must be numeric")
    if (!(nrow(x$GG) == p && ncol(x$GG) == p)) 
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$W)) 
        stop("Component W must be numeric")
    if (!(nrow(x$W) == p && ncol(x$W) == p)) 
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$C0)) 
        stop("Component C0 must be numeric")
    if (!(nrow(x$C0) == p && ncol(x$C0) == p)) 
        stop("Incompatible dimensions of matrices")
    if (!(is.numeric(x$m0) && NCOL(x$m0) == 1 && NROW(x$m0) == 
        p)) 
        stop(paste("Component m0 must be a numeric vector of length", 
            "\n equal to ncol of component FF, or a matrix with one column and", 
            "\n number of rows equal to ncol of component FF"))
    if (!(all.equal(x$C0, t(x$C0)) && all(eigen(x$C0)$values >= 
        0))) 
        stop("C0 is not a valid variance matrix")
    if (any(c(is.na(x$m0), is.na(x$C0)))) 
        stop("Missing values are not allowed in components m0 and C0")
    nm1 <- c("JFF", "JV", "JGG", "JW")
    nm1Ind <- match(nm1, names(x))
    if (all(is.na(nm1Ind))) {
        if (!(all.equal(x$V, t(x$V)) && all(eigen(x$V)$values >= 
            0))) 
            stop("V is not a valid variance matrix")
        if (!(all.equal(x$W, t(x$W)) && all(eigen(x$W)$values >= 
            0))) 
            stop("W is not a valid variance matrix")
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    x[nm1Ind[!is.na(nm1Ind)]] <- lapply(x[nm1Ind[!is.na(nm1Ind)]], 
        function(x) if (!is.null(x)) 
            as.matrix(x))
    if (!is.null(x$JFF)) {
        if (!(is.numeric(x$JFF) && nrow(x$JFF) == m && ncol(x$JFF) == 
            p)) 
            stop("Invalid component JFF")
        JFF <- round(x$JFF)
        if (all(JFF == 0)) 
            JFF <- NULL
    }
    else JFF <- NULL
    if (!is.null(x$JV)) {
        if (!(is.numeric(x$JV) && nrow(x$JV) == m && ncol(x$JV) == 
            m)) 
            stop("Invalid component JV")
        JV <- round(x$JV)
        if (all(JV == 0)) 
            JV <- NULL
    }
    else JV <- NULL
    if (!is.null(x$JGG)) {
        if (!(is.numeric(x$JGG) && nrow(x$JGG) == p && ncol(x$JGG) == 
            p)) 
            stop("Invalid component JGG")
        JGG <- round(x$JGG)
        if (all(JGG == 0)) 
            JGG <- NULL
    }
    else JGG <- NULL
    if (!is.null(x$JW)) {
        if (!(is.numeric(x$JW) && nrow(x$JW) == p && ncol(x$JW) == 
            p)) 
            stop("Invalid component JW")
        JW <- round(x$JW)
        if (all(JW == 0)) 
            JW <- NULL
    }
    else JW <- NULL
    mx <- max(c(JFF, JV, JGG, JW))
    if (mx <= 0) {
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    if (is.null(x$X)) 
        stop("Component X must be provided for time-varying models")
    x$X <- as.matrix(x$X)
    if (!(is.numeric(x$X) && ncol(x$X) >= mx)) 
        stop("Invalid component X")
    mod <- c(x[nmInd], list(JFF = JFF, JV = JV, JGG = JGG, JW = JW, 
        X = x$X))
    class(mod) <- "dlm"
    return(mod)
}
<environment: namespace:dlm>