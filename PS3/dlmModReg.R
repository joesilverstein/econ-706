function (X, addInt = TRUE, dV = 1, dW = rep(0, NCOL(X) + addInt), 
    m0 = rep(0, length(dW)), C0 = 1e+07 * diag(nrow = length(dW))) 
{
    p <- NCOL(X) + addInt
    if (!(length(dV) == 1 && length(dW) == p && length(m0) == 
        p && nrow(C0) == p && ncol(C0) == p)) 
        stop("Inconsistent dimensions of arguments")
    X <- as.matrix(X)
    JFF <- matrix(1:ncol(X), nrow = 1)
    if (addInt) 
        JFF <- cbind(0, JFF)
	mod <- list(m0 = m0, C0 = C0, FF = matrix(1, 1, p), V = as.matrix(dV), 
        GG = diag(nrow = p), W = diag(x = dW, nrow = p), JFF = JFF, 
        JV = NULL, JGG = NULL, JW = NULL, X = X)
    class(mod) <- "dlm"
    return(mod)
}
<environment: namespace:dlm>