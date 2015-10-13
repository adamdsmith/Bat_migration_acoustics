r.squaredGAMM <- function(gammobj, data, offset = NULL) {
    # Modified MuMIn::r.squaredGLMM
    # Only accepts models with a single smooth term for now
    gamobj <- gammobj$gam
    if (is.null(offset)) {
        VarFx <- var(predict(gamobj))
    } else {
        VarFx <- var(predict(gamobj) - data[, offset]) # offset on log scale
    }
    
    gamData <- data.frame(Intercept = 1, Xr = gamobj$model$Xr)
    names(gamData)[1] <- "(Intercept)"
    names(gamData) <- gsub("[.]", "", names(gamData))
    mmRE <- as.matrix(gamData)
    #mmRE <- model.matrix(gammobj$lme$modelStruct$reStruct, 
    #                     data = gamData[ , , drop = FALSE])
    n <- nrow(mmRE)
    sigma2 <- gamobj$sig2
    reStruct <- gammobj$lme$modelStruct$reStruct
    #if ((m <- length(reStruct)) > 1L) {
    #    nams <- names(reStruct)
    #    for (i in seq.int(m)) {
    #        attr(reStruct[[i]], "Dimnames")[[2L]] <- paste(nams[[i]], attr(reStruct[[i]], 
    #                                                                       "Dimnames")[[2L]], sep = ".")
    #}
    varRe <- sum(sapply(reStruct, function(z) {
        sig <- nlme::pdMatrix(z) * sigma2
        mm1 <- mmRE[, rownames(sig), drop = FALSE]
        sum(MuMIn:::matmultdiag(mm1 %*% sig, ty = mm1))/n
    }))
    varTot <- sum(VarFx, varRe)
    res <- c(VarFx, varTot)/(varTot + sigma2)
    names(res) <- c("R2m", "R2c")
    res
}
