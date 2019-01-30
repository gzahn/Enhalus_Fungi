`stressplot.wcmdscale` <-
  function(object, k = 2, pch,  p.col = "blue", l.col = "red", lwd = 2, ...)
  {
    ## Check that original distances can be reconstructed: this
    ## requires that all axes were calculated instead of 'k' first.
    hasdims <- NCOL(object$points)
    if (!is.null(object$negaxes))
      hasdims <- hasdims + NCOL(object$negaxes)
    if (hasdims < length(object$eig))
      stop("observed distances cannot be reconstructed: all axes were not calculated")
    ## Get the ordination distances in k dimensions
    if (k > NCOL(object$points))
      warning(gettextf("max allowed rank is k = %d", NCOL(object$points)))
    k <- min(NCOL(object$points), k)
    w <- sqrt(object$weights)
    u <- diag(w) %*% object$points
    odis <- dist(u[,1:k, drop = FALSE])
    ## Reconstitute the original observed distances
    dis <- dist(u)
    if (!is.null(object$negaxes))
      dis <- sqrt(dis^2 - dist(diag(w) %*% object$negaxes)^2)
    ## Remove additive constant to get original dissimilarities
    if (!is.na(object$ac)) {
      if (object$add == "lingoes")
        dis <- sqrt(dis^2 - 2 * object$ac)
      else if (object$add == "cailliez")
        dis <- dis - object$ac
      else
        stop("unknown Euclidifying adjustment: no idea what to do")
    }
    ##Plot
    if (missing(pch))
      if (length(dis) > 5000)
        pch <- "."
    else
      pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)
  }
