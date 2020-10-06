smacof <- function(D, niter = 100, interc = 1, inicon = NULL, cols = "black", colv = palette()[c(8, 2, 4, 3, 5, 6, 7, 1)], main = "Multidimensional Scaling", k = 2, pch = 16, PLOT = TRUE, VERBOSE = TRUE, ...) {
  dhat <- ddelta <- D <- normdist(D)
  n <- dimdata(D)
  delta <- lowmat(D)
  ycol <- matrix(c(0, 0), nrow = 1, byrow = T)
  
  if (is.null(inicon)) {
    inicon <- cmdscale(D, k)
  }
  yy <- xx <- inicon
  
  dd <- lowmat(dist(xx), ...)
  str0 <- str <- sum((delta - dd)**2)
  for (i in 1:niter)
  {
    if (niter != 0) {
      xx <- (1 / n) * (bmat(delta, dd)) %*% xx
      dd <- dist(xx)
      if (interc == 1) {
        ddd <- dd
        dhat <- lm(dd ~ ddelta)$fitted.values
        delta <- lowmat(normdist(dhat))
        if ((a <- min(dhat)) <= 0) {
          ddd <- .01 + dd + abs(a)
          dhat <- lm(ddd ~ ddelta)$fitted.values
          delta <- lowmat(normdist(dhat))
        }
      }
      dd <- lowmat(dist(xx))
      # str=(sum((delta-dd)**2))**.5
      str <- snorm <- fit(delta, dd)
      if (VERBOSE) cat("it = ", i, ",  crit = ", str, " \n")
      if (abs(str0 - str) < 1e-6) break
      str0 <- str
    }
  }
  xx <- rotat(xx)
  
  
  if (PLOT) {
    lims <- c(
      ll = floor(100 * min(c(xx[, 1], xx[, 2]))) / 100,
      ul = ceiling(100 * max(c(xx[, 1], xx[, 2]))) / 100
    )
    
    # plot(xx,asp=1,type='n',xlab='Dimension 1 ',ylab='Dimension 2')
    # points(xx,asp=1,col="gray",pch=21,bg="gray")
    plot(xx,
         xlab = "Dimension 1", ylab = "Dimension 2",
         main = main, xlim = lims, ylim = lims,
         xaxs = "i", yaxs = "i", pch = pch, col = cols, ...
    )
  }
  
  delta <- normdist(dhat)
  return(list(D = delta, X = xx, interc = interc))
}

plotSmacof <- function(xx, main = "Multidimensional Scaling", pch = 16, cols = "black", ...){
  lims <- c(
    ll = floor(100 * min(c(xx[, 1], xx[, 2]))) / 100,
    ul = ceiling(100 * max(c(xx[, 1], xx[, 2]))) / 100
  )
  
  plot(xx,
       xlab = "Dimension 1", ylab = "Dimension 2",
       main = main, xlim = lims, ylim = lims,
       xaxs = "r", yaxs = "r", pch = pch, col = cols, ...
  )
}

plotNoLims <- function(xx, main = "Plot", xlab = "Dimension 1", ylab = "Dimension 2", pch = 16, cols = "black", ...){
  plot(xx,
       xlab = xlab, ylab = ylab,
       main = main,
       xaxs = "r", yaxs = "r", pch = pch, col = cols, ...
  )
}

bmat <- function(odist, fdist) {
  n <- nrow(odist)
  nnn <- matrix(0, nrow = n, ncol = n)
  bb <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    nnn[i, i] <- 1
  }
  for (i in 1:n)
  {
    sb <- 0
    for (j in 1:n)
    {
      if (fdist[i, j] != 0) {
        bb[i, j] <- odist[i, j] / fdist[i, j]
        sb <- sb - bb[i, j]
      }
    }
    bb[i, i] <- bb[i, i] + sb
  }
  bmat <- -bb
  return(bmat)
}
