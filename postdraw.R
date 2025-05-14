# this function draws parameters
lmer.draw <- function(xobs, yobs, type, clust, rande, fixe, lev, X, Z, Xobs, Zobs){
  est <- lmer.est(xobs, yobs, type, clust, rande, fixe)
  sigma2star <- est$df * est$sigmahat^2/rchisq(1, est$df)
  sigma.star <- sqrt(sigma2star)
  covmat <- sigma2star * chol2inv(est$RX)
  rv <- t(chol(covmat))
  beta.star <- est$betahat + rv %*% rnorm(ncol(rv))
  lambda <- t(est$rancoef) %*% est$rancoef
  df.psi <- nrow(est$rancoef)
  temp.psi.star <- stats::rWishart(1, df.psi, diag(nrow(lambda)))[, , 1]
  temp <- MASS::ginv(lambda)
  ev <- eigen(temp)
  if (sum(ev$values > 0) == length(ev$values)) {
    deco <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
    psi.star <- MASS::ginv(deco %*% temp.psi.star %*% t(deco))
  }
  else {
    try(temp.svd <- svd(lambda))
    if (!inherits(temp.svd, "try-error")) {
      deco <- temp.svd$u %*% diag(sqrt(temp.svd$d), nrow = length(temp.svd$d))
      psi.star <- MASS::ginv(deco %*% temp.psi.star %*% 
                               t(deco))
    }
    else {
      psi.star <- est$psihat
      warning("psi fixed to estimate")
    }
  }
  bi.star <- matrix(NA, nrow = length(lev), ncol = dim(est$rancoef)[2])
  for (jj in lev) {
    if (jj %in% unique(xobs[, clust])) {
      Xi <- Xobs[xobs[, clust] == jj, ]
      Zi <- as.matrix(Zobs[xobs[, clust] == jj, , drop = FALSE])
      yi <- yobs[xobs[, clust] == jj]
      sigma2 <- diag(sigma2star, nrow = nrow(Zi))
      Mi <- psi.star %*% t(Zi) %*% MASS::ginv(Zi %*% psi.star %*% 
                                                t(Zi) + sigma2)
      myi <- Mi %*% (yi - Xi %*% beta.star)
      vyi <- psi.star - Mi %*% Zi %*% psi.star
    }
    else {
      myi <- matrix(0, nrow = nrow(psi.star), ncol = 1)
      vyi <- psi.star
    }
    vyi <- vyi - upper.tri(vyi) * vyi + t(lower.tri(vyi) * 
                                            vyi)
    deco1 <- eigen(vyi)
    if (sum(deco1$values > 0) == length(deco1$values)) {
      A <- deco1$vectors %*% sqrt(diag(deco1$values, nrow = length(deco1$values)))
      bistar <- myi + A %*% rnorm(length(myi))
    }
    else {
      try(deco1 <- svd(vyi))
      if (!inherits(deco1, "try-error")) {
        A <- deco1$u %*% sqrt(diag(deco1$d, nrow = length(deco1$d)))
        bistar <- myi + A %*% rnorm(length(myi))
      }
      else {
        bistar <- myi
        warning("b_", jj, " fixed to estimate")
      }
    }
    bi.star[jj, ] <- bistar
  }
  colnames(bi.star) <- colnames(est$rancoef)
  colnames(psi.star) <- colnames(est$psihat)
  rownames(psi.star) <- rownames(est$psihat)
  rownames(beta.star) <- names(est$betahat)
  parm <- list(est$betahat, est$sigmahat, est$psihat, beta.star, sigma.star, psi.star, est$rancoef, bi.star)
  names(parm) <- c("betahat", "sigmahat", "psihat", "betastar", "sigmastar", "psistar", "bihat", "bistar")
  parm
}



