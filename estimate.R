# This function estimates and extracts the random-effects model parameters
# It performs a similar job as the 'estimice()' function within '.norm.draw()'
lmer.est <- function(xobs, yobs, type, clust, rande, fixe){
  #install.on.demand("lme4", ...)
  fr <- ifelse(length(rande) > 1, paste("+ ( 1 +", paste(rande[-1L], 
                                                         collapse = "+")), "+ ( 1 ")
  randmodel <- paste("yobs ~ ", paste(fixe[-1L], collapse = "+"), 
                     fr, "|", clust, ")")
  suppressMessages(fit <- try(lme4::lmer(formula(randmodel), 
                                         data = data.frame(yobs, xobs)), silent = TRUE))
  if (inherits(fit, "try-error")) {
    warning("lmer does not run. Simplify imputation model")
    return(NA)
  }
  sigma <- function(object, ...) {
    dc <- object@devcomp
    dd <- dc$dims
    if (dd[["useSc"]]) {
      dc$cmp[[if (dd[["REML"]]) {
        "sigmaREML"
      }
      else {
        "sigmaML"
      }]]
    }
    else {
      1
    }
  }
  sigmahat <- sigma(fit)
  df <- nrow(fit@frame) - length(fit@beta)
  betahat <- lme4::fixef(fit)
  RX <- lme4::getME(fit, "RX")
  rancoef <- as.matrix(lme4::ranef(fit)[[1]])
  psi <- function(object, ...){
    dim <- lengths(object@cnms)
    temp <- matrix(unlist(VarCorr(object)), dim, dim)
    colnames(temp) <- rownames(temp) <- object@cnms[[1]] 
    temp
  }
  psihat <- psi(fit)
  return(list(sigmahat = sigmahat, df = df, betahat = betahat, RX = RX, rancoef = rancoef, psihat = psihat))
}


