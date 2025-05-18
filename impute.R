# For now, it is called '2l.lmer.me'
# I am unsure whether the model argument is correctly specified when task == "train". Please check!
mice.impute.2l.lmer.me <- function (y, ry, x, type, wy = NULL, intercept = TRUE, task = "impute", model = NULL, ...){
  check.model.exists(model, task)
  method <- "2l.lmer.me"
  if (is.null(wy)) 
    wy <- !ry
  if (intercept) {
    x <- cbind(1, as.matrix(x))
    type <- c(2, type)
    names(type)[1] <- colnames(x)[1] <- "(Intercept)"
  }
  clust <- names(type[type == -2])
  rande <- names(type[type == 2])
  fixe <- names(type[type > 0])
  lev <- unique(x[, clust])
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]
  if (task == "fill" && check.model.match(model, x, method)){
    for (j in lev) {
      y[wy & x[, clust] == j] <- as.vector(as.matrix(X[wy & x[, clust] == j, , drop = FALSE]) %*% model$betastar + 
                                             as.matrix(Z[wy & x[, clust] == j, , drop = FALSE]) %*% 
                                             as.matrix(model$bistar[j, ]) + 
                                             rnorm(sum(wy & x[, clust] == j)) * model$sigmastar)
    }
    return(y[wy])
  }
  parm <- lmer.draw(xobs, yobs, type, clust, rande, fixe, lev, X, Z, Xobs, Zobs)
  if (task == "train"){
    model <- parm
    model$setup <- list(method = method, task = task)
  }
  for (j in lev) {
    y[wy & x[, clust] == j] <- as.vector(as.matrix(X[wy & x[, clust] == j, , drop = FALSE]) %*% parm$betastar + 
                                           as.matrix(Z[wy & x[, clust] == j, , drop = FALSE]) %*% 
                                           as.matrix(parm$bistar[j, ]) + 
                                           rnorm(sum(wy & x[, clust] == j)) * parm$sigmastar)
  }
  y[wy]
}



