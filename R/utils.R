#' @useDynLib BayesSPsurv
#' @importFrom stats dgamma runif rgamma dnorm model.frame as.formula model.matrix model.response na.omit
#' @importFrom MCMCpack riwish
#' @importFrom coda mcmc
#' @import grDevices
#' @import graphics
#' @importFrom FastGP rcpp_rmvnorm rcpp_log_dmvnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr distinct
#' @importFrom reshape2 acast

.onUnload <- function (libpath) {
  library.dynam.unload("BayesSPsurv", libpath)
}


nameob <- c('capdist', 'numa', 'numb', 'kmdist', 'midist')

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(nameob))
  utils::suppressForeignCheck(c(nameob))
}



formcall <- function(duration,
                     immune,
                     Y0,
                     LY,
                     S = NULL,
                     A = NULL,
                     data,
                     N,
                     burn,
                     thin,
                     w = c(1, 1, 1),
                     m = 10,
                     form,
                     ini.beta =  0,
                     ini.gamma = 0,
                     ini.W = 0,
                     ini.V= 0,
                     prop.varV = NULL,
                     prop.varW = NULL,
                     model = character())
{

  formula1 <- as.formula(duration)
  formula2 <- as.formula(immune)
  variable <- unique(c(all.vars(formula1), all.vars(formula2)))
  mf1 <- model.frame(formula = duration, data = data)
  mf2 <- model.frame(formula = immune,   data = data)
  X <- model.matrix(attr(mf1, "terms"), data = mf1)
  Z <- model.matrix(attr(mf2, "terms"), data = mf2)
  Y <- as.matrix(model.response(mf1))
  C <- as.matrix(model.response(mf2))
  Y0 <- data[,Y0]
  LY <- data[,LY]
  burn <-  burn
  if (is.null(w)) w <- c(1,1,1) else w <- w
  if (is.null(m)) m <- 10 else m <- m
  if (is.null(ini.beta)) ini.beta <- 0 else ini.beta <- ini.beta
  if (is.null(ini.gamma)) ini.gamma <- 0 else ini.gamma <- ini.gamma
  if (is.null(ini.W)) ini.W <- 0 else ini.W <- ini.W
  if (is.null(ini.V)) ini.V <- 0 else ini.V <- ini.V
  form <-  form
  cnx <- colnames(X)
  cnz <- colnames(Z)

  if(model == "SPsurv"){

    dataset <- data.frame(cbind(Y, Y0, C, LY, X, Z))
    dataset <- na.omit(dataset)
    Y  <- as.matrix(dataset[,1])
    Y0 <- as.matrix(dataset[,2])
    C  <- as.matrix(dataset[,3])
    LY <- as.matrix(dataset[,4])
    X  <- as.matrix(dataset[,5:(4+ncol(X))])
    Z  <- as.matrix(dataset[,(5+ncol(X)):ncol(dataset)])
    colnames(X) <- cnx
    colnames(Z) <- cnz
    fm <- list(Y = Y, Y0 = Y0, C = C, LY = LY, X = X, Z = Z, N = N, burn = burn,
               thin = thin, w = w, m = m, ini.beta = ini.beta, ini.gamma = ini.gamma, form = form)

  } else {

    prop.varV <-  prop.varV
    prop.varW <-  prop.varW
    S  <- data[, S]
    dataset <- data.frame(cbind(Y, Y0, C, LY, S, X, Z))
    dataset <- na.omit(dataset)
    Y  <- as.matrix(dataset[,1])
    Y0 <- as.matrix(dataset[,2])
    C  <- as.matrix(dataset[,3])
    LY <- as.matrix(dataset[,4])
    S  <- as.matrix(dataset[,5])
    X  <- as.matrix(dataset[,6:(5+ncol(X))])
    Z  <- as.matrix(dataset[,(6+ncol(X)):ncol(dataset)])
    colnames(X) <- cnx
    colnames(Z) <- cnz
    fm <- list(Y = Y, Y0 = Y0, C = C, LY = LY, X = X, Z = Z, S = S, N = N, burn = burn,
               thin = thin, w = w, m = m, ini.beta = ini.beta, ini.gamma = ini.gamma, ini.W = ini.W, ini.V = ini.V,
               form = form, prop.varV = prop.varV, prop.varW = prop.varW)

  }

  if(!is.null(A)) fm$A <- as.matrix(A)

  return(fm)

}
