# Posterior density functions and Metropolis-Hastings samplers
betas.post = function(betas.p,
                      p, Sigma.b,
                      Y,
                      Y0,
                      X,
                      W,
                      betas,
                      delta,
                      C,
                      LY,
                      rho,
                      form) {
  betas[p] = betas.p
  eXB = exp(-(X %*% betas) + W)
  lprior = rcpp_log_dmvnorm(Sigma.b, rep(0, length(betas)), betas, FALSE)
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}
# @title betas.post
# @description log-posterior distribution of betas with pth element fixed as betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Log Lokelihood)
#
# @return log- posterior density of betas
#
betas.post2 = function(betas.p, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) {
  betas[p] = betas.p
  eXB = exp(-(X %*% betas) + W)
  lprior = rcpp_log_dmvnorm(Sigma.b, rep(0, length(betas)), betas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}
# @title gammas.post
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post = function(gammas.p,
                       p,
                       Sigma.g,
                       Y,
                       Y0,
                       eXB,
                       Z,
                       gammas,
                       C,
                       LY,
                       rho,
                       form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull(Y, Y0,eXB, delta, C, LY, rho) + lprior
  return(lpost)
}


# @title gammas.post2
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post2 = function(gammas.p,
                        p,
                        Sigma.g,
                        Y,
                        Y0,
                        eXB, Z, V, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull(Y, Y0,eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title gammas.post3
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Log Likelihood)
#
# @return log- posterior density of betas
#

gammas.post3 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title gammas.post4
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#

gammas.post4 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) {
  gammas[p] = gammas.p
  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title rho.post
# @description log-posterior distribution of rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#

rho.post = function(Y,
                    Y0,
                    eXB,
                    delta,
                    C,
                    LY,
                    rho,
                    a = 1,
                    b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title rho.post2
# @description log-posterior distribution of rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#


rho.post2 = function(Y, Y0, eXB, delta, C, LY, rho, a = 1, b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W (Weibull or exponential)
#

W.post = function(S,
                  A,
                  lambda,
                  Y,
                  Y0,
                  X,
                  W,
                  betas,
                  delta,
                  C,
                  LY,
                  rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)
  m_j = length(adj)

  if(m_j == 0){
    m_j <- exp(-700)
  } else {
    m_j <- m_j
  }

  wj = S_uniq[which(S_uniq[,1] %in% adj),2]

  if (length(wj) == 0){
    W_j_bar = 0
  } else {
    W_j_bar = mean(wj)
  }

  lprior = lprior + dnorm(S_uniq[s, 2], W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W (Log likelihood)
#
W.post2 = function(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
    adj = which(A[s,] == 1)
    m_j = length(adj)
    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    wj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(wj) == 0){
      W_j_bar = 0
    } else {
      W_j_bar = mean(wj)
    }

    lprior = lprior + dnorm(S_uniq[s, 2], W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}



# @title W.MH.sampling
# @description MH Sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

W.MH.sampling = function(S,
                         A,
                         lambda,
                         Y,
                         Y0,
                         X,
                         W,
                         betas,
                         delta,
                         C,
                         LY,
                         rho,
                         prop.varW) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.varW * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.post(S, A, lambda, Y, Y0,X, W_new[S], betas, delta, C, LY, rho)
  w1[which(w1==-Inf)] <- -740
  w2 = W.post(S, A, lambda, Y,Y0, X, W_old[S], betas, delta, C, LY, rho)
  w2[which(w2==-Inf)] <- -740
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
  	W = W_new[S]
  } else {
  	W = W_old[S]
  }
  return(W)
}


# @title W.MH.sampling2
# @description slice sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
W.MH.sampling2 = function(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho, prop.varW) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.varW * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.post2(S, A, lambda, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  #w1[which(w1==-Inf)] <- -740
  #w1[which(w1==Inf)] <- exp(700)
  w2 = W.post2(S, A, lambda, Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  #w2[which(w2==Inf)] <- exp(700)
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}

# @title W.F.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W
#

W.F.post = function(Sigma.w,
                    S,
                    Y,
                    Y0,
                    X,
                    W,
                    betas,
                    delta,
                    C,
                    LY,
                    rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.w[s,s])
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.F.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of W
#
W.F.post2 = function(Sigma.w, S, Y, Y0, X, W, betas, delta, C, LY, rho) {
  eXB = exp(-(X %*% betas) + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.w[s,s])

  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title W.F.MH.sampling (Cure Model with Frailties)
# @description MH sampling for W
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

W.F.MH.sampling = function(Sigma.w,
                           S,
                           Y,
                           Y0,
                           X,
                           W,
                           betas,
                           delta,
                           C,
                           LY,
                           rho,
                           prop.varW) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.varW * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.F.post(Sigma.w, S, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  w2 = W.F.post(Sigma.w, S, Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}

# @title W.F.MH.sampling2 (Cure Model with Frailties)
# @description MH sampling for W
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
W.F.MH.sampling2 = function(Sigma.w, S, Y, Y0, X, W, betas, delta, C, LY, rho, prop.varW) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.varW * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  w1 = W.F.post2(Sigma.w, S, Y, Y0, X, W_new[S], betas, delta, C, LY, rho)
  w1[which(w1==-Inf)] <- -740
  w2 = W.F.post2(Sigma.w, S,  Y, Y0, X, W_old[S], betas, delta, C, LY, rho)
  w2[which(w2==-Inf)] <- -740
  temp = w1- w2
  alpha = min(0, temp)
  if (u <= alpha) {
    W = W_new[S]
  } else {
    W = W_old[S]
  }
  return(W)
}



# @title V.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas
#

V.post = function(S,
                  A,
                  lambda,
                  Y,
                  Y0,
                  eXB,
                  Z,
                  V,
                  gammas,
                  C,
                  LY,
                  rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)

  m_j = length(adj)
  if(m_j == 0){
    m_j <- exp(-700)
  } else {
    m_j <- m_j
  }

  vj = S_uniq[which(S_uniq[,1] %in% adj),2]
  if (length(vj) == 0){
    V_j_bar = 0
  } else {
    V_j_bar = mean(vj)
  }

  lprior = lprior + dnorm(S_uniq[s, 2], V_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title V.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (log likelihood)
#
V.post2 = function(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C, LY, rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
    adj = which(A[s,] == 1)

    m_j = length(adj)
    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    vj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(vj) == 0){
      V_j_bar = 0
    } else {
      V_j_bar = mean(vj)
    }

    lprior = lprior + dnorm(S_uniq[s, 2], V_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}





# @title V.MH.sampling
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varV proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

V.MH.sampling = function(S,
                         A,
                         lambda,
                         Y,
                         Y0,
                         eXB,
                         Z,
                         V,
                         gammas,
                         C,
                         LY,
                         rho,
                         prop.varV) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.varV * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.post(S, A, lambda, Y,Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  v1[which(v1==-Inf)] <- -740
  v2 = V.post(S, A, lambda, Y, Y0,eXB, Z, V_old[S], gammas, C, LY, rho)
  v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
  	V = V_new[S]
  } else {
  	V = V_old[S]
  }
  return(V)
}

# @title V.MH.sampling2
# @description slice sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varV proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
V.MH.sampling2 = function(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.varV) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.varV * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.post2(S, A, lambda, Y, Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  #v1[which(v1==-Inf)] <- -740
  #V1[which(V1==Inf)] <- exp(700)
  v2 = V.post2(S, A, lambda, Y, Y0, eXB, Z, V_old[S], gammas, C, LY, rho)
  #v2[which(v2==-Inf)] <- -740
  #v2[which(v2==Inf)] <- exp(700)
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}



# @title V.F.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (Weibull and Exponential)
#

V.F.post = function(Sigma.v,
                    S,
                    Y,
                    Y0,
                    eXB,
                    Z,
                    V,
                    gammas,
                    C,
                    LY,
                    rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.v[s,s])
  }
  lpost = llikWeibull(Y,Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}

# @title V.F.post2
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas (log likelihood)
#

V.F.post2 = function(Sigma.v, S, Y, Y0, eXB, Z, V, gammas, C, LY, rho) {

  num = exp(Z %*% gammas + V)
  num[which(is.infinite(num))] <- exp(700)
  denom = (1 + exp(Z %*% gammas + V))
  num[which(is.infinite(denom))] <- exp(700)
  delta = num/denom

  #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:length(unique(S))){
    lprior = lprior + dnorm(S_uniq[s, 2], 0, Sigma.v[s,s])
  }
  lpost = llikLoglog(Y, Y0, eXB, delta, C, LY, rho) + lprior
  return(lpost)
}


# @title V.F.MH.sampling (Cure Model with non-spatial Frailties)
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varV proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#

V.F.MH.sampling = function(Sigma.v,
                           S,
                           Y,
                           Y0,
                           eXB,
                           Z,
                           V,
                           gammas,
                           C,
                           LY,
                           rho,
                           prop.varV) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.varV * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.F.post(Sigma.v,S, Y,Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  v1[which(v1==-Inf)] <- -740
  v2 = V.F.post(Sigma.v,S, Y, Y0,eXB, Z, V_old[S], gammas, C, LY, rho)
  v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}

# @title V.F.MH.sampling2 (Cure Model with non-spatial Frailties)
# @description MH sampling for rcpp_log_dmvnorm
#
# @param S spatial information (e.g. district)
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.varV proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling (log likelihood)
#
V.F.MH.sampling2 = function(Sigma.v, S,  Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.varV) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.varV * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  v1 = V.F.post2(Sigma.v, S, Y, Y0, eXB, Z, V_new[S], gammas, C, LY, rho)
  #v1[which(v1==-Inf)] <- -740
  v2 = V.F.post2(Sigma.v, S, Y, Y0, eXB, Z, V_old[S], gammas, C, LY, rho)
  #v2[which(v2==-Inf)] <- -740
  temp = v1- v2
  alpha = min(0, temp)
  if (u <= alpha) {
    V = V_new[S]
  } else {
    V = V_old[S]
  }
  return(V)
}


# @title lambda.gibbs.sampling2
# @description log-posterior distribution of rho
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param W spatial random effects
# @param V spatial random effects
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#

lambda.gibbs.sampling2 <- function(S,
                                   A,
                                   W,
                                   V,
                                   a = 1,    # prior
                                   b = 1) {  # prior
  S_uniq = unique(cbind(S, W, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  J = nrow(S_uniq)
  sums = 0
  for (j in 1:J) {
    adj = which(A[j,]==1)
    m_j = length(adj)

    if(m_j == 0){
      m_j <- exp(-700)
    } else {
      m_j <- m_j
    }

    wj = S_uniq[which(S_uniq[,1] %in% adj),2]
    if (length(wj) == 0){
      W_j_bar = 0
    } else {
      W_j_bar = mean(wj)
    }

    vj = S_uniq[which(S_uniq[,1] %in% adj),3]
    if (length(vj) == 0){
      V_j_bar = 0
    } else {
      V_j_bar = mean(vj)
    }

    W_j = S_uniq[j, 2]
    V_j = S_uniq[j, 3]

    sums = sums + m_j/2 * ((W_j-W_j_bar)^2 + (V_j-V_j_bar)^2)

  }
  lambda = rgamma(1, J + a, sums + b)
  return(lambda)
}

