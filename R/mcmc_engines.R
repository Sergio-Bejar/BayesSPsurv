# MCMC engine functions
mcmcSP <- function(Y,
                   Y0,
                   C,
                   LY,
                   X,
                   Z,
                   N,
                   burn,
                   thin,
                   w = c(1, 1, 1),
                   m = 10,
                   ini.beta =  0,
                   ini.gamma = 0,
                   form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  rho = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(0, length(Y))
  #V = rep(0, length(Y))
  delta = exp(Z %*% gammas)/ (1 + exp(Z %*% gammas))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)

  pb <- progress::progress_bar$new(total = N)
    for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    # if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas))
    gammas = gammas.slice.sampling(Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      delta.samp[(iter - burn) / thin] = mean(delta)

    }
  }
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp, delta = delta.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma)))
}

# @title mcmcSPlog
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return chain of the variables of interest
#
mcmcSPlog <- function(Y,
                      C,
                      Y0,
                      X,
                      LY,
                      Z,
                      N,
                      burn,
                      thin,
                      w = c(1, 1, 1),
                      m,
                      ini.beta =  0,
                      ini.gamma = 0,
                      ini.W = 0,
                      ini.V= 0,
                      form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  rho = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(0, length(Y))
  #V = rep(0, length(Y))
  delta = exp(Z %*% gammas)/ (1 + exp(Z %*% gammas))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)

  pb <- progress::progress_bar$new(total = N)
  for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    #if (iter %% 1000 == 0) print(iter) #F #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C,  LY, rho, w[1], m, form = form)
    eXB = exp((X %*% betas))
    gammas = gammas.slice.sampling3(Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
    }
  }
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma)))
}



# @title mcmcspatialSP
# @description Markov Chain Monte Carlo (MCMC) routine for Bayesian spatial split population survival model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param LY last observation year
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
# @param form type of parametric model (Exponential or Weibull)
# @param prop.varV proposal variance for Metropolis-Hastings
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcspatialSP <- function(Y,
                          Y0,
                          C,
                          LY,
                          X,
                          Z,
                          S,
                          A,
                          N,
                          burn,
                          thin, w = c(1, 1, 1),
                          m = 10,
                          ini.beta =  0,
                          ini.gamma = 0,
                          ini.W = 0,
                          ini.V= 0,
                          form,
                          prop.varV,
                          prop.varW,
                          id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  rho = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(ini.W, length(Y))
  V = rep(ini.V, length(Y))
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))

  pb <- progress::progress_bar$new(total = N)
  for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    #CAR model
    lambda = lambda.gibbs.sampling2(S, A, W, V)
    W = W.MH.sampling(S, A, lambda, Y, Y0,X, W, betas, delta, C, LY, rho, prop.varW)
    betas = betas.slice.sampling(Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.MH.sampling(S, A, lambda, Y,Y0, eXB, Z, V, gammas, C, LY, rho, prop.varV)
 	  gammas = gammas.slice.sampling2(Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
 	  num = exp(Z %*% gammas + V)
 	  num[which(is.infinite(num))] <- exp(700)
 	  denom = (1 + exp(Z %*% gammas + V))
 	  num[which(is.infinite(denom))] <- exp(700)
 	  delta = num/denom


    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }

  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, delta = delta.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma, ini.W = ini.W, ini.V  = ini.V)))
}


# @title mcmcSpatialLog
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian spatial cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @param prop.varV proposal variance for Metropolis-Hastings
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcSpatialLog <- function(Y,
                           Y0,
                           C,
                           LY,
                           X,
                           Z,
                           S,
                           A,
                           N,
                           burn,
                           thin,
                           w = c(1, 1, 1),
                           m,
                           ini.beta =  0,
                           ini.gamma = 0,
                           ini.W = 0,
                           ini.V= 0,
                           form,
                           prop.varV = .00001,
                           prop.varW = .00001,
                           id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  rho = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(ini.W, length(Y))
  V = rep(ini.V, length(Y))
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  #delta = 1/(1 + exp(-Z %*% gammas + V))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))

  pb <- progress::progress_bar$new(total = N)
  for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    #CAR model
    lambda = lambda.gibbs.sampling2(S, A, W, V)
    W = W.MH.sampling2(S, A, lambda, Y, Y0, X, W, betas, delta, C, LY, rho, prop.varW)
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp((X %*% betas) + W)
    V = V.MH.sampling2(S, A, lambda, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, prop.varV)
    gammas = gammas.slice.sampling4(Sigma.g, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    #delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))

    if (form == "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    # print(delta) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
      # print(100) #### ***** #### ***** #### ***** #### ***** ####
    }
  }
  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp, delta= delta.samp, lambda = lambda.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma, ini.W = ini.W, ini.V  = ini.V)))
}



# @title mcmcfrailtySP
# @description Markov Chain Monte Carlo (MCMC) routine to run Bayesian non-spatial frailties split population survival model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C censoring indicator
# @param LY last observation year
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
# @param form type of parametric model (Exponential or Weibull)
# @param prop.varV proposal variance for Metropolis-Hastings
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcfrailtySP <- function(Y,
                          Y0,
                          C,
                          LY,
                          X,
                          Z,
                          S,
                          N,
                          burn,
                          thin, w = c(1, 1, 1),
                          m,
                          ini.beta =  0,
                          ini.gamma = 0,
                          ini.W = 0,
                          ini.V= 0,
                          form,
                          prop.varV,
                          prop.varW,
                          id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  p3 = length(unique(S))
  p4 = length(unique(S))
  # initial values
  rho = 1
  lambda = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(ini.W, length(Y))
  V = rep(ini.V, length(Y))
  WS = rep(ini.W, p3)
  VS = rep(ini.V, p4)
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  Sigma.w = 10 * p3 * diag(p3)
  Sigma.v = 10 * p4 * diag(p4)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))

  pb <- progress::progress_bar$new(total = N)
  for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
      Sigma.w = riwish(1 + p3, WS %*% t(WS) + p3 * diag(p3))
      Sigma.v = riwish(1 + p4, VS %*% t(VS) + p4 * diag(p4))
    }
    #non-spatial Frailty model
    W = W.F.MH.sampling(Sigma.w, S,Y, Y0,X, W, betas, delta, C, LY, rho, prop.varW)
    betas = betas.slice.sampling(Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.F.MH.sampling(Sigma.v, S, Y,Y0, eXB, Z, V, gammas, C, LY, rho, prop.varV)
    gammas = gammas.slice.sampling2(Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, Y0,eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }


  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, delta = delta.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma, ini.W = ini.W, ini.V  = ini.V)))
}

# @title mcmc Cure with Non-spatial frailties
# @description Markov Chain Monte Carlo (MCMC) to run Bayesian spatial cure model
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param C censoring indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param S spatial information (e.g. district)
# @param N number of MCMC iterations
# @param burn burn-in to be discarded
# @param thin thinning to prevent from autocorrelation
# @param w size of the slice in the slice sampling for (betas, gammas, rho)
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
# @param prop.varV proposal variance for Metropolis-Hastings
# @param prop.varW proposal variance for Metropolis-Hastings
#
# @return chain of the variables of interest
#

mcmcfrailtySPlog <- function(Y,
                             Y0,
                             C,
                             LY,
                             X,
                             Z,
                             S,
                             N,
                             burn,
                             thin,
                             w = c(1, 1, 1),
                             m = 10,
                             ini.beta =  0,
                             ini.gamma = 0,
                             ini.W = 0,
                             ini.V= 0,
                             form,
                             prop.varV = .00001,
                             prop.varW = .00001,
                             id_WV) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  p3 = length(unique(S))
  p4 = length(unique(S))
  # initial values
  rho = 1
  lambda = 1
  betas = rep(ini.beta, p1)
  gammas = rep(ini.gamma, p2)
  W = rep(ini.W, length(Y))
  V = rep(ini.V, length(Y))
  WS = rep(ini.W, p3)
  VS = rep(ini.V, p4)
  delta = exp(Z %*% gammas + V)/ (1 + exp(Z %*% gammas + V))
  prop.varV = prop.varV
  prop.varW = prop.varW
  #delta = 1/(1 + exp(-Z %*% gammas + V))
  Sigma.b = 5 * p1 * diag(p1)
  Sigma.g = 5 * p2 * diag(p2)
  Sigma.w = 5 * p3 * diag(p3)
  Sigma.v = 5 * p4 * diag(p4)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  delta.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))
  V.samp = matrix(NA, nrow = (N - burn) / thin, ncol = length(unique(S)))

  pb <- progress::progress_bar$new(total = N)
  for (iter in 1:N) {
    pb$tick()
    Sys.sleep(1/N)
    #if (iter %% 1000 == 0) print(iter) #### ***** #### ***** #### ***** #### ***** ####
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
      Sigma.w = riwish(1 + p3, WS %*% t(WS) + p3 * diag(p3))
      Sigma.v = riwish(1 + p4, VS %*% t(VS) + p4 * diag(p4))
    }
    #non-spatial Frailty model.
    W = W.F.MH.sampling2(Sigma.w, S,Y, Y0, X, W, betas, delta, C, LY, rho, prop.varW)
    betas = betas.slice.sampling2(Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w[1], m, form = form)
    eXB = exp(-(X %*% betas) + W)
    V = V.F.MH.sampling2(Sigma.v, S, Y, Y0, eXB, Z, V, gammas, C, LY, rho, prop.varV)
    gammas = gammas.slice.sampling4(Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w[2], m, form = form)
    num = exp(Z %*% gammas + V)
    num[which(is.infinite(num))] <- exp(700)
    denom = (1 + exp(Z %*% gammas + V))
    num[which(is.infinite(denom))] <- exp(700)
    delta = num/denom

    if (form == "loglog") {
      rho = rho.slice.sampling2(Y, Y0, eXB, delta, C, LY, rho, w[3], m)
    }
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      delta.samp[(iter - burn) / thin] = mean(delta)
      S_uniq = unique(cbind(S, W, V))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
      V.samp[(iter - burn) / thin, ] = S_uniq[,3]
    }
  }


  colnames(V.samp) <- id_WV #colnames(A)           #ADC
  colnames(W.samp) <- id_WV #colnames(A)           #ADC
  colnames(betas.samp)  <- colnames(X) #adc
  colnames(gammas.samp) <- colnames(Z) #adc
  return(list(betas   = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, delta = delta.samp, W = W.samp, V = V.samp,
              spstats = list(X = X, Z = Z, Y = Y,  Y0 = Y0, C = C, S = S, form = form),
              initial = list(ini.beta = ini.beta, ini.gamma = ini.gamma, ini.W = ini.W, ini.V  = ini.V)))
}




# @title llFun
# @description Log-likelihood function for exchnageable and spatial frailty split population models
#
# @param est parameter estimates from the model
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C immunity indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param W survival stage (spatial) frailty parameter estimates
# @param V split stage (spatial) frailty parameter estimates
# @param data data that contains the variables of interest
#
# @return Log-likelihood of the exchangeable or the spatial split population survial model
#

