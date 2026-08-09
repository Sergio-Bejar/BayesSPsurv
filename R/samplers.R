# Slice sampling functions for betas, gammas, and rho
betas.slice.sampling = function(Sigma.b,
                                Y,
                                Y0,
                                X,
                                W,
                                betas,
                                delta,
                                C,
                                LY,
                                rho,
                                w,
                                m,
                                form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling(betas[p], p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, w, m, form = form)
  }
  return(betas)
}

# @title betas.slice.sampling2
# @description slice sampling for betas
#
# @param Sigma.b variance estimate of betas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
betas.slice.sampling2 = function(Sigma.b, Y, Y0, X, W, betas, delta, C,  LY, rho, w, m, form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling2(betas[p], p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w, m, form = form)
  }
  return(betas)
}


# @title univ.betas.slice.sampling
# @description univariate slice sampling for betas.p
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
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

univ.betas.slice.sampling = function(betas.p,
                                     p,
                                     Sigma.b,
                                     Y,
                                     Y0,
                                     X,
                                     W,
                                     betas,
                                     delta,
                                     C,
                                     LY,
                                     rho,
                                     w,
                                     m,
                                     lower = -Inf,
                                     upper = +Inf,
                                     form) {
  b0 = betas.p
  b.post0 = betas.post(b0, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)

  e.b.post0 = exp(b.post0)
  if (is.infinite(e.b.post0)){e.b.post0 = exp(700)}

  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, e.b.post0))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y,Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y,Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post(b1, p, Sigma.b, Y, Y0,X, W, betas, delta, C, LY, rho, form)

    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}



# @title univ.betas.slice.sampling2
# @description univariate slice sampling for betas.p
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
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#
univ.betas.slice.sampling2 = function(betas.p, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  b0 = betas.p
  b.post0 = betas.post2(b0, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)

  e.b.post0 = exp(b.post0)
  if (is.infinite(e.b.post0)){e.b.post0 = exp(700)}

  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, e.b.post0))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, is.numeric(m)))
    K = (is.numeric(m) - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post2(b1, p, Sigma.b, Y, Y0, X, W, betas, delta, C, LY, rho, form)
    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}




# @title gammas.slice.sampling
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling = function(Sigma.g,
                                 Y,
                                 Y0,
                                 eXB,
                                 Z,
                                 gammas,
                                 C,
                                 LY,
                                 rho,
                                 w,
                                 m,
                                 form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling(gammas[p], p, Sigma.g, Y,Y0, eXB, Z, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}


# @title gammas.slice.sampling2
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling2 = function(Sigma.g,
                                  Y,
                                  Y0,
                                  eXB,
                                  Z,
                                  V,
                                  gammas,
                                  C,
                                  LY,
                                  rho,
                                  w,
                                  m,
                                  form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling2(gammas[p], p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}

# @title gammas.slice.sampling3
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling3 = function(Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling3(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}
# @title gammas.slice.sampling4
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#

gammas.slice.sampling4 = function(Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling4(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, form = form)
  }
  return(gammas)
}



# @title univ.gammas.slice.sampling2
# @description univariate slice sampling for gammas.p
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
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

univ.gammas.slice.sampling2 = function(gammas.p,
                                       p, Sigma.g,
                                       Y,
                                       Y0,
                                       eXB,
                                       Z,
                                       V,
                                       gammas,
                                       C,
                                       LY,
                                       rho,
                                       w,
                                       m,
                                       lower = -Inf,
                                       upper = +Inf,
                                       form) {
  g0 = gammas.p
  g.post0 = gammas.post2(g0, p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y,Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post2(g1, p, Sigma.g, Y, Y0,eXB, Z, V, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title univ.gammas.slice.sampling
# @description univariate slice sampling for gammas.p
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
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#

univ.gammas.slice.sampling = function(gammas.p,
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
                                      w,
                                      m,
                                      lower = -Inf,
                                      upper = +Inf,
                                      form) {
  g0 = gammas.p
  g.post0 = gammas.post(g0, p, Sigma.g, Y, Y0,eXB, Z,gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y,Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post(g1, p, Sigma.g, Y, Y0,eXB, Z, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}


# @title univ.gammas.slice.sampling3
# @description univariate slice sampling for gammas.p
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
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
univ.gammas.slice.sampling3 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post3(g0, p, Sigma.g, Y, Y0, eXB, Z,gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post3(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post3(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post3(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post3(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post3(g1, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title univ.gammas.slice.sampling4
# @description univariate slice sampling for gammas.p
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
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Log Likelihood)
#
# @return One sample update using slice sampling
#
univ.gammas.slice.sampling4 = function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post4(g0, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form)

  e.g.post0 = exp(g.post0)

  if (is.infinite(e.g.post0)){e.g.post0 = exp(700)}

  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, e.g.post0))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post4(L, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post4(R, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, is.numeric(m)))
    K = (is.numeric(m) - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post4(L, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C,  LY, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post4(R, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post4(g1, p, Sigma.g, Y, Y0, eXB, Z, V, gammas, C, LY, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}



# @title rho.slice.sampling
# @description univariate slice sampling for rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling. A vector of values for beta, gamma, rho.
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
#
# @return One sample update using slice sampling
#

rho.slice.sampling = function(Y,
                              Y0,
                              eXB,
                              delta,
                              C,
                              LY,
                              rho,
                              w,
                              m,
                              lower = 0.01,
                              upper = +Inf) {
  l0 = rho
  l.post0 = rho.post(Y, Y0,eXB, delta, C, LY, l0)

  e.l.post0 = exp(l.post0)
  if (is.infinite(e.l.post0)){e.l.post0 = exp(700)}

  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, e.l.post0))}


  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, LY, R))) #browser()
      if (rho.post(Y, Y0,eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, LY, R))) #browser()
      if (rho.post(Y, Y0,eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = rho.post(Y, Y0,eXB, delta, C, LY, l1)

    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}


# @title rho.slice.sampling2
# @description univariate slice sampling for rho
#
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
#
# @return One sample update using slice sampling
#
rho.slice.sampling2 = function(Y, Y0, eXB, delta, C, LY, rho, w, m, lower = 0.01, upper = +Inf) {
  l0 = rho
  l.post0 = rho.post2(Y, Y0, eXB, delta, C, LY, l0)

  e.l.post0 = exp(l.post0)
  if (is.infinite(e.l.post0)){e.l.post0 = exp(700)}

  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, e.l.post0))}


  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post2(Y, Y0, eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (is.na(rho.post2(Y, Y0, eXB, delta, C, LY, R))) #browser()
        if (rho.post2(Y, Y0, eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (rho.post2(Y, Y0, eXB, delta, C, LY, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (is.na(rho.post2(Y, Y0, eXB, delta, C, LY, R))) #browser()
        if (rho.post2(Y, Y0, eXB, delta, C, LY, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = rho.post2(Y, Y0, eXB, delta, C, LY, l1)

    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
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
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#


