# Internal helper functions (likelihood, JointCount, data_plots)
llFun <- function(est,
                  Y,
                  Y0,
                  C,
                  X,
                  Z,
                  W,
                  V,
                  data,
                  form){		#Note the extra variable Y0 passed to the time varying MF Weibull

  n     <- nrow(data)
  llik  <- matrix(0, nrow = n, ncol = 1)
  gamma <- est[1:ncol(Z)]
  beta  <- est[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  p     <- est[length(est)]
  XB    <- X %*% beta
  ZG    <- Z %*% gamma
  phi   <- exp(-ZG + V)/(1+exp(-ZG + V))
  eXB   <- exp(-XB + W)
  if(form == "loglog"){
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*(1+((eXB*Y0)^(p)))/(1+exp(-(eXB*Y)^(p)))^2)+(1-C)*(log(phi+(1-phi)*(1+((eXB*Y0)^p)))^2/(1+((eXB*Y)^p)))^2)
  } else {
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*exp(-(eXB*Y))^p/exp(-(eXB*Y0))^p))+(1-C)*(log(phi+(1-phi)*((exp(-eXB*Y))^p)/((exp(-eXB*Y0))^p)))
  }
  one   <- nrow(llik)
  llik  <- subset(llik, is.finite(llik))
  two   <- nrow(llik)
  llik  <- subset(llik, llik[,1] > -1000)
  three <- nrow(llik)
  llik  <- -1*sum(llik)
  list(llik = llik, one = one, two = two, three = three)

}

# @title rllFun
# @description Log-likelihood function for pooled split population models
#
# @param est parameter estimates from the model
# @param Y the time (duration) dependent variable for the survival stage (t)
# @param Y0 the elapsed time since inception until the beginning of time period (t-1)
# @param C immunity indicator
# @param X covariates for betas
# @param Z covariates for gammas
# @param data data that contains the variables of interest
#
# @return Log-likelihood of the poold slit population survival model
#

rllFun <- function(est,
                   Y,
                   Y0,
                   C,
                   X,
                   Z,
                   data,
                   form){		#Note the extra variable Y0 passed to the time varying MF Weibull

  n     <- nrow(data)
  llik  <- matrix(0, nrow = n, ncol = 1)
  gamma <- est[1:ncol(Z)]
  beta  <- est[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  p     <- est[length(est)]
  XB    <- X %*% beta
  ZG    <- Z %*% gamma
  phi   <- exp(-ZG )/(1+exp(-ZG ))
  eXB   <- exp(-XB )
  if(form == "loglog"){
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*(1+((eXB*Y0)^(p)))/(1+exp(-(eXB*Y)^(p)))^2)+(1-C)*(log(phi+(1-phi)*(1+((eXB*Y0)^p)))^2/(1+((eXB*Y)^p)))^2)
  } else {
    llik <- C*(log((1-phi)*eXB*p*((eXB*Y)^(p-1))*exp(-(eXB*Y))^p/exp(-(eXB*Y0))^p))+(1-C)*(log(phi+(1-phi)*((exp(-eXB*Y))^p)/((exp(-eXB*Y0))^p)))
  }
  one   <- nrow(llik)
  llik  <- subset(llik, is.finite(llik))
  two   <- nrow(llik)
  llik  <- subset(llik, llik[,1] > -1000)
  three <- nrow(llik)
  llik  <- -1*sum(llik)
  list(llik = llik, one = one, two = two, three = three)

}


JointCount<- function(X, W) {
  PB<- mean (X)
  PW<- 1- mean(X)
  PBB<- PB*PB
  PWW<- PW*PW
  PBW<- 2*PB*(1-PB)
  n<-length(X)
  nB<-sum(X)
  EBB<- 0.5*sum(W)*PBB
  EWW<- 0.5*sum(W)*PWW
  EBW<- 0.5*sum(W)*PBW
  s.2j<-rep(NA,ncol(W))
  s.2i<-rep(NA,nrow(W))
  s.3j<-rep(NA,ncol(W))
  s.3i<-rep(NA,nrow(W))
  s_2<-matrix(NA,nrow(W),ncol(W))
  bb<-matrix(NA,nrow(W),ncol(W))
  ww<-matrix(NA,nrow(W),ncol(W))
  bw<-matrix(NA,nrow(W),ncol(W))
  for (i in 1:nrow(W)){
    for (j in 1:ncol(W)){
      s_2[i,j]<-(W[i,j]+W[j,i])^2
      bb[i,j]<-W[i,j]*X[i]*X[j]
      ww[i,j]<-W[i,j]*(1-X[i])*(1-X[j])
      bw[i,j]<-W[i,j]*((X[i]-X[j])^2)
    }
    s.3j[i]<-(sum(W[i,])+ sum(W[,i]))^2
  }
  BB<-0.5*sum(bb)
  WW<-0.5*sum(ww)
  BW<-0.5*sum(bw)
  s2<-sum(s_2)*0.5
  s3<-sum(s.3j)
  s1<-sum(W)
  varBW<-0.25*(((2*s2*nB*(n-nB))/(n*(n-1)))+(((s3-s1)*nB*(n-nB))/(n*(n-1)))+((4*(s1^2+s2-s3)*nB*(nB-1)*(n-nB)*(n-nB-1))/(n*(n-1)*(n-2)*(n-3)) ))-(EBW^2)
  ZBW<-(BW-EBW)/sqrt(varBW)
  out<-data.frame(BW,EBW,sqrt(varBW),ZBW)
  colnames(out)<-c("Obs","Exp","sd","Z")
  return(out)
}

data_plots <- function(data, var_id = character(), var_time = character(), n = 0, threshold){

  dat  <- lapply(split(data, data[,var_time]), function(x){x[!duplicated(x[, var_id]), ]})
  dat  <- dat[sapply(dat, nrow) >= n]
  mats <- lapply(dat, function(x){BayesSPsurv::spatial_SA(x, var_ccode = var_id, threshold = threshold)$matrixA})
  w2 <- list()
  for(i in 1:length(dat)){w2[[i]] <- dat[[i]][order(dat[[i]][, var_id]),]}
  list(mats = mats, w2 = w2)

}





