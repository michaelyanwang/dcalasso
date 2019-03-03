#' Generate simulation data to test adaptive lasso
#'
#' \code{SIM.FUN} generates continuous-time survival response data that are associated with design matrix.
#' The design matrix comes from a correlated multivariate normal. The default signals (beta0) are sparse.
#'
#' @param nn sample size
#' @param p.x number of covariates
#' @param cor correlation of covariates
#' @param family the family of response data taking c('binary','count','Cox')
#' @param beta0 the coefficients for the design, including intercept
#' @return For survival data, it returns a matrix with the first column U, second column delta (0,1), and rest = design matrix.
#' @author Yan Wang, Tianxi, Chuan Hong
#' @export
#' @examples SIM.FUN(nn = 1e6, p.x = 50, family = 'binary')
SIM.FUN = function(nn,p.x=50, cor = 0.2, family = c('binary','count','Cox'),beta0=NULL){
  if(is.null(beta0)){
    bb=c(rep(0.8,3),rep(0.4,3),rep(0.2,3))
    beta0=c(1,bb,rep(0,p.x-length(bb)))
  }
  Sig0.X=cor+(1-cor)*diag(p.x)

  xx = mvrnorm(nn,mu=rep(0,p.x),Sigma=Sig0.X);
  if (family=='Cox'){
    ## Weibull ##
    h = c(cbind(xx)%*%beta0[-1])
    t = rweibull(nn, shape = 2, scale = (0.5*exp(h))^(-0.5))
    c = rexp(nn, rate = exp(0.5))
    delta = t<=c
    u = ifelse(t<=c, t, c)
    return(cbind(u,delta,xx))
  }else{
    stop('Unrecognized family argument.')
  }
}


#' Generate simulation data to test time-dependent Cox model with adaptive lasso
#'
#' \code{SIM.FUN.TVC} generates time-dependent survival response with four time-intervals 0-1, 1-2, 2-3, 3-4 for each subject data. All subjects are administratively censored at 4, if T>4. T comes from a Weibull distribution with shape of 2.
#' The design matrix comes from a correlated multivariate normal. The default signals (beta0) are sparse.
#'
#' @param p.ti number of time-invariant covariates
#' @param p.tv number of time-varying covariates
#' @param n.subject number of subjects
#' @param cor correlation between time-varying and each interval's time-varying covairates
#' @param beta.ti coefficient for time-invariant covariates
#' @param beta.tv coefficients for time-varying covariates
#' @return a matrix with the first column starting time, second column ending time, third column event (0,1), and rest = design matrix + ID for subject.
#' @author Yan Wang, Tianxi Cai, Chuan Hong
#' @references Section 3.3 in Austin, P.C., 2012. Generating survival times to simulate Cox proportional hazards models with time‐varying covariates. Statistics in medicine, 31(29), pp.3946-3958.
#' @export
#' @examples SIM.FUN.tvc()
SIM.FUN.TVC = function(p.ti=50, p.tv=50, n.subject = 1e6, cor = 0.2, beta0.ti = NULL, beta0.tv = NULL){

  if ((is.null(beta0.ti))&(is.null(beta0.tv))&(p.ti>9)&(p.tv>9)){
    beta0.ti = c(c(rep(0.08,3),rep(0.04,3),rep(0.02,3)),rep(0,p.ti-9))
    beta0.tv = c(c(rep(0.08,3),rep(0.04,3),rep(0.02,3)),rep(0,p.tv-9))
  }

  lambda = 0.05
  # Consider intervals D1 = [0,t1), D2 = [t1,t2), D3 = [t2,t3), D4 = [t3,+infty)
  interval.max = 4
  nu = 2; t = 0:(interval.max-1)

  # Coefficient vector
  beta0 = c(beta0.ti, rep(beta0.tv, each = interval.max))

  # Design matrix has column number TI+TV*MaxInterval
  p.x = p.ti + p.tv*interval.max
  Sig0.X=cor+(1-cor)*diag(p.x)

  xx0 = mvrnorm(n.subject, mu=rep(0,p.x), Sigma=Sig0.X)
  xx0 = cbind(xx0, 1:n.subject)
  gc()

  u = runif(n.subject, 0, 1)
  # Calculate R1-R4, H1-H4
  H = R = Hinv = flag = matrix(0, nrow = n.subject, ncol = interval.max)
  for (i in 1:interval.max){
    ind = c(1:p.ti, p.ti+(0:(p.tv-1))*interval.max+i-1)
    H[,i] = lambda*exp(xx0[,ind]%*%c(beta0.ti,beta0.tv))
    if (i>1){
      R[,i] = R[,(i-1)] + H[,(i-1)]*(t[i]^nu-t[i-1]^nu)
      flag[,i] = -log(u)<R[,i]
    }
    Hinv[,i] = ((-log(u)-R[,i])/H[,i]+t[i]^nu)^(1/nu)
  }
  rm(list=c('H','R'));gc()
  h.ind = interval.max-rowSums(flag)
  rm(list=c('flag'));gc()
  for (i in 1:interval.max){
    ind.col = c(1:p.ti, p.ti+(0:(p.tv-1))*interval.max+i-1, dim(xx0)[2])
    ind.elig = h.ind>=i
    t0 = t[i]
    t1 = Hinv[ind.elig,i]
    t1[t1>t[i]+1] = t[i]+1
    status = t1==Hinv[ind.elig,i]
    if (i==1){
      xx = cbind(t0,t1,status,xx0[ind.elig, ind.col])
    }else{
      xx = rbind(xx, cbind(t0,t1,status,xx0[ind.elig, ind.col]))
    }
  }
  xx[which(xx[,2]-xx[,1]<=4*sqrt(.Machine$double.eps)),2] = xx[which(xx[,2]-xx[,1]<=4*sqrt(.Machine$double.eps)),1] + 4*sqrt(.Machine$double.eps)
  rm(list=c('xx0','Hinv','t0','ind.elig','status','u','h.ind'));gc()
  return(xx)
}
