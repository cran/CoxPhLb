#' @title Test the Proportional Hazards Assumption of Cox Model with Right-Censored Length-Biased Data
#' @description Tests the proportional hazards assumption for a Cox model fit (coxphlb).
#' @import survival
#' @usage coxphlb.phtest(fit, data, spec.p = NULL, n.sim = 1000,
#' seed.n = round(runif(1,1,1e09)), digits = 3L)
#' @param fit The result of fitting a Cox model, using the \code{coxphlb} function.
#' @param data A data frame containing the variables in the model.
#' @param spec.p An integer specifying which covariate to be tested. Default is  NULL. If NULL, global test is conducted. If specified, the per-variable test is conducted.
#' @param n.sim The number of resampling. Default is 1000.
#' @param seed.n An integer specifying seed number.
#' @param digits An integer controlling the number of digits to print.
#' @return A list containing the following components:
#' \item{p.value}{A p-value.}
#' \item{}{The list is returned as an object of the \code{coxphlb.phtest} class. Objects of this class have methods for the function \code{print}. The object also contains the following: \code{spec.p}; \code{n.sim}; \code{stat.mat.t}, the test statistic; \code{sim.mat.t}, samples from the null distribution; \code{yy}, the observed ordered failure times; \code{varnames}, the variable that is tested; \code{result}, the table output.}
#' @details The proportional hazards assumption is checked by constructing test statistics based on asymptotically mean-zero processes. The asymptotic distribution of the test statistics is approximated via resampling. This function computes the p-value by comparing the test statistics with n.sim number of resamples. If the p-value is small (e.g., <0.05), it is likely that the assumption is violated. The test can be done either per variable or globally. The global test checks if the proportional hazards assumption is valid for the overall covariates.
#' @seealso \code{\link{coxphlb}}, \code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest.plot}}
#' @references Lee, C.H., Ning, J., and Shen, Y. Model diagnostics for proportional hazards model with length-biased data. \emph{Lifetime Data Analysis} 25(1), 79-96.
#' @examples
#' \dontrun{
#' # Fit a Cox model
#' fit.ee <- coxphlb(Surv(a, y, delta) ~ x1 + x2, data = ExampleData1,
#'                  method = "EE")
#'
#' # Check the Proportional Hazards Assumption
#' ptest1 <- coxphlb.phtest(fit.ee, data = ExampleData1, spec.p = 2,
#'                         seed.n = 1234)
#' print(ptest1)			# display the results
#'
#' # Run a Global Test
#' ptest2 <- coxphlb.phtest(fit.ee, data = ExampleData1, spec.p = NULL,
#'                         seed.n=1234)
#' print(ptest2)			# display the results
#' }
#' @export

coxphlb.phtest <- function(fit, data, spec.p = NULL, n.sim = 1000, seed.n = round(runif(1,1,1e09)), digits = 3L) {
  if (!inherits(fit, 'coxphlb')) stop ("Argument must be the result of coxphlb")

  ## construct data structure
  if (!is.data.frame(data)) stop ("The data must be a data frame!")
  varnames = fit$varnames
  nobservation = length(data[[varnames[1]]])
  rawdata = data
  data = cbind(rawdata[[varnames[2]]],rawdata[[varnames[3]]],rawdata[[varnames[2]]]-rawdata[[varnames[1]]])
  for (j in 4:length(varnames)) data = cbind(data,rawdata[[varnames[j]]])

  data = as.matrix(data)
  ## sort data by failure time
  data = data[order(data[,1]),]
  n = nrow(data)
  p = ncol(data)-3

  ### construct a matrix among failured times
  fdata = data[data[,2]==1,]
  yy = fdata[,1]
  m = length(yy)
  zz = as.matrix(fdata[,-c(1:3)])

  weight.fun = function(data,fdata) {
    #data: original data sorted by failure times
    #fdata: uncensored data sorted by failure times
    if (!sum(!data[,2])) wc = data[,1]
    else {
      yf = fdata[,1]
      mf = length(yf)
      cen.surv = survfit(Surv(data[,3],1-data[,2])~1) #survival function of residual censoring
      cen.time = c(0, summary(cen.surv)$time) #ordered residual censoring times
      cen.prob = c(1, summary(cen.surv)$surv) #survival probability at residual censoring times

      wc = rep(1,m) #integral of Sc
      cen.surv1 = summary(cen.surv,yy)
      Sc = cen.surv1$surv
      if (m>length(Sc)) {Sc[(length(Sc)+1):m]=min(cen.prob)}

      for (i in 1:m) {
        x0 = c(cen.time[cen.time<yy[i]],yy[i])
        x1 = c(0, cen.time[cen.time<yy[i]])
        wc[i] = sum(c(cen.prob[cen.time<yy[i]],Sc[i])*(x0-x1))
      }
    }
    return(list(wc=wc,w.hat=1/wc))
  }

  weight.out = weight.fun(data,fdata)
  wc = weight.out$wc
  w.hat = weight.out$w.hat
  beta.hat = as.vector(fit$coefficients)

  ##subcomponents
  exp.bz = exp(as.numeric(zz%*%as.matrix(beta.hat)))
  S0 = wc*rev(cumsum(rev(w.hat*exp.bz)))
  S0.1 = 1/S0
  Lam.0 = cumsum(wc*S0.1)

  S1 = wc*apply(as.matrix(w.hat*zz*exp.bz),2,function(x)return(rev(cumsum(rev(x)))))
  E1 = S1/S0
  S2 = array(NA,c(m,p,p))
  for (j in 1:p) S2[,,j] = wc*apply(as.matrix(w.hat*zz[,j]*zz*exp.bz),2,function(x)return(rev(cumsum(rev(x)))))
  E1.2 = t(E1)%*%E1
  E2 = matrix(NA,nrow=p,ncol=p)
  for (j in 1:p){
    for (k in 1:p){
      E2[j,k] =sum(S2[,j,k]/S0)
    }
  }
  Gam = (E2-E1.2)/n

  U0.1 = zz-as.vector(E1)
  U0.2 = zz*w.hat*exp.bz*Lam.0-w.hat*exp.bz*apply(wc*E1*S0.1,2,cumsum)
  U0 = U0.1-U0.2


  ###-------general process given data
  stat.mat.t = apply(as.matrix(zz-as.vector(E1)),2,cumsum)/sqrt(n)
  if (is.null(spec.p)) {
    ##global test statistic for proportional hazards assumption
    tstat.t = max(apply(abs(stat.mat.t),1,sum)) } else {
      ##test statistic of specified covariate
      tstat.t = max(abs(stat.mat.t[,spec.p])) }


  ##null distribution by resampling
  b1 = b2 = b3 = array(0,c(m,p,m))
  H1t = array(0,c(p,p,m))
  H2t = array(0,c(p,p,m))
  for (j in 1:m){
    b1[1:j,,j] =zz[1:j,]*(1-w.hat[1:j]*exp.bz[1:j]*Lam.0[1:j])
    if (j<m) b1[(j+1):m,,j] = zz[(j+1):m,]*(-w.hat[(j+1):m]*exp.bz[(j+1):m]*Lam.0[j])
    sub.Lam.0 = c(Lam.0[1:j],rep(Lam.0[j],m-j))
    H1t[,,j] = t(zz)%*%(zz*w.hat*exp.bz*sub.Lam.0)/n
    if (j==1) H2t[,,j] = (as.matrix(S1)[j,]/n)%*%t(as.matrix(S1/(S0)^2)[j,])
    if (j>1) H2t[,,j] = t(as.matrix(S1)[1:j,]/n)%*%(as.matrix(S1/(S0)^2)[1:j,])
    b2[,,j] = t((H1t[,,j]-H2t[,,j])%*%solve(Gam)%*%t(U0))
  }
  b3[1,,1] = as.matrix(S1)[1,]/S0[1]-as.matrix(S1)[1,]*w.hat[1]*exp.bz[1]*wc[1]*S0.1[1]/S0[1]
  b3[2:m,,1] = -matrix(rep(as.matrix(S1)[1,]*wc[1]*S0.1[1]/S0[1],m-1),nrow=m-1,byrow=TRUE)*w.hat[2:m]*exp.bz[2:m]
  for (j in 2:m){
    b3[1:j,,j] = as.matrix(S1)[1:j,]/S0[1:j]-apply(as.matrix(as.matrix(S1)[1:j,]*wc[1:j]*S0.1[1:j]/S0[1:j]),2,cumsum)*w.hat[1:j]*exp.bz[1:j]
    if (j<m) b3[(j+1):m,,j] = -matrix(rep(apply(as.matrix(as.matrix(S1)[1:j,]*wc[1:j]*S0.1[1:j]/S0[1:j]),2,sum),m-j),nrow=m-j,byrow=TRUE)*w.hat[(j+1):m]*exp.bz[(j+1):m]
  }
  b0 = b1-b2-b3

  set.seed(seed.n)
  sim.mat.t = array(NA,c(n.sim,m,p))
  for (j in 1:n.sim){
    rv = rnorm(m,0,1)
    for (k in 1:m){
      sim.mat.t[j,k,] = apply(as.matrix(b0[,,k])*rv,2,sum)/sqrt(n)
    }
  }
  crit.t = rep(NA,n.sim)
  for (j in 1:n.sim){
    if (is.null(spec.p)) {
      ##global test
      crit.t[j] = max(apply(abs(as.matrix(sim.mat.t[j,,])),1,sum)) } else {
        crit.t[j] = max(abs(sim.mat.t[j,,spec.p])) }
  }

  p.val.t = sum(crit.t>tstat.t)/n.sim

  tab = cbind(p.value = format.pval(p.val.t, digits, eps = 1e-03))
  varnames = fit$varnames[-(1:3)]
  if (is.null(spec.p)) {rownames(tab) = "GLOBAL"} else {
  rownames(tab) = varnames[spec.p]}
  tab = as.table(tab)
  print(tab)
  out=list(n.sim = n.sim, stat.mat.t = stat.mat.t, sim.mat.t = sim.mat.t, spec.p = spec.p, varnames = varnames[spec.p], yy = yy,
           p.value = c(p.val.t), result = tab)
  class(out) <- "coxphlb.phtest"
  invisible(out)
}
