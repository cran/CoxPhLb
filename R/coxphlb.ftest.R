#' @title Test the Functional Form of Covariates in Cox Model with Right-Censored Length-Biased Data
#' @description Tests the functional form of covariates assumed for a Cox model fit (coxph.lb).
#' @import survival
#' @usage coxphlb.ftest(y, a, delta, z, fit, spec.p = c(1), n.sim = 1000, z0 = NULL, plot.p = FALSE, n.plot = 20, seed.n = round(runif(1,1,1e09)), print.out = FALSE, digits = 3L)
#' @param y A vector of failure time.
#' @param a A vector of left-truncation time.
#' @param delta A vector of censoring indicator, 0=censored, 1=uncensored.
#' @param z A vector or a matrix of covariates.
#' @param fit A vector of coefficients.
#' @param spec.p An integer specifying which covariate to be tested. Default is 1. If set to 1, the first column of the covariate matrix is tested.
#' @param n.sim The number of resampling. Default is 1000.
#' @param z0 A vector of grid points to use for the specified covariate. The default is a vector of 100 equally distributed numeric values within the range of the specified covariate.
#' @param plot.p A logical variable. If TRUE, the function is followed by a plot. Default is FALSE.
#' @param n.plot The number of randomly selected realizations. Default is 20.
#' @param seed.n An integer specifying seed number.
#' @param print.out A logical variable. If TRUE, the function prints the outputs. Default is FALSE.
#' @param digits An integer controlling the number of digits to print.
#' @return
#' \item{p.value}{A p-value.}
#' @details The functional form of a continuous covariate is checked by constructing test statistics based on asymptotically mean-zero processes. The asymptotic distribution of the test statistics is approximated via resampling. This function computes the p-value by comparing the test statistics with n.sim number of resamples. If the p-value is small (e.g., <0.05), it is likely that the assumption is violated. It may return a plot with the test statistics in a black line and 20 processes randomly sampled from the pool of resamples in grey lines. When the test statistics lie within the randomly sampled lines, it suggests that the model assumption is valid. The test should be done per variable for continuous covariates.
#' @seealso \code{\link{coxph.lb}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test}}, \code{\link{station.test.plot}}
#' @references Lee, C.H., Ning, J., and Shen, Y. Model diagnostics for proportional hazards model with length-biased data. \emph{Manuscript available upon request}.
#' @examples
#' # Fit a Cox model
#' fit.ee <- coxph.lb(ExampleData1$y, ExampleData1$a, ExampleData1$delta,  
#'                      cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), method="EE")
#'
#' # Check the Functional Form of the Cox Model
#' ftest <- coxphlb.ftest(ExampleData1$y, ExampleData1$a, ExampleData1$delta, 
#'                          cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), fit=as.vector(fit.ee$coef),
#'                          spec.p=c(2), plot.p=TRUE, n.plot=50, seed.n=1234, print.out=TRUE)
#' print(ftest)			# display the results
#' @export

coxphlb.ftest <- function(y, a, delta, z, fit, spec.p = c(1), n.sim = 1000, z0 = NULL, plot.p = FALSE, n.plot = 20, seed.n = round(runif(1,1,1e09)), print.out = FALSE, digits = 3L) {
  ## construct data structure
  data = cbind(y, delta, y-a, z)
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

  if (is.null(z0)) z0 = seq(min(zz[,spec.p]),max(zz[,spec.p]),length=100)
  n.z0 = length(z0)

  weight.fun = function(data,fdata) {
  #data: original data sorted by failure times
  #fdata: uncensored data sorted by failure times
  if (!sum(!delta)) wc = data[,1]
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
  beta.hat = fit

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
  ##test statistic for functional form
  z.mat = matrix(rep(as.matrix(zz)[,spec.p],n.z0),ncol=n.z0)
  z.ind = t(t(z.mat)<=z0)
  Mi.hat = 1 - w.hat*exp.bz*Lam.0
  stat.mat.z = apply(z.ind*Mi.hat/sqrt(n),2,sum)
  tstat.z = max(abs(stat.mat.z))


  ##null distribution by resampling
  a1 = z.ind*Mi.hat
  H1z = matrix(NA,nrow=n.z0,ncol=p)
  for (j in 1:p) H1z[,j] = apply(z.ind*w.hat*exp.bz*as.matrix(zz)[,j]*Lam.0,2,sum)/n
  Sz0 = wc*apply(z.ind*w.hat*exp.bz,2,function(x)return(rev(cumsum(rev(x)))))/n
  #dAt = (S1/(S0)^2)
  H2z = matrix(NA,nrow=n.z0,ncol=p)
  for (j in 1:p) H2z[,j] =apply(Sz0*S1[,j]/(S0)^2,2,sum)
  a2 = t((H1z-H2z)%*%solve(Gam)%*%t(U0))
  a3 = Sz0/(S0/n)-w.hat*exp.bz*apply(wc*Sz0*S0.1/(S0/n),2,cumsum)
  a0 = a1-a2-a3

  set.seed(seed.n)
  sim.mat.z = matrix(NA,nrow=n.sim,ncol=n.z0)
  for (j in 1:n.sim){
    rv = rnorm(m,0,1)
    sim.mat.z[j,] = apply(a0*rv,2,sum)/sqrt(n)
  }
  crit.z = apply(sim.mat.z,1,function(x)return(max(abs(x))))

  p.val.z = sum(crit.z>tstat.z)/n.sim

  if (plot.p == TRUE) {
  samp.ind = sample(1:n.sim,n.plot,replace=FALSE)
  plot(z0,stat.mat.z,type="l",ylim=c(-max(abs(sim.mat.z)),max(abs(sim.mat.z))),xlab="z",ylab="Partial-Mean Processes")
  for (j in samp.ind) {
    lines(z0,sim.mat.z[j,],col="grey")
  }
  lines(z0,stat.mat.z,lwd=1.5)
  }

  if (print.out==TRUE) {
  cat("\tp.value = ", format.pval(p.val.z, digits, eps = 1e-03), "\n", sep = "")}
  out=list(p.value=c(p.val.z))
}

