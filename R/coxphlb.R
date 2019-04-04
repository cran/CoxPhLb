#' @title Fit Cox Model to Right-Censored Length-Biased Data
#' @description Fits a Cox model to right-censored length-biased data.
#' @import survival
#' @usage coxphlb(formula, data, method = c("Bootstrap","EE"),
#' boot.iter = 500, seed.n = round(runif(1,1,1e09)), digits = 3L)
#' @param formula A formula object with the form \code{response ~ predictors}. For \code{response}, use \code{Surv} object.
#' @param data A data frame containing the variables in the model.
#' @param method A character string specifying the method for variance estimation. The bootstrap resampling method ("Bootstrap") is used as the default. The estimating equation method ("EE")  uses the asymptotic variance estimation.
#' @param boot.iter The number of bootstrap iterations. Default is 500.
#' @param seed.n An integer specifying seed number.
#' @param digits An integer controlling the number of digits to print.
#' @return A list containing the following components:
#' \item{coefficients}{The vector of coefficients.}
#' \item{var}{The variance matrix of the coefficients.}
#' \item{std.err}{The standard error of the coefficients.}
#' \item{z.score}{z scores for the coefficients.}
#' \item{p.value}{p-values for the coefficients.}
#' \item{lower.95}{Lower 95\% confidence intervals of the coefficients.}
#' \item{upper.95}{Upper 95\% confidence intervals of the coefficients.}
#' \item{method}{The approach used to obtain the standard error of the coefficients.}
#' \item{}{The list is returned as an object of the \code{coxphlb} class to represent a fitted proportional hazards model. Objects of this class have methods for the functions \code{coef}, \code{print}, \code{summary}, and \code{vcov}. The object also contains the following: \code{formula}; \code{varnames}, the variables used in the model; \code{result}, the table output.}
#' @details This function uses the weighted estimating equation proposed by Qin and Shen (2010). It returns coefficient estimates and the corresponding variance estimates based on either the asymptotic variance or the bootstrap resampling method. It also tests the null hypothesis that the coefficients are equal to 0.
#' @seealso \code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test}}, \code{\link{station.test.plot}}
#' @references Qin J. and Shen Y. (2010). Statistical Methods for Analyzing Right-Censored Length-Biased Data under Cox Model. \emph{Biometrics} 66(2), 382-392.
#' @examples
#' \dontrun{
#' # Fit a Cox model using model based variance estimation
#' fit.ee <- coxphlb(Surv(a, y, delta) ~ x1 + x2, data = ExampleData1,
#'                  method = "EE")
#' summary(fit.ee)			# display the results
#'
#' # Fit a Cox model using bootstrap resampling method
#' fit.bs <- coxphlb(Surv(a, y, delta) ~ x1 + x2, data = ExampleData1,
#'                  method = "Bootstrap", seed.n = 1234)
#' summary(fit.bs)			# display the results
#' }
#' @export

coxphlb <- function(formula, data, method = c("Bootstrap","EE"), boot.iter = 500, seed.n = round(runif(1,1,1e09)), digits = 3L) {
  ## construct data structure
  if (!is.data.frame(data)) stop ("The data must be a data frame!")
  varnames = all.vars(as.formula(formula))
  nobservation = length(data[[varnames[1]]])
  rawdata = data
  data = cbind(rawdata[[varnames[2]]],rawdata[[varnames[3]]],rawdata[[varnames[2]]]-rawdata[[varnames[1]]])
  for (j in 4:length(varnames)) data = cbind(data,rawdata[[varnames[j]]])

  data = as.matrix(data)
  ## sort data by failure time
  data = data[order(data[,1]),]
  n = nrow(data)
  covdim = ncol(data)-3

  ### construct a matrix among failured times
  fdata = data[data[,2]==1,]
  yy = fdata[,1]
  m = length(yy)

  weightfun <- function(data, fdata) {

    ## PH model solving from estimating equation II in Biometrics paper.
    ## weight:  int_0^y S_c(t) dt
    if (!sum(!data[,2])) what = data[,1]
    else {
      yf = fdata[,1]
      mf = length(yf)
      H00 = survfit(Surv(data[,3],1-data[,2])~1) ## using residual cen
      rescen = c(0, summary(H00)$time) ## ordered residual cen times
      cenprob = c(1, summary(H00)$surv) ## surv prob at res cen times

      ## predict H00 on all residual failures (sorted by original t)

      Sc = rep(1,mf)

      for (k in 1:mf) {
        Sc[k] = min(cenprob[rescen<=fdata[k,3]])
      }

      what = rep(1,mf)

      # use integral of S_c
      ## estimate KM survival curve for censoring variable
      ## predict surv prob of H00 on all failures (sorted)
      H10 = summary(H00,yf)
      if (length(H10$time)) {
        H1 = rep(0,mf)
        H1[1:length(H10$time)] = H10$surv
        if (mf > length(H10$time))
        {H1[(length(H10$time)+1):mf] = min(cenprob)}
      } else {H1 = rep(min(cenprob), mf)}

      rescen = rescen[-1]
      cenprob = cenprob[-1]
      for (i in 1:mf) {
        x0 = c(rescen[rescen < yf[i]], yf[i])
        x1 = c(0, rescen[rescen < yf[i]])

        what[i]= sum(c(cenprob[rescen < yf[i]], H1[i])*(x0-x1))
      }

    }

    what = (1/what)

    return(what)
  }

  what = weightfun(data,fdata)
  whatmat = diag(what)

  bootvar <- function(data,iter) {

    ## use bootstrap method to est variance for all EE II under PH model

    estmat=NULL

    for (k in 1:iter) {
      bsmat = data[sample(seq(1,n),n,replace=TRUE),]
      if (sum(bsmat[,2])) {
        bsmat = bsmat[order(bsmat[,1]),]
        bsfdata = bsmat[bsmat[,2]==1,]

        bsw = weightfun(bsmat,bsfdata)

        fit=coxph(Surv(bsfdata[,1],rep(1,nrow(bsfdata)))~bsfdata[,-c(1:3)]+offset(log(bsw)))

        estmat=rbind(estmat, c(fit$coef)) }
    }

    estmat = matrix(c(estmat),nrow=nrow(estmat),ncol=ncol(estmat))

    varest = apply(estmat,2,var,na.rm=TRUE)

    return(varest)
  }

  phfun <- function(par0,r){

    ## PH model solving from estimating equation II in Biometrics paper.
    ## weight:  int_0^y S_c(t) dt
    ## par0: covariate coefficient for two covariates (z1,z2)
    ## r: # of max iteration

    estfun = matrix(0,nrow=covdim,ncol=1)

    for (j in 1:(m-1)) {
      estfun = estfun+
        as.matrix(fdata[j,-c(1:3)]) - t(fdata[j:m,-c(1:3)])%*%(whatmat[j:m,j:m])%*%t(exp(par0%*%t(fdata[j:m,-c(1:3)])))/
        sum(whatmat[j:m,j:m]%*%t(exp(par0%*%t(fdata[j:m,-c(1:3)]))))	}

    estfun = estfun-r
    return(estfun)
  }

  # find the derivative for the estimating function
  dg <- function(par0){
    ## using numeric approximation with two covariates

    h <- 0.00000001

    te<-matrix(0,covdim,covdim)
    for(i in 1:covdim)
    {
      s1 <- par0
      s2 <- par0
      s1[i] <- s1[i]+h
      s2[i] <- s2[i]-h

      te[,i] <- (phfun(s1,0)-phfun(s2,0))/(2*h)
    }
    return(te)
  }

  eesol <- function(beta0,iter=1) {
    ## mimic JQ's sampl: modified Newtow-Raphson method
    ## initial value for par estimate = beta_0
    ## let maximum step of iteration
    #iter=30

    rot = beta0
    r = phfun(rot,0)
    step = r/iter

    for (i in 1:iter) { # this is another layer of iteration in addition to Newton-Raphson, may be unnecessary if function is smooth enough
      check = 1
      r = r-step
      while (check>0.000001) {
        ddt = det(dg(rot))
        if (identical(all.equal(ddt,0),TRUE))
        {rot = rep(NA,covdim);break} else
          tm = solve(dg(rot),phfun(rot,r))
        rot = rot-t(tm)

        check = sqrt(mean(tm^2))
      }
    }
    check2 = phfun(rot,0)
    list(rot=rot, check=check, fitback=check2)
  }

  eevar <- function(par0){
    wc = 1/what
    S0 = rev(cumsum(rev((what)*exp(par0%*%t(fdata[,-c(1:3)])))))*wc
    #lambda0=cumsum(1/S0)  # cumulated hazards
    lambda0 = 1/S0

    if (covdim==1) {
      S1 = wc*rev(cumsum(rev(fdata[,-c(1:3)]*as.vector((what)*exp(par0%*%t(fdata[,-c(1:3)]))))))

      S2 = wc*rev(cumsum(rev(fdata[,-c(1:3)]^2*as.vector((what)*exp(par0%*%t(fdata[,-c(1:3)]))))))

      temp=S2/S0 - (S1/S0)^2

      ## information matrix
      gamma = matrix(sum(temp),nrow=covdim,ncol=covdim)
    } else {
      S1 = wc*apply(fdata[,-c(1:3)]*as.vector((what)*exp(par0%*%t(fdata[,-c(1:3)]))),2,function(x){return(rev(cumsum(rev(x))))}) ## yy times to each corresponding row

      S2 = wc*apply(t(apply(fdata[,-c(1:3)],1,function(x){return(x%o%x)}))*as.vector((what)*exp(par0%*%t(fdata[,-c(1:3)]))),2,function(x){return(rev(cumsum(rev(x))))}) ## this should be an array

      temp=S2/S0 - t(apply((S1/S0),1,function(x){return(x%o%x)}))

      ## information matrix
      gamma = matrix(apply(temp,2,sum),nrow=covdim,ncol=covdim)
    }

    sigma = matrix(0,nrow=covdim,ncol=covdim)

    S0 = S0/wc

    for (j in 1:(m-1)) {
      temp = as.matrix(fdata[j,-c(1:3)]) - t(fdata[j:m,-c(1:3)])%*%(whatmat[j:m,j:m])%*%t(exp(par0%*%t(fdata[j:m,-c(1:3)])))/sum(whatmat[j:m,j:m]%*%t(exp(par0%*%t(fdata[j:m,-c(1:3)]))))

      dMj = 1 - whatmat[j,j]*exp(par0%*%fdata[j,-c(1:3)])/S0[j]
      temp = as.vector(dMj)*temp
      sigma = sigma + as.vector(dMj)*temp%*%t(temp)/n
      #sigma=sigma+ temp%*%t(temp)/n
    }

    estvar = solve(gamma)%*%sigma%*%solve(gamma)*n

    return(list(estvar=estvar, inf=solve(gamma)))
  }

  par_est = eesol(rep(0,covdim),iter=1)$rot

  if(length(method)>1) method="Bootstrap"

  if (method=="Bootstrap") {set.seed(seed.n); varcov = NA; estvar = bootvar(data,iter=boot.iter)} else if (method=="EE") {varcov = eevar(par_est)$estvar; estvar = diag(varcov)}

  std.err = sqrt(estvar)
  z.score = par_est/std.err
  p.value = (1-pnorm(abs(z.score)))*2
  lower.95 = par_est - 1.96*std.err
  upper.95 = par_est + 1.96*std.err


  tab = cbind(coef = round(as.vector(par_est),digits),
              variance = round(estvar,digits),
              std.err = round(std.err,digits),
              z.score = round(as.vector(z.score),digits-1),
              p.value = format.pval(as.vector(p.value), digits, eps = 1e-03),
              lower.95 = round(as.vector(lower.95),digits),
              upper.95 = round(as.vector(upper.95),digits))
  rownames(tab) = varnames[-(1:3)]
  tab = as.table(tab)
  form = paste(formula)
    cat("Call:\n","coxphlb(formula = ",form[2]," ",form[1]," ",form[3],", method = ", method,")\n", sep = "")
    print(tab)

  out=list(formula = formula, coefficients = par_est, var = varcov, std.err = std.err, z.score = z.score, p.value = p.value,
           lower.95 = lower.95, upper.95 = upper.95, method = method, varnames = varnames, result=tab)
  class(out) <- "coxphlb"
  invisible(out)
}
