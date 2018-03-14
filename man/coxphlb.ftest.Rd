\name{coxphlb.ftest}
\alias{coxphlb.ftest}
\title{Test the Functional Form of Covariates in Cox Model with Right-Censored Length-Biased Data}
\usage{
coxphlb.ftest(y, a, delta, z, fit, spec.p = c(1), n.sim = 1000, z0 = NULL,
      plot.p = FALSE, n.plot = 20, seed.n = round(runif(1,1,1e09)),
      print.out = FALSE, digits = 3L)
}
\arguments{
\item{y}{A vector of failure time.}

\item{a}{A vector of left-truncation time.}

\item{delta}{A vector of censoring indicator, 0=censored, 1=uncensored.}

\item{z}{A vector or a matrix of covariates.}

\item{fit}{A vector of coefficients.}

\item{spec.p}{An integer specifying which covariate to be tested. Default is 1. If set to 1, the first column of the covariate matrix is tested.}

\item{n.sim}{The number of resampling. Default is 1000.}

\item{z0}{A vector of grid points to use for the specified covariate. The default is a vector of 100 equally distributed numeric values within the range of the specified covariate.}

\item{plot.p}{A logical variable. If TRUE, the function is followed by a plot. Default is FALSE.}

\item{n.plot}{The number of randomly selected realizations. Default is 20.}

\item{seed.n}{An integer specifying seed number.}

\item{print.out}{A logical variable. If TRUE, the function prints the outputs. Default is FALSE.}

\item{digits}{An integer controlling the number of digits to print.}
}
\value{
\item{p.value}{A p-value.}
}
\description{
Tests the functional form of covariates assumed for a Cox model fit (coxph.lb).
}
\details{
The functional form of a continuous covariate is checked by constructing test statistics based on asymptotically mean-zero processes. The asymptotic distribution of the test statistics is approximated via resampling. This function computes the p-value by comparing the test statistics with n.sim number of resamples. If the p-value is small (e.g., <0.05), it is likely that the assumption is violated. It may return a plot with the test statistics in a black line and 20 processes randomly sampled from the pool of resamples in grey lines. When the test statistics lie within the randomly sampled lines, it suggests that the model assumption is valid. The test should be done per variable for continuous covariates.
}
\examples{
# Fit a Cox model
fit.ee <- coxph.lb(ExampleData1$y, ExampleData1$a, ExampleData1$delta,
          cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), method="EE")

# Check the Functional Form of the Cox Model
ftest <- coxphlb.ftest(ExampleData1$y, ExampleData1$a, ExampleData1$delta,
         cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), fit=as.vector(fit.ee$coef),
         spec.p=c(2), plot.p=TRUE, n.plot=50, seed.n=1234, print.out=TRUE)
print(ftest)			# display the results
}
\references{
Lee, C.H., Ning, J., and Shen, Y. (2018). Model diagnostics for proportional hazards model with length-biased data. \emph{Lifetime Data Analysis} doi:10.1007/s10985-018-9422-y.
}
\seealso{
\code{\link{coxph.lb}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test}}, \code{\link{station.test.plot}}
}
