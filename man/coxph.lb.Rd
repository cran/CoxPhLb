\name{coxph.lb}
\alias{coxph.lb}
\title{Fit Cox Model to Right-Censored Length-Biased Data}
\usage{
coxph.lb(y, a, delta, z, method = c("Bootstrap","EE"), boot.iter = 500,
      seed.n = round(runif(1,1,1e09)), print.out = FALSE, digits = 3L)
}
\arguments{
\item{y}{A vector of failure time.}

\item{a}{A vector of left-truncation time.}

\item{delta}{A vector of censoring indicator, 0=censored, 1=uncensored.}

\item{z}{A vector or a matrix of covariates.}

\item{method}{A character string specifying the method for variance estimation. The bootstrap resampling method ("Bootstrap") is used as the default. The estimating equation method ("EE")  uses the asymptotic variance estimation.}

\item{boot.iter}{The number of bootstrap iterations. Default is 500.}

\item{seed.n}{An integer specifying seed number.}

\item{print.out}{A logical variable. If TRUE, the function prints the outputs. Default is FALSE.}

\item{digits}{An integer controlling the number of digits to print.}
}
\value{
A list containing the following components:
\item{coef}{The vector of coefficients.}
\item{variance}{The variance of the coefficients.}
\item{std.err}{The standard error of the coefficients.}
\item{z.score}{z scores for the coefficients.}
\item{p.value}{p-values for the coefficients.}
\item{lower.95}{Lower 95\% confidence intervals of the coefficients.}
\item{upper.95}{Upper 95\% confidence intervals of the coefficients.}
\item{method}{The approach used to obtain the standard error of the coefficients.}
}
\description{
Fits a Cox model to right-censored length-biased data.
}
\details{
This function uses the weighted estimating equation proposed by Qin and Shen (2010). It returns coefficient estimates and the corresponding variance estimates based on either the asymptotic variance or the bootstrap resampling method. It also tests the null hypothesis that the coefficients are equal to 0.
}
\examples{
\dontrun{
# Fit a Cox model using model based variance estimation
fit.ee <- coxph.lb(ExampleData1$y, ExampleData1$a, ExampleData1$delta,
          cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), method="EE", print.out = TRUE)
print(fit.ee)			# display the results

# Fit a Cox model using bootstrap resampling method
fit.bs <- coxph.lb(ExampleData1$y, ExampleData1$a, ExampleData1$delta,
          cbind(x1=ExampleData1$x1, x2=ExampleData1$x2), method="Bootstrap",
          seed.n=1234, print.out=TRUE)
print(fit.bs)			# display the results
}
}
\references{
Qin J. and Shen Y. (2010). Statistical Methods for Analyzing Right-Censored Length-Biased Data under Cox Model. \emph{Biometrics} 66(2), 382-392.
}
\seealso{
\code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test}}, \code{\link{station.test.plot}}
}
