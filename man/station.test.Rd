\name{station.test}
\alias{station.test}
\title{Test the Stationarity Assumption}
\usage{
station.test(a, v, delta, print.out = FALSE, digits = 3L)
}
\arguments{
\item{a}{A vector of backward recurrence time (i.e., left-truncation time).}

\item{v}{A vector of forward recurrence time (i.e., failure time minus left-truncation time).}

\item{delta}{A vector of censoring indicator, 0=censored, 1=uncensored.}

\item{print.out}{A logical variable. If TRUE, the function prints the outputs. Default is FALSE.}

\item{digits}{An integer controlling the number of digits to print.}
}
\value{
A list containing the following components:
\item{test.statistic}{A test statistic.}
\item{p.value}{A p-value based on two-sided test.}
}
\description{
Tests the null hypothesis that the incidence process is stationary.
}
\details{
The stationarity assumption is checked by computing the test statistic and the corresponding p-value. A large p-value suggests strong evidence of stationarity. When the p-value is small (e.g., <0.05), it is likely that the stationarity assumption is violated.
}
\examples{
# Check the Stationarity Assumption
stest1 <- station.test(ExampleData1$a, ExampleData1$y-ExampleData1$a, ExampleData1$delta,
         print.out=TRUE)
print(stest1) 			# display the results

stest2 <- station.test(ExampleData2$a, ExampleData2$y-ExampleData2$a, ExampleData2$delta,
         print.out=TRUE)
print(stest2)   		# display the results
}
\references{
Addona, V. and Wolfson, D. B. (2006). A formal test for the stationarity of the incidence rate using data from a prevalent cohort study with follow-up. \emph{Lifetime data analysis}, 12(3), 267-284.
}
\seealso{
\code{\link{coxph.lb}}, \code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test.plot}}
}
