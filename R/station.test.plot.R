#' @title Graphical Test of the Stationarity Assumption
#' @description Returns a plot of two Kaplan-Meier curves for forward recurrence time and backward recurrence time.
#' @import survival
#' @usage station.test.plot(a, v, delta)
#' @param a A vector of backward recurrence time (i.e., left-truncation time).
#' @param v A vector of forward recurrence time (i.e., failure time minus left-truncation time).
#' @param delta A vector of censoring indicator, 0=censored, 1=uncensored.
#' @return
#' \item{A plot of Kaplan-Meier curves for backward and forward recurrence time.}{ }
#' @details The stationarity assumption can be checked by comparing the Kaplan-Meier curves. More overlap of the two survival curves suggests stronger evidence of stationarity.
#' @seealso \code{\link{coxph.lb}}, \code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test}}
#' @references Asgharian, M., Wolfson, D. B., and Zhang, X. (2006). Checking stationarity of the incidence rate using prevalent cohort survival data. \emph{Statistics in medicine}, 25(10), 1751-1767.
#' @examples
#' # Check the Stationarity Assumption Graphically
#' station.test.plot(ExampleData1$a, ExampleData1$y-ExampleData1$a, ExampleData1$delta)			# plot curves
#'
#' station.test.plot(ExampleData2$a, ExampleData2$y-ExampleData2$a, ExampleData2$delta)			# plot curves
#' @export

station.test.plot <- function(a, v, delta)
{
  len<-length(delta)

  plot(survfit(Surv(v,delta)~1), conf.int = FALSE, lty=1, ylab="Estimated survival function")
  lines(survfit(Surv(a,rep(1,len))~1), conf.int = FALSE, lty=2)
  legend("topright", c("forward recurrence time","backward recurrence time"), lty=c(1,2))
}
