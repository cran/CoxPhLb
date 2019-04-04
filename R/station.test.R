#' @title Test the Stationarity Assumption
#' @description Tests the null hypothesis that the incidence process is stationary.
#' @import survival
#' @usage station.test(a, v, delta, digits = 3L)
#' @param a A vector of backward recurrence time (i.e., left-truncation time).
#' @param v A vector of forward recurrence time (i.e., failure time minus left-truncation time).
#' @param delta A vector of censoring indicator, 0=censored, 1=uncensored.
#' @param digits An integer controlling the number of digits to print.
#' @return A list containing the following components:
#' \item{test.statistic}{A test statistic.}
#' \item{p.value}{A p-value based on two-sided test.}
#' \item{}{The list is returned as an object of the \code{station.test} class. Objects of this class have methods for the function \code{print}. The object also contains the following: \code{result}, the table output.}
#' @details The stationarity assumption is checked by computing the test statistic and the corresponding p-value. A large p-value suggests strong evidence of stationarity. When the p-value is small (e.g., <0.05), it is likely that the stationarity assumption is violated.
#' @seealso \code{\link{coxphlb}}, \code{\link{coxphlb.ftest}}, \code{\link{coxphlb.phtest}}, \code{\link{station.test.plot}}
#' @references Addona, V. and Wolfson, D. B. (2006). A formal test for the stationarity of the incidence rate using data from a prevalent cohort study with follow-up. \emph{Lifetime data analysis}, 12(3), 267-284.
#' @examples
#' # Check the Stationarity Assumption
#' stest1 <- station.test(ExampleData1$a, ExampleData1$y-ExampleData1$a,
#'                        ExampleData1$delta)
#' print(stest1) 			# display the results
#'
#' stest2 <- station.test(ExampleData2$a, ExampleData2$y-ExampleData2$a,
#'                        ExampleData2$delta)
#' print(stest2)   		# display the results
#' @export

station.test <- function(a, v, delta, digits = 3L)
{
  len<-length(delta)

  sum<-0

  for(i in 1:len)
    for(j in 1:len)
    {
      sum<- sum + ifelse(a[i]>v[j] & delta[j]==1, 1, 0) - ifelse(a[i]<v[j], 1, 0)
    }
  # Calculate test statistic

  wtest <- sum/len^2
  # Calculate the variance of the test statistic
  var <- 0
  for( j in 1:len)
  {
    var <- var + sum(ifelse(v>a[j], 1, 0))^2+
      sum(ifelse(a>v[j]& delta[j]==1, 1, 0))^2+
      2*sum(ifelse(v<=a[j] & delta==1, 1, 0))*sum(ifelse(a>v[j], 1, 0))*delta[j]+
      2*sum(ifelse(v>a[j], 1, 0))*sum(ifelse(a<=v[j], 1, 0))+
      2*sum(ifelse(a<=v[j], 1, 0))*sum(ifelse(v<=a[j] & delta==1, 1, 0))-
      2*sum(ifelse(a>v[j], 1, 0))*sum(ifelse(v>a[j], 1, 0))*delta[j]

  }
  # standardized test statistics
  test.statistic <- sqrt(len)*wtest/sqrt(var/len^3)
  p.value <- 1 + pnorm(-abs(test.statistic)) - pnorm(abs(test.statistic))

  tab = cbind(test.statistic = format(round(test.statistic, digits)), p.value = format.pval(p.value, digits, eps = 1e-03))
  rownames(tab) = " "
  tab = as.table(tab)
  print(tab)
  out=list(test.statistic = test.statistic, p.value = p.value, result = tab)
  class(out) <- "station.test"
  invisible(out)
}
