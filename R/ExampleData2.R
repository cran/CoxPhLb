#' @title Example Data 2
#' @description A simulated right-censored left-truncated data set. The example data set contains failure time, left-truncation time, censoring indicator, and the covariates.
#' @name ExampleData2
#' @docType data
#' @usage ExampleData2
#' @format A data frame with 200 observations.
#' \describe{
#'   \item{y}{failure time}
#'   \item{a}{left-truncation time}
#'   \item{delta}{censoring status (0=censored, 1=uncensored)}
#'   \item{x1}{first covariate (binary variable, 0 or 1)}
#'   \item{x2}{second covariate (continuous variable, range from 0.0 to 1.0)}
#' }
#' @note The stationarity assumption is violated for this example data.
#' @keywords ExampleData2
NULL
