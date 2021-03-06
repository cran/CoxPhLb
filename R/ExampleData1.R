#' @title Example Data 1
#' @description A simulated right-censored length-biased data set. The example data set contains failure time, left-truncation time, censoring indicator, and the covariates.
#' @name ExampleData1
#' @docType data
#' @usage ExampleData1
#' @format A data frame with 200 observations.
#' \describe{
#'   \item{y}{failure time}
#'   \item{a}{left-truncation time}
#'   \item{delta}{censoring status (0=censored, 1=uncensored)}
#'   \item{x1}{first covariate (binary variable, 0 or 1)}
#'   \item{x2}{second covariate (continuous variable, range from 0.0 to 1.0)}
#' }
#' @note This example data satisfy the stationarity assumption.
#' @keywords ExampleData1
NULL
