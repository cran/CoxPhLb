#' @title Graphical Test of the Proportional Hazards Assumption of Covariates in Cox Model with Right-Censored Length-Biased Data
#' @description Returns a plot of the cumulative sums of mean zero stochastic processes.
#' @import survival
#' @usage coxphlb.phtest.plot(x, n.plot = 20, seed.n = round(runif(1,1,1e09)))
#' @param x The result of the \code{coxphlb.phtest} function.
#' @param n.plot The number of randomly selected realizations. Default is 20.
#' @param seed.n An integer specifying seed number.
#' @details The function returns a plot with the test statistics in a black line and 20 processes randomly sampled from the pool of resamples in grey lines. When the test statistics lie within the randomly sampled lines, it suggests that the model assumption is valid. A plot cannot be generated for the global test.
#' @seealso \code{\link{coxphlb}}, \code{\link{coxphlb.phtest}}
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
#' coxphlb.phtest.plot(ptest1, n.plot = 50, seed.n = 1234)			# display the plot
#' }
#' @export

coxphlb.phtest.plot <- function(x, n.plot = 20, seed.n = round(runif(1,1,1e09))) {
  if (!inherits(x, 'coxphlb.phtest')) stop ("Argument must be the result of coxphlb.phtest")

  n.sim = x$n.sim
  stat.mat.t = x$stat.mat.t
  sim.mat.t = x$sim.mat.t
  varnames = x$varnames
  spec.p = x$spec.p
  yy = x$yy

  if (is.null(spec.p)) stop ("To print plot, spec.p should be specified in coxphlb.phtest")

  set.seed(seed.n)
  samp.ind = sample(1:n.sim, n.plot, replace=FALSE)
    tmp.p = spec.p
    plot(yy, stat.mat.t[,tmp.p], type="l", ylim=c(-max(abs(sim.mat.t[,,tmp.p])), max(abs(sim.mat.t[,,tmp.p]))), xlab="t", ylab="Pseudo-Score Processes", main=varnames)
    for (j in samp.ind) {
      lines(yy, sim.mat.t[j,,tmp.p], lwd = 1, col="grey")
    }
    lines(yy, stat.mat.t[,tmp.p], lwd=1.5)
}
