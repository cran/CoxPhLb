#' @title Graphical Test of the Functional Form of Covariates in Cox Model with Right-Censored Length-Biased Data
#' @description Returns a plot of the cumulative sums of mean zero stochastic processes.
#' @import survival
#' @usage coxphlb.ftest.plot(x, n.plot = 20, seed.n = round(runif(1,1,1e09)))
#' @param x The result of the \code{coxphlb.ftest} function.
#' @param n.plot The number of randomly selected realizations. Default is 20.
#' @param seed.n An integer specifying seed number.
#' @details The function returns a plot with the test statistics in a black line and 20 processes randomly sampled from the pool of resamples in grey lines. When the test statistics lie within the randomly sampled lines, it suggests that the model assumption is valid.
#' @seealso \code{\link{coxphlb}}, \code{\link{coxphlb.ftest}}
#' @references Lee, C.H., Ning, J., and Shen, Y. Model diagnostics for proportional hazards model with length-biased data. \emph{Lifetime Data Analysis} 25(1), 79-96.
#' @examples
#' \dontrun{
#' # Fit a Cox model
#' fit.ee <- coxphlb(Surv(a, y, delta) ~ x1 + x2, data = ExampleData1,
#'                  method = "EE")
#'
#' # Check the Functional Form of the Cox Model
#' ftest <- coxphlb.ftest(fit.ee, data = ExampleData1, spec.p = 2,
#'                       seed.n = 1234)
#' coxphlb.ftest.plot(ftest, n.plot = 50, seed.n = 1234)			# display the plot
#' }
#' @export

coxphlb.ftest.plot <- function(x, n.plot = 20, seed.n = round(runif(1,1,1e09))) {
  if (!inherits(x, 'coxphlb.ftest')) stop ("Argument must be the result of coxphlb.ftest")

  n.sim = x$n.sim
  z0 = x$z0
  stat.mat.z = x$stat.mat.z
  sim.mat.z = x$sim.mat.z
  varnames = x$varnames

  set.seed(seed.n)
  samp.ind = sample(1:n.sim, n.plot, replace=FALSE)
    plot(z0, stat.mat.z, type="l", ylim=c(-max(abs(sim.mat.z)), max(abs(sim.mat.z))), xlab="z", ylab="Partial-Mean Processes", main=varnames)
    for (j in samp.ind) {
      lines(z0, sim.mat.z[j,], lwd = 1, col = "grey")
    }
    lines(z0, stat.mat.z, lwd = 1.5)
}
