#' @export
summary.coxphlb = function(object, ...) {
  x = object
  if (!inherits(x, 'coxphlb')) stop ("Argument must be the result of coxphlb")
  form = paste(x$formula)
  cat("Call:\n","coxphlb(formula = ",form[2]," ",form[1]," ",form[3],", method = ", x$method,")\n", sep = "")
  print(x$result)
}

#' @export
print.coxphlb = function(x, ...) {
  if (!inherits(x, 'coxphlb')) stop ("Argument must be the result of coxphlb")
  print(x$result)
}

#' @export
vcov.coxphlb = function(object, ...) {
  x = object
  if (!inherits(x, 'coxphlb')) stop ("Argument must be the result of coxphlb")
  out = x$var
  colnames(out) = x$varnames[-(1:3)]
  rownames(out) = x$varnames[-(1:3)]
  out
}

#' @export
coef.coxphlb = function(object, ...) {
  x = object
  if (!inherits(x, 'coxphlb')) stop ("Argument must be the result of coxphlb")
  out = as.vector(x$coefficient)
  names(out) = x$varnames[-(1:3)]
  out
}

#' @export
print.coxphlb.ftest = function(x, ...) {
  if (!inherits(x, 'coxphlb.ftest')) stop ("Argument must be the result of coxphlb.ftest")
  print(x$result)
}

#' @export
print.coxphlb.phtest = function(x, ...) {
  if (!inherits(x, 'coxphlb.phtest')) stop ("Argument must be the result of coxphlb.phtest")
  print(x$result)
}

#' @export
print.station.test = function(x, ...) {
  if (!inherits(x, 'station.test')) stop ("Argument must be the result of station.test")
  print(x$result)
}
