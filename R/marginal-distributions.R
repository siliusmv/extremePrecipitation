
#' x: The data to be transformed
#' a: shape-parameter of the gamma distribution
#' b: rate-parameter of the gamma distribution
#' u: threshold for the GP distribution
#' u_prob: The probability corresponding to the threshold u
#' xi: Tail parameter of the GP distribution
#' s: Scale parameter of the GP distribution
#' alpha: probability deciding which quantile the scale parameter is equal to
#' @export
transform = function(x, a, b, u, u_prob, xi, s, alpha) {
  res = rep(NA_real_, length(x))
  small_index = which(x <= u)
  if (any(small_index)) {
    res[small_index] = pgamma(x[small_index], shape = a, rate = b)
  }
  big_index = which(x > u)
  if (any(big_index)) {
    res[big_index] = u_prob + (1 - u_prob) * pgp(x[big_index], s, xi, u = u, alpha = alpha)
  }
  qlaplace(res)
}

#' x: The data to be transformed
#' a: shape-parameter of the gamma distribution
#' b: rate-parameter of the gamma distribution
#' u: threshold for the GP distribution
#' u_prob: The probability corresponding to the threshold u
#' xi: Tail parameter of the GP distribution
#' s: Scale parameter of the GP distribution
#' alpha: probability deciding which quantile the scale parameter is equal to
#' @export
inv_transform = function(x, a, b, u, u_prob, xi, s, alpha) {
  p = plaplace(x)
  res = rep(NA_real_, length(p))
  small_index = which(p <= u_prob)
  if (any(small_index)) {
    res[small_index] = qgamma(p[small_index], shape = a, rate = b)
  }
  big_index = which(p > u_prob)
  if (any(big_index)) {
    res[big_index] = qgp((p[big_index] - u_prob) / (1 - u_prob), s, xi, u = u, alpha = alpha)
  }
  res
}

#' @export
pgp = function(x, s, xi, u = 0, alpha = .5) {
  if (length(xi) == 1) xi = rep(xi, length(x))
  stopifnot(length(xi) == length(x))

  z = pmax((x - u) / s, 0)
  res = 1 - ifelse(
    test = (xi == 0),
    yes = exp(z * log(1 - alpha)),
    no = pmax(1 + z * ((1 - alpha)^-xi - 1), 0) ^ (-1 / xi))

  res
}

#' @export
qgp = function(p, s, xi, u = 0, alpha = .5) {
  if (length(xi) == 1) xi = rep(xi, length(p))
  stopifnot(length(xi) == length(p))
  stopifnot(all(p >= 0 & p <= 1, na.rm = TRUE))

  res = u + s * ifelse(
    test = (xi == 0),
    yes = log(1 - p) / log(1 - alpha),
    no = ((1 - p)^-xi - 1) / ((1 - alpha)^-xi - 1))

  res
}

#' @export
plaplace = function(x) {
  res = x
  below_zero = which(x <= 0)
  above_zero = which(x > 0)
  if (any(below_zero)) res[below_zero] = exp(x[below_zero]) / 2
  if (any(above_zero)) res[above_zero] = 1 - exp(-x[above_zero]) / 2
  res
}

#' @export
qlaplace = function(p) {
  stopifnot(all(p >= 0 & p <= 1, na.rm = TRUE))
  res = p
  below_median = which(p <= .5 & p >= 0)
  above_median = which(p > .5 & p <= 1)
  if (any(below_median)) res[below_median] = log(2 * p[below_median])
  if (any(above_median)) res[above_median] = -log(2 * (1 - p[above_median]))
  res
}
