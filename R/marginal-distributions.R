
#' Transform a precipitation observations from the precipitation scale
#' to the Laplace scale using the probability integral transform
#'
#' Input variables:
#' x: The data to be transformed
#' a: shape-parameter of the gamma distribution
#' b: rate-parameter of the gamma distribution
#' u: threshold for the GP distribution
#' p_u: The probability corresponding to the threshold u
#' xi: Tail parameter of the GP distribution
#' s: Scale parameter of the GP distribution
#' alpha: probability deciding which quantile the scale parameter is equal to
#' @export
transform = function(x, a, b, u, p_u, xi, s, alpha) {
  res = rep(NA_real_, length(x))
  small_index = which(x <= u)
  if (any(small_index)) {
    res[small_index] = pgamma(x[small_index], shape = a, rate = b)
  }
  big_index = which(x > u)
  if (any(big_index)) {
    res[big_index] = p_u + (1 - p_u) * pgp(x[big_index], s, xi, u = u, alpha = alpha)
  }
  qlaplace(res)
}

#' Transform a random variable on the Laplace scale back to the precipitation scale
#' using the probability integral transform
#'
#' x: The data to be transformed
#' a: shape-parameter of the gamma distribution
#' b: rate-parameter of the gamma distribution
#' u: threshold for the GP distribution
#' p_u: The probability corresponding to the threshold u
#' xi: Tail parameter of the GP distribution
#' s: Scale parameter of the GP distribution
#' alpha: probability deciding which quantile the scale parameter is equal to
#' @export
inv_transform = function(x, a, b, u, p_u, xi, s, alpha) {
  p = plaplace(x)
  res = rep(NA_real_, length(p))
  small_index = which(p <= p_u)
  if (any(small_index)) {
    res[small_index] = qgamma(p[small_index], shape = a, rate = b)
  }
  big_index = which(p > p_u)
  if (any(big_index)) {
    res[big_index] = qgp((p[big_index] - p_u) / (1 - p_u), s, xi, u = u, alpha = alpha)
  }
  res
}

#' Cumulative distribution function for the generalised Pareto (GP) distribution
#'
#' Input variables:
#' x: The input x to the function F(x)
#' s: The scale parameter of the GP distribution
#' xi: The shape parameter of the GP distribution
#' u: The location parameter of the GP distribution
#' alpha: a probability such that the scale parameter s is equal to the
#'   alpha-quantile of the GP distribution
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

#' Quantile function for the generalised Pareto (GP) distribution
#'
#' Input variables:
#' p: The input p to the function Q(p)
#' s: The scale parameter of the GP distribution
#' xi: The shape parameter of the GP distribution
#' u: The location parameter of the GP distribution
#' alpha: a probability such that the scale parameter s is equal to the
#'   alpha-quantile of the GP distribution
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

#' Cumulative distribution function for the Laplace distribution
#' @export
plaplace = function(x) {
  res = x
  below_zero = which(x <= 0)
  above_zero = which(x > 0)
  if (any(below_zero)) res[below_zero] = exp(x[below_zero]) / 2
  if (any(above_zero)) res[above_zero] = 1 - exp(-x[above_zero]) / 2
  res
}

#' Quantile function for the Laplace distribution
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
