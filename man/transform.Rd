% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal-distributions.R
\name{transform}
\alias{transform}
\title{Transform a precipitation observations from the precipitation scale
to the Laplace scale using the probability integral transform}
\usage{
transform(x, a, b, u, p_u, xi, s, alpha)
}
\description{
Input variables:
x: The data to be transformed
a: shape-parameter of the gamma distribution
b: rate-parameter of the gamma distribution
u: threshold for the GP distribution
p_u: The probability corresponding to the threshold u
xi: Tail parameter of the GP distribution
s: Scale parameter of the GP distribution
alpha: probability deciding which quantile the scale parameter is equal to
}
