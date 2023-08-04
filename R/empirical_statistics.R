
#' This function takes in a data object that has been created by the
#' extract_extreme_fields() function. Then, a sliding window approach
#' is performed over the distances to the conditioning sites and the
#' values of the threshold exceedances at the conditioning sites,
#' for computing the empirical mean and standard deviation of all observations
#' with a distance d ± Δd to a conditioning site where a threshold exceedance log(y0) ± Δy0 was
#' observed, for a range of different values of d and y0.
#'
#' The input variables are:
#' data: This is a data object created by the extract_extreme_fields() function
#' dist_centers: The values of d used for computing empirical moments
#' log_y0_centers: The values of log(y0) used for computing empirical moments
#' dist_radius: The value of Δd
#' log_y0_radius: The value of Δy0
#' add_zero_dist: Boolean variable. Should we add rows for d = 0 in the output?
#'
#' The output of the function is a data.frame with one row for each combination of
#' the given values of d and log(y0). The data.frame has columns:
#' mean: Empirical mean of all observations with distance d ± Δd to a conditioning site where a
#'   threshold exceedance log(y0) ± Δy0 was observed.
#' sd: Empirical standard deviation of all observations with distance d ± Δd to a conditioning site
#'   where a threshold exceedance log(y0) ± Δy0 was observed.
#' log_y0: Value of log(y0).
#' dist: Value of d.
#' n: Number of observations distance d ± Δd to a conditioning site where a
#'   threshold exceedance log(y0) ± Δy0 was observed.
#' @export
empirical_moments = function(data,
                             dist_centers,
                             log_y0_centers,
                             dist_radius = .5,
                             log_y0_radius = .05,
                             add_zero_dist = TRUE) {
  # Preallocate matrices for computing the sum, squared sum and number of observations
  # inside each of the sliding windows
  y_sum = matrix(0, length(dist_centers), length(log_y0_centers))
  y_sq_sum = matrix(0, length(dist_centers), length(log_y0_centers))
  y_n = matrix(0, length(dist_centers), length(log_y0_centers))

  # Loop over all the sliding windows and print the progress
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) { # Loop over all conditioning sites in the data
    dist = data$dist_to_s0[[j]]
    log_y0 = log(data$y0[[j]])
    for (d in seq_along(dist_centers)) { # Loop over all distances of interest
      # Find which observations have a distance d ± Δd to the conditioning site
      d_index = which(dist > dist_centers[d] - dist_radius & dist < dist_centers[d] + dist_radius)
      if (!any(d_index)) next
      for (i in seq_along(log_y0_centers)) { # Loop over all threshold exceedances of interest
        # Find which threshold exceedanecs have a value of exp(log(y0) ± Δy0)
        y0_index = which(
          log_y0 > log_y0_centers[i] - log_y0_radius &
          log_y0 < log_y0_centers[i] + log_y0_radius)
        if (!any(y0_index)) next

        # Add information to the preallocated result matrices
        y_sum[d, i] = y_sum[d, i] + sum(data$y[[j]][d_index, y0_index], na.rm = TRUE)
        y_sq_sum[d, i] = y_sq_sum[d, i] + sum(data$y[[j]][d_index, y0_index]^2, na.rm = TRUE)
        y_n[d, i] = y_n[d, i] + sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()

  # Compute the mean and standard deviation for all combinations of d and log(y0)
  y_mean = y_sum / y_n
  y_sq_mean = y_sq_sum / y_n
  y_sd = sqrt(y_sq_mean - y_mean^2)

  # Create a data.frame with the results
  res = data.frame(
    mean = as.numeric(y_mean),
    sd = as.numeric(y_sd),
    log_y0 = rep(log_y0_centers, each = length(dist_centers)),
    dist = rep(dist_centers, length(log_y0_centers)),
    n = as.numeric(y_n))

  # Possibly add rows for d = 0 into the data.frame
  if (add_zero_dist) {
    zero_dist_df = data.frame(
      log_y0 = unique(res$log_y0),
      dist = 0,
      mean = unique(exp(res$log_y0)),
      sd = 0,
      n = 1)
    res = rbind(res, zero_dist_df)
  }

  # Return the results
  res$y0 = exp(res$log_y0)
  res
}

#' This function takes in a data object that has been created by the
#' extract_extreme_fields() function. Then, a sliding window approach
#' is performed over the distances to the conditioning sites,
#' for computing the empirical exceedance probabilities given a threshold τ, for all observations
#' with a distance d ± Δd to a conditioning site, for a range of different values of d and τ.
#'
#' The input variables are:
#' data: This is a data object created by the extract_extreme_fields() function
#' threshold: A vector of thresholds used for computing empirical exceedance probabilities
#' dist_centers: The values of d used for computing empirical exceedance probabilities
#' dist_radius: The value of Δd
#' add_zero_dist: Boolean variable. Should we add rows for d = 0 in the output?
#'
#' The output of the function is a data.frame with one row for each combination of
#' the given values of d and τ. The data.frame has columns:
#' chi: Empirical exceedance probability, for the threshold τ, for all observations with
#'   distance d ± Δd to a conditioning site of interest.
#' n: Number of observations with distance d ± Δd to a conditioning site that exceeds a
#'   threshold of τ
#' Threshold: Value of τ.
#' dist: Value of d.
#' @export
empirical_chi = function(data,
                         thresholds,
                         dist_centers,
                         dist_radius = .5,
                         add_zero_dist = TRUE) {
  # Preallocate matrices for computing the number of observations with distance d ± Δd to a
  # conditioning site that exceeds the threshold τ, and the number of observations that also
  # exceed the same threshold
  n_big = matrix(0, length(dist_centers), length(thresholds))
  n_all = matrix(0, length(dist_centers), length(thresholds))

  # Loop over all the sliding windows and print the progress
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) { # Loop over all conditioning sites in the data
    dist = data$dist_to_s0[[j]]
    for (d in seq_along(dist_centers)) { # Loop over all distances of interest
      # Find which observations have a distance d ± Δd to the conditioning site
      d_index = which(dist > dist_centers[d] - dist_radius & dist < dist_centers[d] + dist_radius)
      if (!any(d_index)) next
      for (i in seq_along(thresholds)) { # Loop over all thresholds of interest
        y0_index = which(data$y0[[j]] > thresholds[i])
        if (!any(y0_index)) next

        # Find the number of observations with the correct distance that exceed the threshold τ,
        # and the total number of observations with the correct distance
        n_big[d, i] = n_big[d, i] +
          sum(data$y[[j]][d_index, y0_index] > thresholds[i], na.rm = TRUE)
        n_all[d, i] = n_all[d, i] +
          sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()

  # Compute the empirical exceedance probabilities
  chi = n_big / n_all

  # Possibly add information for when d = 0
  if (add_zero_dist) {
    chi = rbind(1, chi)
    n_all = rbind(1, n_all)
    dist_centers = c(0, dist_centers)
  }

  # Return the results in a data.frame
  data.frame(
    chi = as.numeric(chi),
    n = as.numeric(n_all),
    threshold = rep(thresholds, each = length(dist_centers)),
    dist = rep(dist_centers, length(thresholds)))
}
