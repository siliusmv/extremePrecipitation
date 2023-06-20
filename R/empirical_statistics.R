#' @export
empirical_moments = function(data,
                             dist_centers,
                             log_y0_centers,
                             dist_radius = .5,
                             log_y0_radius = .05,
                             add_zero_dist = TRUE) {
  y_sum = matrix(0, length(dist_centers), length(log_y0_centers))
  y_sq_sum = matrix(0, length(dist_centers), length(log_y0_centers))
  y_n = matrix(0, length(dist_centers), length(log_y0_centers))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = data$dist_to_s0[[j]]
    log_y0 = log(data$y0[[j]])
    for (d in seq_along(dist_centers)) {
      d_index = which(dist > dist_centers[d] - dist_radius & dist < dist_centers[d] + dist_radius)
      if (!any(d_index)) next
      for (i in seq_along(log_y0_centers)) {
        y0_index = which(
          log_y0 > log_y0_centers[i] - log_y0_radius &
          log_y0 < log_y0_centers[i] + log_y0_radius)
        if (!any(y0_index)) next
        y_sum[d, i] = y_sum[d, i] + sum(data$y[[j]][d_index, y0_index], na.rm = TRUE)
        y_sq_sum[d, i] = y_sq_sum[d, i] + sum(data$y[[j]][d_index, y0_index]^2, na.rm = TRUE)
        y_n[d, i] = y_n[d, i] + sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  y_mean = y_sum / y_n
  y_sq_mean = y_sq_sum / y_n
  y_sd = sqrt(y_sq_mean - y_mean^2)
  res = data.frame(
    mean = as.numeric(y_mean),
    sd = as.numeric(y_sd),
    log_y0 = rep(log_y0_centers, each = length(dist_centers)),
    dist = rep(dist_centers, length(log_y0_centers)),
    n = as.numeric(y_n))
  if (add_zero_dist) {
    zero_dist_df = data.frame(
      log_y0 = unique(res$log_y0),
      dist = 0,
      mean = unique(exp(res$log_y0)),
      sd = 0,
      n = 1)
    res = rbind(res, zero_dist_df)
  }
  res$y0 = exp(res$log_y0)
  res
}

#' @export
empirical_chi = function(data,
                         thresholds,
                         dist_centers,
                         dist_radius = .5,
                         add_zero_dist = TRUE) {
  n_big = matrix(0, length(dist_centers), length(thresholds))
  n_all = matrix(0, length(dist_centers), length(thresholds))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = data$dist_to_s0[[j]]
    for (d in seq_along(dist_centers)) {
      d_index = which(dist > dist_centers[d] - dist_radius & dist < dist_centers[d] + dist_radius)
      if (!any(d_index)) next
      for (i in seq_along(thresholds)) {
        y0_index = which(data$y0[[j]] > thresholds[i])
        if (!any(y0_index)) next
        n_big[d, i] = n_big[d, i] +
          sum(data$y[[j]][d_index, y0_index] > thresholds[i], na.rm = TRUE)
        n_all[d, i] = n_all[d, i] +
          sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  chi = n_big / n_all
  if (add_zero_dist) {
    chi = rbind(1, chi)
    n_all = rbind(1, n_all)
    dist_centers = c(0, dist_centers)
  }
  data.frame(
    chi = as.numeric(chi),
    n = as.numeric(n_all),
    threshold = rep(thresholds, each = length(dist_centers)),
    dist = rep(dist_centers, length(thresholds)))
}
