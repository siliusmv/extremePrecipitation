#' Search through data after threshold exceedances at a set of conditioning sites.
#' Extract relevant information and observations for each of the threshold exceedances.
#'
#' The input variables are:
#' data: A (d x n)-dimensional matrix of observations.
#' coords: A (d x 2)-dimensional matrix with the coordinates of the data.
#' s0_index: a vector of indices describing which of the d coordinates are chosen
#'   as conditioning sites.
#' threshold: The threshold t for the conditional extremes model. This can either be a single
#'   number, or a vector with one element for each of the n columns of the data
#' remove_y0_from_y: a Boolean stating whether y0 should be removed from y or not
#' obs_index: A vector of indices showing which coordinates we care about.
#' n, r: Used for only extracting a subset of the observations in a circle around
#'   the conditioning site. See the docs for extract_thinned_out_circles() for more info.
#' If obs_index is non-null, we will not use the values of n and r.
#'
#' The output of the function is a list where each element is a list or vector
#' of length equal, containing these variables:
#' s0_index: This is a subset of all the s0_indices given as input such that each
#'   conditioning site has at least one threshold exceedance.
#' s0: A list where element nr. i is a 2-dimensional vector with the coordinates of
#'   conditioning site nr. i.
#' obs_index: A list where element nr. i is a vector of indices for all the
#'   locations in coords we use for performing inference using conditioning site nr. i.
#' dist_to_s0: A list where element nr. i is a vector with all the distances between
#'   conditioning site nr. i and the locations used for inference with conditioning
#'   site nr. i.
#' y0: A list where element nr. i is a vector of all threshold exceedances for
#'   conditioning site nr. i.
#' y: A list where element nr. i is a matrix consisting of one column per threshold
#'   exceedance in y0[[i]]. Each column contains observations from the locations
#'   described in obs_index[[i]].
#' time_index: A list where element nr. i contains the times of all threshold
#'   exceedances from conditioning site nr. i.
#' n: A vector where element nr. i is the number of threshold exceedances at
#'   conditioning site nr. i.
#' @export
extract_extreme_fields = function(data,
                                  coords,
                                  s0_index,
                                  threshold,
                                  obs_index = NULL,
                                  n = NULL,
                                  r = NULL,
                                  verbose = FALSE,
                                  remove_y0_from_y = TRUE,
                                  n_cores = 1) {
  # threshold must have length 1 or lenght equal to the number of columns in data
  stopifnot(length(threshold) %in% c(1, ncol(data)))
  # Either obs_index or (n, r) must be non-null
  stopifnot(!is.null(obs_index) || (!is.null(n) && !is.null(r)))
  # s0_indices and obs_indices must be between 1 and nrow(coords)
  if (!is.null(obs_index)) stopifnot(1 <= min(obs_index) && max(obs_index) <= nrow(coords))
  stopifnot(1 <= min(s0_index) && max(s0_index) <= nrow(coords))
  n_s0 = length(s0_index)

  # Look through the conditioning sites
  res = parallel::mclapply(
    X = seq_len(n_s0),
    mc.cores = n_cores,
    FUN = function(i) {
      res = list(s0_index = s0_index[i])

      # Extract all variables at s0 nr. i, and locate threshold exceedances
      y0 = data[s0_index[i], ]
      res$time_index = which(y0 > threshold)
      # If there are no threshold exceedances, skip to next conditioning site
      if (length(res$time_index) == 0) return(NULL)

      # Extract all threshold exceedances
      res$y0 = y0[res$time_index]

      # Locate the thinned out observations locations
      if (!is.null(obs_index)) {
        res$obs_index = obs_index
      } else {
        res$obs_index = extract_thinned_out_circles(
          coords = coords,
          center = coords[s0_index[i], , drop = FALSE],
          n = n,
          r = r,
          index_only = TRUE)
      }
      if (remove_y0_from_y) res$obs_index = res$obs_index[res$obs_index != s0_index[i]]

      # Extract all observations of interest for the threshold exceedances
      res$y = data[res$obs_index, res$time_index, drop = FALSE]

      # Fill in information about s0 and the dist_to_s0
      res$s0 = coords[s0_index[i], , drop = FALSE]
      res$dist_to_s0 = as.numeric(dist_euclid(
        x = coords[s0_index[i], , drop = FALSE],
        y = coords[res$obs_index, , drop = FALSE]))

      if (verbose) message(i, " / ", n_s0)

      res
    })

  # Only keep information about the conditioning sites with at least one threshold exceedance
  good_index = which(sapply(res, length) > 0)
  res = purrr::transpose(res[good_index])

  res$s0_index = unlist(res$s0_index)

  # Find the number of threshold exceedances for each conditioning site
  res$n = sapply(res$y0, length)
  res$n_loc = sapply(res$dist_to_s0, length)

  res
}

#' @export
get_s0_index = function(coords, delta_s0 = 1) {
  stopifnot(all(as.integer(delta_s0) == delta_s0))
  x_coords = coords[, 1] |> unique() |> sort()
  y_coords = coords[, 2] |> unique() |> sort()
  s0_locs = expand.grid(
    x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
    y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
    as.matrix()
  s0_index = lapply(
    X = 1:nrow(s0_locs),
    FUN = function(i) which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2]))
  s0_index = s0_index[sapply(s0_index, length) > 0] |>
    unlist() |>
    unname()
  s0_index
}

#' Given a matrix of coordinates on a regular grid, extract a subset of the coordinates
#' that consists of circles centred around some center, with varying degrees of
#' densities. The input variables are:
#' coords: A matrix with two columns and multiple rows, containing (x, y) coordinates
#'   in a regular grid and with a distance-preserving projection.
#' center: One (x, y) tuple of coordinates, explaining the center of the circles we
#'   wish to extract.
#' n: A vector of densities for the different circles we wish to extract. If the ith
#'   value of n is k, then we only keep every kth of the original coords, both in
#'   the x-direction and in the y-direction.
#' r: A vector of radii for the different circles we wish to extract. If e.g.
#'   n = (1, 3) and r = (2, 5), then we will extract every coordinate that is closer
#'   to the center than a radius of 2. Then, we will only keep every 3rd coordinate,
#'   both in x-direction and y-direction, for all coordinates that have a distance
#'   between 2 and 5 from the center. Finally, all coords that are further away from
#'   the center than 5 will be dropped.
#' index_only: A bool describing if the function should return the extracted
#'   coordinates, or only their indices inside the coords object.
#' @export
extract_thinned_out_circles = function(coords,
                                       center,
                                       n,
                                       r = Inf,
                                       index_only = FALSE) {
  stopifnot(ncol(coords) == 2)
  stopifnot(length(n) == length(r))
  stopifnot(is.matrix(coords))

  # Locate all unique x-values and y-values
  x_vals = sort(unique(coords[, 1]))
  y_vals = sort(unique(coords[, 2]))

  # Compute distances and find the coord that is closest to the center
  dist_to_center = as.numeric(dist_euclid(center, coords))
  closest_to_center = coords[which.min(dist_to_center), ]
  x_center_index = which(x_vals == closest_to_center[1])
  y_center_index = which(y_vals == closest_to_center[2])

  # This will be the list of all indices we wish to extract from coords
  chosen_indices = NULL

  # Locate all the indices we wish to extract
  r = c(0, r)
  for (i in seq_along(n)) {
    # Indicies for all coords that have the correct distance from center
    dist_index = which(r[i] <= dist_to_center & dist_to_center <= r[i + 1])

    # Indices for all coords on a regular grid of resolution n[i] x n[i]
    thinned_x_index = c(seq(x_center_index, 1, by = -n[i]),
                        seq(x_center_index, length(x_vals), by = n[i]))
    thinned_y_index = c(seq(y_center_index, 1, by = -n[i]),
                        seq(y_center_index, length(y_vals), by = n[i]))
    thinned_x_vals = x_vals[thinned_x_index]
    thinned_y_vals = y_vals[thinned_y_index]
    thinned_index = which(coords[, 1] %in% thinned_x_vals & coords[, 2] %in% thinned_y_vals)

    # The indices of interest are the intersection of the two index sets above
    chosen_indices = c(chosen_indices, intersect(dist_index, thinned_index))
  }
  # Some times, duplicates may appear. Remove these
  chosen_indices = unique(chosen_indices)

  if (index_only) {
    res = chosen_indices
  } else {
    res = coords[chosen_indices, ]
  }
  res
}
