#' @export
threshold_occurrence = function(samples, nonzero_prob) {
  threshold = quantile(samples, probs = 1 - nonzero_prob, na.rm = TRUE)
  samples[samples < threshold] = 0
  samples
}

#' @export
simulate_intensity_posterior = function(samples,
                                        threshold,
                                        coords,
                                        s0_index,
                                        spde,
                                        get_a_func,
                                        get_b_func,
                                        get_Q,
                                        get_tau,
                                        verbose = FALSE,
                                        n_cores = 1,
                                        n_per_sample = 1) {
  # Simulate all threshold exceedances
  n_samples = nrow(samples)
  y0 = rexp(n_samples * n_per_sample) + threshold

  s0 = coords[s0_index, , drop = FALSE]
  coords = coords[-s0_index, , drop = FALSE]
  A = inla.spde.make.A(spde$mesh, coords)

  # Compute distances to the conditioning site
  dist_to_s0 = as.numeric(dist_euclid(s0, coords))
  dist_to_s0_from_mesh = as.numeric(dist_euclid(s0, spde$mesh$loc[, -3]))

  # Remove elements from the projection matrix that correspond to the conditioning site
  ii = which(dist_to_s0_from_mesh == 0)
  dist_to_s0_from_mesh = dist_to_s0_from_mesh[-ii]
  A = A[, -ii]

  # Create a data-object of same type as that returned by the function extract_extreme_fields(),
  # for saving the simulated realisations
  data = list(
    y0 = list(y0),
    dist_to_s0 = list(dist_to_s0),
    coords = coords,
    s0 = s0,
    n = n_samples * n_per_sample)

  # Simulate from the fitted conditional extremes model
  data$y =
    local({
      RNGkind("L'Ecuyer-CMRG")
      parallel::mclapply(
        X = seq_len(n_samples),
        mc.cores = n_cores,
        FUN = function(i) {
          # These are the parameters from model-parameter sample nr. i
          params = samples[i, ]

          # Create and constrain the precision matrix
          Q = get_Q(params, spde)
          Q = Q[-ii, -ii]

          # Simulate from the spatial conditional distribution
          res = rconditional(
            y0 = list(y0[(i - 1) * n_per_sample + 1:n_per_sample]),
            a_func = get_a_func(params),
            b_func = get_b_func(params),
            Q = Q,
            tau = get_tau(params),
            dist_to_s0 = list(dist_to_s0),
            dist_to_s0_from_mesh = list(dist_to_s0_from_mesh),
            A = list(A))[[1]]

          if (verbose && (i %% 10 == 0)) {
            message("Finished sample nr. ", i, " / ", n_samples)
          }

          res
        })
    })

  data$y = list(do.call(cbind, data$y))
  data
}

#' Function for simulating from the conditional extremes model given the values of y(s0).
#' The arguments of the function are:
#' y0: A list of vectors containing the values of y(s0) for each of the possible
#'   conditioning sites s0. The vectors can have length 0.
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors may be of different size.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must be of the same size.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix may contain a different number of
#'   rows, but all must have the same number of columns, which is equal to the number
#'   nodes in the triangular mesh.
#' @export
rconditional = function(y0,
                        a_func,
                        b_func,
                        Q,
                        tau,
                        dist_to_s0,
                        dist_to_s0_from_mesh,
                        A) {
  n = sapply(y0, length)
  n_loc = length(y0)
  # Check that the input has correct lengths
  lengths = c(length(A), length(dist_to_s0), length(dist_to_s0_from_mesh))
  stopifnot(all(lengths == n_loc))
  # Check that the A matrices have the same number of rows as the lengths of dist_to_s0
  stopifnot(all(sapply(dist_to_s0, length) == sapply(A, nrow)))
  # Check that all the dist_to_s0_from_mesh dist_to_s0ances have length equal to the
  # number of nodes in the triangular mesh
  stopifnot(all(sapply(dist_to_s0_from_mesh, length) == nrow(Q)))

  # Sample from the SPDE approximation
  z = rnorm_spde(n = sum(n), Q = Q)

  # Preallocate the result
  res = vector("list", n_loc)

  # Sample from the conditional extremes model, given the values of a, b and z
  for (i in 1:n_loc) {
    if (n[i] == 0) next
    index = sum(n[seq_len(i - 1)]) + seq_len(n[i])
    res[[i]] = a_func(y0[[i]], dist_to_s0[[i]]) +
      as.matrix(A[[i]] %*% (b_func(y0[[i]], dist_to_s0_from_mesh[[i]]) * z[, index])) +
      rnorm(length(y0[[i]]) * length(dist_to_s0[[i]]), sd = tau^-.5)
  }

  res
}

#' This is a wrapper for the inla.qsample() function for sampling from a Gaussian
#' random field with sparse precision matrix Q. The wrapper allows us to get
#' reproduceable results if the seed outside the function scope is known.
#' It also suppresses warnings and removes some unnecessary column and row names
#' @export
rnorm_spde = function(n, Q, ...) {
  res = suppressWarnings(INLA::inla.qsample(n, Q, seed = round(runif(1, 0, 1e6)), ...))
  colnames(res) = NULL
  rownames(res) = NULL
  res
}

#' @export
simulate_spat_probit_posterior = function(hyperpar_samples,
                                          fixed_samples,
                                          coords,
                                          spde,
                                          s0_index,
                                          constr_index,
                                          create_X,
                                          create_Q,
                                          replace_zeros_at_s0 = FALSE,
                                          verbose = FALSE,
                                          n_cores = 1,
                                          n_per_sample = 1) {
  n_samples = nrow(hyperpar_samples)
  stopifnot(nrow(fixed_samples) == n_samples)

  s0 = coords[s0_index, , drop = FALSE]
  A = inla.spde.make.A(spde$mesh, coords)
  if (!is.null(constr_index)) A = A[, -constr_index]

  dist_to_s0 = as.numeric(dist_euclid(s0, coords))
  X = create_X(dist_to_s0, df = FALSE)

  simulations = local({
    RNGkind("L'Ecuyer-CMRG")
    parallel::mclapply(
      X = seq_len(n_samples),
      mc.cores = n_cores,
      FUN = function(i) {
        hyperpar = hyperpar_samples[i, ]
        fixed_par = fixed_samples[i, ]

        Q = create_Q(hyperpar, spde)
        if (!is.null(constr_index)) Q = Q[-constr_index, -constr_index]

        mean_trend = as.numeric(X %*% fixed_par)

        simulations = matrix(NA_real_, nrow = nrow(A), ncol = n_per_sample)
        bad_index = seq_len(n_per_sample)
        while (TRUE) {
          simulations[, bad_index] = as.matrix(A %*% rnorm_spde(length(bad_index), Q))
          p = pnorm(simulations[, bad_index] + mean_trend)
          simulations[, bad_index] = as.integer(rbinom(length(p), 1, p))
          bad_index = which(simulations[s0_index, ] == 0)
          if (!replace_zeros_at_s0 || length(bad_index) == 0) break
        }

        if (verbose && (i %% 10 == 0)) {
          message("Finished sample nr. ", i, " / ", n_samples)
        }

        simulations
      })
  })

  list(
    simulations = do.call(cbind, simulations),
    coords = coords,
    dist = dist_to_s0)
}

#' @export
simulate_probit_posterior = function(fixed_samples,
                                     coords,
                                     s0_index,
                                     create_X,
                                     replace_zeros_at_s0 = FALSE,
                                     verbose = FALSE,
                                     n_cores = 1,
                                     n_per_sample = 1) {
  n_samples = nrow(fixed_samples)

  s0 = coords[s0_index, , drop = FALSE]

  dist_to_s0 = as.numeric(dist_euclid(s0, coords))
  X = create_X(dist_to_s0, df = FALSE)

  simulations = local({
    RNGkind("L'Ecuyer-CMRG")
    parallel::mclapply(
      X = seq_len(n_samples),
      mc.cores = n_cores,
      FUN = function(i) {
        fixed_par = fixed_samples[i, ]

        mean_trend = as.numeric(X %*% fixed_par)

        simulations = matrix(NA_real_, nrow = nrow(coords), ncol = n_per_sample)
        bad_index = seq_len(n_per_sample)
        while (TRUE) {
          p = pnorm(rep(mean_trend, length(bad_index)))
          simulations[, bad_index] = as.integer(rbinom(length(p), 1, p))
          bad_index = which(simulations[s0_index, ] == 0)
          if (!replace_zeros_at_s0 || length(bad_index) == 0) break
        }

        if (verbose && (i %% 10 == 0)) {
          message("Finished sample nr. ", i, " / ", n_samples)
        }

        simulations
      })
  })

  list(
    simulations = do.call(cbind, simulations),
    coords = coords,
    dist = dist_to_s0)
}

