
#' Function for the threshold occurrence model. This functions takes
#' a matrix `samples` as input, together with a probability `nonzero_prob`.
#' Then, all data points in `samples` that are smaller than the `nonzero_prob`-quantile
#' of `samples` gets set equal to zero.
#' @export
threshold_occurrence = function(samples, nonzero_prob) {
  threshold = quantile(samples, probs = 1 - nonzero_prob, na.rm = TRUE)
  samples[samples < threshold] = 0
  samples
}

#' Simulate spatial realisations of the precipitation intensity process.
#'
#' The input variables are:
#' samples: an (n x m)-dimensional matrix of hyperparameter samples from the fitted
#'   spatial conditional extremes model, where n is the total number of samples and
#'   m is the number of hyperparameters in our chosen model
#' threshold: The threshold of the spatial conditional extremes model
#' coords: a (d x 2)-dimensional matrix of the coordinates for all the d locations in which
#'   we want to simulate data
#' s0_index: an integer between 1 and d, describing which of the d elements in coords that
#'   correspond to the conditioning site
#' spde: an INLA SPDE object which we use for sampling from the SPDE approximation
#' get_a_func: a function that takes in an m-dimensional vector of hyperparameters and returns
#'   a function that computes the standardising function a() of the spatial conditional
#'   extremes model
#' get_b_func: a function that takes in an m-dimensional vector of hyperparameters and returns
#'   a function that computes the standardising function b() of the spatial conditional
#'   extremes model
#' get_Q: a function that takes in an m-dimensional vector of hyperparameters and returns
#'   a function that computes the precision matrix of the SPDE approximation
#' get_tau: a function that takes in an m-dimensional vector of hyperparameters and returns
#'   a function that computes the precision of the nugget effect
#' verbose: a Boolean: Should we be verbose or not?
#' n_cores: an integer describing the number of parallel cores that should be used for
#'   running the function
#' n_per_sample: number of spatial realisations of the intensity process that we will
#'   simulate for each of the n hyperparameter samples
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

#' Simulate spatial realisations of the precipitation occurrence process using
#' the spatial probit model
#'
#' The input variables are:
#' hyperpar_samples: an (n x m)-dimensional matrix of hyperparameter samples from the fitted
#'   spatial probit model, where n is the total number of samples and
#'   m is the number of hyperparameters in our chosen model
#' fixed_samples: an (n x k)-dimensional matrix of samples for the parameters of the fixed
#'   effects of the spatial probit model, i.e., those describing the spline function in
#'   the mean structure of the latent field
#' coords: a (d x 2)-dimensional matrix of the coordinates for all the d locations in which
#'   we want to simulate data
#' spde: an INLA SPDE object which we use for sampling from the SPDE approximation
#' s0_index: an integer between 1 and d, describing which of the d elements in coords that
#'   correspond to the conditioning site
#' constr_index: an index describing which mesh nodes that should be constrained to zero
#'   in the SPDE approximation
#' create_X: a function that takes in the distance to the conditioning site and outputs
#'   a design matrix with all the basis functions for the spline function
#' create_Q: a function that takes in an m-dimensional vector of hyperparameters and an INLA spde
#'   object, and returns the precision matrix of the SPDE approximation
#' replace_zeros_at_s0: a Boolean: Should we enforce I(s0) = 1 by resampling from the model
#'   whenever we get a realisation that gives I(s0) = 0?
#' verbose: a Boolean: Should we be verbose or not?
#' n_cores: an integer describing the number of parallel cores that should be used for
#'   running the function
#' n_per_sample: number of spatial realisations of the occurrence process that we will
#'   simulate for each of the n hyperparameter samples
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

  # Compute distances to the conditioning site and build the design matrix X
  dist_to_s0 = as.numeric(dist_euclid(s0, coords))
  X = create_X(dist_to_s0, df = FALSE)

  # Simulate from the fitted occurrence model
  simulations = local({
    RNGkind("L'Ecuyer-CMRG")
    parallel::mclapply(
      X = seq_len(n_samples),
      mc.cores = n_cores,
      FUN = function(i) {
        # These are the parameters from model-parameter sample nr. i
        hyperpar = hyperpar_samples[i, ]
        fixed_par = fixed_samples[i, ]

        # Create and maybe constrain the precision matrix
        Q = create_Q(hyperpar, spde)
        if (!is.null(constr_index)) Q = Q[-constr_index, -constr_index]

        # Compute the spline function in the mean
        mean_trend = as.numeric(X %*% fixed_par)

        # Preallocate the matrix of occurrence samples
        simulations = matrix(NA_real_, nrow = nrow(A), ncol = n_per_sample)
        bad_index = seq_len(n_per_sample)
        while (TRUE) {
          # Simulate from the SPDE model, and transform linear predictor using the probit link
          simulations[, bad_index] = as.matrix(A %*% rnorm_spde(length(bad_index), Q))
          p = pnorm(simulations[, bad_index] + mean_trend)
          # Simulate occurrences using the Bernoulli distribution
          simulations[, bad_index] = as.integer(rbinom(length(p), 1, p))
          bad_index = which(simulations[s0_index, ] == 0)
          # Should we replace zeros at s0 or not?
          if (!replace_zeros_at_s0 || length(bad_index) == 0) break
        }

        if (verbose && (i %% 10 == 0)) {
          message("Finished sample nr. ", i, " / ", n_samples)
        }

        simulations
      })
  })

  # Return the results
  list(
    simulations = do.call(cbind, simulations),
    coords = coords,
    dist = dist_to_s0)
}

#' Simulate spatial realisations of the precipitation occurrence process using
#' the non-spatial probit model
#'
#' The input variables are:
#' fixed_samples: an (n x k)-dimensional matrix of samples for the parameters of the fixed
#'   effects of the spatial probit model, i.e., those describing the spline function in
#'   the mean structure of the latent field
#' coords: a (d x 2)-dimensional matrix of the coordinates for all the d locations in which
#'   we want to simulate data
#' s0_index: an integer between 1 and d, describing which of the d elements in coords that
#'   correspond to the conditioning site
#' create_X: a function that takes in the distance to the conditioning site and outputs
#'   a design matrix with all the basis functions for the spline function
#' replace_zeros_at_s0: a Boolean: Should we enforce I(s0) = 1 by resampling from the model
#'   whenever we get a realisation that gives I(s0) = 0?
#' verbose: a Boolean: Should we be verbose or not?
#' n_cores: an integer describing the number of parallel cores that should be used for
#'   running the function
#' n_per_sample: number of spatial realisations of the occurrence process that we will
#'   simulate for each of the n hyperparameter samples
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

  # Compute distances to the conditioning site and build the design matrix X
  dist_to_s0 = as.numeric(dist_euclid(s0, coords))
  X = create_X(dist_to_s0, df = FALSE)

  # Simulate from the fitted occurrence model
  simulations = local({
    RNGkind("L'Ecuyer-CMRG")
    parallel::mclapply(
      X = seq_len(n_samples),
      mc.cores = n_cores,
      FUN = function(i) {
        fixed_par = fixed_samples[i, ]

        # Compute the spline function in the mean
        mean_trend = as.numeric(X %*% fixed_par)

        # Preallocate the matrix of occurrence samples
        simulations = matrix(NA_real_, nrow = nrow(coords), ncol = n_per_sample)
        bad_index = seq_len(n_per_sample)
        while (TRUE) {
          # transform linear predictor using the probit link
          p = pnorm(rep(mean_trend, length(bad_index)))
          # Simulate occurrences using the Bernoulli distribution
          simulations[, bad_index] = as.integer(rbinom(length(p), 1, p))
          bad_index = which(simulations[s0_index, ] == 0)
          # Should we replace zeros at s0 or not?
          if (!replace_zeros_at_s0 || length(bad_index) == 0) break
        }

        if (verbose && (i %% 10 == 0)) {
          message("Finished sample nr. ", i, " / ", n_samples)
        }

        simulations
      })
  })

  # Return the results
  list(
    simulations = do.call(cbind, simulations),
    coords = coords,
    dist = dist_to_s0)
}
