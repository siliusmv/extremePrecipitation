#' @export
a_model = function(y0,
                   threshold,
                   dist_to_s0,
                   init,
                   priors,
                   is_fixed = rep(FALSE, n_theta),
                   debug = FALSE) {
  n_theta = 5
  stopifnot(length(y0) == length(dist_to_s0))
  stopifnot(all(y0 >= threshold))
  stopifnot(length(is_fixed) == length(init))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == n_theta)
  stopifnot(all(c(
    "lambda0", "lambda_lambda", "kappa0", "lambda_kappa", "kappa_kappa")[!is_fixed]
    %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))

  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "a_model"
  args$shlib = file.path(cgeneric_dir(), "a.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = length(y0)
  args$is_fixed = as.integer(is_fixed)
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$threshold = threshold
  args$init = init

  if (!is_fixed[1]) args$lambda0_prior = priors$lambda0
  if (!is_fixed[2]) args$lambda_lambda_prior = priors$lambda_lambda
  if (!is_fixed[3]) args$kappa0_prior = priors$kappa0
  if (!is_fixed[4]) args$lambda_kappa_prior = priors$lambda_kappa
  if (!is_fixed[5]) args$kappa_kappa_prior = priors$kappa_kappa

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' @export
spde_b_model = function(n,
                        y0,
                        spde,
                        init,
                        priors,
                        dist_to_s0,
                        is_fixed = rep(FALSE, n_theta),
                        debug = FALSE) {
  n_theta = 5
  if (!all(class(spde) == "list")) spde = list(spde)
  stopifnot(length(n) == length(spde))
  stopifnot(length(n) == length(dist_to_s0))
  stopifnot(sum(n) == length(y0))
  stopifnot(all(sapply(spde, `[[`, "n.spde") == sapply(dist_to_s0, length)))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == n_theta)
  stopifnot(length(is_fixed) == length(init))
  stopifnot(all(c(
    "rho", "sigma", "beta0", "lambda", "kappa")[!is_fixed] %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))
  stopifnot(all(sapply(dist_to_s0, min) == 0))

  # Locate and remove the s0_locs from dist_to_s0
  s0_location = sapply(dist_to_s0, function(x) which(x == 0))
  stopifnot(is.integer(s0_location))
  for (i in seq_along(dist_to_s0)) dist_to_s0[[i]] = dist_to_s0[[i]][-s0_location[i]]

  # Start building the argument list for the cgeneric model
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "spde_b_model"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = sum(sapply(dist_to_s0, length) * n)
  args$s0_index = rep(seq_along(n), n) - 1L # 0-indexation in C
  args$is_fixed = as.integer(is_fixed)
  args$s0_location = s0_location - 1L # 0-indexation in C
  args$init = init
  args$y0 = y0
  if (!is_fixed[1]) args$rho_prior = priors$rho
  if (!is_fixed[2]) args$sigma_prior = priors$sigma
  if (!is_fixed[3]) args$beta0_prior = priors$beta0
  if (!is_fixed[4]) args$lambda_prior = priors$lambda
  if (!is_fixed[5]) args$kappa_prior = priors$kappa

  for (i in seq_along(n)) {
    args[[paste0("dist_to_s0_", i)]] = dist_to_s0[[i]]
    args[[paste0("B0_", i)]] = spde[[i]]$param.inla$B0
    args[[paste0("B1_", i)]] = spde[[i]]$param.inla$B1
    args[[paste0("B2_", i)]] = spde[[i]]$param.inla$B2
    args[[paste0("M0_", i)]] = spde[[i]]$param.inla$M0
    args[[paste0("M1_", i)]] = spde[[i]]$param.inla$M1
    args[[paste0("M2_", i)]] = spde[[i]]$param.inla$M2
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' @export
spde_b_model2 = function(n,
                         y0,
                         spde,
                         init,
                         priors,
                         dist_to_s0,
                         is_fixed = rep(FALSE, n_theta),
                         debug = FALSE) {
  n_theta = 5
  if (!all(class(spde) == "list")) spde = list(spde)
  stopifnot(length(n) == length(spde))
  stopifnot(length(n) == length(dist_to_s0))
  stopifnot(sum(n) == length(y0))
  stopifnot(all(sapply(spde, `[[`, "n.spde") == sapply(dist_to_s0, length)))
  stopifnot(is.logical(is_fixed))
  stopifnot(length(init) == n_theta)
  stopifnot(length(is_fixed) == length(init))
  stopifnot(all(c(
    "rho", "sigma", "beta0", "lambda", "kappa")[!is_fixed] %in% names(priors)))
  stopifnot(all(sapply(priors, length) == 2))
  stopifnot(all(sapply(dist_to_s0, min) == 0))

  # Locate and remove the s0_locs from dist_to_s0
  s0_location = sapply(dist_to_s0, function(x) which(x == 0))
  stopifnot(is.integer(s0_location))
  for (i in seq_along(dist_to_s0)) dist_to_s0[[i]] = dist_to_s0[[i]][-s0_location[i]]

  # Start building the argument list for the cgeneric model
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "spde_b_model2"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  # Put all the arguments into args in the order defined in the c-function
  args$n = sum(sapply(dist_to_s0, length) * n)
  args$s0_index = rep(seq_along(n), n) - 1L # 0-indexation in C
  args$is_fixed = as.integer(is_fixed)
  args$s0_location = s0_location - 1L # 0-indexation in C
  args$init = init
  args$y0 = y0
  if (!is_fixed[1]) args$rho_prior = priors$rho
  if (!is_fixed[2]) args$sigma_prior = priors$sigma
  if (!is_fixed[3]) args$beta0_prior = priors$beta0
  if (!is_fixed[4]) args$lambda_prior = priors$lambda
  if (!is_fixed[5]) args$kappa_prior = priors$kappa

  for (i in seq_along(n)) {
    args[[paste0("dist_to_s0_", i)]] = dist_to_s0[[i]]
    args[[paste0("B0_", i)]] = spde[[i]]$param.inla$B0
    args[[paste0("B1_", i)]] = spde[[i]]$param.inla$B1
    args[[paste0("B2_", i)]] = spde[[i]]$param.inla$B2
    args[[paste0("M0_", i)]] = spde[[i]]$param.inla$M0
    args[[paste0("M1_", i)]] = spde[[i]]$param.inla$M1
    args[[paste0("M2_", i)]] = spde[[i]]$param.inla$M2
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}



#' spdes: a list of M spdes
#' priors: a list of PC priors for the range and standard deviation parameters
#' init: initial values for the range and standard deviation parameters
#' mesh_index: a vector of length N containing integers between 1 and M, stating
#'   which mesh design is used for which time point.
#' constr_index: either NULL or a list of length M, with constraining indices
#'   for each mesh. If an element of the list is a negative number, then we
#'   don't do any constraining at the corresponding mesh
#' is_fixed: a vector of bools stating if the parameters should be estimated or
#'   fixed equal to their initial values
#' debug: a bool stating if INLA should print debug information
#' @export
multimesh_spde_model = function(spdes,
                                priors,
                                init,
                                constr_index = NULL,
                                mesh_index = seq_along(spdes),
                                is_fixed = rep(FALSE, 2),
                                debug = FALSE) {
  stopifnot(length(constr_index) %in% c(0, length(spdes)))
  stopifnot(length(init) == 2, length(is_fixed) == 2, length(priors) == 2)
  stopifnot(c("rho", "sigma") %in% names(priors))
  stopifnot(all(mesh_index >= 1 & mesh_index <= length(spdes)))
  stopifnot(all(mesh_index == as.integer(mesh_index)))
  stopifnot(all(sapply(priors, length) == 2))
  stopifnot(is.logical(is_fixed))
  if (length(constr_index) > 0) {
    stopifnot(all(unlist(constr_index) == as.integer(unlist(constr_index))))
  }
  for (i in seq_along(constr_index)) {
    stopifnot(all(constr_index[[i]] > 0 & constr_index[[i]] <= spdes[[i]]$n.spde))
  }

  args = list(debug = debug)
  # The name and location of the required c-function for defining the cgeneric model
  args$model = "spde_model"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  n_per_mesh = sapply(spdes, `[[`, "n.spde")
  if (length(constr_index) > 0) {
    n_per_mesh = n_per_mesh - sapply(constr_index, function(x) sum(x > 0))
  }

  # Put all the arguments into args in the order defined in the c-function
  args$n = sum(n_per_mesh[mesh_index])
  args$mesh_index = as.integer(mesh_index) - 1L # C uses 0-indexing
  args$is_fixed = as.integer(is_fixed)
  args$init = init
  args$rho_prior = priors$rho
  args$sigma_prior = priors$sigma

  for (i in seq_along(constr_index)) {
    args[[paste0("constr_", i)]] = as.integer(constr_index[[i]]) - 1L # C uses 0-indexing
  }

  for (i in seq_along(spdes)) {
    args[[paste0("B0_", i)]] = spdes[[i]]$param.inla$B0
    args[[paste0("B1_", i)]] = spdes[[i]]$param.inla$B1
    args[[paste0("B2_", i)]] = spdes[[i]]$param.inla$B2
    args[[paste0("M0_", i)]] = spdes[[i]]$param.inla$M0
    args[[paste0("M1_", i)]] = spdes[[i]]$param.inla$M1
    args[[paste0("M2_", i)]] = spdes[[i]]$param.inla$M2
  }

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}
