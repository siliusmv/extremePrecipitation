devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(INLA)
library(inlabru)
library(patchwork)
library(lubridate)
library(RhpcBLASctl)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

threshold = qlaplace(.9)

gp_filename = file.path(results_dir(), "gp_model.rds")
gamma_filename = file.path(results_dir(), "gamma_model.rds")
gp_res = readRDS(gp_filename)
gamma_res = readRDS(gamma_filename)

filename = file.path(results_dir(), "intensity_process.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
heights = radar$coords$height
n_loc = nrow(coords)
n_time = length(radar$times)
radar$data = t(radar$data) # We will mostly extract data for all locations at different time points

# Remove the Rissa radar and its closest neighbours
rissa = radar$rissa
dist_to_rissa = as.numeric(st_distance(rissa, radar$coords))
bad_radius = 5
bad_index = which(dist_to_rissa <= bad_radius)
if (length(bad_index) > 0) {
  radar$coords = radar$coords[-bad_index, ]
  radar$data = radar$data[-bad_index, ]
  coords = coords[-bad_index, ]
  heights = radar$coords$height
  n_loc = nrow(coords)
  n_time = nrow(radar$data)
}

# ==============================================================================
# Transform the marginal distributions
# ==============================================================================

# First remove the zeros (X = [Y | Y > 0]). Then, set
# P(X <= x) = P_gamma(X <= x) if x <= u
# P(X <= x) = p + (1 - p) P_gp(X <= x | X > u)

zero_threshold = gamma_res$zero_threshold
radar$data[radar$data <= zero_threshold] = NA
# Shift the data to correctly align with the fitted marginals
radar$data = radar$data - zero_threshold

pb = progress_bar(n_time)
for (i in seq_len(nrow(gamma_res$b))) {
  index = which(radar$day == gamma_res$u$day[i] & radar$year == gamma_res$u$year[i])
  if (any(index)) {
    radar$data[, index] = transform(
      x = radar$data[, index],
      a = gamma_res$a,
      b = gamma_res$b$b[i],
      u = gamma_res$u$mean[i],
      u_prob = gamma_res$prob,
      xi = gp_res$xi,
      s = gp_res$s$mean[i],
      alpha = gp_res$prob)
    pb$tick(length(index))
  }
}
pb$terminate()

if (FALSE) {
  hist(as.numeric(radar$data), breaks = seq(-20, 20, by = .1))
  summary(as.numeric(radar$data))
  # This looks almost perfect. And we don't really care that much
  # about the bump all the way to the left...
}

# ==============================================================================
# Examine empirical conditional moments in the data
# ==============================================================================

delta_s0 = 1
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = 1,
  n = 1,
  r = Inf)

summary(data$n)
sum(data$n)

dist_radius = .5
log_y0_radius = .025

log_y0_centers = seq(
  from = ceiling(log(threshold) / (log_y0_radius * 2)) * (log_y0_radius * 2),
  to = 2.2,
  by = log_y0_radius * 2)
dist_centers = seq(1, 70, by = dist_radius * 2)

moments = empirical_moments(
  data = data,
  dist_centers = dist_centers,
  log_y0_centers = log_y0_centers,
  dist_radius = dist_radius,
  log_y0_radius = log_y0_radius)
chi = empirical_chi(
  data = data,
  thresholds = exp(log_y0_centers),
  dist_centers = dist_centers,
  dist_radius = dist_radius)

rm(data)
gc()

moments$alpha = moments$mean / moments$y0
moments = moments |>
  dplyr::group_by(dist) |>
  dplyr::mutate(beta = log(sd / sd[1]) / log(y0 / y0[1])) |>
  dplyr::ungroup()
moments$sigma = moments$sd / (moments$y0 ^ moments$beta)
moments$sigma[moments$dist == 0] = 0
attr(moments, "dist_radius") = dist_radius
attr(moments, "log_y0_radius") = log_y0_radius
attr(chi, "dist_radius") = dist_radius

local({
  tmp = readRDS(filename)
  tmp$threshold = threshold
  tmp$moments = moments
  tmp$chi = chi
  saveRDS(tmp, filename)
})

plot = local({
  moments = readRDS(filename)$moments
  chi = readRDS(filename)$chi
  df1 = moments |>
    tidyr::pivot_longer(c(mean, sd, alpha, beta, sigma)) |>
    dplyr::select(-log_y0)
  df2 = chi |>
    dplyr::rename(value = chi, y0 = threshold) |>
    dplyr::mutate(name = "chi")
  plot_df = rbind(df1, df2) |>
    dplyr::select(-n) |>
    dplyr::filter(dist <= 70, log(y0) > min(log(y0)) + .1) |>
    dplyr::mutate(name = factor(
      name,
      levels = c("mean", "sd", "alpha", "beta", "sigma", "chi"),
      labels = paste0(
        "$\\hat \\",
        c("mu", "zeta", "alpha", "beta", "sigma", "chi"),
        "(d; y_0",
        c(")$", ", y_1)$")[c(1, 1, 1, 2, 2, 1)])))
  lims_df = plot_df |>
    dplyr::group_by(name) |>
    dplyr::reframe(
      dist = c(min(dist), max(dist)),
      value = c(0, max(value, na.rm = TRUE)))
  ggplot(plot_df) +
    geom_line(aes(x = dist, y = value, group = y0, col = y0)) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "$d$", y = "Value", col = "$y_0$") +
    geom_blank(data = lims_df, aes(x = dist, y = value)) +
    theme_light() +
    theme(
      text = element_text(size = 18),
      strip.text = element_text(colour = "black", size = 20),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
})

plot_tikz(
  plot,
  file = file.path(image_dir(), "empirical_conditional.pdf"),
  width = 10,
  height = 6)

local({
  tmp = readRDS(filename)
  tmp$empirical_lims = plot$layers[[2]]$data
  saveRDS(tmp, filename)
})

# ==============================================================================
# Minimise MSE to estimate parameters of α for each value of y0
# ==============================================================================

delta_s0 = 4
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = 1,
  n = 2,
  r = Inf)
summary(data$n)
sum(data$n)

f = function(θ, y, d, y0) {
  λ = exp(θ[1])
  κ = exp(θ[2])
  a = y0 * exp(-(d / λ)^κ)
  sum((a - y)^2)
}
a_pars = list()
max_y0 = 11
max_dist = 90
y0_rounding = .2
y0_steps = seq(round(threshold / y0_rounding) * y0_rounding, max_y0, by = y0_rounding)
pb = progress_bar(length(data$n))
df = lapply(
  X = seq_along(data$n),
  FUN = function(i) {
    ii = which(data$dist_to_s0[[i]] < max_dist)
    if (length(ii) == 0) return(NULL)
    pb$tick()
    data.frame(
      y = as.numeric(data$y[[i]][ii, ]),
      dist = rep(data$dist_to_s0[[i]][ii], data$n[i]),
      y0 = rep(data$y0[[i]], each = length(ii))) |>
      dplyr::filter(!is.na(y))
  }) |>
  dplyr::bind_rows()
pb$terminate()
pb = progress_bar(length(y0_steps))
for (i in seq_along(y0_steps)) {
  tmp_data = dplyr::filter(df, abs(y0 - y0_steps[i]) < y0_rounding / 2)
  if (is.null(tmp_data)) {
    a_pars[[i]] = NULL
  } else {
    a_pars[[i]] = optim(
      par = log(c(30, .5)),
      fn = f,
      y = tmp_data$y,
      d = tmp_data$dist,
      control = list(maxit = 2000),
      y0 = tmp_data$y0)
    a_pars[[i]]$y0 = y0_steps[i]
  }
  pb$tick()
}
pb$terminate()

local({
  tmp = readRDS(filename)
  tmp$a_pars = a_pars
  saveRDS(tmp, filename)
})

rm(df)
gc()

# ==============================================================================
# Estimate MSE-parameters for α with the λ = λ0 - λ1 y0, κ = κ0 - κ1 y0 model
# ==============================================================================

f2 = function(θ, y, d, y0) {
  λ0 = exp(θ[1])
  λ_λ = exp(θ[2])
  κ0 = exp(θ[3])
  λ_κ = exp(θ[4])
  κ_κ = exp(θ[5])
  λ = λ0 * exp(-(y0 - threshold) / λ_λ)
  κ = κ0 * exp(-((y0 - threshold) / λ_κ)^κ_κ)
  a = y0 * exp(-(d / λ)^κ)
  sum((a - y)^2)
}
df = lapply(
  X = seq_along(data$n),
  FUN = function(i) {
    ii = which(data$dist_to_s0[[i]] < 70)
    if (length(ii) == 0) return(NULL)
    data.frame(
      y = as.numeric(data$y[[i]][ii, ]),
      dist = rep(data$dist_to_s0[[i]][ii], data$n[i]),
      y0 = rep(data$y0[[i]], each = length(ii))) |>
      dplyr::filter(!is.na(y))
  }) |>
  dplyr::bind_rows()
a_pars2 = optim(
  par = log(c(45, 4, .8, 5, 2)),
  fn = f2,
  y = df$y,
  d = df$dist,
  control = list(trace = 6),
  y0 = df$y0)

rm(df)
gc()

local({
  tmp = readRDS(filename)
  tmp$a_pars2 = a_pars2
  saveRDS(tmp, filename)
})

tmp = readRDS(filename)
moments = tmp$moments
a_pars2 = tmp$a_pars2
a_pars = tmp$a_pars
rm(tmp)

plot = local({
  max_lambda = 50
  tmp = lapply(a_pars, \(x) exp(x$par)) |>
    do.call(what = rbind) |>
    as.data.frame() |>
    dplyr::rename(lambda = V1, kappa = V2) |>
    dplyr::mutate(y0 = sapply(a_pars, `[[`, "y0")) |>
    dplyr::filter(lambda <= max_lambda, y0 > threshold) |>
    tidyr::pivot_longer(-y0)
  tmp2 = data.frame(y0 = sapply(a_pars, `[[`, "y0")) |>
    dplyr::mutate(
      lambda = exp(a_pars2$par[1]) * exp(-(y0 - threshold) / exp(a_pars2$par[2])),
      kappa = exp(a_pars2$par[3]) * exp(-((y0 - threshold) / exp(a_pars2$par[4]))^exp(a_pars2$par[5]))) |>
    tidyr::pivot_longer(-y0) |>
    dplyr::filter(y0 > threshold)
  tmp$name = factor(
    tmp$name,
    levels = c("lambda", "kappa"),
    labels = c("$\\hat \\lambda_a(y_0)$", "$\\hat \\kappa_a(y_0)$"))
  tmp2$name = factor(
    tmp2$name,
    levels = c("lambda", "kappa"),
    labels = c("$\\hat \\lambda_a(y_0)$", "$\\hat \\kappa_a(y_0)$"))
  plot1 = ggplot() +
    geom_point(data = tmp, aes(x = y0, y = value)) +
    geom_line(data = tmp2, aes(x = y0, y = value)) +
    facet_wrap(~name, scales = "free") +
    labs(x = "$y_0$", y = "Value")
  tmp3 = expand.grid(
    dist = unique(moments$dist),
    y0 = unique(moments$y0)) |>
    as.data.frame() |>
    dplyr::filter(y0 > threshold) |>
    dplyr::mutate(
      lambda = exp(a_pars2$par[1]) * exp(-(y0 - threshold) / exp(a_pars2$par[2])),
      kappa = exp(a_pars2$par[3]) * exp(-((y0 - threshold) / exp(a_pars2$par[4]))^exp(a_pars2$par[5])),
      alpha = exp(-(dist / lambda)^kappa),
      a = alpha * y0) |>
    tidyr::pivot_longer(c(a, alpha)) |>
    dplyr::filter(dist <= 70) |>
    dplyr::mutate(name = factor(
      name,
      levels = c("a", "alpha"),
      labels = c("$\\hat a(d; y_0)$", "$\\hat \\alpha(d; y_0)$")))
  plot2 = tmp3 |>
    ggplot() +
    geom_line(aes(x = dist, y = value, group = y0, col = y0)) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "$d$", col = "$y_0$", y = "Value")
  plot = patchwork::wrap_plots(plot1, plot2, design = "AAAA#BBBB") *
    theme_light() *
    theme(
      text = element_text(size = 18),
      strip.text = element_text(colour = "black", size = 20),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
  plot
})

plot_tikz(
  plot = plot,
  file = file.path(image_dir(), "a_mse_fit.pdf"),
  width = 11,
  height = 4)

# ==============================================================================
# Create and save functions describing how we create our model
# ==============================================================================
get_tau = function(theta) {
  exp(theta[1])
}
get_Q = function(theta, spde) {
  log_rho = theta[2]
  log_sigma = theta[3]
  inla.spde2.precision(spde, c(log_rho, log_sigma))
}
get_a_func = local({
  function(theta) {
    lambda0 = exp(tail(theta, 5)[1])
    lambda_lambda = exp(tail(theta, 4)[1])
    kappa0 = exp(tail(theta, 3)[1])
    lambda_kappa = exp(tail(theta, 2)[1])
    kappa_kappa = exp(tail(theta, 1))
    function(y, dist) {
      lambda = lambda0 * exp(-rep((y - threshold), each = length(dist)) / lambda_lambda)
      kappa = kappa0 * exp(-(rep((y - threshold), each = length(dist)) / lambda_kappa)^kappa_kappa)
      alpha = exp(-(rep(dist, length(y)) / lambda)^kappa)
      a = rep(y, each = length(dist)) * alpha
      matrix(a, nrow = length(dist), ncol = length(y))
    }
  }
})
environment(get_a_func)$threshold = threshold
get_b_func = function(theta) {
  beta0 = exp(theta[4])
  lambda_b = exp(theta[5])
  kappa_b = exp(theta[6])
  function(y, dist) {
    beta = beta0 * exp(-(rep(dist, length(y)) / lambda_b)^kappa_b)
    b = rep(y, each = length(dist)) ^ beta
    matrix(b, nrow = length(dist), ncol = length(y))
  }
}
get_b_func = function(theta) {
  beta0 = exp(theta[4]) / (1 + exp(theta[4]))
  lambda_b = exp(theta[5])
  kappa_b = exp(theta[6])
  function(y, dist) {
    beta = beta0 * exp(-(rep(dist, length(y)) / lambda_b)^kappa_b)
    b = rep(y, each = length(dist)) ^ beta
    matrix(b, nrow = length(dist), ncol = length(y))
  }
}


local({
  tmp = readRDS(filename)
  tmp$get_b_func = get_b_func
  tmp$get_Q = get_Q
  tmp$get_tau = get_tau
  tmp$get_a_func = get_a_func
  saveRDS(tmp, filename)
})

# ==============================================================================
# Prepare data for inference with R-INLA
# ==============================================================================
s0 = radar$s0 |>
  st_coordinates()
s0_index = sapply(
  X = seq_len(nrow(s0)),
  FUN = function(i) {
    which(abs(coords[, 1] - s0[i, 1]) < .01 & abs(coords[, 2] - s0[i, 2]) < .01)
  }) |>
  as.numeric()

data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n_cores = 1,
  n = 1,
  r = Inf)

# Create a triangulated mesh for every of our chosen conditioning sites.
# The function INLA::inla.mesh.2d() creates triangulated meshes, but it will
# fail to create a mesh and run indefinitely for certain combinations of input arguments.
# When creating many different meshes, it is really hard to ensure that INLA::inla.mesh.2d()
# will work at a first try for all of them. Therefore, we need to create the meshes
# using a while loop that slightly changes the input arguments if the mesh creation fails
multimesh_data = parallel::mclapply(
  X = seq_along(data$obs_index),
  mc.cores = 1,
  mc.preschedule = FALSE,
  FUN = function(i) {
    convex = 40 # The original convexity used for creating the mesh boundary
    while (TRUE) {
      mesh_loc = extract_thinned_out_circles(
        coords = coords,
        center = data$s0[[i]],
        n = 4,
        r = Inf)
      # The original mesh arguments
      args = list(
        #loc = coords[c(data$obs_index[[i]], data$s0_index[[i]]), ],
        loc = rbind(mesh_loc, data$s0[[i]]),
        boundary = inla.nonconvex.hull(coords[data$obs_index[[i]], ], convex = convex),
        max.edge = 50)
      # Try to create a mesh using a function that automatically fails if the mesh
      # creation takes more than `timeout` seconds
      mesh = create_mesh(args)
      if (!is.null(mesh)) break # Break out of the loop if this succeeded
      convex = convex + 5 # If not, increase the convexity range and try again
    }
    # Create an SPDE object, a projection matrix and the distance from all mesh nodes to
    # the conditioning site for the created mesh
    spde = inla.spde2.pcmatern(
      mesh,
      prior.range = c(60, .95),
      prior.sigma = c(4, .05),
      alpha = 1.5)
    dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
    A = inla.spde.make.A(mesh, coords[data$obs_index[[i]], ])
    message("Completed mesh nr. ", i, " of ", length(data$obs_index))
    list(spde = spde, dist = dist_to_s0_from_mesh, A = A)
  })
multimesh_data = purrr::transpose(multimesh_data)

# ==============================================================================
# Perform full inference using INLA
# ==============================================================================

start_total_time = proc.time()
n_s0 = length(data$n)
for (my_index in seq_len(n_s0)) {
  local({
    start_iter_time = proc.time()

    data = lapply(data, `[`, my_index)
    multimesh_data = lapply(multimesh_data, `[`, my_index)
    df = lapply(
      X = seq_along(data$n),
      FUN = function(i) {
        data.frame(
          y = as.numeric(data$y[[i]]),
          y0 = rep(data$y0[[i]], each = nrow(data$y[[i]])),
          dist = rep(data$dist_to_s0[[i]], data$n[i]))
      }) |>
      dplyr::bind_rows()
    na_index = which(is.na(df$y))
    df = df[-na_index, ]

    beta_par = c(.65, 8.5, .5)
    beta_par[-1] = log(beta_par[-1])
    beta_par[1] = log(beta_par[1]) - log(1 - beta_par[1])

    spde_priors = list(
      rho = c(60, .95),
      sigma = c(4, .05),
      beta0 = c(beta_par[1], 5),
      lambda = c(beta_par[2], 5),
      kappa = c(beta_par[3], 5))
    b_model = spde_b_model(
      n = data$n,
      y0 = unlist(data$y0),
      spde = multimesh_data$spde,
      #init = c(log(40), log(1.3), beta_par),
      init = c(log(40), log(1.2), beta_par),
      priors = spde_priors,
      dist_to_s0 = multimesh_data$dist)

    a_par = readRDS(filename)$a_pars2
    #a_par = list(par = c(3.6, 1.4, -.1, 2.3, 2.7))

    a_priors = list(
      lambda0 = c(a_par$par[1], 5),
      lambda_lambda = c(a_par$par[2], 5),
      kappa0 = c(a_par$par[3], 5),
      lambda_kappa = c(a_par$par[4], 5),
      kappa_kappa = c(a_par$par[5], 5))
    a_model = a_model(
      y0 = df$y0,
      dist_to_s0 = df$dist,
      threshold = threshold,
      init = a_par$par,
      priors = a_priors)

    A_spatial = local({
      A = multimesh_data$A
      for (i in seq_along(data$n)) {
        ii = which(multimesh_data$dist[[i]] == 0)
        A[[i]] = A[[i]][, -ii]
        A[[i]] = Matrix::bdiag(rep(list(A[[i]]), data$n[i]))
      }
      A = Matrix::bdiag(A)
      A[-na_index, ]
    })
    identical(dim(A_spatial), c(nrow(df), b_model$f$n))

    stack = inla.stack(
      data = list(y = df$y),
      A = list(
        spatial = A_spatial,
        a = 1
      ),
      effects = list(
        spatial = seq_len(ncol(A_spatial)),
        a = seq_along(df$y)
      ))

    formula = y ~ -1 +
      f(spatial, model = b_model) +
      f(a, model = a_model)

    message("Start running INLA for iter nr. ", my_index, " / ", n_s0)
    fit = inla(
      formula = formula,
      family = "gaussian",
      data = inla.stack.data(stack),
      control.compute = list(config = TRUE),
      control.family = list(hyper = list(prec = list(
        initial = 2,
        prior = "pc.prec",
        param = c(1, .01)))),
      control.inla = list(int.strategy = "eb"),
      inla.mode = "experimental",
      control.predictor = list(A = inla.stack.A(stack)),
      verbose = TRUE,
      num.threads = 1)

    set.seed(1)
    samples = inla.hyperpar.sample(
      n = 1e3,
      result = fit,
      improve.marginals = TRUE,
      intern = TRUE)

    res = list(
      samples = samples,
      mode = fit$mode$theta,
      s0 = radar$s0[my_index],
      cpu = fit$cpu,
      mlik = fit$mlik,
      hyperpar = fit$summary.hyperpar,
      time_index = data$time_index,
      obs_index = data$obs_index)

    tmp = readRDS(filename)
    if (is.null(tmp[["inla"]])) tmp[["inla"]] = list()
    tmp[["inla"]][[my_index]] = res
    saveRDS(tmp, filename)

    end_iter_time = proc.time()
    iter_time = (end_iter_time - start_iter_time)[3]
    total_time = (end_iter_time - start_total_time)[3]

    message(
      "Done with iter nr. ", my_index, " / ", n_s0, ".\n",
      "Time spent in this iteration: ", lubridate::seconds_to_period(iter_time), ".\n",
      "Time spent in total: ", lubridate::seconds_to_period(total_time), ".\n")

  })
  gc()
}

# ==============================================================================
# Compute conditional moments and other properties of the model fits
# ==============================================================================

moments = readRDS(filename)$moments

s0 = radar$s0 |>
  st_coordinates()
s0_index = sapply(
  X = seq_len(nrow(s0)),
  FUN = function(i) {
    which(abs(coords[, 1] - s0[i, 1]) < .01 & abs(coords[, 2] - s0[i, 2]) < .01)
  }) |>
  as.numeric()

pb = progress_bar(length(s0_index))
model_properties = list()
fits = readRDS(filename)
for (i in seq_along(s0_index)) {
#for (i in 1:2) {
  fit = fits$inla[[i]]
  matern_corr = function(dist, rho, nu = 1.5) {
    kappa = sqrt(8 * nu) / rho
    res = 2 ^ (1 - nu) / gamma(nu) * (kappa * dist) ^ nu * besselK(kappa * dist, nu)
    res[dist == 0] = 1
    res
  }
  n_samples = nrow(fit$samples)
  model_properties[[i]] = local({
    d = c(0, unique(moments$dist))
    y0 = unique(moments$y0)
    res = list(
      arrays = rep(list(array(0, dim = c(length(d), length(y0), n_samples))), 5),
      mats = rep(list(matrix(0, nrow = length(d), ncol = n_samples)), 1))
    names(res$arrays) = c("a", "alpha", "beta", "zeta", "chi")
    names(res$mats) = "sd"
    for (j in seq_len(nrow(fit$samples))) {
      res$mats$sd[, j] = exp(fit$samples[j, 3]) *
        sqrt(1 - matern_corr(d, exp(fit$samples[j, 2]), .5)^2)
      res$arrays$a[, , j] = fits$get_a_func(fit$samples[j, ])(y0, d)
      res$arrays$alpha[, , j] = res$arrays$a[, , j] / rep(y0, each = length(d))
      res$arrays$beta[, , j] = log(fits$get_b_func(fit$samples[j, ])(exp(1), d))
      res$arrays$zeta[, , j] = sqrt(
        fits$get_b_func(fit$samples[j, ])(y0, d)^2 * rep(res$mats$sd[, j]^2, length(y0))
        + 1 / fits$get_tau(fit$samples[j, ]))
      # compute χ via integration
      delta_y0 = .05
      for (k in 0:100) {
        centers = y0 + (k - .5) * delta_y0
        a = fits$get_a_func(fit$samples[j, ])(centers, d)
        b = fits$get_b_func(fit$samples[j, ])(centers, d)
        zeta = sqrt(
          b^2 * rep(res$mats$sd[, j]^2, length(y0))
          + 1 / fits$get_tau(fit$samples[j, ]))
        heights = dexp(centers - y0) *
          pnorm(rep(y0, each = length(d)), mean = a, sd = zeta, lower.tail = FALSE)
        areas = heights * delta_y0
        res$arrays$chi[, , j] = res$arrays$chi[, , j] + areas
      }
      res$arrays$chi[1, , j] = 1
    }

    for (name in names(res$arrays)) {
      res$arrays[[name]] = data.frame(
        y0 = rep(y0, each = length(d)),
        d = rep(d, length(y0)),
        mean = as.numeric(apply(res$arrays[[name]], c(1, 2), mean)),
        lower = as.numeric(apply(res$arrays[[name]], c(1, 2), quantile, probs = .025)),
        upper = as.numeric(apply(res$arrays[[name]], c(1, 2), quantile, probs = .975))) |>
        dplyr::mutate(name = !!name)
    }
    for (name in names(res$mats)) {
      res$mats[[name]] = data.frame(
        d = d,
        mean = as.numeric(apply(res$mats[[name]], 1, mean)),
        lower = as.numeric(apply(res$mats[[name]], 1, quantile, probs = .025)),
        upper = as.numeric(apply(res$mats[[name]], 1, quantile, probs = .975))) |>
        dplyr::mutate(name = !!name, y0 = NA)
    }

    res = rbind(do.call(rbind, res$arrays), do.call(rbind, res$mats))
    row.names(res) = NULL
    res
  })
  gc()
  pb$tick()
}
pb$terminate()

local({
  tmp = readRDS(filename)
  tmp$model_properties = model_properties
  saveRDS(tmp, filename)
})

# ==============================================================================
# Plot the computed properties
# ==============================================================================

model_properties = readRDS(filename)$model_properties

plot_data = list()
for (i in seq_along(model_properties)) {
  plot_data[[i]] = model_properties[[i]] |>
    dplyr::filter(d <= 70, is.na(y0) | log(y0) > min(log(y0), na.rm = TRUE) + .1) |>
    dplyr::mutate(
      y0 = ifelse(name == "beta", NA, y0),
      lower = ifelse(name %in% c("beta", "sd"), lower, NA),
      upper = ifelse(name %in% c("beta", "sd"), upper, NA)) |>
    dplyr::mutate(name = factor(
      name,
      levels = c("a", "zeta", "alpha", "beta", "sd", "chi"),
      labels = paste0(
        "$\\",
        c("mu", "zeta", "alpha", "beta", "sigma", "chi"),
        "(d",
        c(")$", "; y_0)$")[c(2, 2, 2, 1, 1, 2)])))
}

plot_lims = list()
empirical_lims = readRDS(filename)$empirical_lims
for (i in seq_along(plot_data)) {
  plot_lims[[i]] = plot_data[[i]] |>
    dplyr::group_by(name) |>
    dplyr::reframe(
      d = c(min(d), max(d)),
      value = c(0, max(c(upper, mean), na.rm = TRUE)))
  for (j in seq_along(levels(plot_lims[[i]]$name))) {
    plot_lims[[i]]$value[j * 2] = max(plot_lims[[i]]$value[j * 2], empirical_lims$value[j * 2])
    plot_lims[[i]]$d[j * 2] = max(plot_lims[[i]]$d[j * 2], empirical_lims$dist[j * 2])
    plot_lims[[i]]$value[j * 2 - 1] = min(plot_lims[[i]]$value[j * 2 - 1], empirical_lims$value[j * 2 - 1])
    plot_lims[[i]]$d[j * 2 - 1] = min(plot_lims[[i]]$d[j * 2 - 1], empirical_lims$dist[j * 2 - 1])
  }
}
plot_lims2 = plot_lims
for (i in seq_along(plot_lims2)) {
  plot_lims2[[i]] = dplyr::mutate(
    plot_lims2[[i]],
    name = factor(name, levels = levels(plot_lims[[i]]$name), labels = levels(empirical_lims$name)))
}

obs_plot = local({
  moments = readRDS(filename)$moments
  chi = readRDS(filename)$chi
  df1 = moments |>
    tidyr::pivot_longer(c(mean, sd, alpha, beta, sigma)) |>
    dplyr::select(-log_y0)
  df2 = chi |>
    dplyr::rename(value = chi, y0 = threshold) |>
    dplyr::mutate(name = "chi")
  rbind(df1, df2) |>
    dplyr::select(-n) |>
    dplyr::filter(dist <= 70, y0 > min(y0) + .2) |>
    dplyr::mutate(name = factor(
      name,
      levels = c("mean", "sd", "alpha", "beta", "sigma", "chi"),
      labels = paste0(
        "$\\hat \\",
        c("mu", "zeta", "alpha", "beta", "sigma", "chi"),
        "(d; y_0",
        c(")$", ", y_1)$")[c(1, 1, 1, 2, 2, 1)])
    )) |>
    ggplot() +
    geom_line(aes(x = dist, y = value, group = y0, col = y0)) +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    labs(x = "$d$", y = "Value", col = "$y_0$") +
    theme_light() +
    theme(
      text = element_text(size = 18),
      strip.text = element_text(colour = "black", size = 20),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
})

plots = list()
for (i in seq_along(model_properties)) {
  plots[[i]] = local({
    tmp = ggplot(plot_data[[i]]) +
      geom_line(aes(x = d, y = mean, group = y0, col = y0)) +
      geom_ribbon(aes(x = d, ymin = lower, ymax = upper, group = y0, fill = y0), alpha = .35) +
      scale_color_continuous(na.value = "black") +
      scale_fill_continuous(na.value = "black") +
      facet_wrap(~name, scales = "free_y", nrow = 3) +
      labs(x = "$d$", y = "Value", col = "$y_0$", fill = "$y_0$") +
      theme_light() +
      theme(
        text = element_text(size = 18),
        strip.text = element_text(colour = "black", size = 20),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
      geom_blank(data = plot_lims[[i]], aes(x = d, y = value))
    obs_plot = obs_plot +
      geom_blank(data = plot_lims2[[i]], aes(x = d, y = value))
    patchwork::wrap_plots(
      obs_plot, tmp, ncol = 2, guides = "collect", design = "AAAAAAA#BBBBBBB") +
      patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")
  })
}

plot_tikz(
  plots,
  file = file.path(image_dir(), "intensity_results.pdf"),
  width = 12,
  height = 6)

obs_plot2 = obs_plot +
  facet_wrap(~name, scales = "free_y", nrow = 2)

plot = local({
  plot_lims[[1]] = plot_lims[[1]]
  plot_lims2[[1]] = plot_lims2[[1]]
  for (i in seq_along(plot_lims)[-1]) {
    for (j in seq_len(nrow(plot_lims[[1]]) / 2)) {
      plot_lims[[1]]$value[2 * j] = max(
        plot_lims[[1]]$value[2 * j],
        plot_lims[[i]]$value[2 * j])
      plot_lims[[1]]$value[2 * j - 1] = min(
        plot_lims[[1]]$value[2 * j - 1],
        plot_lims[[i]]$value[2 * j - 1])
      plot_lims2[[1]]$value[2 * j] = max(
        plot_lims2[[1]]$value[2 * j],
        plot_lims2[[i]]$value[2 * j])
      plot_lims2[[1]]$value[2 * j - 1] = min(
        plot_lims2[[1]]$value[2 * j - 1],
        plot_lims2[[i]]$value[2 * j - 1])
    }
  }
  for (i in seq_along(plots)) {
    plots[[i]] = plots[[i]][[2]] +
      labs(title = paste("Model fit nr.", i)) +
      facet_wrap(~name, scales = "free_y", nrow = 2) +
      geom_blank(data = plot_lims[[1]], aes(x = d, y = value))
  }
  obs_plot2 = obs_plot2 +
    labs(title = "Observations") +
    geom_blank(data = plot_lims2[[1]], aes(x = d, y = value))
  patchwork::wrap_plots(c(list(obs_plot2), plots), ncol = 2, guides = "collect")
})

plot_tikz(
  plot,
  file = file.path(image_dir(), "intensity_results2.pdf"),
  width = 14,
  height = 14)
