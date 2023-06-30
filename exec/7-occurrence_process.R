devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(INLA)
library(lubridate)
library(patchwork)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

gp_filename = file.path(results_dir(), "gp_model.rds")
gamma_filename = file.path(results_dir(), "gamma_model.rds")
intensity_filename = file.path(results_dir(), "intensity_process.rds")

filename = file.path(results_dir(), "occurrence_process.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

n_cores = 5

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# ==============================================================================
# Load the data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
gp_res = readRDS(gp_filename)
gamma_res = readRDS(gamma_filename)
intensity_res = readRDS(intensity_filename)

threshold = intensity_res$threshold
zero_threshold = gamma_res$zero_threshold
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
  n_loc = nrow(coords)
  n_time = nrow(radar$data)
}

# Locate neighbour_indices
neighbour_radius = 1
neighbour_indices = lapply(
  X = 1:nrow(coords),
  FUN = function(i) {
    as.numeric(which(
      0 < dist_euclid(coords[i, ], coords) &
      dist_euclid(coords[i, ], coords) <= neighbour_radius))
  })

# ==============================================================================
# Transform the threshold to the precipitation scale
# ==============================================================================
precipitation_thresholds = rep(NA_real_, length(radar$times))
for (i in seq_along(gamma_res$b$b)) {
  index = which(gamma_res$b$day[i] == radar$day & gamma_res$b$year[i] == radar$year)
  precipitation_thresholds[index] = inv_transform(
      x = threshold,
      a = gamma_res$a,
      b = gamma_res$b$b[i],
      u = gamma_res$u$mean[i],
      u_prob = gamma_res$prob,
      xi = gp_res$xi,
      s = gp_res$s$mean[i],
      alpha = gp_res$prob)
}
precipitation_thresholds = precipitation_thresholds + zero_threshold

# ==============================================================================
# Extract interesting data
# ==============================================================================

delta_s0 = 4
s0_index = get_s0_index(coords, delta_s0)

data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_index,
  threshold = precipitation_thresholds,
  n_cores = 1,
  n = c(1, 2, 3),
  remove_y0_from_y = FALSE,
  r = c(16, 48, Inf))

data$neighbour_mean = parallel::mclapply(
  X = seq_along(data$n),
  mc.cores = 10,
  mc.preschedule = FALSE,
  FUN = function(i) {
    res = matrix(NA_real_, nrow = data$n_loc[i], ncol = data$n[i])
    for (l in seq_len(data$n_loc[i])) {
        res[l, ] = apply(radar$data[
          neighbour_indices[[data$obs_index[[i]][l]]],
          data$time_index[[i]]], 2, mean)
    }
    message(i, " / ", length(data$n))
    res
  })

data$y0_laplace = list()
for (i in seq_along(data$n)) {
  data$y0_laplace[[i]] = rep(NA_real_, data$n[i])
  for (j in seq_len(data$n[i])) {
    index = which(
      gamma_res$b$day == radar$day[data$time_index[[i]][j]] &
      gamma_res$b$year == radar$year[data$time_index[[i]][j]])
    data$y0_laplace[[i]][j] = transform(
      x = data$y0[[i]][j] - zero_threshold,
      a = gamma_res$a,
      b = gamma_res$b$b[index],
      u = gamma_res$u$mean[index],
      u_prob = gamma_res$prob,
      xi = gp_res$xi,
      s = gp_res$s$mean[index],
      alpha = gp_res$prob)
  }
}

# ==============================================================================
# Perform some exploratory analysis
# ==============================================================================

df = data.frame(
  I = as.numeric(unlist(data$y) > zero_threshold),
  dist = unlist(rep(data$dist_to_s0, data$n)),
  neighbours = as.numeric(unlist(data$neighbour_mean)),
  y0 = unlist(lapply(
    X = seq_along(data$y0_laplace),
    FUN = function(i) rep(data$y0_laplace[[i]], each = data$n_loc[i]))))
na_index = which(is.na(df$I))
if (length(na_index) > 0) df = df[-na_index, ]

plot1 = df |>
  dplyr::mutate(
    y0 = floor(y0 / .2) * .2,
    dist = round(dist)) |>
  dplyr::group_by(dist) |>
  dplyr::summarise(p = mean(I)) |>
  dplyr::filter(dist < 65) |>
  ggplot() +
  geom_line(aes(x = dist, y = p)) +
  labs(x = "$d$", y = "$\\hat p(d)$") +
  theme_light() +
  theme(axis.title.y = element_text(vjust = .5, angle = 0))

pseudo_log = function(x) asinh(x)

plot2 = df |>
  dplyr::mutate(
    neighbours = round(neighbours, 1)) |>
  dplyr::group_by(neighbours) |>
  dplyr::summarise(p = mean(I)) |>
  ggplot() +
  geom_line(aes(x = pseudo_log(neighbours), y = p)) +
  labs(x = "$\\bar y$", y = "$\\hat p(\\bar y)$", col = "$y_0$") +
  scale_x_continuous(breaks = pseudo_log(c(0, exp(0:5))), labels = c("0", "1", paste0("e$^", 1:5, "$"))) +
  theme_light() +
  theme(axis.title.y = element_text(vjust = .5, angle = 0))

s0 = radar$s0 |>
  st_coordinates()
s0_index = sapply(
  X = seq_len(nrow(s0)),
  FUN = function(i) {
    which(abs(coords[, 1] - s0[i, 1]) < .01 & abs(coords[, 2] - s0[i, 2]) < .01)
  }) |>
  as.numeric()

my_s0_index = s0_index[4]
time_index = which(radar$data[my_s0_index, ] > threshold)
y0 = radar$data[my_s0_index, time_index]
plot_data = list()
set.seed(123)
for (i in 1:6) {
  j = sample(seq_along(y0)[y0 > i + 1 & y0 < i + 2], 1)
  plot_data[[i]] = data.frame(
    I = radar$data[, time_index[j]] > zero_threshold,
    y = ifelse(radar$data[, time_index[j]] > zero_threshold, radar$data[, time_index[j]], NA),
    x_coord = coords[, 1],
    y_coord = coords[, 2],
    y0 = y0[j],
    i = i)
}
plot_data = do.call(rbind, plot_data)

plot3 = plot_data |>
  ggplot() +
  geom_raster(aes(x = x_coord, y = y_coord, fill = pseudo_log(y))) +
  scale_fill_viridis_c(breaks = pseudo_log(c(0, exp(0:5))), labels = c("0", "1", paste0("e$^", 1:5, "$"))) +
  geom_sf(data = radar$s0[4], size = 1, col = "red") +
  facet_wrap(~i) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Easting", y = "Northing", fill = "mm/h")
plot3 = latex_friendly_map_plot(plot3)

plot = patchwork::wrap_plots(plot1, plot3, plot2, nrow = 1, widths = c(.2, .6, .2)) +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

plot_tikz(
  plot,
  file = file.path(image_dir(), "zeros_exploratory.pdf"),
  width = 11,
  height = 4)

rm(df, data)
gc()

# ==============================================================================
# Perform probit modelling with INLA
# ==============================================================================

s0 = radar$s0 |>
  st_coordinates()
s0_index = sapply(
  X = seq_len(nrow(s0)),
  FUN = function(i) {
    which(abs(coords[, 1] - s0[i, 1]) < .01 & abs(coords[, 2] - s0[i, 2]) < .01)
  }) |>
  as.numeric()

start_total_time = proc.time()
for (my_index in seq_along(s0_index)) {
  start_iter_time = proc.time()

  data = extract_extreme_fields(
    data = radar$data,
    coords = coords,
    s0_index = s0_index[my_index],
    threshold = precipitation_thresholds,
    n_cores = 1,
    remove_y0_from_y = FALSE,
    n = 1,
    r = Inf)

  df = data.frame(
    I = as.numeric(unlist(data$y) > zero_threshold),
    dist = unlist(rep(data$dist_to_s0, data$n)),
    y0 = unlist(lapply(seq_along(data$y0), function(i) rep(data$y0[[i]], each = data$n_loc[i]))))
  na_index = which(is.na(df$y))
  if (length(na_index) > 0) df = df[-na_index, ]

  spde = local({
    mesh_coords = coords[get_s0_index(coords, 4), ]
    convex1 = -.1 # The original convexity used for creating the inner mesh boundary
    convex2 = -1 # The original convexity used for creating the outer mesh boundary
    while (TRUE) {
      # The original mesh arguments
      mesh_args = list(
        loc = rbind(data$s0[[1]], mesh_coords),
        boundary = list(
          inla.nonconvex.hull(mesh_coords, convex = convex1),
          inla.nonconvex.hull(mesh_coords, convex = convex2)),
        max.edge = c(10, 100))
      mesh = create_mesh(mesh_args)
      if (!is.null(mesh)) break # Break out of the loop if this succeeded
      convex2 = convex2 - .1 # If not, increase the convexity range and try again
    }
    inla.spde2.pcmatern(mesh, prior.range = c(40, .5), prior.sigma = c(1, .5), alpha = 1.5)
  })

  dist_to_s0_from_mesh = as.numeric(dist_euclid(data$s0[[1]], spde$mesh$loc[, -3]))
  constr_index = which(dist_to_s0_from_mesh == 0)

  spde_priors = list(
    rho = c(70, .95),
    sigma = c(5, .99))
  spde_model = multimesh_spde_model(
    spdes = list(spde),
    priors = spde_priors,
    init = log(c(35, 10)),
    mesh_index = rep(seq_along(data$n), data$n),
    constr_index = constr_index)

  A_spatial = local({
    A = inla.spde.make.A(spde$mesh, coords[data$obs_index[[1]], ])
    if (!is.null(constr_index)) A = A[, -constr_index]
    A = Matrix::bdiag(rep(list(A), data$n))
    if (length(na_index) > 0) A = A[-na_index, ]
    A
  })
  identical(nrow(df), nrow(A_spatial))
  identical(spde_model$f$cgeneric$data$ints$n, ncol(A_spatial))

  create_X = function(dist, df = TRUE) {
    res = data.frame(
      intercept = 1,
      x1 = pmin(dist, 5),
      x2 = pmin(pmax(dist - 5, 0), 5),
      x3 = pmin(pmax(dist - 10, 0), 5),
      x4 = pmin(pmax(dist - 15, 0), 15),
      x5 = pmin(pmax(dist - 30, 0), 15),
      x6 = pmax(dist - 45, 0))
    if (!df) res = as.matrix(res)
    res
  }

  my_df = create_X(df$dist)

  stack = inla.stack(
    data = list(I = df$I),
    A = list(spatial = A_spatial, 1),
    effects = list(spatial = seq_len(ncol(A_spatial)), my_df))

  formula = I ~
    -1 + intercept +
    x1 + x2 + x3 + x4 + x5 + x6 +
    f(spatial, model = spde_model)

  create_Q = function(theta, spde) {
    log_rho = theta[1]
    log_sigma = theta[2]
    inla.spde2.precision(spde, c(log_rho, log_sigma))
  }

  message("Start the spatial probit inference")
  fit = inla(
    formula = formula,
    family = "binomial",
    control.family = list(control.link = list(model = "probit")),
    data = inla.stack.data(stack),
    control.predictor = list(A = inla.stack.A(stack)),
    inla.mode = "experimental",
    control.fixed = list(prec = .01),
    control.inla = list(int.strategy = "eb"),
    control.compute = list(config = TRUE),
    num.threads = 1)
  message("Finished the spatial probit inference")

  set.seed(1)
  n_samples = 1e3
  hyperpar_samples = inla.hyperpar.sample(
    n = n_samples,
    result = fit,
    improve.marginals = TRUE,
    intern = TRUE)
  fixed_samples = local({
    selection = rep(list(1), nrow(fit$summary.fixed))
    names(selection) = row.names(fit$summary.fixed)
    inla.posterior.sample(
      n = n_samples,
      result = fit,
      seed = 1L,
      selection = selection)
  })
  fixed_samples = sapply(fixed_samples, `[[`, "latent") |>
    t()
  colnames(fixed_samples) = row.names(fit$summary.fixed)

  res = list(
    create_X = create_X,
    create_Q = create_Q,
    s0 = radar$s0[my_index],
    time_index = data$time_index,
    obs_index = data$obs_index)

  res$spatial = list(
    fixed_samples = fixed_samples,
    hyperpar_samples = hyperpar_samples,
    constr = !is.null(constr_index),
    cpu = fit$cpu,
    mlik = fit$mlik,
    hyperpar = fit$summary.hyperpar,
    summary.hyperpar = fit$summary.hyperpar,
    summary.fixed = fit$summary.fixed)

  formula = I ~
    -1 + intercept +
    x1 + x2 + x3 + x4 + x5 + x6

  my_df$I = df$I

  message("Start the non-spatial probit inference")
  fit = inla(
    formula = formula,
    family = "binomial",
    data = my_df,
    control.fixed = list(prec = .01),
    control.family = list(control.link = list(model = "probit")),
    inla.mode = "experimental",
    control.inla = list(int.strategy = "eb"),
    control.compute = list(config = TRUE),
    num.threads = 1)
  message("Finished the non-spatial probit inference")

  set.seed(1)
  fixed_samples = local({
    selection = rep(list(1), nrow(fit$summary.fixed))
    names(selection) = row.names(fit$summary.fixed)
    inla.posterior.sample(
      n = n_samples,
      result = fit,
      seed = 1L,
      selection = selection)
  })
  fixed_samples = sapply(fixed_samples, `[[`, "latent") |>
    t()
  colnames(fixed_samples) = row.names(fit$summary.fixed)

  res$non_spatial = list(
    fixed_samples = fixed_samples,
    cpu = fit$cpu,
    mlik = fit$mlik,
    hyperpar = fit$summary.hyperpar,
    summary.fixed = fit$summary.fixed)

  tmp = readRDS(filename)
  if (is.null(tmp$inla)) tmp$inla = list()
  tmp$inla[[my_index]] = res
  saveRDS(tmp, filename)

  end_iter_time = proc.time()
  iter_time = (end_iter_time - start_iter_time)[3]
  total_time = (end_iter_time - start_total_time)[3]

  message(
    "Done with iter nr. ", my_index, " / ", length(s0_index), ".\n",
    "Time spent in this iteration: ", lubridate::seconds_to_period(iter_time), ".\n",
    "Time spent in total: ", lubridate::seconds_to_period(total_time), ".\n")

  # Clean up the largest object you made
  rm(data, df, spde, my_df, stack, fit, tmp, res, hyperpar_samples, fixed_samples)
  gc()
}
