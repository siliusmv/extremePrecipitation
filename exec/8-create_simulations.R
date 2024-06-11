library(extremePrecipitation)
library(sf)
library(INLA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Filenames for all of the fitted model results
gp_filename = file.path(results_dir(), "gp_model.rds")
gamma_filename = file.path(results_dir(), "gamma_model.rds")
intensity_filename = file.path(results_dir(), "intensity_process.rds")
occurrence_filename = file.path(results_dir(), "occurrence_process.rds")

# ==============================================================================
# Load the data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
gp_res = readRDS(gp_filename)
gamma_res = readRDS(gamma_filename)
intensity_res = readRDS(intensity_filename)
occurrence_res = readRDS(occurrence_filename)
catchment = st_read(file.path(downloads_dir(), "stordalselva.geojson")) |>
  st_transform(st_crs(radar$coords))

zero_threshold = gamma_res$zero_threshold
threshold = intensity_res$threshold
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)
n_time = length(radar$times)
radar$data = t(radar$data) # We will mostly extract data for all locations at different time points
radar$data[radar$data <= zero_threshold] = 0

# Find out which of the locations in `coords` are inside the Stordalselva catchment
inside_catchment_index = radar$coords |>
  st_within(catchment, sparse = FALSE) |>
  which()
coords_inside_catchment = coords[inside_catchment_index, ]

# Remove locations close to the Rissa radar
rissa = radar$rissa
dist_to_rissa = as.numeric(st_distance(rissa, radar$coords))
bad_radius = 5
bad_index = which(dist_to_rissa <= bad_radius)
if (length(bad_index) > 0) {
  radar$coords = radar$coords[-bad_index, ]
  radar$data = radar$data[-bad_index, ]
  coords = coords[-bad_index, ]
  n_loc = nrow(coords)
}

# ==============================================================================
# Extract the extremes over the catchment
# ==============================================================================

# Transform the conditional extremes threshold to the precipitation scale for
# each location in `coords`
precipitation_thresholds = rep(NA_real_, length(radar$times))
for (i in seq_along(gamma_res$b$b)) {
  index = which(gamma_res$b$day[i] == radar$day & gamma_res$b$year[i] == radar$year)
  precipitation_thresholds[index] = inv_transform(
      x = threshold,
      a = gamma_res$a,
      b = gamma_res$b$b[i],
      u = gamma_res$u$mean[i],
      p_u = gamma_res$prob,
      xi = gp_res$xi,
      s = gp_res$s$mean[i],
      alpha = gp_res$prob)
}
precipitation_thresholds = precipitation_thresholds + zero_threshold

# Compute the indices for the rows of `coords` that represent the five chosen
# conditioning sites we will use for inference
s0_indices = sapply(
  X = radar$s0,
  FUN = \(x) which(abs(x[1] - coords[, 1]) <= .1 & abs(x[2] - coords[, 2]) <= .1))

# Extract all data where a threshold exceedance is observed at one of the five conditioning sites
obs_data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_indices,
  threshold = precipitation_thresholds,
  remove_y0_from_y = FALSE,
  obs_index = inside_catchment_index)

# ==============================================================================
# Create the simulations
# ==============================================================================

# Define an spde model in order to simulate from the fitted models
mesh = inla.mesh.2d(
  loc = coords_inside_catchment,
  boundary = list(
    inla.nonconvex.hull(coords_inside_catchment, convex = -.1),
    inla.nonconvex.hull(coords_inside_catchment, convex = -.8)),
  max.edge = c(20, 100))
spde = inla.spde2.pcmatern(mesh, prior.range = c(30, .5), prior.sigma = c(2, .5))

n_samples = 1000
n_per_sample = 1
n_cores = 10

# Simulate precipitation data from location nr. i, for i = 1, 2, ..., 5
set.seed(123)
simulations = list()
pb = progress_bar(length(radar$s0))
for (i in seq_along(radar$s0)) {
  simulations[[i]] = local({
    intensity_index = sapply(intensity_res$inla, \(x) identical(x$s0, radar$s0[i])) |>
      which()
    tmp = intensity_res[["inla"]][[intensity_index]]

    # Find the index for conditioning site nr. i
    s0_index = which(
      abs(coords_inside_catchment[, 1] - tmp$s0[[1]][1]) < .1 &
      abs(coords_inside_catchment[, 2] - tmp$s0[[1]][2]) < .1)

    # Simulate intensity samples from model fit nr. i
    intensity_samples = simulate_intensity_posterior(
      samples = tmp$samples[seq_len(min(n_samples, nrow(tmp$samples))), ],
      threshold = threshold,
      coords = coords_inside_catchment,
      s0_index = s0_index,
      spde = spde,
      get_a_func = intensity_res$get_a_func,
      get_b_func = intensity_res$get_b_func,
      get_Q = intensity_res$get_Q,
      get_tau = intensity_res$get_tau,
      n_cores = n_cores,
      n_per_sample = n_per_sample)
    intensity_samples$y = intensity_samples$y[[1]]
    intensity_samples$y0 = intensity_samples$y0[[1]]
    intensity_samples$dist_to_s0 = intensity_samples$dist_to_s0[[1]]

    # Add the value of the threshold exceedance y0 into the simulated intensities
    intensity_samples$y = rbind(
      intensity_samples$y[seq_len(s0_index - 1), ],
      intensity_samples$y0,
      intensity_samples$y[-seq_len(s0_index - 1), ])
    intensity_samples$dist_to_s0 = c(
      intensity_samples$dist_to_s0[seq_len(s0_index - 1)],
      0,
      intensity_samples$dist_to_s0[-seq_len(s0_index - 1)])
    intensity_samples$coords = rbind(
      intensity_samples$coords[seq_len(s0_index - 1), ],
      coords_inside_catchment[s0_index, ],
      intensity_samples$coords[-seq_len(s0_index - 1), ])

    # Draw time points uniformly in order to backtransform the data to
    # the precipitation scale
    time_index = sample.int(n_time, intensity_samples$n, replace = FALSE)
    day = radar$day[time_index]
    year = radar$year[time_index]
    intensity_samples$time_index = time_index

    # Backtransform the intensity samples to the precipitation scale
    all_indices = NULL
    for (j in seq_along(gamma_res$b$b)) {
      index = which(
        day == gamma_res$u$day[j] &
        year == gamma_res$u$year[j])
      if (any(index)) {
        intensity_samples$y[, index] = inv_transform(
          x = intensity_samples$y[, index],
          a = gamma_res$a,
          b = gamma_res$b$b[j],
          u = gamma_res$u$mean[j],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[j],
          alpha = gp_res$prob)
        all_indices = c(all_indices, index)
      }
    }
    intensity_samples$y = intensity_samples$y + zero_threshold

    occurrence_index = sapply(occurrence_res$inla, \(x) identical(x$s0, radar$s0[i])) |>
      which()
    tmp = occurrence_res$inla[[occurrence_index]]

    # Spatial probit model
    # ------------------------------------------------------------------------------

    # Find out if the SPDE model of the spatial probit model is constrained or not
    if (tmp$spatial$constr) {
      constr_index = which(
        abs(spde$mesh$loc[, 1] - tmp$s0[[1]][1]) < .1 &
        abs(spde$mesh$loc[, 2] - tmp$s0[[1]][2]) < .1)
      stopifnot(length(constr_index) == 1)
    } else {
      constr_index = NULL
    }

    # Simulate precipitation occurrences from the spatial probit model
    occurrence_samples = simulate_spat_probit_posterior(
      hyperpar_samples = tmp$spatial$hyperpar_samples[
        seq_len(min(n_samples, nrow(tmp$spatial$hyperpar_samples))), ],
      fixed_samples = tmp$spatial$fixed_samples[
        seq_len(min(n_samples, nrow(tmp$spatial$fixed_samples))), ],
      n_per_sample = n_per_sample,
      n_cores = n_cores,
      coords = coords_inside_catchment,
      replace_zeros_at_s0 = !tmp$spatial$constr,
      spde = spde,
      s0_index = s0_index,
      constr_index = constr_index,
      create_X = tmp$create_X,
      create_Q = tmp$create_Q)

    # Multiply intensity samples with occurrence samples
    stopifnot(all(dim(intensity_samples$y) == dim(occurrence_samples$simulations)))
    intensity_samples$y_spatial_probit = intensity_samples$y * occurrence_samples$simulations

    # Simulate precipitation occurrences from the non-spatial probit model
    occurrence_samples = simulate_probit_posterior(
      fixed_samples = tmp$non_spatial$fixed_samples[
        seq_len(min(n_samples, nrow(tmp$non_spatial$fixed_samples))), ],
      n_per_sample = n_per_sample,
      n_cores = n_cores,
      coords = coords_inside_catchment,
      s0_index = s0_index,
      create_X = tmp$create_X)

    # Multiply intensity samples with occurrence samples
    stopifnot(all(dim(intensity_samples$y) == dim(occurrence_samples$simulations)))
    intensity_samples$y_probit = intensity_samples$y * occurrence_samples$simulations

    # Set small precipitation intensity values equal to zero,
    # using thethreshold occurrence model
    intensity_samples$y_threshold = threshold_occurrence(
      samples = intensity_samples$y,
      nonzero_prob = mean(obs_data$y[[i]] > zero_threshold))

    intensity_samples
  })
  pb$tick()
}
pb$terminate()

simulations = purrr::transpose(simulations)
simulations$n = unlist(simulations$n)

simulations$names = c("y_spatial_probit", "y_probit", "y_threshold", "y")
simulations$pretty_names = factor(
  simulations$names,
  levels = simulations$names,
  labels = c("Spatial probit", "Probit", "Threshold", "Nonzero"))

# ==============================================================================
# Compare sums of aggregated precipitation using QQ plots
# ==============================================================================

# Define the radii of the balls used for computing aggregated precipitation
dists = c(seq(5, 20, by = 5), Inf)

# Create QQ plots of simulated aggregated precipitation for each ball radius
plot_data = list()
for (i in seq_along(dists)) {
  plot_data[[i]] = local({
    sums = list()
    sums$obs = lapply(
      X = seq_along(obs_data$n),
      FUN = function(j) {
        index = which(obs_data$dist_to_s0[[j]] <= dists[i])
        apply(obs_data$y[[j]][index, ], 2, sum)
      })
    for (name in simulations$names) {
      sums[[name]] = lapply(
        X = seq_along(simulations$n),
        FUN = function(j) {
          index = which(simulations$dist_to_s0[[j]] <= dists[i])
          apply(simulations[[name]][[j]][index, ], 2, sum)
        })
    }
    probs = head(seq(0, 1, by = .05), -1)[-1]
    df = list()
    for (j in seq_along(simulations$n)) {
      for (l in seq_along(simulations$names)) {
        df[[length(df) + 1]] = data.frame(
          sim = quantile(sums[[simulations$names[l]]][[j]], probs),
          obs = quantile(sums$obs[[j]], probs),
          prob = probs,
          j = j,
          tag = paste0("Conditioning site nr. ", j),
          dist = dists[i],
          name = simulations$names[l],
          pretty_name = simulations$pretty_names[l])
      }
    }
    do.call(rbind, df)
  })
}
plot_data = do.call(rbind, plot_data)

# Merge all of the QQ plots together into one single plot
plot = local({
  res = list()
  for (i in seq_along(dists)) {
    df = plot_data |>
      dplyr::bind_rows() |>
      dplyr::filter(dist == dists[i])
    res[[i]] = df |>
      ggplot() +
      geom_point(aes(x = sim, y = obs, col = pretty_name)) +
      geom_abline(slope = 1, intercept = 0) +
      facet_wrap(~tag, nrow = 1) +
      #facet_grid(i ~ j) +
      labs(x = "Simulation quantiles", y = "Observation quantiles",
           col = "Occurrence model",
           title = paste(
             "$d =",
             ifelse(is.finite(dists[i]), dists[i], "\\infty"),
             "$")
           ) +
      coord_equal(
        xlim = c(0, max(c(df$sim, df$obs))),
        ylim = c(0, max(c(df$sim, df$obs)))) +
      theme_light() +
      theme(
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
  }
  patchwork::wrap_plots(res, ncol = 1, guides = "collect")
})

plot_tikz(
  plot,
  file = file.path(image_dir(), "precipitation_sum_qq.pdf"),
  width = 6.5 * 1.5,
  height = 8 * 1.5)

# Create QQ plots divided over multiple pages, one page for each ball radius
plots = local({
  res = list()
  for (j in unique(plot_data$j)) {
    df = plot_data |>
      dplyr::bind_rows() |>
      dplyr::filter(j == !!j) |>
      dplyr::mutate(
        dist = factor(
          dist,
          levels = sort(unique(dist)),
          labels = paste0("$d = ", ifelse(
            is.finite(sort(unique(dist))), sort(unique(dist)), "\\infty"),
            "$"))) |>
      dplyr::group_by(dist) |>
      dplyr::mutate(maxval = max(c(sim, obs)))
    res[[j]] = df |>
      ggplot() +
      geom_point(aes(x = sim, y = obs, col = pretty_name)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_blank(aes(x = maxval, y = maxval)) +
      facet_wrap(~dist, nrow = 1, scales = "free") +
      labs(x = "Simulation quantiles", y = "Observation quantiles",
           col = "Occurrence model") +
      theme_light() +
      theme(
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        aspect.ratio = 1,
        legend.position = "top")
  }
  res
})

plot_tikz(
  plots,
  file = file.path(image_dir(), "precipitation_sum_qq_individual.pdf"),
  width = 11 * .9,
  height = 3 * .9)

# ==============================================================================
# Plot simulated and observed precipitation
# ==============================================================================

# Only plot simulated and observed data where the observation at the
# conditioning site exceeds `my_threshold`
my_threshold = 5

# Sample random precipitation realisations from the observed and the
# simulated data sets, and format the realisations in a way that is easy
# to plot with ggplot
set.seed(123)
plot_data = local({
  plot_data = list()
  for (i in seq_along(simulations$n)) {
    # simulations
    index = which(simulations$y0[[i]] > my_threshold)
    stopifnot(length(index) > 0)
    simulations$y0[[i]] = simulations$y0[[i]][index]
    for (name in simulations$names) {
      simulations[[name]][[i]] = simulations[[name]][[i]][, index, drop = FALSE]
    }
    simulations$n[i] = length(index)
    # obs_data
    index = which(obs_data$y0[[i]] > my_threshold)
    stopifnot(length(index) > 0)
    obs_data$y0[[i]] = obs_data$y0[[i]][index]
    obs_data$y[[i]] = obs_data$y[[i]][, index, drop = FALSE]
    obs_data$n[i] = length(index)
  }
  s0_choice = rep(1:5, each = 2)
  for (i in seq_along(s0_choice)) {
    j = sample.int(obs_data$n[s0_choice[i]], 1)
    plot_data[[i]] = data.frame(
      Observations = obs_data$y[[s0_choice[i]]][, j],
      x_coord = coords_inside_catchment[, 1],
      y_coord = coords_inside_catchment[, 2],
      i = i,
      s0_choice = s0_choice[i])
    j = sample.int(simulations$n[s0_choice[i]], 1)
    for (ii in seq_along(simulations$names)) {
      plot_data[[i]][[levels(simulations$pretty_names)[ii]]] =
        simulations[[simulations$names[ii]]][[s0_choice[i]]][, j]
    }
    plot_data[[i]] = plot_data[[i]] |>
      tidyr::pivot_longer(
        tidyselect::all_of(c("Observations", levels(simulations$pretty_names)))
      ) |>
      dplyr::mutate(
        value = ifelse(value <= zero_threshold, NA_real_, value),
        name = factor(name, levels = c("Observations", levels(simulations$pretty_names))))
  }
  attr(plot_data, "s0_choice") = s0_choice
  plot_data
})

# Create an sf object of the locations of the five conditioning sites to add to the
# plots of simulated and observed precipitation realisations
s0_choice = attr(plot_data, "s0_choice")
s0_data = st_as_sf(radar$s0)[rep(s0_choice, each = length(levels(plot_data[[1]]$name))), ] |>
  dplyr::mutate(
    name = rep(levels(plot_data[[1]]$name), length(plot_data)),
    i = rep(seq_along(plot_data), each = length(levels(plot_data[[1]]$name))),
    name = factor(name, levels = unique(name)))

# Plot everything
plot = plot_data |>
  dplyr::bind_rows() |>
  ggplot() +
  geom_raster(aes(x = x_coord, y = y_coord, fill = sqrt(value))) +
  geom_sf(data = catchment, fill = NA) +
  geom_sf(data = s0_data, col = "red", size = 1) +
  scale_fill_viridis_c(
    breaks = asinh(c(0, exp(0:5))),
    labels = c("0", "1", paste0("e$^", 1:5, "$"))) +
  facet_grid(i ~ name) +
  theme_light() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Easting", y = "Northing", fill = "mm/h") +
  scale_x_continuous(
    breaks = c(10, 10.4, 10.8),
    labels = paste0(c(10, 10.4, 10.8), "$^\\circ$E"),
    expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(63.9, 64, 64.1),
    labels = paste0(c(63.9, 64, 64.1), "$^\\circ$E"),
    expand = c(0, 0))

plot_tikz(
  plot,
  file = file.path(image_dir(), "precipitation_simulations.pdf"),
  width = 10,
  height = 12)

# ==============================================================================
# Create simulations over a larger domain than the Stordalselva catchment
# (i.e., basicly repeat what we did previously, but this time for the entire
# `coords` object, instead of for only the locations inside the catchment)
# ==============================================================================

# Define an spde model in order to simulate from the fitted models
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(
    inla.nonconvex.hull(coords, convex = -.1),
    inla.nonconvex.hull(coords, convex = -.8)),
  max.edge = c(20, 100))
spde = inla.spde2.pcmatern(mesh, prior.range = c(30, .5), prior.sigma = c(2, .5))

n_samples = 200
n_per_sample = 10
n_cores = 10

# Extract all observed precipitation data for the entire spatial domain
# from when we have a threshold exceedance at one of the five chosen
# conditioning sites
obs_data = extract_extreme_fields(
  data = radar$data,
  coords = coords,
  s0_index = s0_indices,
  threshold = precipitation_thresholds,
  n = 1,
  r = Inf,
  remove_y0_from_y = FALSE)

# Simulate precipitation data from location nr. i, for i = 1, 2, ..., 5
set.seed(123)
simulations = list()
pb = progress_bar(length(radar$s0))
for (i in seq_along(radar$s0)) {
  simulations[[i]] = local({
    intensity_index = sapply(intensity_res$inla, \(x) identical(x$s0, radar$s0[i])) |>
      which()
    tmp = intensity_res[["inla"]][[intensity_index]]

    # Find the index for conditioning site nr. i
    s0_index = which(
      abs(coords[, 1] - tmp$s0[[1]][1]) < .1 &
      abs(coords[, 2] - tmp$s0[[1]][2]) < .1)

    # Simulate intensity samples from model fit nr. i
    intensity_samples = simulate_intensity_posterior(
      samples = tmp$samples[seq_len(min(n_samples, nrow(tmp$samples))), ],
      threshold = threshold,
      coords = coords,
      s0_index = s0_index,
      spde = spde,
      get_a_func = intensity_res$get_a_func,
      get_b_func = intensity_res$get_b_func,
      get_Q = intensity_res$get_Q,
      get_tau = intensity_res$get_tau,
      n_cores = n_cores,
      n_per_sample = n_per_sample)
    intensity_samples$y = intensity_samples$y[[1]]
    intensity_samples$y0 = intensity_samples$y0[[1]]
    intensity_samples$dist_to_s0 = intensity_samples$dist_to_s0[[1]]

    # Add the value of the threshold exceedance y0 into the simulated intensities
    intensity_samples$y = rbind(
      intensity_samples$y[seq_len(s0_index - 1), ],
      intensity_samples$y0,
      intensity_samples$y[-seq_len(s0_index - 1), ])
    intensity_samples$dist_to_s0 = c(
      intensity_samples$dist_to_s0[seq_len(s0_index - 1)],
      0,
      intensity_samples$dist_to_s0[-seq_len(s0_index - 1)])
    intensity_samples$coords = rbind(
      intensity_samples$coords[seq_len(s0_index - 1), ],
      coords[s0_index, ],
      intensity_samples$coords[-seq_len(s0_index - 1), ])

    # Draw time points uniformly in order to backtransform the data to
    # the precipitation scale
    time_index = sample.int(n_time, intensity_samples$n, replace = FALSE)
    day = radar$day[time_index]
    year = radar$year[time_index]
    intensity_samples$time_index = time_index

    # Backtransform the intensity samples to the precipitation scale
    all_indices = NULL
    for (j in seq_along(gamma_res$b$b)) {
      index = which(
        day == gamma_res$u$day[j] &
        year == gamma_res$u$year[j])
      if (any(index)) {
        intensity_samples$y[, index] = inv_transform(
          x = intensity_samples$y[, index],
          a = gamma_res$a,
          b = gamma_res$b$b[j],
          u = gamma_res$u$mean[j],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[j],
          alpha = gp_res$prob)
        all_indices = c(all_indices, index)
      }
    }
    intensity_samples$y = intensity_samples$y + zero_threshold

    occurrence_index = sapply(occurrence_res$inla, \(x) identical(x$s0, radar$s0[i])) |>
      which()

    tmp = occurrence_res$inla[[occurrence_index]]

    # Spatial probit model
    # ------------------------------------------------------------------------------

    # Find out if the SPDE model of the spatial probit model is constrained or not
    if (tmp$spatial$constr) {
      constr_index = which(
        abs(spde$mesh$loc[, 1] - tmp$s0[[1]][1]) < .1 &
        abs(spde$mesh$loc[, 2] - tmp$s0[[1]][2]) < .1)
      stopifnot(length(constr_index) == 1)
    } else {
      constr_index = NULL
    }

    # Simulate precipitation occurrences from the spatial probit model
    occurrence_samples = simulate_spat_probit_posterior(
      hyperpar_samples = tmp$spatial$hyperpar_samples[
        seq_len(min(n_samples, nrow(tmp$spatial$hyperpar_samples))), ],
      fixed_samples = tmp$spatial$fixed_samples[
        seq_len(min(n_samples, nrow(tmp$spatial$fixed_samples))), ],
      n_per_sample = n_per_sample,
      n_cores = n_cores,
      coords = coords,
      replace_zeros_at_s0 = !tmp$spatial$constr,
      spde = spde,
      s0_index = s0_index,
      constr_index = constr_index,
      create_X = tmp$create_X,
      create_Q = tmp$create_Q)

    # Multiply intensity samples with occurrence samples
    stopifnot(all(dim(intensity_samples$y) == dim(occurrence_samples$simulations)))
    intensity_samples$y_spatial_probit = intensity_samples$y * occurrence_samples$simulations

    # Simulate precipitation occurrences from the non-spatial probit model
    occurrence_samples = simulate_probit_posterior(
      fixed_samples = tmp$non_spatial$fixed_samples[
        seq_len(min(n_samples, nrow(tmp$non_spatial$fixed_samples))), ],
      n_per_sample = n_per_sample,
      n_cores = n_cores,
      coords = coords,
      s0_index = s0_index,
      create_X = tmp$create_X)

    # Multiply intensity samples with occurrence samples
    stopifnot(all(dim(intensity_samples$y) == dim(occurrence_samples$simulations)))
    intensity_samples$y_probit = intensity_samples$y * occurrence_samples$simulations

    # Set small precipitation intensity values equal to zero,
    # using thethreshold occurrence model
    intensity_samples$y_threshold = threshold_occurrence(
      samples = intensity_samples$y,
      nonzero_prob = mean(obs_data$y[[i]] > zero_threshold))

    intensity_samples
  })
  pb$tick()
}
pb$terminate()

simulations = purrr::transpose(simulations)
simulations$n = unlist(simulations$n)

simulations$names = c("y_spatial_probit", "y_probit", "y_threshold", "y")
simulations$pretty_names = factor(
  simulations$names,
  levels = simulations$names,
  labels = c("Spatial probit", "Probit", "Threshold", "Nonzero"))

# ==============================================================================
# Examine occurrence properties for the simulated data,
# just as how we did it in exec/7-occurrence_process.R
# ==============================================================================

# Compute the indices of all neighbours to location nr. i, for all i = 1, 2, ..., nrow(coords)
neighbour_radius = 1
neighbour_indices = lapply(
  X = seq_len(nrow(coords)),
  FUN = function(i) {
    dd = dist_euclid(coords[i, ], coords)
    as.numeric(which(0 < dd & dd <= neighbour_radius))
  })

# Compute the mean of the <=4 neighbours to a certain location at a certain time,
# for all locations and times in the simulated and observed data sets
neighbour_means = list()
for (name in simulations$names) {
  neighbour_means[[name]] = parallel::mclapply(
    X = seq_along(simulations$n),
    mc.cores = length(simulations$n),
    FUN = function(i) {
      res = matrix(NA_real_, nrow = nrow(coords), ncol = simulations$n[i])
      for (l in seq_len(nrow(coords))) {
        res[l, ] = apply(simulations[[name]][[i]][neighbour_indices[[l]], , drop = FALSE], 2, mean)
      }
      res
    })
  message(name)
}
neighbour_means[["obs"]] = parallel::mclapply(
  X = seq_along(obs_data$n),
  mc.cores = length(obs_data$n),
  FUN = function(i) {
    res = matrix(NA_real_, nrow = nrow(coords), ncol = obs_data$n[i])
    for (l in seq_len(nrow(coords))) {
      res[l, ] = apply(obs_data$y[[i]][neighbour_indices[[l]], , drop = FALSE], 2, mean)
    }
    res
  })

# Prepare the data for plotting with ggplot
df = local({
  res = list()
  for (j in seq_along(simulations$names)) {
    res[[j]] = lapply(
      X = seq_along(simulations$n),
      FUN = function(i) {
        data.frame(
          I = as.numeric(simulations[[simulations$name[j]]][[i]] > zero_threshold),
          neighbour_mean = as.numeric(neighbour_means[[simulations$name[j]]][[i]]),
          d = simulations$dist_to_s0[[i]],
          name = simulations$pretty_name[j])
      }) |>
      dplyr::bind_rows()
  }
  res[["obs"]] = lapply(
    X = seq_along(obs_data$n),
    FUN = function(i) {
      data.frame(
        I = as.numeric(obs_data$y[[i]] > zero_threshold),
        neighbour_mean = as.numeric(neighbour_means[["obs"]][[i]]),
        d = obs_data$dist_to_s0[[i]],
        name = "Observations")
    }) |>
    dplyr::bind_rows()
  dplyr::bind_rows(res)
})

# Plot the empirical precipitation occurrence probability as a function of the
# distance to the conditioning site
plot1 = df |>
  dplyr::mutate(d = round(d)) |>
  dplyr::group_by(d, name) |>
  dplyr::summarise(p = mean(I)) |>
  dplyr::mutate(name = factor(
    name, levels = c("Observations", levels(simulations$pretty_names)))) |>
  ggplot() +
  geom_line(aes(x = d, y = p, group = name, col = name, size = name, linetype = name)) +
  labs(x = "$d$", y = "$\\hat p(d)$", col = "", linetype = "", size = "") +
  scale_linetype_manual(values = c("solid", rep("dashed", length(simulations$names)))) +
  scale_size_manual(values = c(1.5, rep(.8, length(simulations$names)))) +
  scale_color_manual(values = c("black", scales::hue_pal()(length(simulations$names)))) +
  theme_light() +
  theme(
    axis.title.y = element_text(vjust = .5, angle = 0),
    text = element_text(size = 18),
    legend.text = element_text(size = 15))

# Plot the empirical precipitation occurrence probability as a function of the
# mean of the <=4 neighbouring locations
plot2 = df |>
  dplyr::mutate(m = ifelse(neighbour_mean < 1, round(neighbour_mean, 1), round(neighbour_mean))) |>
  dplyr::group_by(m, name) |>
  dplyr::summarise(p = mean(I)) |>
  dplyr::mutate(name = factor(
    name, levels = c("Observations", levels(simulations$pretty_names)))) |>
  ggplot() +
  geom_line(aes(x = asinh(m), y = p, group = name, col = name, size = name, linetype = name)) +
  labs(x = "$\\bar y$", y = "$\\hat p(\\bar y)$", col = "", linetype = "", size = "") +
  scale_linetype_manual(values = c("solid", rep("dashed", length(simulations$names)))) +
  scale_size_manual(values = c(1.5, rep(.8, length(simulations$names)))) +
  scale_color_manual(values = c("black", scales::hue_pal()(length(simulations$names)))) +
  scale_x_continuous(
    breaks = asinh(c(0, exp(0:3))),
    labels = c("0", "1", paste0("e$^", 1:3, "$")),
    limits = c(0, asinh(exp(3)))) +
  theme_light() +
  theme(
    axis.title.y = element_text(vjust = .5, angle = 0),
    text = element_text(size = 18),
    legend.text = element_text(size = 12))

# Extract some random realisations of precipitation occurrence given a threshold
# exceedance at conditioning site nr. 4
my_s0_index = s0_indices[4]
time_index = which(radar$data[my_s0_index, ] > precipitation_thresholds)
y0 = radar$data[my_s0_index, time_index]
plot_data = list()
set.seed(101)
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

latex_friendly_map_plot = function(x) {
  info = ggplot2::ggplot_build(x)$layout$panel_params[[1]]$graticule
  east_ticks = info$degree[info$type == "E" & info$y_start == 0]
  north_ticks = info$degree[info$type == "N" & info$x_start == 0]
  x +
    ggplot2::scale_x_continuous(breaks = east_ticks, labels = paste0(east_ticks, "$^\\circ$E"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = north_ticks, labels = paste0(north_ticks, "$^\\circ$N"), expand = c(0, 0))
}

# Plot the random occurrence realisations
plot3 = plot_data |>
  ggplot() +
  geom_raster(aes(x = x_coord, y = y_coord, fill = asinh(y))) +
  scale_fill_viridis_c(
    breaks = asinh(c(0, exp(0:5))),
    labels = c("0", "1", paste0("e$^", 1:5, "$"))) +
  geom_sf(data = radar$s0[4], size = 1, col = "red") +
  facet_wrap(~i) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 18),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 12),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Easting", y = "Northing", fill = "mm/h")
plot3 = latex_friendly_map_plot(plot3)

# Merge all of the occurrence plots together
plot = patchwork::wrap_plots(
  plot3, plot1 + guides(col = "none", size = "none", linetype = "none"), plot2,
  nrow = 1,
  widths = c(.6, .2, .2)) +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

plot_tikz(
  plot,
  file = file.path(image_dir(), "zero_simulation_properties.pdf"),
  width = 14,
  height = 4)

# ==============================================================================
# Compute χ_p for the observed and simulated data
# ==============================================================================

# This is easier if we transform everything to the Laplace scale, so we
# do that for both the observed and the simulated data

# Transform observed data
obs_data_laplace = local({
  index_count = rep(0, length(obs_data$n))
  pb = progress_bar(sum(obs_data$n))
  for (i in seq_along(gamma_res$b$b)) {
    for (j in seq_along(obs_data$n)) {
      index = which(
        radar$day[obs_data$time_index[[j]]] == gamma_res$u$day[i] &
        radar$year[obs_data$time_index[[j]]] == gamma_res$u$year[i])
      if (any(index)) {
        obs_data$y[[j]][, index] = transform(
          x = obs_data$y[[j]][, index] - zero_threshold,
          a = gamma_res$a,
          b = gamma_res$b$b[i],
          u = gamma_res$u$mean[i],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[i],
          alpha = gp_res$prob)
        obs_data$y0[[j]][index] = transform(
          x = obs_data$y0[[j]][index] - zero_threshold,
          a = gamma_res$a,
          b = gamma_res$b$b[i],
          u = gamma_res$u$mean[i],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[i],
          alpha = gp_res$prob)
        index_count[[j]] = index_count[[j]] + length(index)
        pb$tick(length(index))
      }
    }
  }
  pb$terminate()
  stopifnot(max(abs(index_count - obs_data$n)) < .1)
  obs_data
})

# Transformed simulated data
simulations_laplace = local({
  index_count = rep(0, length(simulations$n))
  pb = progress_bar(sum(simulations$n))
  for (i in seq_along(gamma_res$b$b)) {
    for (j in seq_along(simulations$n)) {
      index = which(
        radar$day[simulations$time_index[[j]]] == gamma_res$u$day[i] &
        radar$year[simulations$time_index[[j]]] == gamma_res$u$year[i])
      if (any(index)) {
        for (name in simulations$names) {
        simulations[[name]][[j]][, index] = transform(
          x = simulations[[name]][[j]][, index] - zero_threshold,
          a = gamma_res$a,
          b = gamma_res$b$b[i],
          u = gamma_res$u$mean[i],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[i],
          alpha = gp_res$prob)
        }
        simulations$y0[[j]][index] = transform(
          x = simulations$y0[[j]][index] - zero_threshold,
          a = gamma_res$a,
          b = gamma_res$b$b[i],
          u = gamma_res$u$mean[i],
          p_u = gamma_res$prob,
          xi = gp_res$xi,
          s = gp_res$s$mean[i],
          alpha = gp_res$prob)
        index_count[[j]] = index_count[[j]] + length(index)
        pb$tick(length(index))
      }
    }
  }
  pb$terminate()
  stopifnot(max(abs(index_count - simulations$n)) < .1)
  simulations
})

# Define different thresholds for computing χ
thresholds = seq(ceiling(intensity_res$threshold / .1) * .1, 6, by = .1)
chi = list()
for (i in seq_along(obs_data_laplace$n)) {
  chi[[i]] = list()
  chi[[i]][["Observations"]] = empirical_chi(
    data = lapply(obs_data_laplace, `[`, i),
    thresholds = thresholds,
    dist_centers = 1:70,
    dist_radius = .5) |>
    dplyr::mutate(name = "Observations")
  for (j in seq_along(simulations$names)) {
    chi[[i]][[simulations$names[j]]] = local({
      simulations_laplace$y = simulations_laplace[[simulations$names[j]]]
      empirical_chi(
        data = lapply(simulations_laplace, `[`, i),
        thresholds = thresholds,
        dist_centers = 1:70,
        dist_radius = .5) |>
        dplyr::mutate(name = levels(simulations$pretty_names)[j])
    })
  }
}

# Plot the estimators for χ
plots = list()
for (i in seq_along(chi)) {
  plots[[i]] = chi[[i]] |>
    dplyr::bind_rows() |>
    dplyr::filter(dist <= 60) |>
    dplyr::mutate(
      name = factor(
        name, levels = c("Observations", levels(simulations$pretty_names))),
      p = plaplace(threshold),
      minus_log_one_minus_p = -log2(1 - p)) |>
    ggplot() +
    geom_line(aes(
      x = dist,
      y = chi,
      group = minus_log_one_minus_p,
      col = minus_log_one_minus_p)) +
    facet_wrap(~name, nrow = 1) +
    theme_light() +
    labs(x = "$d$", y = "$\\hat \\chi_p(d)$", col = "$p$",
         title = paste("Conditioning site nr.", i)) +
    scale_color_continuous(
      breaks = 1:10,
      labels = paste0("$1 - 2^{-", 1:10, "}$")) +
    scale_y_continuous(breaks = seq(0, 1, by = .25), labels = seq(0, 1, by = .25)) +
    theme(
      text = element_text(size = 13),
      axis.title.y = element_text(angle = 0, vjust = .5),
      strip.text = element_text(colour = "black", size = 15),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
}

plot_tikz(
  plots,
  file = file.path(image_dir(), "precipitation_chi.pdf"),
  width = 14,
  height = 3.5)
