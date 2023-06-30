devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(INLA)
library(inlabru)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

gamma_filename = file.path(results_dir(), "gamma_model.rds")
filename = file.path(results_dir(), "gp_model.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

# ==============================================================================
# Prepare the data for modelling with the GP
# ==============================================================================

gamma_res = readRDS(gamma_filename)
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
heights = radar$coords$height
rissa = radar$rissa
n_time = nrow(radar$data)
n_loc = nrow(coords)

df = data.frame(
  y = as.numeric(radar$data),
  day = rep(radar$day, n_loc),
  year = rep(radar$year, n_loc),
  week = rep(radar$week, n_loc),
  month = rep(radar$month, n_loc),
  x_coord = rep(as.numeric(coords[, 1]), n_time),
  y_coord = rep(as.numeric(coords[, 2]), n_time)) |>
  dplyr::mutate(
    y = y - gamma_res$zero_threshold # Shift the data so their distributions align with gamma_res
    ) |>
  dplyr::filter(y > 0) |> # Speed things up by removing to small entries
  dplyr::left_join(dplyr::select(gamma_res$u, day, year, mean), by = c("day", "year")) |>
  dplyr::rename(u = mean) |>
  dplyr::filter(y >= u) |>
  dplyr::mutate(y = y - u)

# ==============================================================================
# Do some exploratory analysis
# ==============================================================================

spatial_plots = local({
  months = 6:8
  month_names = c("June", "July", "August")
  month_names = factor(month_names, levels = month_names)
  plot_data = list()
  for (i in seq_along(months)) {
    plot_data[[i]] = df |>
      dplyr::filter(month == months[i]) |>
      dplyr::group_by(x_coord, y_coord) |>
      dplyr::summarise(
        median = median(y),
        mean = mean(y),
        sd = sd(y),
        upper = as.numeric(quantile(y, .95))) |>
      dplyr::mutate(
        month = months[i],
        month_name = month_names[i]) |>
      tidyr::pivot_longer(c(median, mean, sd, upper))
  }
  plot_data = dplyr::bind_rows(plot_data)
  res = list()
  for (name in unique(plot_data$name)) {
    res[[name]] = plot_data |>
      dplyr::filter(name == !!name) |>
      ggplot() +
      geom_sf(data = rissa, size = .001, alpha = .001) +
      geom_raster(aes(x = x_coord, y = y_coord, fill = value)) +
      scale_fill_viridis_c() +
      facet_wrap(~month_name, nrow = 1) +
      labs(x = "Easting", y = "Northing", fill = paste("mm/h")) +
      theme_light() +
      theme(
        text = element_text(size = 11),
        strip.text = element_text(colour = "black", size = 13),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
    res[[name]] = latex_friendly_map_plot(res[[name]])
  }
  res
})

temporal_plots = local({
  plot_data = list()
  for (y in unique(df$year)) {
    for (d in unique(df$day)) {
      plot_data[[length(plot_data) + 1]] = df |>
        dplyr::filter(year == !!y, abs(day - d) <= 4) |>
        dplyr::summarise(
          median = median(y),
          mean = mean(y),
          sd = sd(y),
          upper = as.numeric(quantile(y, .95))) |>
        dplyr::mutate(year = y, day = !!d, week = day / 7 + 1) |>
        tidyr::pivot_longer(c(median, mean, sd, upper))
    }
  }
  plot_data = dplyr::bind_rows(plot_data)
  res = list()
  for (name in unique(plot_data$name)) {
    res[[name]] = plot_data |>
      dplyr::filter(name == !!name) |>
      ggplot() +
      geom_line(aes(x = week, y = value)) +
      facet_wrap(~year, nrow = 2) +
      labs(x = "Week nr.", y = "mm/h") +
      theme_light() +
      theme(
        text = element_text(size = 11),
        strip.text = element_text(colour = "black", size = 13),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
  }
  res
})

# If we have already executed this script and performed inference for the marginal models,
# then read the results and add the fitted splined to the exploratory plot
if (!is.null(readRDS(filename)$s)) {
  temporal_plots$median = local({
    res = readRDS(filename)
    res$s$week = res$s$day / 7 + 1
    res$s$year = factor(res$s$year)
    temporal_plots$median +
      geom_line(
        data = res$s,
        mapping = aes(x = week, y = mean),
        col = "red")
  })
}

pretty_names = list(
  mean = "Mean",
  sd = "Standard deviation",
  median = "Median",
  upper = "$95\\%$ quantile")

plots = list()
for (name in names(spatial_plots)) {
  plots[[name]] = patchwork::wrap_plots(
    spatial_plots[[name]],
    temporal_plots[[name]],
    nrow = 2,
    heights = c(.55, .45)) +
    patchwork::plot_annotation(title = pretty_names[[name]])
}

for (name in names(spatial_plots)) {
  spatial_plots[[name]] = spatial_plots[[name]] +
    labs(title = pretty_names[[name]])
}

plot_tikz(
  plots,
  file = file.path(image_dir(), "threshold_exceedances.pdf"),
  width = 9,
  height = 5)

plot_tikz(
  spatial_plots,
  file = file.path(image_dir(), "threshold_exceedances_1.pdf"),
  width = 9,
  height = 2.8)

for (name in names(temporal_plots)) {
  temporal_plots[[name]] = temporal_plots[[name]] + labs(title = "B)")
}

plot_tikz(
  temporal_plots,
  file = file.path(image_dir(), "threshold_exceedances_2.pdf"),
  width = 9,
  height = 2.5)

# ==============================================================================
# Remove the Rissa radar and its closest neighbours
# ==============================================================================

dist_to_rissa = as.numeric(st_distance(rissa, radar$coords))
bad_radius = 5
bad_index = which(dist_to_rissa <= bad_radius)
if (length(bad_index) > 0) {
  radar$coords = radar$coords[-bad_index, ]
  radar$data = radar$data[, -bad_index]
  coords = coords[-bad_index, ]
  heights = radar$coords$height
  n_loc = nrow(coords)
  n_time = nrow(radar$data)
}

df = data.frame(
  y = as.numeric(radar$data),
  day = rep(radar$day, n_loc),
  week = rep(radar$week, n_loc),
  year = rep(radar$year, n_loc),
  month = rep(radar$month, n_loc),
  x_coord = rep(as.numeric(coords[, 1]), n_time),
  y_coord = rep(as.numeric(coords[, 2]), n_time))
df$y = df$y - gamma_res$zero_threshold # Shift the data so their distributions align with gamma_res
df = dplyr::filter(df, y > 0) |> # Speed things up by removing to small entries
  dplyr::left_join(dplyr::select(gamma_res$u, day, year, mean), by = c("day", "year")) |>
  dplyr::rename(u = mean) |>
  dplyr::filter(y >= u) |>
  dplyr::mutate(y = y - u)

# ==============================================================================
# Start the modelling
# ==============================================================================

day_model = local({
  day_mesh = inla.mesh.1d(
    loc = seq(min(df$day), max(df$day), by = 14),
    degree = 1)
  inla.spde2.pcmatern(day_mesh, prior.range = c(28, .95), prior.sigma = c(3, .05))
})

components = y ~
  -1 +
  day(day, model = day_model, group = year, control.group = list(model = "iid"))

prob = .5

fit = bru(
  components,
  like(
    formula = y ~ .,
    family = "gp",
    data = df,
    control.family = list(
      control.link = list(quantile = prob),
      hyper = list(theta = list(prior = "pc.gevtail", param = c(7, 0, .5))))),
  options = list(
    verbose = TRUE,
    num.threads = 4,
    inla.mode = "experimental"))

summary(fit)

# Save the necessary results of the model fit
res = local({
  set.seed(1)
  all_days = sort(unique(radar$day))
  all_years = sort(unique(radar$year))
  pred_df = data.frame(
    day = rep(all_days, length(all_years)),
    year = rep(all_years, each = length(all_days)))
  predictions = predict(
    object = fit,
    data = pred_df,
    formula = ~ exp(day),
    n.samples = 500,
    seed = 1)
  list(
    prob = prob,
    s = predictions,
    u = gamma_res$u,
    fixed = fit$summary.fixed,
    hyperpar = fit$summary.hyperpar,
    xi = fit$summary.hyperpar[1, 1],
    coords = coords,
    cpu = fit$cpu.used,
    mlik = fit$mlik)
})

saveRDS(res, filename)

# ==============================================================================
# Evaluate the results
# ==============================================================================

res = readRDS(filename)

df = dplyr::left_join(df, res$s, by = c("day", "year")) |>
  dplyr::rename(s = mean)
df$pit = pgp(df$y, df$s, res$xi, alpha = res$prob)
df$y_standardised = df$y / df$s

probs = 1 - 2^-seq(0, 9, by = .2)
qq_data = local({
  df1 = df |>
    dplyr::reframe(
      empirical = quantile(y_standardised, probs),
      model = qgp(probs, s = 1, xi = res$xi, alpha = res$prob)) |>
    dplyr::mutate(month = 0)
  df2 = df |>
    dplyr::group_by(month) |>
    dplyr::reframe(
      empirical = quantile(y_standardised, probs),
      model = qgp(probs, s = 1, xi = res$xi, alpha = res$prob))
  rbind(df1, df2)
})
qq_data$month = factor(
  qq_data$month,
  levels = c(6:8, 0),
  labels = c("June", "July", "August", "All"))
qq_data$n_larger = sapply(qq_data$empirical, \(x) sum(df$y_standardised > x))

qq_plot = qq_data |>
  ggplot() +
  geom_point(aes(x = model, y = empirical)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~month, nrow = 1) +
  labs(x = "Model quantiles", y = "Observation quantiles") +
  coord_equal() +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(-.5, 17), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, by = 5), limits = c(-.5, 17), expand = c(0, 0)) +
  theme_light() +
  theme(
    text = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 13),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))

plot_tikz(
  qq_plot,
  file = file.path(image_dir(), "gp_qq.pdf"),
  width = 9,
  height = 2.8)
