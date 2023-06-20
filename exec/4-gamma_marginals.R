devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(INLA)
library(inlabru)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

zero_threshold = .1
filename = file.path(results_dir(), "gamma_model.rds")

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
heights = radar$coords$height
n_loc = nrow(coords)
n_time = nrow(radar$data)
rissa = radar$rissa

# ==============================================================================
# Do some plotting
# ==============================================================================

# Plot nonzero proportions in time and space with different zero thresholds
# ------------------------------------------------------------------------------
zero_thresholds = c(0, .05, .1, .2)

plot1 = local({
  small_proportion = sapply(
    X = zero_thresholds,
    FUN = function(x) {
      apply(radar$data, 2, function(y) mean(y <= x, na.rm = TRUE))
    })
  as.data.frame(coords) |>
    cbind(small_proportion) |>
    tidyr::pivot_longer(-c(X, Y)) |>
    dplyr::mutate(name = paste("$\\tau_0 =", zero_thresholds[as.numeric(name)], "$ mm/h")) |>
    ggplot() +
    geom_sf(data = rissa) +
    geom_raster(aes(x = X, y = Y, fill = value)) +
    facet_wrap(~name, nrow = 1) +
    scale_fill_viridis_c() +
    theme_light() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      text = element_text(size = 12),
      strip.text = element_text(colour = "black", size = 15),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    labs(
      x = "Easting",
      y = "Northing",
      fill = "$P\\left(X(\\bm s, \\cdot) \\leqslant \\tau_0\\right)$"
    )
})
plot2 = local({
  small_sum = sapply(
    X = zero_thresholds,
    FUN = function(x) {
      apply(radar$data, 1, function(y) sum(y <= x, na.rm = TRUE))
    })
  non_na_sum = apply(radar$data, 1, function(y) sum(!is.na(y)))
  months = list(6, 7, 8)
  month_names = factor(
    c("June", "July", "August"),
    levels = c("June", "July", "August"))
  plot_data = list()
  for (i in seq_along(months)) {
    plot_data[[i]] = data.frame(
      year = radar$year,
      month = radar$month,
      non_na_sum = non_na_sum) |>
      cbind(small_sum) |>
      dplyr::filter(month %in% months[[i]]) |>
      tidyr::pivot_longer(-c(year, month, non_na_sum)) |>
      dplyr::mutate(name = paste("$\\tau_0 =", zero_thresholds[as.numeric(name)], "$ mm/h")) |>
      dplyr::group_by(year, name) |>
      dplyr::summarise(mean = sum(value) / sum(non_na_sum)) |>
      dplyr::mutate(month_name = month_names[i])
  }
  plot_data |>
    dplyr::bind_rows() |>
    ggplot() +
    geom_line(aes(
      x = year, y = mean, group = month_name, col = month_name)) +
    facet_wrap(~name, nrow = 1) +
    scale_fill_viridis_c() +
    theme_light() +
    scale_x_continuous(breaks = seq(min(radar$year), max(radar$year), by = 5), expand = c(0, 0)) +
    theme(
      text = element_text(size = 12),
      strip.text = element_text(colour = "black", size = 15),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
    labs(
      x = "Year",
      y = "$P\\left(X(\\cdot, t) \\leqslant \\tau_0\\right)$",
      col = "Month")
})

plot1 = latex_friendly_map_plot(plot1)
plot = patchwork::wrap_plots(plot1, plot2, nrow = 2, heights = c(.5, .5))

plot_tikz(
  plot,
  file = file.path(image_dir(), "small_proportions.pdf"),
  width = 11,
  height = 6)

# Plot statistics of the marginal distributions in time and space
# ------------------------------------------------------------------------------
info = local({
  res = apply(
    radar$data,
    2,
    function(x) {
      sapply(
        X = unique(radar$month),
        FUN = function(m) {
          index = which(x > zero_threshold & radar$month == m)
          c(mean = mean(x[index], na.rm = TRUE),
            sd = sd(x[index], na.rm = TRUE),
            median = median(x[index], na.rm = TRUE),
            upper = as.numeric(quantile(x[index], .95, na.rm = TRUE)))
        }) |>
        as.numeric()
    })
  res = cbind(coords, t(res)) |>
    as.data.frame() |>
    tidyr::pivot_longer(-c(X, Y)) |>
    dplyr::mutate(
      month = rep(rep(c("June", "July", "August"), each = 4), n_loc),
      name = rep(rep(c("mean", "sd", "median", "upper"), 3), n_loc),
      height = rep(heights, each = 12)) |>
    dplyr::mutate(month = factor(month, levels = c("June", "July", "August")))
  res
})

spatial_plots = list()
for (name in unique(info$name)) {
  spatial_plots[[name]] = info |>
    dplyr::filter(name == !!name) |>
    ggplot() +
    geom_sf(data = rissa) +
    geom_raster(aes(x = X, y = Y, fill = value)) +
    facet_wrap(~month, nrow = 1) +
    scale_fill_viridis_c() +
    labs(x = "Easting", y = "Northing", fill = "mm/h") +
    theme_light() +
    theme(
      text = element_text(size = 11),
      strip.text = element_text(colour = "black", size = 13),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
  spatial_plots[[name]] = latex_friendly_map_plot(spatial_plots[[name]])
}

info2 = local({
  res = list()
  for (y in unique(radar$year)) {
    res[[length(res) + 1]] = sapply(
      X = unique(radar$day),
      FUN = function(d) {
        index = which(abs(radar$day - d) <= 4 & radar$year == y)
        x = as.numeric(radar$data[index, ])
        x = x[x > zero_threshold]
        c(mean = mean(x, na.rm = TRUE),
          sd = sd(x, na.rm = TRUE),
          upper = as.numeric(quantile(x, .95, na.rm = TRUE)),
          median = median(x, na.rm = TRUE))
      }) |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        day = unique(radar$day),
        year = y) |>
      tidyr::pivot_longer(-c(day, year)) |>
      dplyr::mutate(week = day / 7 + 1)
  }
  dplyr::bind_rows(res)
})

temporal_plots = list()
for (name in unique(info2$name)) {
  temporal_plots[[name]] = info2 |>
    dplyr::filter(name == !!name) |>
    dplyr::mutate(year = factor(year)) |>
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

# If we have already executed this script and performed inference for the marginal models,
# then read the results and add the fitted splined to the exploratory plot
if (!is.null(readRDS(filename)$u)) {
  temporal_plots$upper = local({
    res = readRDS(filename)
    res$u$week = res$u$day / 7 + 1
    res$u$year = factor(res$u$year)
    temporal_plots$upper +
      geom_line(
        data = res$u,
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
  file = file.path(image_dir(), "marginal_bulk.pdf"),
  width = 9,
  height = 5)

plot_tikz(
  spatial_plots,
  file = file.path(image_dir(), "marginal_bulk_1.pdf"),
  width = 9,
  height = 2.8)

for (name in names(temporal_plots)) {
  temporal_plots[[name]] = temporal_plots[[name]] + labs(title = "A)")
}

plot_tikz(
  temporal_plots,
  file = file.path(image_dir(), "marginal_bulk_2.pdf"),
  width = 9,
  height = 2.5)

# ==============================================================================
# Remove neighbours to Rissa
# ==============================================================================

# Remove the Rissa radar and its closest neighbours
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

# ==============================================================================
# Start modelling the data
# ==============================================================================

#delta_obs = 2
delta_obs = 1
obs_index = get_s0_index(coords, delta_obs)

df = data.frame(
  y = as.numeric(radar$data[, obs_index]),
  time = rep(radar$times, length(obs_index)),
  day = rep(radar$day, length(obs_index)),
  year = rep(radar$year, length(obs_index)),
  sqrt_height = rep(sqrt(heights[obs_index]), n_time))
df$y[which(df$y <= zero_threshold)] = NA_real_
na_index = which(is.na(df$y))
df = df[-na_index, ]
df$y = df$y - zero_threshold

day_model = local({
  day_mesh = inla.mesh.1d(
    loc = seq(min(df$day), max(df$day), by = 12),
    degree = 1)
  inla.spde2.pcmatern(day_mesh, prior.range = c(28, .95), prior.sigma = c(3, .05))
})

components = y ~
  -1 +
  day(day, model = day_model, group = year, control.group = list(model = "iid"))

prob = .95

fit = bru(
  components,
  like(
    formula = y ~ .,
    family = "gamma",
    data = df,
    control.family = list(
      control.link = list(model = "quantile", quantile = prob),
      hyper = list(prec = list(prior = "loggamma", param = c(1, .5))))),
  options = list(
    control.inla = list(int.strategy = "eb"),
    num.threads = 1,
    verbose = TRUE,
    bru_verbose = TRUE,
    inla.mode = "experimental"))

cpu1 = fit$cpu

fit = bru_rerun(fit)

summary(fit)

# Save the necessary results of the model fit
res = local({
  all_days = sort(unique(radar$day))
  all_years = sort(unique(radar$year))
  pred_df = data.frame(
    day = rep(all_days, length(all_years)),
    year = rep(all_years, each = length(all_days)))
  predictions = predict(
    object = fit,
    data = pred_df,
    formula = ~ exp(day),
    seed = 1)
  list(
    prob = prob,
    u = predictions,
    fixed = fit$summary.fixed,
    hyperpar = fit$summary.hyperpar,
    a = fit$summary.hyperpar[1, 1],
    b = dplyr::mutate(pred_df, b = qgamma(prob, fit$summary.hyperpar[1, 1], 1) / predictions$mean),
    zero_threshold = zero_threshold,
    delta_obs = delta_obs,
    coords = coords[obs_index, ],
    cpu = fit$cpu.used,
    mlik = fit$mlik)
})

saveRDS(res, filename)


# ==============================================================================
# Plot the results
# ==============================================================================
res = readRDS(filename)

all_data = data.frame(
  y = as.numeric(radar$data),
  sqrt_height = rep(sqrt(heights), n_time))
all_data$y[which(all_data$y <= zero_threshold)] = NA_real_
all_data$day = rep(radar$day, n_loc)
all_data$year = rep(radar$year, n_loc)
all_data$month = rep(radar$month, n_loc)
na_index = which(is.na(all_data$y))
all_data = all_data[-na_index, ]
all_data$intercept = 1
all_data$y = all_data$y - zero_threshold
all_data = dplyr::left_join(all_data, res$u, by = c("day", "year"))
all_data$b = qgamma(res$prob, res$a, 1) / all_data$mean
all_data$pit = pgamma(all_data$y, shape = res$a, rate = all_data$b)
# Standardise the data so we can create QQ-plots
all_data$y_standardised = all_data$y * all_data$b

# Create QQ plots
probs = head(seq(0, 1, length.out = 50), -1)
qq_data = local({
  df1 = all_data |>
    dplyr::reframe(
      empirical = quantile(y_standardised, probs),
      model = qgamma(probs, res$a)) |>
    dplyr::mutate(month = 0)
  df2 = all_data |>
    dplyr::group_by(month) |>
    dplyr::reframe(
      empirical = quantile(y_standardised, probs),
      model = qgamma(probs, res$a))
  rbind(df1, df2)
})
qq_data$month = factor(
  qq_data$month,
  levels = c(6:8, 0),
  labels = c("June", "July", "August", "All"))

qq_plot = qq_data |>
  ggplot() +
  geom_point(aes(x = model, y = empirical)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~month, nrow = 1) +
  labs(x = "Model quantiles", y = "Observation quantiles") +
  coord_equal() +
  scale_x_continuous(breaks = seq(0, 2.5, by = .5), limits = c(-.1, 2.8), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2.5, by = .5), limits = c(-.1, 2.8), expand = c(0, 0)) +
  theme_light() +
  theme(
    text = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 13),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))

plot_tikz(
  qq_plot,
  file = file.path(image_dir(), "gamma_qq.pdf"),
  width = 9,
  height = 2.8)
