library(extremePrecipitation)
library(sf)
library(ggplot2)
library(dplyr)
library(INLA)
library(inlabru)

# Activate pardiso if available
INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

# Choose a filename for the gamma_model results
filename = file.path(results_dir(), "gamma_model.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

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

# Plot nonzero proportions in time with different zero thresholds
# ------------------------------------------------------------------------------
zero_thresholds = c(0, .01, .1, .2)

plot = local({
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

plot_tikz(
  plot,
  file = file.path(image_dir(), "small_proportions.pdf"),
  width = 11,
  height = 3)

# Plot statistics of the marginal nonzero precipitation distributions
# ------------------------------------------------------------------------------
zero_threshold = .1 # Our chosen zero threshold

info = local({
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

plots = list()
for (name in unique(info$name)) {
  plots[[name]] = info |>
    dplyr::filter(name == !!name) |>
    dplyr::mutate(year = factor(year)) |>
    ggplot() +
    geom_line(aes(x = week, y = value)) +
    facet_wrap(~year, nrow = 2) +
    labs(x = "Week nr.", y = "mm/h", title = "A)") +
    theme_light() +
    theme(
      text = element_text(size = 11),
      strip.text = element_text(colour = "black", size = 13),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
}

# If we have already executed this script and performed inference for the marginal models,
# then read the results and add the fitted splined to the exploratory plot
if (!is.null(readRDS(filename)$u)) {
  plots$upper = local({
    res = readRDS(filename)
    res$u$week = res$u$day / 7 + 1
    res$u$year = factor(res$u$year)
    plots$upper +
      geom_line(
        data = res$u,
        mapping = aes(x = week, y = mean),
        col = "red")
  })
}

plot_tikz(
  plots,
  file = file.path(image_dir(), "marginal_bulk.pdf"),
  width = 9,
  height = 2.5)

# ==============================================================================
# Remove neighbours to the Rissa radar
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

zero_threshold = .1 # Our chosen zero threshold

# Get the indices of coordinates from a 2x2 subgrid of the data
delta_obs = 2
obs_index = get_s0_index(coords, delta_obs)

# Create a data.frame with all positive precipitation observations on
# the 2x2 subgrid
df = data.frame(
  y = as.numeric(radar$data[, obs_index]),
  time = rep(radar$times, length(obs_index)),
  day = rep(radar$day, length(obs_index)),
  year = rep(radar$year, length(obs_index)))
df$y[which(df$y <= zero_threshold)] = NA_real_
na_index = which(is.na(df$y))
df = df[-na_index, ]
df$y = df$y - zero_threshold

# Create a model component for inlabru of the temporal Gaussian smoothing spline
day_model = local({
  day_mesh = inla.mesh.1d(
    loc = seq(min(df$day), max(df$day), by = 12),
    degree = 1)
  inla.spde2.pcmatern(day_mesh, prior.range = c(28, .95), prior.sigma = c(3, .05))
})

# Define the component formula for inlabru
components = y ~
  -1 +
  day(day, model = day_model, group = year, control.group = list(model = "iid"))

# The probability p_u such that the threshold u is the p_u-quantile
u_prob = .95

# Perform inference with inlabru
fit = bru(
  components,
  like(
    formula = y ~ .,
    family = "gamma",
    data = df,
    control.family = list(
      control.link = list(model = "quantile", quantile = u_prob),
      hyper = list(prec = list(prior = "loggamma", param = c(1, .5))))),
  options = list(
    control.inla = list(int.strategy = "eb"),
    num.threads = 1,
    verbose = TRUE,
    bru_verbose = TRUE,
    inla.mode = "experimental"))

# Save the necessary information from the model fit
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
    n.samples = 500L,
    seed = 1)
  list(
    prob = u_prob,
    u = predictions,
    fixed = fit$summary.fixed,
    hyperpar = fit$summary.hyperpar,
    a = fit$summary.hyperpar[1, 1],
    b = dplyr::mutate(
      pred_df,
      b = qgamma(u_prob, fit$summary.hyperpar[1, 1], 1) / predictions$mean
    ),
    zero_threshold = zero_threshold,
    delta_obs = delta_obs,
    coords = coords[obs_index, ],
    cpu = fit$cpu,
    mlik = fit$mlik)
})

saveRDS(res, filename)

# ==============================================================================
# Plot the results
# ==============================================================================

# Load the results
res = readRDS(filename)

# Create a data.frame used for plotting the results
all_data = data.frame(y = as.numeric(radar$data))
all_data$y[which(all_data$y <= zero_threshold)] = NA_real_
all_data$day = rep(radar$day, n_loc)
all_data$year = rep(radar$year, n_loc)
all_data$month = rep(radar$month, n_loc)
na_index = which(is.na(all_data$y))
all_data = all_data[-na_index, ]
all_data$y = all_data$y - zero_threshold
all_data = dplyr::left_join(all_data, res$u, by = c("day", "year"))
# Estimate the rate parameters of the fitted gamma distributions
all_data$b = qgamma(res$prob, res$a, 1) / all_data$mean
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
