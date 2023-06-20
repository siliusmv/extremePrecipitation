library(ggplot2)
library(dplyr)
library(rayshader)
library(stars)
devtools::load_all()

# ==============================================================================
# Load the data
# ==============================================================================
heights = stars::read_stars(file.path(downloads_dir(), "dem.tif"))
radar_filename = file.path(downloads_dir(), "radar.rds")
radar = readRDS(radar_filename)

poly = st_read(file.path(downloads_dir(), "stordalselva.geojson"))
poly = st_transform(poly, st_crs(heights))
boundary = st_boundary(poly)

bbox = radar$coords |>
  st_transform(st_crs(heights)) |>
  st_bbox()

s0 = list(
  st_point(c(277000, 7114000)),
  st_point(c(294000, 7113000)),
  st_point(c(268000, 7101000)),
  st_point(c(285000, 7099000)),
  st_point(c(255000, 7092000))) |>
  st_as_sfc()
st_crs(s0) = st_crs(heights)

radar$s0 = s0 |>
  st_transform(st_crs(radar$coords))

saveRDS(radar, radar_filename)

# ==============================================================================
# Plot a map of Norway with a box containing the spatial domain
# ==============================================================================

poly_bbox = st_as_sfc(bbox) |>
  st_transform(4326)

plot = ggplot() +
  geom_sf(data = poly_bbox, fill = NA, linewidth = 1.2) +
  add_norway_map(
    crs = st_crs(poly_bbox),
    bbox = c(4, 58, 32, 71),
    linewidth = .5) +
  theme_light() +
  theme(text = element_text(size = 15))

pdf(file.path(image_dir(), "norway_map.pdf"), height = 5, width = 5)
plot
dev.off()

# ==============================================================================
# Plot a fancy height map
# ==============================================================================

tmp = st_crop(heights, bbox) |>
  st_as_stars()
tmp2 = tmp[[1]]

rissa = radar$rissa |>
  st_transform(st_crs(tmp))

height_plot = tmp2 |>
  height_shade() |>
  add_overlay(sphere_shade(tmp2, texture = "desert", zscale = 10), alphalayer = .5) |>
  add_shadow(lamb_shade(tmp2, zscale = 10), .8) |>
  add_shadow(texture_shade(tmp2), .8) |>
  add_water(detect_water(tmp2, min_area = 10), color = "imhof1")
tryCatch(dev.off(), error = function(e) NULL)

legend_cols = local({
  max_height = max(tmp2)
  heights = seq(0, max_height, length.out = 500) |>
    rep(5) |>
    matrix(ncol = 5)
  rgb_cols = heights |>
    height_shade() |>
    add_overlay(sphere_shade(heights, texture = "desert", zscale = 10), alphalayer = .5) |>
    add_shadow(lamb_shade(heights, zscale = 10), .8) |>
    add_shadow(texture_shade(heights), .8)
  cols = rgb(rgb_cols[1, , 1], rgb_cols[1, , 2], rgb_cols[1, , 3])
  tryCatch(dev.off(), error = function(e) NULL)
  data.frame(height = heights[, 1], col = cols)
})

df = data.frame(
  red = as.numeric(height_plot[, , 1]),
  green = as.numeric(height_plot[, , 2]),
  blue = as.numeric(height_plot[, , 3]),
  y = rep(st_get_dimension_values(tmp, "y"), dim(tmp)[1]),
  x = rep(st_get_dimension_values(tmp, "x"), each = dim(tmp)[2]))

plot = ggplot(df) +
  geom_point(
    data = data.frame(
      x = st_get_dimension_values(tmp, "x")[10],
      y = st_get_dimension_values(tmp, "y")[10],
      z = legend_cols$height),
    aes(x = x, y = y, col = z)) +
  geom_raster(aes(x = x, y = y), fill = rgb(red = df$red, green = df$green, blue = df$blue)) +
  geom_sf(data = boundary) +
  theme_light() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_line(colour = "black", linewidth = .001, color = "gray", linetype = "dashed"),
    panel.ontop = TRUE,
    panel.background = element_blank()) +
  geom_sf(data = rissa, shape = 17, size = 1.5) +
  geom_sf(data = s0, size = 1.5) +
  geom_sf_text(data = dplyr::mutate(st_as_sf(s0), label = 1:5), aes(label = label),
               nudge_y = 2000, size = 4) +
  geom_sf_text(data = dplyr::mutate(st_as_sf(rissa), name = "Rissa"), aes(label = name),
               nudge_y = 3000, size = 6) +
  scale_color_gradientn(colours = legend_cols$col) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Easting", y = "Northing", col = "Altitude [m]")

ggsave(
  filename = "data_domain.jpg",
  width = 9.5 * .75,
  path = image_dir(),
  height = 6 * .75,
  dpi = 200,
  plot = plot)
