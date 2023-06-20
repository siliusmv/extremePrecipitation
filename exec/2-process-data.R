devtools::load_all()
library(dplyr)
library(sf)
library(lubridate)
library(raster)

# Load radar data
radar_file = file.path(downloads_dir(), "radar.nc")
data = raster::brick(radar_file)

# Coordinates of the data
coords = raster::rasterToPoints(data[[1]], spatial = TRUE) |>
  sf::st_as_sf() |>
  {\(x) cbind(x, sf::st_coordinates(x))}() |>
  dplyr::select(geometry)

# Transform from [m] to [km] to remove annoying zeros
myproj = st_crs(coords)[["input"]] |>
  sub("units=m", "units=km", x = _) |>
  st_crs()
coords = st_transform(coords, myproj)

# Add covariates
height_raster = raster::raster(file.path(downloads_dir(), "dem.tif"))
transformed_coords = st_coordinates(st_transform(coords, st_crs(height_raster)))
coords$height = raster::extract(height_raster, transformed_coords)
coords$height = ifelse(is.na(coords$height), 0, coords$height)
coords$height[coords$height < 0] = 0 # a few locations get negative heights

# Extract the times of all observations
time_series = data[1, 1]
times = colnames(time_series) |>
  sub("X", "", x = _) |>
  lubridate::ymd_hms()

# Check if the time series is monotonely increasing
time_diff = tail(times, -1) - head(times, -1)
min(time_diff)

rissa = st_point(c(10.203845, 63.690527)) |>
  st_sfc(crs = 4326) |>
  st_transform(st_crs(coords))

# Extract the data into a matrix
data = raster::as.matrix(data)
colnames(data) = NULL

# Reorder so time is the first dimension. This makes it faster to extract for any given location
data = t(data)

# Only keep observations from inside this box
tmp_coords = st_coordinates(coords)
loc_index = which(
  tmp_coords[, 1] >= 230
  & tmp_coords[, 1] <= 320
  & tmp_coords[, 2] >= 7070
  & tmp_coords[, 2] <= 7140)
data = data[, loc_index]
coords = coords[loc_index, ]

# Save the processed data
radar_data = list(
  data = data,
  coords = coords,
  times = times,
  rissa = rissa)

radar_data$year = lubridate::year(radar_data$times)
radar_data$month = lubridate::month(radar_data$times)
radar_data$week = lubridate::week(radar_data$times)
radar_data$day = lubridate::yday(radar_data$times) -
  ifelse(lubridate::leap_year(radar_data$times), 1, 0)

saveRDS(radar_data, file = file.path(downloads_dir(), "radar.rds"))
