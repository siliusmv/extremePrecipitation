library(extremePrecipitation)
library(dplyr)
library(sf)
library(lubridate)
library(raster)
library(ncdf4)
library(stars)

# Load radar data
radar_file = file.path(downloads_dir(), "radar.nc")
data = raster::brick(radar_file)

# Coordinates of the data
coords = raster::rasterToPoints(data[[1]], spatial = TRUE) |>
  sf::st_as_sf() |>
  (\(x) cbind(x, sf::st_coordinates(x)))() |>
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
  sub("(\\d{4}.\\d{2}.\\d{2})$", "\\1.00.00.00", x = _) |>
  lubridate::ymd_hms()

# Check if the time series is monotonely increasing
time_diff = tail(times, -1) - head(times, -1)
min(time_diff)

# This is the location of the Rissa radar
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

# Define the list of our five chosen conditioning sites
s0 = list(
  st_point(c(277000, 7114000)),
  st_point(c(294000, 7113000)),
  st_point(c(268000, 7101000)),
  st_point(c(285000, 7099000)),
  st_point(c(255000, 7092000))) |>
  st_as_sfc()
st_crs(s0) = st_crs(height_raster)

# Transform the conditioning sites to be on the same
# projection as the radar data
s0 = st_transform(s0, st_crs(coords))

# Save the processed data
radar_data = list(
  data = data,
  coords = coords,
  times = times,
  rissa = rissa,
  s0 = s0,
  year = lubridate::year(times),
  month = lubridate::month(times),
  week = lubridate::week(times),
  day = lubridate::yday(times) - ifelse(lubridate::leap_year(times), 1, 0)
)

saveRDS(radar_data, file = file.path(downloads_dir(), "radar.rds"))
