library(stars)
library(extremePrecipitation)

# ==============================================================================
# Download radar data.
# This script requires that you have installed the CDO program on your computer.
# See https://code.mpimet.mpg.de/projects/cdo for info and installation of the program.
# ==============================================================================

# Find all valid urls containing radar data
dates = seq(as.Date("2010-01-01"), as.Date("2022-12-31"), by = "1 month")
dates = dates[lubridate::month(dates) %in% c(6, 7, 8)]
urls = unlist(lapply(dates, get_radar_url))

# Arguments to CDO for downloading the data
bbox = c(xmin = 9, ymin = 63.4, xmax = 11.5, ymax = 64.6)
args = c(
  paste0("-sellonlatbox,", bbox[1], ",", bbox[3], ",", bbox[2], ",", bbox[4]),
  # This is the name of the variable that contains hourly mean precipitation estimates
  "-selname,lwe_precipitation_rate")

filename = file.path(downloads_dir(), "radar.nc")

# This will take some time, and may need to be restarted several times if you have
# a bad internet connection
download(urls, args, filename)

# ============================================================================================
# Download DEM
# ============================================================================================
filename = file.path(downloads_dir(), "dem.tif")
url = "https://hoydedata.no/LaserServices/REST/DownloadFile.ashx?id=56"
zipfile = file.path(downloads_dir(), "dem.zip")
download.file(url, zipfile)
dem_dir = file.path(downloads_dir(), "dem")
unzip(zipfile, exdir = dem_dir)
file.remove(zipfile)

files = list.files(dem_dir, full.names = TRUE)
tif_files = grep("*.tif$", files, value = TRUE)

# Combine all the tif_files into one big raster file
stars::st_mosaic(tif_files, dst = filename)
