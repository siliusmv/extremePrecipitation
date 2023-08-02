
#' Given a vector of urls from https://thredds.met.no,
#' download data using the arguments args,
#' then concatenate all the data and save it into the file filename.
#' Underways, all the data are downloaded in a temporary directory.
#' Since it may take a long time to download all the data, and connection
#' issues may arise, the function only downloads data from a certain url
#' if the tempfile for that url does not already exist
#' @export
download = function(urls, args, filename) {
  # Create a temporary directory for storing the data from each url
  temp_dir = sub(".[^.]+$", "", filename)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  # Create one temporary filename for each url
  tempfiles = file.path(temp_dir, paste0(seq_along(urls), ".nc"))

  # Download the data
  for (i in seq_along(urls)) {
    # Check if the data already has been downloaded
    file_already_exists = file.exists(tempfiles[i])
    if (!file_already_exists) {
      execute_shell_script(
        command = "cdo",
        args = c(args, urls[i], shQuote(tempfiles[i])))
    }
  }

  # Concatenate all the tempfiles and save them into filename
  execute_shell_script(
    command = "cdo",
    args = c("-cat", tempfiles, shQuote(filename)))
}

#' For a given date on a format that can be understood by lubridate,
#' return a vector of urls containing all the radar data archives
#' for that month
#' @export
get_radar_url = function(date) {
  catalog_url = paste0(
    "https://thredds.met.no/thredds/catalog/remotesensingradaraccr/",
    lubridate::year(date), "/",
    sprintf("%02d", lubridate::month(date)),
    "/catalog.html")
  files = locate_nc_files_in_metno_page(catalog_url)
  files = sort(files)
  get_full_url_of_metno_file(catalog_url, files)
}

#' Given the url to a thredds.met.no webpage,
#' locate all files in the webpage that satisfy certain requirements
locate_nc_files_in_metno_page = function(url) {
  page = xml2::read_html(url)
  links = rvest::html_nodes(page, "a")
  links = rvest::html_text(links)

  # Only extract the nc-files
  links = grep(".nc$", links, value = TRUE)

  # These files are bad and unwanted
  links = grep(".nc.nc$", links, invert = TRUE, value = TRUE)
  links = grep("nohur", links, invert = TRUE, value = TRUE)

  links
}

get_full_url_of_metno_file = function(path, filename) {
  if (length(filename) == 0 || length(path) == 0) {
    warning("invalid metno_url")
    return(NULL)
  }
  path = sub("catalog/", "dodsC/", path)
  path = sub("catalog.html", "", path)
  paste0(path, filename)
}
