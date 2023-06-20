#' Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
#' using the tikzDevice package
#' @export
plot_tikz = function(plot = NULL, expression = NULL, file = "Rplots.pdf", ...) {
  expression = substitute(expression)
  if (is.null(plot) && is.null(expression)) {
    stop("Either `plot` or `expression` must be non-NULL")
  }

  # Ensure that you are on an operating system that you have tested
  operating_system = Sys.info()[["sysname"]]
  if (operating_system == "Windows") {
    proceed = readline(paste("This function was written on a Mac,",
                             "I have no idea if it will work on Windows.",
                             "Proceed? (y/n) "))
    if (proceed != "y") return()
  }

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(suppressWarnings(file.remove(tmp)), add = TRUE)

  # Extract default tex usepackages and add the bm and amsmath packages
  opt = options()
  on.exit(options(opt)) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)
  tex_packages = options()$tikzLatexPackages
  if (!any(grepl("usepackage\\{bm\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{bm}\n")
  }
  if (!any(grepl("usepackage\\{amsmath\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{amsmath}\n")
  }
  if (!any(grepl("usepackage\\{amssymb\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{amssymb}\n")
  }

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit(dev.off(current_device))

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(expression)
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf using lualatex
  system2("lualatex", shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  file.remove(files_to_clean)
}

#' @export
latex_friendly_map_plot = function(x) {
  info = ggplot2::ggplot_build(x)$layout$panel_params[[1]]$graticule
  east_ticks = info$degree[info$type == "E" & info$y_start == 0]
  north_ticks = info$degree[info$type == "N" & info$x_start == 0]
  x +
    ggplot2::scale_x_continuous(breaks = east_ticks, labels = paste0(east_ticks, "$^\\circ$E"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = north_ticks, labels = paste0(north_ticks, "$^\\circ$N"), expand = c(0, 0))
}

#' @export
add_norway_map = function(crs = NULL, bbox = NULL, ...) {
  map = rnaturalearth::ne_countries(scale = 50, country = "Norway", returnclass = "sf")
  c(ggplot2::layer_sf(
    geom = ggplot2::GeomSf, data = map, mapping = ggplot2::aes(),
    stat = "sf", position = "identity", show.legend = NA,
    inherit.aes = TRUE,
    params = list(na.rm = FALSE, fill = NA, ...)),
    ggplot2::coord_sf(default = TRUE, crs = crs,
                      xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)]))
}
