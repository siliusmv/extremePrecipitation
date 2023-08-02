
#' @export
data_dir = function() file.path(here::here(), "raw_data")
#' @export
downloads_dir = function() file.path(data_dir(), "downloads")
#' @export
image_dir = function() file.path(data_dir(), "images")
#' @export
results_dir = function() file.path(data_dir(), "results")
#' @export
cgeneric_dir = function() file.path(here::here(), "cgeneric")

#' Compute the Euclidean distance between x and y, where
#' x is either a row vector, or a matrix of row vectors,
#' and y is a matrix of row vectors
#' @export
dist_euclid = function(x, y) {
  if (!is.matrix(y)) y = matrix(y, nrow = 1)
  if (is.matrix(x)) {
    res = dist_euclid_mat_mat(x, y)
    if (!is.null(dim(res)) && (nrow(res) == 1 || ncol(res) == 1)) res = as.numeric(res)
  } else {
    res = dist_euclid_vec_mat(x, y)
  }
  res
}

#' Execute a shell script, and provide feedback if it crashes
#' @export
execute_shell_script = function(command, args = character(), ...) {
  output = system2(command, args, ...)
  success = (output == 0)
  if (!success) {
    formatted_args = paste(args, collapse = " ")
    stop("shell script ", command,
         " gave error code ", output,
         " with arguments ", formatted_args)
  }
  0
}

#' Call a makefile used for compiling and linking cgeneric scripts,
#' in order to use them with R-INLA
#' @export
make_cgeneric = function(cmd) {
  current_path = getwd()
  on.exit(setwd(current_path))
  setwd(cgeneric_dir())
  execute_shell_script("make", cmd)
}

#' Create a progress bar for tracking long-running processes.
#' This is a thin wrapper around the progress package
#' @export
progress_bar = function(n) {
  pb = progress::progress_bar$new(
    format = ":percent [:bar] time elapsed: :elapsedfull, eta: :eta",
    total = n, width = 70, clear = FALSE)
  # pb$tick() throws an error if we tick too many times, which can potentially stop
  # a script that is otherwise working fine. We don't want that to happen,
  # so we change the tick function slightly
  res = list(
    terminate = pb$terminate,
    tick = function(...) tryCatch(pb$tick(...), error = function(e) NULL))
  res
}

#' Necessary for the dist_euclid() function
dist_euclid_vec_mat = function(x, y) {
  stopifnot(length(x) == ncol(y))
  res = rep(0, nrow(y))
  for (i in seq_len(ncol(y))) {
    res = res + (y[, i] - x[i])^2
  }
  sqrt(res)
}

#' Necessary for the dist_euclid() function
dist_euclid_mat_mat = function(x, y) {
  apply(x, 1, dist_euclid_vec_mat, y = y)
}
