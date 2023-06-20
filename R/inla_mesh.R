#' This function is a wrapper around INLA::inla.mesh.2d().
#' For some input combinations to inla.mesh.2d(), the function will never
#' manage to create a triangulated mesh, and will instead contnue running indefinitely,
#' meaning that the function must be actively stopped by the user.
#' The function create_mesh() ensures that this does not happen by automatically
#' killing the subprocess started by inla.mesh.2d() and returning NULL after
#' a certain time has passed.
#' NB: Even though the subprocess started by inla.mesh.2d() is killed, the subprocesses
#' that are again started by this subprocess are not killed. This must be performed
#' manually by the user, by, e.g., calling `system("pkill fmesher")`, if you are
#' working on a *nix system.
#' create_mesh() has two arguments:
#' - args is a names list of arguments to use in the function inla.mesh.2d()
#' - timeout is the time, in seconds, to wait before killing the subprocess started by
#'   inla.mesh.2d()
#' @export
create_mesh = function(args, timeout = 5) {
  # Fork the process to create a child process that can create the mesh
  p = parallel:::mcfork()

  if (inherits(p, "masterProcess")) { # This is only TRUE for the child process
    mesh = do.call(INLA::inla.mesh.2d, args) # create the mesh
    mesh_file = tempfile() # Create a tempfile for saving the mesh
    saveRDS(mesh, mesh_file)
    parallel:::mcexit(0L, mesh_file) # Send the filename to the parent process, and exit
  }

  # Everything below is read only by the parent process
  # Wait `timeout` seconds for a message from the child process
  stdout = parallel:::readChildren(timeout = timeout)

  if (is.raw(stdout)) { # This means that we received a message
    child_pid = attr(stdout, "pid") # Read the process ID of the created child
    mesh_file = unserialize(stdout) # Read the filename of the created mesh

    # Try to kill the child process.
    # For some reason, mcexit() does not always seem to close the process,
    # so we do it once more to be on the safe side
    for (i in 1:5) tryCatch(parallel:::rmChild(child_pid), error = function(e) NULL)

    # Read the created mesh and return it
    mesh = readRDS(mesh_file)
    file.remove(mesh_file)
    return(mesh)
  }

  # Kill the child process and return NULL since it never sent the name of the mesh file
  for (child in parallel:::children()) {
    for (i in 1:5) tryCatch(parallel:::rmChild(child), error = function(e) NULL)
    tryCatch(parallel:::rmChild(child), error = function(e) NULL)
  }
  NULL
}
