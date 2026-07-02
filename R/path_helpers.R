repo_root <- function() {
  root <- getOption("mnp_spatial.repo_root", unset = NA_character_)
  if (!is.na(root)) return(normalizePath(root, mustWork = TRUE))
  normalizePath(getwd(), mustWork = TRUE)
}

repo_path <- function(...) {
  file.path(repo_root(), ...)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
