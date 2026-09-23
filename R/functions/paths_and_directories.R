# Ensures output directories exist (creates parent directories as needed).
#
# Created: 1/25/2026
# Author: X. Takver

ensure_dir <- function(path) {
  # Create directory only if it does not already exist
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  # Return the path invisibly for piping/assignment convenience
  invisible(path)
}
