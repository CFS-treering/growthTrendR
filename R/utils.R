#' Check optional dependencies
#' @keywords internal
check_optional_deps <- function() {
  required_suggests <- c("ggforce", "ggeffects", "geosphere", "scales", "purrr", "future", "furrr", "parallel", "MuMIn", "MASS", "rmarkdown", "curl",
                         "rstudioapi", "gstat", "sf", "sp", "spdep", "magick", "pryr")  # Remove "other_package"

  missing <- required_suggests[!sapply(required_suggests, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Required packages missing: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0("'", missing, "'", collapse = ", "), "))"
    )
  }
}
