# Source ".Rprofile_local" if it exists. Used for machine-specific configurations,
# e.g., setting a local conda environment. This file should not be committed to version control.

if (file.exists(".Rprofile_local")) {
    source(".Rprofile_local")
}


# If in an interactive R session, silently load essential development packages.
# This helps streamline the development workflow without flooding the console with messages.

if (interactive()) {
    require(devtools)
    require(usethis)
    require(pkgdown)
    require(testthat)
}
