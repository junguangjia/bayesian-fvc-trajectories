#!/usr/bin/env Rscript
# Run from the repository root. Installs only into this project's libraries.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!file.exists("DESCRIPTION")) stop("Run from the repository root.")
if (!requireNamespace("renv", quietly = TRUE)) {
  dir.create(".bootstrap", showWarnings = FALSE)
  install.packages("renv", lib = ".bootstrap")
  .libPaths(c(normalizePath(".bootstrap"), .libPaths()))
}
if (file.exists("renv.lock")) {
  renv::restore(library = renv::paths$library(), prompt = FALSE)
} else {
  sources <- .libPaths()
  renv::init(bare = TRUE, restart = FALSE)
  renv::settings$snapshot.type("explicit")
  renv::hydrate(packages = c("rstan", "StanHeaders", "posterior", "digest", "jsonlite", "rmarkdown", "knitr"),
                sources = sources, prompt = FALSE)
  renv::snapshot(prompt = FALSE)
}
cat("Project library ready. A C++ toolchain and pdflatex are required for sampling and PDF reports.\n")
