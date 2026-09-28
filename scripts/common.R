# Entrypoints are run from the repository root; no machine-specific paths.
if (!file.exists("DESCRIPTION") || !dir.exists("stan")) stop("Run from the repository root.")
required <- c("rstan", "posterior", "digest", "jsonlite", "rmarkdown", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing project dependencies; run Rscript scripts/setup.R.")
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
arg_value <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) stop(paste("Missing value for", flag))
  args[i + 1]
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
