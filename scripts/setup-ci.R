#!/usr/bin/env Rscript
# Linux CI: compile the pinned Stan/TBB stack together instead of combining
# prebuilt packages linked against potentially different TBB ABIs.
if (!identical(Sys.info()[['sysname']], 'Linux')) stop('This bootstrap is for Linux CI.')
if (!requireNamespace('renv', quietly = TRUE)) stop('The project renv bootstrap did not load.')
lock <- renv::lockfile_read('renv.lock')
binary_repo <- Sys.getenv('RSPM')
if (!nzchar(binary_repo)) stop('setup-r must supply its Linux binary repository.')
Sys.setenv(RENV_CONFIG_REPOS_OVERRIDE = binary_repo,
           RENV_CONFIG_INSTALL_JOBS = '1')
options(renv.config.repos.override = binary_repo,
        renv.config.install.jobs = 1L)
stan_stack <- c('RcppParallel', 'StanHeaders', 'rstan')
stopifnot(all(stan_stack %in% names(lock$Packages)))
renv::restore(packages = setdiff(names(lock$Packages), stan_stack),
              prompt = FALSE, retry = FALSE)
Sys.setenv(RENV_CONFIG_REPOS_OVERRIDE = 'https://cloud.r-project.org',
           RENV_CONFIG_PPM_ENABLED = 'FALSE')
Sys.unsetenv(c('TBB_ROOT', 'TBB_INC', 'TBB_LIB'))
options(pkgType = 'source', repos = c(CRAN = 'https://cloud.r-project.org'),
        renv.config.ppm.enabled = FALSE,
        renv.config.repos.override = 'https://cloud.r-project.org')
renv::restore(packages = stan_stack, rebuild = stan_stack,
              prompt = FALSE, retry = FALSE)
for (package in names(lock$Packages)) {
  stopifnot(packageVersion(package) == numeric_version(lock$Packages[[package]]$Version))
}
stopifnot(requireNamespace('rstan', quietly = TRUE))
cat('Pinned project environment and source-built Stan stack are ready.\n')
