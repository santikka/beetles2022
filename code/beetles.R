# Session info of the system used to run the simulation:

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
#
# Matrix products: default
#
# locale:
# [1] LC_COLLATE=Finnish_Finland.1252  LC_CTYPE=Finnish_Finland.1252
# [3] LC_MONETARY=Finnish_Finland.1252 LC_NUMERIC=C
# [5] LC_TIME=Finnish_Finland.1252
#
# attached base packages:
# [1] splines   parallel  stats   graphics  grDevices utils  datasets  methods
# [9] base
#
# other attached packages:
# [1] forcats_0.5.1    stringr_1.4.0    tidyr_1.2.0      dplyr_1.0.8
# [5] coda_0.19-4      nimble_0.12.2    tibble_3.1.6     ggnewscale_0.4.7
# [9] ggrepel_0.9.1    ggplot2_3.3.5
#
# loaded via a namespace (and not attached):
#  [1] igraph_1.2.11    Rcpp_1.0.8.3     magrittr_2.0.2   tidyselect_1.1.2
#  [5] munsell_0.5.0    lattice_0.20-45  colorspace_2.0-3 R6_2.5.1
#  [9] rlang_1.0.2      fansi_1.0.2      tools_4.1.2      grid_4.1.2
# [13] gtable_0.3.0     utf8_1.2.2       cli_3.2.0        withr_2.5.0
# [17] ellipsis_0.3.2   lifecycle_1.0.1  crayon_1.5.0     purrr_0.3.4
# [21] vctrs_0.3.8      glue_1.6.2       stringi_1.7.6    compiler_4.1.2
# [25] pillar_1.7.0     generics_0.1.2   scales_1.1.1     pkgconfig_2.0.3

rm(list = ls())

library(coda)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(forcats)
library(moments)
library(nimble)
library(parallel)
library(splines)
library(stringr)
library(tidyr)
library(tibble)


# Data and Models ---------------------------------------------------------

data_path <- "./data/"
source("code/functions.R")
source("code/data.R")
source("code/mcmc.R")


# MCMC Simulation ---------------------------------------------------------

# Run the simulations
# source("code/runs.R")
#
# Or load existing samples
# load("data/samples.RData")


# Diagnostics -------------------------------------------------------------

# Transformed family effects f are filtered out from
# the PSRF and ESS calculations because the samples
# also contain  the raw effects f_raw.

# 24h Survival
mcmc_psrf(surv_fit_24h, filter = "f[")
mcmc_ess(surv_fit_24h, filter = "f[") |> min()
mcmc_summary(surv_samples_24h, filter = "f_raw") |> print(n = Inf)

# Larva to adult survival
mcmc_psrf(surv_fit_la, filter = "f[")
mcmc_ess(surv_fit_la, filter = "f[") |> min()
mcmc_summary(surv_samples_la, filter = "f_raw") |> print(n = Inf)

# Body mass (females)
mcmc_psrf(bm_fit_f, filter = "f[")
mcmc_ess(bm_fit_f, filter = "f[") |> min()
mcmc_summary(bm_samples_f, filter = "f_raw") |> print(n = Inf)

# Body mass (males)
mcmc_psrf(bm_fit_m, filter = "f[")
mcmc_ess(bm_fit_m, filter = "f[") |> min()
mcmc_summary(bm_samples_m, filter = "f_raw") |> print(n = Inf)

# Reproduction
mcmc_psrf(repro_fit, filter = "phi[5") # last phi is removed due to sum constraint
mcmc_ess(repro_fit) |> min()
mcmc_summary(repro_samples) |> print(n = Inf)


# Figures and Tables ------------------------------------------------------

source("code/results.R")
source("code/figtab.R")
