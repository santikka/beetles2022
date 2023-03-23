#############
# MCMC Runs #
#############

set.seed(1000)
n_cores <- detectCores()
n_iter <- 2600000
n_burnin <- 100000
n_chains <- 12
thin <- 2500

# 24-hour Survival
{
  cl <- makeCluster(n_cores)
  surv_fit_24h <- parLapply(
    cl,
    X = 1:n_chains,
    fun = surv_mcmc,
    mcmc_const = surv_const_24h,
    mcmc_data = surv_data_24h,
    n_iter = n_iter,
    n_burnin = n_burnin,
    thin = thin
  )
  stopCluster(cl)
}

# Larva-to-adult survival
{
  cl <- makeCluster(n_cores)
  surv_fit_la <- parLapply(
    cl,
    X = 1:n_chains,
    fun = surv_mcmc,
    mcmc_const = surv_const_la,
    mcmc_data = surv_data_la,
    n_iter = n_iter,
    n_burnin = n_burnin,
    thin = thin
  )
  stopCluster(cl)
}

# Body mass (females)
{
  cl <- makeCluster(n_cores)
  bm_fit_f <- parLapply(
    cl,
    X = 1:n_chains,
    fun = bm_mcmc,
    mcmc_const = bm_const_f,
    mcmc_data = bm_data_f,
    n_iter = n_iter,
    n_burnin = n_burnin,
    thin = thin
  )
  stopCluster(cl)
}

# Body mass (males)
{
  cl <- makeCluster(n_cores)
  bm_fit_m <- parLapply(
    cl,
    X = 1:n_chains,
    fun = bm_mcmc,
    mcmc_const = bm_const_m,
    mcmc_data = bm_data_m,
    n_iter = n_iter,
    n_burnin = n_burnin,
    thin = thin
  )
  stopCluster(cl)
}

# Reproduction
{
  cl <- makeCluster(n_cores)
  repro_fit <- parLapply(
    cl,
    X = 1:n_chains,
    fun = repro_mcmc,
    mcmc_const = repro_const,
    mcmc_data = repro_data,
    n_iter = n_iter,
    n_burnin = n_burnin,
    thin = thin
  )
  stopCluster(cl)
}

# Combine the posterior samples from different chains
surv_samples_24h <- do.call(rbind, surv_fit_24h)
surv_samples_la <- do.call(rbind, surv_fit_la)
bm_samples_f <- do.call(rbind, bm_fit_f)
bm_samples_m <- do.call(rbind, bm_fit_m)
repro_samples <- do.call(rbind, repro_fit)

# Create predictive samples
surv_nimble_24h <- surv_mcmc(seed = 0, surv_const_24h, surv_data_24h, run = FALSE)
surv_nimble_la <- surv_mcmc(seed = 0, surv_const_la, surv_data_la, run = FALSE)
bm_nimble_f <- bm_mcmc(seed = 0, bm_const_f, bm_data_f, run = FALSE)
bm_nimble_m <- bm_mcmc(seed = 0, bm_const_m, bm_data_m, run = FALSE)
repro_nimble <- repro_mcmc(seed = 0, repro_const, repro_data, run = FALSE)
surv_pp_samples_24h <- get_predictive_sample(surv_nimble_24h, surv_samples_24h)
surv_pp_samples_la <- get_predictive_sample(surv_nimble_la, surv_samples_la)
bm_pp_samples_f <- get_predictive_sample(bm_nimble_f, bm_samples_f)
bm_pp_samples_m <- get_predictive_sample(bm_nimble_m, bm_samples_m)
repro_pp_samples <- get_predictive_sample(repro_nimble, repro_samples)

# Save all samples into a .RData file
# NOTE: file size is may be large! (~ 1 Gb with the default MCMC sample sizes)
save(
  surv_fit_24h, surv_samples_24h, surv_pp_samples_24h,
  surv_fit_la, surv_samples_la, surv_pp_samples_la,
  bm_fit_f, bm_samples_f, bm_pp_samples_f,
  bm_fit_m, bm_samples_m, bm_pp_samples_m,
  repro_fit, repro_samples, repro_pp_samples,
  file = "samples.RData"
)
