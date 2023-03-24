#################
# Nimble Models #
#################

surv_mcmc <- function(seed, mcmc_const, mcmc_data,
                      n_iter = 10000, n_burnin = 5000, thin = 1, run = TRUE) {
  library(nimble)
  code <- nimbleCode({
    # Likelihood
    trt_spline[1:N_c] <- sp_basis_c[1:N_c, 1:S_c] %*% sp_c[1:S_c]
    trt_spline[(N_c + 1):N] <- sp_basis_t[1:N_t, 1:S_t] %*% sp_t[1:S_t]
    for (i in 1:N) {
      logit(p[i]) <- beta[1] * pa_trt[i] + beta[2] * age4[i] + beta[3] * age5[i] + trt_spline[i] + f[fam[i]]
      surv[i] ~ dbern(p[i])
    }
    # Priors
    for (j in 1:3) {
      beta[j] ~ dnorm(0, 0.5)
    }
    # B-spline parameters
    log_sigma_sp_c ~ dnorm(sigma_sp_init_c, 3)
    log_sigma_sp_t ~ dnorm(sigma_sp_init_t, 3)
    sp_c[1] ~ dnorm(sp_init_c[1], 2)
    sp_t[1] ~ dnorm(sp_init_t[1], 2)
    for (j in 2:S_c) {
      sp_c[j] ~ dnorm(sp_c[j - 1], sd = exp(log_sigma_sp_c))
    }
    for (j in 2:S_t) {
      sp_t[j] ~ dnorm(sp_t[j - 1], sd = exp(log_sigma_sp_t))
    }
    # Family effects
    rho_f ~ dunif(0, 1)
    tau_f[1] ~ dgamma(1, 1)
    tau_f[2] ~ dgamma(1, 1)
    for (r in 1:(fam_max_size - 1)) {
      for (c in (r + 1):fam_max_size) {
        R[r, c] <- rho_f
        R[c, r] <- rho_f
        tau_diag[r, c] <- 0
        tau_diag[c, r] <- 0
      }
      tau_diag[r + 1, r + 1] <- sqrt(tau_f[2])
      R[r, r] <- 1
      zeros[r] <- 0
    }
    tau_diag[1, 1] <- sqrt(tau_f[1])
    R[fam_max_size, fam_max_size] <- 1
    zeros[fam_max_size] <- 0
    for (j in 1:N_fam_tg) {
      P[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], j] <-
        tau_diag[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]]] %*%
        R[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]]] %*%
        tau_diag[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]]]
      f_raw[1:fam_size[fam_tg[j]], fam_tg[j]] ~ dmnorm(zeros[1:fam_size[fam_tg[j]]],
        prec = P[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], j]
      )
      for (l in 1:fam_size[fam_tg[j]]) {
        f[fam_mat[fam_tg[j], l]] <- f_raw[l, fam_tg[j]]
      }
    }
    for (j in 1:N_fam_g) {
      f_raw[1, fam_g[j]] ~ dnorm(0, tau_f[1])
      f[fam_mat[fam_g[j], 1]] <- f_raw[1, fam_g[j]]
    }
  })
  fam_dims_vec <- c(mcmc_const$fam_max_size, mcmc_const$N_fam1)
  fam_dims_mat <- c(mcmc_const$fam_max_size, fam_dims_vec)
  const_dims_mat <- rep(mcmc_const$fam_max_size, 2)
  f_raw_init <- array(0, dim = fam_dims_vec)
  sapply(1:mcmc_const$N_fam1, function(i) {
    f_raw_init[1:mcmc_const$fam_size[i], i] <- rnorm(mcmc_const$fam_size[i], 0, 0.01)
  })
  mod <- nimbleModel(
    code = code,
    constants = mcmc_const,
    data = mcmc_data,
    dimensions = list(
      f_raw = fam_dims_vec,
      P = fam_dims_mat,
      R = const_dims_mat,
      tau_diag = const_dims_mat,
      zeros = mcmc_const$fam_max_size
    ),
    inits = list(
      beta = rnorm(3, 0, 0.01),
      f_raw = f_raw_init,
      rho_f = runif(1, 0, 0.1),
      log_sigma_sp_c = rnorm(1, mcmc_const$sigma_sp_init_c, 0.01),
      log_sigma_sp_t = rnorm(1, mcmc_const$sigma_sp_init_t, 0.01),
      sp_c = rnorm(mcmc_const$S_c, mcmc_const$sp_init_c, 0.01),
      sp_t = rnorm(mcmc_const$S_t, mcmc_const$sp_init_t, 0.01),
      tau_f = 1 + rexp(2)
    )
  )
  conf <- configureMCMC(
    mod,
    monitors = c(
      "beta", "f", "f_raw", "rho_f",
      "log_sigma_sp_c", "log_sigma_sp_t",
      "sp_c", "sp_t", "tau_f"
    )
  )
  conf$removeSampler("sp_c[]")
  conf$removeSampler("sp_t[]")
  conf$removeSampler("log_sigma_sp_c")
  conf$removeSampler("log_sigma_sp_t")
  conf$addSampler(
    target = c("log_sigma_sp_c", paste0("sp_c[1:", mcmc_const$S_c, "]")),
    type = "RW_block",
    control = mcmc_const$control[[1]]
  )
  conf$addSampler(
    target = c("log_sigma_sp_t", paste0("sp_t[1:", mcmc_const$S_t, "]")),
    type = "RW_block",
    control = mcmc_const$control[[2]]
  )
  Rmcmc <- buildMCMC(conf)
  Cmod <- compileNimble(mod)
  if (run) {
    Cmcmc <- compileNimble(Rmcmc)
    runMCMC(
      Cmcmc,
      niter = n_iter,
      nburnin = n_burnin,
      nchains = 1,
      thin = thin,
      setSeed = seed
    )
  } else {
    list(model = mod, Cmodel = Cmod, mcmc = Rmcmc)
  }
}

bm_mcmc <- function(seed = 0, mcmc_const, mcmc_data,
                    n_iter = 10000, n_burnin = 5000, thin = 1, run = TRUE) {
  library(nimble)
  code <- nimbleCode({
    # Likelihood
    for (i in 1:N) {
      mu[i, 1] <- beta[1, 1] + beta[2, 1] * trt[i] + beta[3, 1] * pa_trt[i] +
        beta[4, 1] * ti[i] + beta[5, 1] * age4[i] + beta[6, 1] * age5[i] + f[fam[i], 1]
      mu[i, 2] <- beta[1, 2] + beta[2, 2] * trt[i] + beta[3, 2] * pa_trt[i] +
        beta[4, 2] * ti[i] + beta[5, 2] * age4[i] + beta[6, 2] * age5[i] + beta[7, 2] * bm[i, 1] + f[fam[i], 2]
      bm[i, 1] ~ dnorm(mu[i, 1], tau_bm[gen[i], 1])
      bm[i, 2] ~ dnorm(mu[i, 2], tau_bm[gen[i], 2])
    }
    # Priors
    beta[1, 1] ~ dnorm(11, 1)
    beta[1, 2] ~ dnorm(5.5, 1)
    beta[7, 1] <- 0
    beta[7, 2] ~ dnorm(0.7, 1)
    for (k in 1:2) {
      for (j in 2:6) {
        beta[j, k] ~ dnorm(0, 1)
      }
      tau_bm[1, k] ~ dgamma(1, 0.5)
      tau_bm[2, k] ~ dgamma(1, 0.5)
      # Family effects
      rho_f[k] ~ dunif(0, 1)
      tau_f[1, k] ~ dgamma(1, 0.1)
      tau_f[2, k] ~ dgamma(1, 0.1)
      for (r in 1:(fam_max_size - 1)) {
        for (c in (r + 1):fam_max_size) {
          R[r, c, k] <- rho_f[k]
          R[c, r, k] <- rho_f[k]
          tau_diag[r, c, k] <- 0
          tau_diag[c, r, k] <- 0
        }
        tau_diag[r + 1, r + 1, k] <- sqrt(tau_f[2, k])
        R[r, r, k] <- 1
        zeros[r, k] <- 0
      }
      tau_diag[1, 1, k] <- sqrt(tau_f[1, k])
      R[fam_max_size, fam_max_size, k] <- 1
      zeros[fam_max_size, k] <- 0
      for (j in 1:N_fam_tg) {
        P[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], j, k] <-
          tau_diag[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], k] %*%
          R[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], k] %*%
          tau_diag[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], k]
        f_raw[1:fam_size[fam_tg[j]], fam_tg[j], k] ~ dmnorm(zeros[1:fam_size[fam_tg[j]], k],
          prec = P[1:fam_size[fam_tg[j]], 1:fam_size[fam_tg[j]], j, k]
        )
        for (l in 1:fam_size[fam_tg[j]]) {
          f[fam_mat[fam_tg[j], l], k] <- f_raw[l, fam_tg[j], k]
        }
      }
      for (j in 1:N_fam_g) {
        f_raw[1, fam_g[j], k] ~ dnorm(0, tau_f[1, k])
        f[fam_mat[fam_g[j], 1], k] <- f_raw[1, fam_g[j], k]
      }
    }
  })
  bm_init <- as.matrix(mcmc_data$bm)
  bm_obs <- !is.na(bm_init)
  bm_init[!bm_obs] <- rnorm(sum(!bm_obs), 125, 5)
  bm_init[bm_obs] <- NA
  fam_dims_vec <- c(mcmc_const$fam_max_size, mcmc_const$N_fam1, 2)
  fam_dims_mat <- c(mcmc_const$fam_max_size, fam_dims_vec)
  const_dims_mat <- c(rep(mcmc_const$fam_max_size, 2), 2)
  f_raw_init <- array(0, dim = fam_dims_vec)
  sapply(1:mcmc_const$N_fam1, function(i) {
    f_raw_init[1:mcmc_const$fam_size[i], i, 1:2] <- rnorm(mcmc_const$fam_size[i] * 2)
  })
  mod <- nimbleModel(
    code = code,
    constants = mcmc_const,
    data = mcmc_data,
    dimensions = list(
      zeros = const_dims_mat[-1],
      tau_diag = const_dims_mat,
      R = const_dims_mat,
      P = fam_dims_mat,
      f_raw = fam_dims_vec
    ),
    inits = list(
      beta = rbind(
        rnorm(2, c(11, 5.5), sd = 0.25),
        matrix(rnorm(10, 0, 0.25), ncol = 2),
        c(0, rnorm(1, 0.7, sd = 0.05))
      ),
      bm = bm_init,
      f_raw = f_raw_init,
      rho_f = rbeta(2, 2, 2),
      tau_bm = matrix(rexp(4), ncol = 2),
      tau_f = cbind(
        rnorm(2, 12, 0.25),
        rnorm(2, 8, 0.25)
      )
    )
  )
  conf <- configureMCMC(mod, monitors = c("beta", "f", "f_raw", "rho_f", "tau_bm", "tau_f"))
  conf$removeSampler("beta[]")
  conf$addSampler(target = c("beta[1, 2]", "beta[7, 2]"), type = "AF_slice")
  conf$addSampler(
    target = paste0("beta[", 1:6, ", 1]"),
    type = "RW_block",
    control = mcmc_const$control[[1]]
  )
  conf$addSampler(
    target = paste0("beta[", 2:6, ", 2]"),
    type = "RW_block",
    control = mcmc_const$control[[2]]
  )
  Rmcmc <- buildMCMC(conf)
  Cmod <- compileNimble(mod)
  if (run) {
    Cmcmc <- compileNimble(Rmcmc)
    runMCMC(
      Cmcmc,
      niter = n_iter,
      nburnin = n_burnin,
      nchains = 1,
      thin = thin,
      setSeed = seed
    )
  } else {
    list(model = mod, Cmodel = Cmod, mcmc = Rmcmc)
  }
}

repro_mcmc <- function(seed, mcmc_const, mcmc_data,
                       n_iter = 10000, n_burnin = 5000, thin = 1, run = TRUE) {
  library(nimble)
  dgennorm <- nimbleFunction(
    run = function(x = double(0), location = double(0), scale = double(0),
                   shape = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      log_prob <- log(shape) - log(2 * scale) - lgamma(1 / shape) -
        (abs(x - location) / scale)^shape
      if (log) {
        return(log_prob)
      } else {
        return(exp(log_prob))
      }
    }
  )
  pgennorm <- nimbleFunction(
    run = function(q = double(0), location = double(0), scale = double(0),
                   shape = double(0), lower.tail = integer(0, default = 1),
                   log.p = integer(0, default = 0)) {
      returnType(double(0))
      x <- (abs(q - location) / scale)^shape
      # unnormalized incomplete lower gamma:
      # \gamma(a, b) = pgamma(b, a, lower.tail = TRUE) * gamma(a)
      p <- 0.5 + 0.5 * (-1)^(q < location) * pgamma(x, 1 / shape, lower.tail = TRUE)
      if (!lower.tail) {
        if (log.p) {
          return(log(1 - p))
        } else {
          return(1 - p)
        }
      } else {
        if (log.p) {
          return(log(p))
        } else {
          return(p)
        }
      }
    }
  )
  qgennorm <- nimbleFunction(
    run = function(p = double(0), location = double(0), scale = double(0),
                   shape = double(0), lower.tail = integer(0, default = 1),
                   log.p = integer(0, default = 0)) {
      returnType(double(0))
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      return((-1)^(p < 0.5) * (scale^shape * qgamma(2 * abs(p - 0.5), 1 / shape))^(1 / shape) + location)
    }
  )
  rgennorm <- nimbleFunction(
    run = function(n = integer(0), location = double(0),
                   scale = double(0), shape = double(0)) {
      returnType(double(0))
      u <- runif(1, 0, 1)
      return((-1)^(u < 0.5) * (scale^shape * qgamma(2 * abs(u - 0.5), 1 / shape))^(1 / shape) + location)
    }
  )
  assign("dgennorm", dgennorm, envir = .GlobalEnv)
  assign("rgennorm", rgennorm, envir = .GlobalEnv)
  assign("pgennorm", pgennorm, envir = .GlobalEnv)
  assign("qgennorm", qgennorm, envir = .GlobalEnv)
  code <- nimbleCode({
    for (i in 1:N) {
      # Number of batches
      batches[i] ~ dcat(phi[1:B, gen[i]])
      # Number of eggs
      mu[i, 1] <- beta[1, 1] + beta[2, 1] * trt[i] + beta[3, 1] * pa_trt[i] +
        beta[4, 1] * ti[i] + beta[5, 1] * age4[i] + beta[6, 1] * age5[i] +
        beta[7, 1] * gen2[i] + beta[8, 1] * bm[i] + beta[9, 1] * (batches[i] - 1)
      log(scale_rate[i]) <- kappa[1] + kappa[2] * gen2[i]
      log(shape_rate[i]) <- lambda[1] + lambda[2] * gen2[i]
      rate[i] ~ dgennorm(mu[i, 1], scale_rate[i], shape_rate[i])
      total_trunc[i] <- max(trunc(batches[i] * rate[i]), batches[i])
      tt[i] <- log(total_trunc[i]) / log(2) - log_total_mean
      # Hatching
      logit(mu[i, 2]) <- beta[1, 2] + beta[2, 2] * trt[i] + beta[3, 2] * pa_trt[i] +
        beta[4, 2] * ti[i] + beta[5, 2] * age4[i] + beta[6, 2] * age5[i] + beta[7, 2] * gen2[i] +
        beta[8, 2] * bm[i] + beta[9, 2] * tt[i]
      log(nu[i]) <- theta[1] + theta[2] * trt[i] + theta[3] * pa_trt[i] +
        theta[4] * ti[i] + theta[5] * age4[i] + theta[6] * age5[i] + theta[7] * gen2[i] +
        theta[8] * bm[i] + theta[9] * tt[i]
      a[i] <- mu[i, 2] * nu[i]
      b[i] <- (1 - mu[i, 2]) * nu[i]
      p[i] ~ dbeta(a[i], b[i])
      hatched[i] ~ dbinom(p[i], total_trunc[i])
    }
    # Priors
    for (j in 1:B) {
      phi_init[j] <- 2^j
    }
    phi[1:B, 1] ~ ddirch(phi_init[1:B])
    phi[1:B, 2] ~ ddirch(phi_init[1:B])
    beta[1, 1] ~ dnorm(20, 0.1)
    beta[1, 2] ~ dnorm(0, 0.1)
    for (k in 1:2) {
      for (j in 2:9) {
        beta[j, k] ~ dnorm(0, 0.1)
      }
    }
    kappa[1] ~ dnorm(2, 0.1)
    kappa[2] ~ dnorm(0, 0.1)
    lambda[1] ~ dnorm(1, 0.1)
    lambda[2] ~ dnorm(0, 0.1)
    for (j in 1:9) {
      theta[j] ~ dnorm(0, 0.1)
    }
  })
  # Initial values for missing data
  batches_init <- mcmc_data$batches
  batches_obs <- !is.na(batches_init)
  batches_init[!batches_obs] <- sample(mcmc_const$B, sum(!batches_obs), TRUE)
  rate_init <- mcmc_data$rate
  rate_obs <- !is.na(rate_init)
  rate_init[!rate_obs] <- dnorm(sum(!rate_obs), 25, 0.01)
  hatched_init <- mcmc_data$hatched
  hatched_obs <- !is.na(hatched_init)
  hatched_init[!hatched_obs] <- rbinom(sum(!hatched_obs), trunc(rate_init[!hatched_obs] * batches_init[!hatched_obs]) + 5, 0.5)
  group_init <- rep(1, mcmc_const$N)
  group_init[hatched_init == 0] <- 0
  batches_init[batches_obs] <- NA
  rate_init[rate_obs] <- NA
  hatched_init[hatched_obs] <- NA
  mod <- nimbleModel(
    code = code,
    constants = mcmc_const,
    data = mcmc_data,
    inits = list(
      batches = batches_init,
      beta = cbind(
        c(rnorm(9, c(20, rep(0, 8)), 0.01)),
        rnorm(9, 0, 0.05)
      ),
      hatched = hatched_init,
      lambda = rnorm(2, c(1, 0), 0.01),
      kappa = rnorm(2, c(1, 0), 0.01),
      p = runif(mcmc_const$N, 0.45, 0.55),
      phi = cbind(
        rdirch(1, rep(1, mcmc_const$B)),
        rdirch(1, rep(1, mcmc_const$B))
      ),
      rate = rate_init,
      theta = rnorm(9, c(1.8, rep(0, 8)), 0.01)
    )
  )
  conf <- configureMCMC(mod, monitors = c("beta", "kappa", "lambda", "phi", "theta"))
  conf$removeSampler("beta[]")
  conf$removeSampler("kappa[]")
  conf$removeSampler("lambda[]")
  conf$removeSampler("theta[]")
  conf$addSampler(
    target = paste0("beta[", 1:9, ", 1]"),
    type = "RW_block",
    control = mcmc_const$control[[1]]
  )
  conf$addSampler(
    target = c("kappa[1]", "kappa[2]", "lambda[1]", "lambda[2]"),
    type = "RW_block",
    control = mcmc_const$control[[2]]
  )
  conf$addSampler(
    target = paste0("beta[", 1:9, ", 2]"),
    type = "RW_block",
    control = mcmc_const$control[[3]]
  )
  conf$addSampler(
    target = paste0("theta[", 1:9, "]"),
    type = "RW_block",
    control = mcmc_const$control[[4]]
  )
  Rmcmc <- buildMCMC(conf)
  Cmod <- compileNimble(mod)
  if (run) {
    Cmcmc <- compileNimble(Rmcmc)
    runMCMC(
      Cmcmc,
      niter = n_iter,
      nburnin = n_burnin,
      nchains = 1,
      thin = thin,
      setSeed = seed
    )
  } else {
    list(model = mod, Cmodel = Cmod, mcmc = Rmcmc)
  }
}
