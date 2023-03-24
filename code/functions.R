#############
# Functions #
#############

generate_consts <- function(x, s, cov_inits = list(),
                            sp_init_c = NULL, sp_init_t = NULL,
                            sigma_sp_init_c = NULL, sigma_sp_init_t = NULL) {
  mother <- match(x$Mother_ID, x$ID, incomparables = NA)
  if (missing(s)) {
    N <- nrow(x)
    N_c <- sum(!x$treatment_interaction)
  } else {
    s_ix <- which(x$sex == s)
    N <- length(s_ix)
    N_c <- sum(!x[s_ix, ]$treatment_interaction)
  }
  N_t <- N - N_c
  N_fam1 <- length(unique(x$family1)) - 1
  N_fam2 <- length(unique(x$family2)) - 1
  N_fam <- N_fam1 + N_fam2
  gen <- x$generation
  ex <- x$experiment
  age <- as.integer(x$treat_date - x$hatching_date)
  age4 <- 1L * (age == 4)
  age5 <- 1L * (age == 5)
  fam <- x$family1
  fam[is.na(fam)] <- N_fam1 + x$family2[!is.na(x$family2)]
  fam2_ord <- order(unique(na.exclude(x$family2)))
  mo <- x$family1[unique(na.exclude(mother))[fam2_ord]]
  fam_ix <- sapply(1:N_fam1, function(x) c(x, N_fam1 + which(mo == x)))
  fam_size <- sapply(fam_ix, length)
  fam_max_size <- max(fam_size)
  fam_seq <- seq_len(fam_max_size)
  fam_mat <- t(sapply(fam_ix, "[", i = fam_seq))
  fam_mat[is.na(fam_mat)] <- 0
  trt <- x$treatment
  trt_c <- trt - mean(trt, na.rm = TRUE)
  pa_trt <- x$parental_treatment
  ti <- x$treatment_interaction
  sp_basis_c <- bs(
    x = trt[!ti],
    knots = seq(0.0, 8.0, by = 0.5),
    degree = 3,
    Boundary.knots = c(0.0, 8.0)
  )
  sp_basis_t <- bs(
    x = trt[ti],
    knots = seq(0.0, 6.0, by = 0.5),
    degree = 3,
    Boundary.knots = c(0.0, 6.0)
  )
  S_t <- ncol(sp_basis_t)
  S_c <- ncol(sp_basis_c)
  fam_tg <- which(fam_size > 1)
  fam_g <- which(fam_size == 1)
  N_fam_tg <- length(fam_tg)
  N_fam_g <- length(fam_g)
  fam_arr <- sapply(
    fam,
    function(x) {
      which(fam_mat == x, arr.ind = TRUE)
    }
  )
  n_cov <- length(cov_inits)
  control <- vector(mode = "list", length = n_cov)
  for (i in seq_len(n_cov)) {
    if (is.null(cov_inits[[i]])) {
      control[[i]] <- list()
    } else {
      control[[i]] <- list(propCov = cov_inits[[i]])
    }
  }
  if (is.null(sp_init_c)) {
    sp_init_c <- rep(1, S_c)
  }
  if (is.null(sp_init_t)) {
    sp_init_t <- rep(1, S_t)
  }
  if (is.null(sigma_sp_init_c)) {
    sigma_sp_init_c <- 0
  }
  if (is.null(sigma_sp_init_t)) {
    sigma_sp_init_t <- 0
  }
  named_list(
    N, N_c, N_t, N_fam, N_fam1, N_fam2,
    fam, fam_mat, fam_size, fam_max_size,
    N_fam_tg, N_fam_g, fam_tg, fam_g,
    mo, gen, ex, age4, age5, trt, pa_trt, ti,
    S_t, S_c, sp_init_c, sp_init_t,
    sigma_sp_init_c, sigma_sp_init_t,
    sp_basis_c, sp_basis_t, control
  )
}

generate_repro_consts <- function(x, cov_inits = list()) {
  N <- nrow(x)
  B <- length(unique(x$batches)[-1])
  gen <- x$generation
  gen2 <- 1L * (x$generation == 2)
  age4 <- 1L * (x$age == 4)
  age5 <- 1L * (x$age == 5)
  trt <- x$treatment
  pa_trt <- x$parental_treatment
  ti <- x$treatment_interaction
  bm <- (x$weight_10d_mg - mean(x$weight_10d_mg)) / 100
  log_total_mean <- mean(log(x$total) / log(2), na.rm = TRUE)
  n_cov <- length(cov_inits)
  control <- vector(mode = "list", length = n_cov)
  for (i in seq_len(n_cov)) {
    if (is.null(cov_inits[[i]])) {
      control[[i]] <- list()
    } else {
      control[[i]] <- list(propCov = cov_inits[[i]])
    }
  }
  named_list(
    N, B, gen, gen2, age4, age5, trt, pa_trt, ti, bm, log_total_mean, control
  )
}

is_varying <- function(x) {
  if (any(is.na(x))) {
    return(FALSE)
  }
  v <- var(x)
  if (v > 0) {
    return(TRUE)
  }
  return(FALSE)
}

mcmc_summary <- function(x, filter, transforms, effects, digits = 3) {
  if (!missing(filter)) {
    pars <- colnames(x)
    drop_cols <- integer(0)
    for (i in seq_along(filter)) {
      drop_cols <- c(drop_cols, which(startsWith(pars, filter[i])))
    }
    drop_cols <- unique(drop_cols)
    x <- x[, -drop_cols]
  }
  pars <- colnames(x)
  npar <- ncol(x)
  pr_lims <- numeric(npar)
  names(pr_lims) <- colnames(x)
  if (!missing(transforms)) {
    for (i in seq_along(transforms)) {
      cols <- transforms[[i]]$cols
      cols <- cols[cols %in% pars]
      x[, cols] <- transforms[[i]]$fun(x[, cols])
      pr_lims[cols] <- transforms[[i]]$fun(pr_lims[cols])
    }
  }
  means <- apply(x, 2, mean) |> round(digits)
  sds <- apply(x, 2, sd) |> round(digits)
  medians <- apply(x, 2, median) |> round(digits)
  credibles <- apply(x, 2, quantile, c(0.025, 0.975)) |> round(digits)
  probs <- rep(NA_real_, npar)
  if (!missing(effects)) {
    probs <- numeric(npar)
    effects_ix <- which(pars %in% effects)
    for (i in seq_len(npar)) {
      if (i %in% effects_ix) {
        if (means[i] < pr_lims[i]) {
          probs[i] <- mean(x[, i] <= pr_lims[i])
        } else {
          probs[i] <- mean(x[, i] >= pr_lims[i])
        }
      } else {
        probs[i] <- NA_real_
      }
    }
    probs <- round(probs, digits)
  }
  stoch <- which(!is.na(sds) & (sds > 0))
  tibble(
    Parameter = pars[stoch],
    Mean = means[stoch],
    Median = medians[stoch],
    SD = sds[stoch],
    `2.5%` = credibles[1, stoch],
    `97.5%` = credibles[2, stoch],
    EP = probs[stoch]
  )
}

mcmc_clean <- function(x, filter) {
  if (!missing(filter)) {
    keep <- !startsWith(colnames(x), filter)
  } else {
    keep <- rep(TRUE, ncol(x))
  }
  varying <- apply(x, 2, is_varying)
  mcmc(x[, keep & varying])
}

mcmc_psrf <- function(x, filter, autoburnin = FALSE) {
  x |>
    lapply(mcmc_clean, filter = filter) |>
    mcmc.list() |>
    gelman.diag(autoburnin = autoburnin)
}

mcmc_ess <- function(x, filter) {
  if (is.list(x)) {
    x |>
      lapply(mcmc_clean, filter = filter) |>
      mcmc.list() |>
      coda::effectiveSize()
  } else {
    x |>
      mcmc_clean() |>
      coda::effectiveSize()
  }
}

read_spss2 <- function(path, suffix) {
  x <- haven::read_spss(paste0(path, suffix)) |>
    haven::zap_labels() |>
    haven::zap_formats()
  x[] <- lapply(x, \(y) {
    attr(y, "display_width") <- NULL
    y
  })
  x
}

named_list <- function(...) {
  y <- list(...)
  names(y) <- as.character(substitute(list(...)))[-1]
  y
}

row_na <- function(x) {
  cols <- ncol(x)
  t(apply(x, 1, function(y) {
    if (any(is.na(y))) {
      rep(NA, cols)
    } else {
      y
    }
  }))
}

get_matches <- function(pattern, x) {
  reg <- regexec(pattern, x)
  matches <- regmatches(x, reg)
  sapply(matches, "[", 2)
}

# Based on: https://r-nimble.org/posterior-predictive-sampling-and-other-post-mcmc-use-of-samples-in-nimble
pp_sampler_NF <- nimbleFunction(
  setup = function(model, mcmc) {
    data_nodes <- model$getNodeNames(dataOnly = TRUE)
    parent_nodes <- model$getParents(data_nodes, stochOnly = TRUE)
    dep_nodes <- model$getDependencies(parent_nodes, self = TRUE, includeData = FALSE)
    sampled_nodes <- mcmc$mvSamples$getVarNames()
    dep_expand <- dep_nodes |>
      model$expandNodeNames(returnScalarComponents = TRUE)
    data_expand <- data_nodes |>
      model$expandNodeNames(returnScalarComponents = TRUE)
    sampled_expand <- sampled_nodes |>
      model$expandNodeNames(returnScalarComponents = TRUE)
    sim_nodes <- setdiff(unique(c(dep_expand, data_expand)), sampled_expand)
    n <- data_expand |>
      length()
  },
  run = function(samples = double(2)) {
    n_samp <- dim(samples)[1]
    pp_samples <- matrix(nrow = n_samp, ncol = n)
    for (i in 1:n_samp) {
      values(model, sampled_expand) <<- samples[i, ]
      model$simulate(sim_nodes, includeData = TRUE)
      pp_samples[i, ] <- values(model, data_expand)
    }
    returnType(double(2))
    return(pp_samples)
  }
)

get_predictive_sample <- function(x, samples) {
  sampler <- pp_sampler_NF(x$model, x$mcmc)
  Csampler <- compileNimble(sampler, project = x$model)
  vars <- x$mcmc$mvSamples$getVarNames() |>
    x$model$expandNodeNames(returnScalarComponents = TRUE)
  out <- Csampler$run(samples[, vars])
  colnames(out) <- x$model$getNodeNames(dataOnly = TRUE) |>
    x$model$expandNodeNames(returnScalarComponents = TRUE)
  out[, str_order(colnames(out), numeric = TRUE)]
}

compute_stats <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    minimum = min(x, na.rm = TRUE),
    maximum = max(x, na.rm = TRUE),
    variance = var(x, na.rm = TRUE),
    skewness = skewness(x, na.rm = TRUE),
    kurtosis = kurtosis(x, na.rm = TRUE)
  )
}

bayes_p_values <- function(obs, pred) {
  obs_stats <- compute_stats(obs)
  apply(pred, 1, compute_stats) |>
    t() |>
    apply(1, function(x) x < obs_stats) |>
    t() |>
    colMeans()
}

family_effects <- function(samples, idx, fun = exp, slice_vec) {
  if (missing(idx)) {
    pattern <- "^f\\[([0-9]+)\\]$"
  } else {
    pattern <- paste0("^f\\[([0-9]+), ", idx, "\\]$")
  }
  out <- samples |>
    as.data.frame() |>
    select(matches(pattern, perl = TRUE)) |>
    fun() |>
    gather(key = "family", value = "effect") |>
    group_by(family) |>
    summarise(
      mean = mean(effect),
      low = quantile(effect, 0.025),
      high = quantile(effect, 0.975),
      logmean = mean(log(effect))
    ) |>
    mutate(
      family = gsub(pattern, "\\1", family, perl = TRUE) |> as.integer(),
      generation = as.factor(2 - 1 * (family <= 25)),
      experiment = factor(2 - 1 * (family <= 25) + 1 * (family >= 75)),
      extinct = family %in% mis_families
    )
  if (missing(slice_vec)) {
    out |>
      arrange(experiment, desc(mean)) |>
      mutate(rank = factor(1:n()))
  } else {
    out |>
      arrange(family) |>
      slice(as.integer(slice_vec[1:n()])) |>
      mutate(rank = factor(1:n()))
  }
}

family_effects_plot <- function(family, group, fun = c("identity", "log"),
                                breaks, labels, y_label) {
  fun <- match.arg(fun)
  if (identical(fun, "log")) {
    label_fun <- function(x) {
      out <- as.character(x)
      y <- sapply(x < 1, isTRUE)
      out[y] <- paste0("1/", 1 / x[y])
      out
    }
  } else {
    label_fun <- waiver()
    breaks <- waiver()
  }
  x_seq <- seq(1, max(family$family), by = 2)
  ggplot(family, aes(x = rank, y = mean, colour = {{ group }})) +
    scale_y_continuous(
      trans = fun,
      labels = label_fun,
      breaks = breaks
    ) +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    scale_x_discrete(breaks = x_seq) +
    scale_colour_discrete(name = "Generation", labels = labels) +
    labs(x = "Family index", y = y_label)
}

summarise_curve <- function(x, axis_vals = seq(0, 6, by = 0.1), prob = FALSE,
                            group, transfun = expit, centerfun = median) {
  out <- x |>
    transfun() |>
    as.data.frame() |>
    `colnames<-`(axis_vals) |>
    mutate(id = row_number()) |>
    pivot_longer(cols = !id, names_to = c("treatment")) |>
    mutate(treatment = as.numeric(treatment), id = id) |>
    group_by(treatment)
  if (prob) {
    out |> summarise(
      prob = mean(value > transfun(0)),
      group = group,
    )
  } else {
    out |> summarise(
      center = centerfun(value),
      group = group,
      low = quantile(value, 0.025),
      high = quantile(value, 0.975)
    )
  }
}

survival_curve <- function(samples, const, fun = expit, prob = FALSE) {
  trt_plot <- seq(0, 6, by = 0.1)
  sp_basis_c <- t(bs(trt_plot,
    knots = attr(const$sp_basis_c, "knots"), degree = 3,
    Boundary.knots = attr(const$sp_basis_c, "Boundary.knots")
  ))
  sp_basis_t <- t(bs(trt_plot,
    knots = attr(const$sp_basis_t, "knots"), degree = 3,
    Boundary.knots = attr(const$sp_basis_t, "Boundary.knots")
  ))
  surv_spline_c <- samples[, paste0("sp_c[", 1:const$S_c, "]")] %*% sp_basis_c
  surv_spline_t <- samples[, paste0("sp_t[", 1:const$S_t, "]")] %*% sp_basis_t
  surv_transgen <- samples[, "beta[1]"]
  if (prob) {
    surv_control <- surv_spline_c
    surv_parental_a <- -surv_spline_c + surv_spline_t + surv_transgen * 1.00
    surv_parental_b <- -surv_spline_c + surv_spline_t + surv_transgen * 8.00
  } else {
    surv_control <- surv_spline_c
    surv_parental_a <- surv_spline_t + surv_transgen * 1.00
    surv_parental_b <- surv_spline_t + surv_transgen * 8.00
  }
  bind_rows(
    summarise_curve(surv_control, group = "1", transfun = fun, prob = prob),
    summarise_curve(surv_parental_a, group = "2", transfun = fun, prob = prob),
    summarise_curve(surv_parental_b, group = "3", transfun = fun, prob = prob)
  )
}

curve_plot <- function(curve, fun = c("identity", "logit"),
                       breaks = "", y_label = "") {
  fun <- match.arg(fun)
  if (fun == "identity") {
    breaks <- waiver()
  }
  curve |> ggplot(aes(x = treatment, y = center, colour = group)) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 6),
      labels = function(x) format(x, nsmall = 1)
    ) +
    scale_colour_manual(
      name = "Transgenerational\ntreatment",
      values = c("black", "red", "blue"),
      labels = c(
        "Control",
        "Treated (1.0 \u00b5g/ml)",
        "Treated (8.0 \u00b5g/ml)"
      )
    ) +
    scale_y_continuous(
      trans = fun,
      breaks = breaks
    ) +
    geom_line() +
    labs(
      x = "Within-generational treatment (\u00b5g/ml)",
      y = y_label
    ) +
    geom_line(aes(x = treatment, y = low), lty = 2, alpha = 0.5) +
    geom_line(aes(x = treatment, y = high), lty = 2, alpha = 0.5)
}

benefit_plot <- function(benefit) {
  benefit |>
    filter(group != 1) |>
    ggplot(aes(x = treatment, y = prob, colour = group)) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 6),
      labels = function(x) format(x, nsmall = 1)
    ) +
    scale_colour_manual(
      name = "Transgenerational\ntreatment",
      values = c("red", "blue"),
      labels = c("Treated (1.0 \u00b5g/ml)", "Treated (8.0 \u00b5g/ml)")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_line() +
    labs(
      x = "Within-generational treatment (\u00b5g/ml)",
      y = "Probability of increased survival"
    )
}

survival_predictive <- function(pp_samples, data) {
  pp_samples |>
    t() |>
    as.data.frame() |>
    mutate(
      treatment = data$treatment,
      parental_treatment = data$parental_treatment,
      trt_combn = interaction(treatment, parental_treatment, drop = TRUE, sep = "/")
    ) |>
    group_by(trt_combn) |>
    summarise(across(starts_with("V"), ~ mean(.x))) |>
    ungroup() |>
    pivot_longer(cols = starts_with("V"), names_to = "iter", values_to = "surv")
}

survival_predictive_plot <- function(pred, surv_var, data, y_label = "") {
  pred |>
    group_by(trt_combn) |>
    summarise(
      mean = mean(surv),
      low = quantile(surv, 0.25),
      high = quantile(surv, 0.75)
    ) |>
    mutate(datamean = data |>
      mutate(trt_combn = interaction(treatment, parental_treatment, drop = TRUE, sep = "/")) |>
      group_by(trt_combn) |>
      summarise(mean = mean({{ surv_var }})) |>
      pull(mean)) |>
    ggplot(aes(x = trt_combn, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    geom_point(aes(x = trt_combn, y = datamean, colour = "red")) +
    guides(colour = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1, vjust = .5)) +
    labs(
      x = "Within-generational treatment / Transgenerational treatment (\u00b5g/ml)",
      y = y_label
    )
}

generational_differences <- function(samples, transfun = identity, comparefun = `-`, idx) {
  pattern <- character(2)
  if (missing(idx)) {
    pattern[1] <- "^f\\[(2[6-9]|[3-9][0-9])\\]$"
    pattern[2] <- "^f\\[([1-9]|1[0-9]|2[0-5])\\]$"
  } else {
    pattern[1] <- paste0("^f\\[(2[6-9]|[3-9][0-9]), ", idx, "\\]$")
    pattern[2] <- paste0("^f\\[([1-9]|1[0-9]|2[0-5]), ", idx, "\\]$")
  }
  comparefun(
    samples |>
      as.data.frame() |>
      select(matches(pattern[1])) |>
      transfun() |>
      rowMeans(),
    samples |>
      as.data.frame() |>
      select(matches(pattern[2])) |>
      transfun() |>
      rowMeans()
  ) |>
    as.data.frame() |>
    `colnames<-`("effect") |>
    summarise(
      mean = mean(effect),
      sign = sign(mean),
      EP = if (sign == 1) {
        mean(effect > transfun(0))
      } else {
        mean(effect < transfun(0))
      }
    )
}

body_mass_curve <- function(samples, idx, data) {
  trt_plot <- seq(0, 6, by = 0.1)
  bm_icpt <- samples[, paste0("beta[1, ", idx, "]")]
  bm_trt <- samples[, paste0("beta[2, ", idx, "]")] %*% t(trt_plot)
  bm_transgen <- samples[, paste0("beta[3, ", idx, "]")]
  bm_interact <- samples[, paste0("beta[4, ", idx, "]")]
  bm_ad <- 0
  if (idx == 2) {
    bm_ad <- samples[, "beta[7, 2]"] * mean(data$bm |> pull(var = 1), na.rm = TRUE)
  }
  bind_rows(
    summarise_curve(
      bm_icpt + bm_trt + bm_ad,
      group = "1",
      transfun = \(x) x^2
    ),
    summarise_curve(
      bm_icpt + bm_trt + bm_transgen * 1.00 + bm_interact + bm_ad,
      group = "2",
      transfun = \(x) x^2
    ),
    summarise_curve(
      bm_icpt + bm_trt + bm_transgen * 8.00 + bm_interact + bm_ad,
      group = "3",
      transfun = \(x) x^2
    )
  )
}

repro_curve <- function(samples, idx, transfun) {
  trt_plot <- seq(0, 6, by = 0.1)
  repro_icpt <- samples[, paste0("beta[1, ", idx, "]")]
  repro_trt <- samples[, paste0("beta[2, ", idx, "]")] %*% t(trt_plot)
  repro_transgen <- samples[, paste0("beta[3, ", idx, "]")]
  repro_interact <- samples[, paste0("beta[4, ", idx, "]")]
  bind_rows(
    summarise_curve(
      repro_icpt + repro_trt,
      group = "1",
      transfun = transfun
    ),
    summarise_curve(
      repro_icpt + repro_trt + repro_transgen * 1.00 + repro_interact,
      group = "2",
      transfun = transfun
    ),
    summarise_curve(
      repro_icpt + repro_trt + repro_transgen * 8.00 + repro_interact,
      group = "3",
      transfun = transfun
    )
  )
}

# To conserve memory, family effects in the linear predictors can be omitted
# from the MCMC samples. Afterwards, they can be computed from the raw effects
# using this function (raw effects contain the same information, but use a different indexing)
restore_families <- function(mcmc_const, mcmc_samples, nvars = 2) {
  f <- matrix(0, nrow = nrow(mcmc_samples), ncol = mcmc_const$N_fam * nvars)
  if (nvars == 2) {
    colnames(f) <- c(
      paste0("f[", 1:mcmc_const$N_fam, ", ", 1, "]"),
      paste0("f[", 1:mcmc_const$N_fam, ", ", 2, "]")
    )
    for (k in 1:nvars) {
      for (j in 1:mcmc_const$N_fam_tg) {
        for (l in 1:mcmc_const$fam_size[mcmc_const$fam_tg[j]]) {
          ix_raw <- paste0("f_raw[", l, ", ", mcmc_const$fam_tg[j], ", ", k, "]")
          ix <- paste0("f[", mcmc_const$fam_mat[mcmc_const$fam_tg[j], l], ", ", k, "]")
          f[, ix] <- mcmc_samples[, ix_raw]
        }
      }
      for (j in 1:mcmc_const$N_fam_g) {
        ix_raw <- paste0("f_raw[", 1, ", ", mcmc_const$fam_g[j], ", ", k, "]")
        ix <- paste0("f[", mcmc_const$fam_mat[mcmc_const$fam_g[j], 1], ", ", k, "]")
        f[, ix] <- mcmc_samples[, ix_raw]
      }
    }
  } else if (nvars == 1) {
    colnames(f) <- paste0("f[", 1:mcmc_const$N_fam, "]")
    for (j in 1:mcmc_const$N_fam_tg) {
      for (l in 1:mcmc_const$fam_size[mcmc_const$fam_tg[j]]) {
        ix_raw <- paste0("f_raw[", l, ", ", mcmc_const$fam_tg[j], "]")
        ix <- paste0("f[", mcmc_const$fam_mat[mcmc_const$fam_tg[j], l], "]")
        f[, ix] <- mcmc_samples[, ix_raw]
      }
    }
    for (j in 1:mcmc_const$N_fam_g) {
      ix_raw <- paste0("f_raw[", 1, ", ", mcmc_const$fam_g[j], "]")
      ix <- paste0("f[", mcmc_const$fam_mat[mcmc_const$fam_g[j], 1], "]")
      f[, ix] <- mcmc_samples[, ix_raw]
    }
  } else {
    stop("Invalid nvars value")
  }
  cbind(mcmc_samples, f)
}

surv_table <- function(x, idx, suffix) {
  effects_surv <- paste0("beta[", 1:3, "]")
  prec_family <- paste0("tau_f[", 1:2, "]")
  transform_surv <- list(
    list(cols = effects_surv, fun = exp),
    list(cols = prec_family, fun = function(x) 1 / sqrt(x))
  )
  mcmc_summary(x,
    filter = c("f", "sp", "log_sigma"),
    transforms = transform_surv, effects = effects_surv
  ) |>
    select(!Median) |>
    mutate(Parameter = recode(
      Parameter,
      "tau_f[1]" = "\\sigma_{f1}^{(idx)}",
      "tau_f[2]" = "\\sigma_{f2}^{(idx)}",
      "rho_f"    = "\\rho_{f}^{(idx)}",
      "beta[1]"  = "\\exp(\\beta_{1}^{(idx)})",
      "beta[2]"  = "\\exp(\\beta_{2}^{(idx)})",
      "beta[3]"  = "\\exp(\\beta_{3}^{(idx)})"
    )) |>
    mutate(Parameter = gsub("idx", idx, Parameter)) |>
    slice(str_order(Parameter, numeric = TRUE)) |>
    mutate(Description = c(
      paste0("Transgen. treatment ", suffix),
      paste0("Age is 4 days ", suffix),
      paste0("Age is 5 days ", suffix),
      paste0("Family RE correlation ", suffix),
      paste0("1st gen. family RE sd. ", suffix),
      paste0("2nd gen. family RE sd. ", suffix)
    )) |>
    relocate(Description)
}

bm_table <- function(x, idx) {
  effects_bm <- paste0("beta[", rep(2:7, 2), ", ", rep(1:2, each = 6), "]")
  spline_sigma <- c("log_sigma_sp_c", "log_sigma_sp_t")
  prec_bm <- paste0("tau_bm[", rep(1:2, 2), ", ", rep(1:2, each = 2), "]")
  prec_family <- paste0("tau_f[", rep(1:2, 2), ", ", rep(1:2, each = 2), "]")
  transform_bm <- list(
    list(cols = spline_sigma, fun = exp),
    list(cols = c(prec_bm, prec_family), fun = function(x) 1 / sqrt(x))
  )
  mcmc_summary(x,
    filter = c("f"),
    transforms = transform_bm, effects = effects_bm
  ) |>
    select(!Median) |>
    mutate(Parameter = recode(
      Parameter,
      "beta[1, 1]"   = "\\beta_{01}^{(idx)}",
      "beta[2, 1]"   = "\\beta_{11}^{(idx)}",
      "beta[3, 1]"   = "\\beta_{21}^{(idx)}",
      "beta[4, 1]"   = "\\beta_{31}^{(idx)}",
      "beta[5, 1]"   = "\\beta_{41}^{(idx)}",
      "beta[6, 1]"   = "\\beta_{51}^{(idx)}",
      "beta[1, 2]"   = "\\beta_{02}^{(idx)}",
      "beta[2, 2]"   = "\\beta_{12}^{(idx)}",
      "beta[3, 2]"   = "\\beta_{22}^{(idx)}",
      "beta[4, 2]"   = "\\beta_{32}^{(idx)}",
      "beta[5, 2]"   = "\\beta_{42}^{(idx)}",
      "beta[6, 2]"   = "\\beta_{52}^{(idx)}",
      "beta[7, 2]"   = "\\beta_{62}^{(idx)}",
      "rho_f[1]"     = "\\rho_{f1}^{(idx)}",
      "rho_f[2]"     = "\\rho_{f2}^{(idx)}",
      "tau_bm[1, 1]" = "\\sigma_{m11}^{(idx)}",
      "tau_bm[1, 2]" = "\\sigma_{m12}^{(idx)}",
      "tau_bm[2, 1]" = "\\sigma_{m21}^{(idx)}",
      "tau_bm[2, 2]" = "\\sigma_{m22}^{(idx)}",
      "tau_f[1, 1]"  = "\\sigma_{f11}^{(idx)}",
      "tau_f[1, 2]"  = "\\sigma_{f12}^{(idx)}",
      "tau_f[2, 1]"  = "\\sigma_{f21}^{(idx)}",
      "tau_f[2, 2]"  = "\\sigma_{f22}^{(idx)}"
    )) |>
    mutate(Parameter = gsub("idx", idx, Parameter)) |>
    slice(c(1:6, 14, 20:21, 16:17, 7:13, 15, 22:23, 18:19)) |>
    mutate(Description = c(
      "Intercept (a)",
      "Within-gen. treatment (a)",
      "Transgen. treatment (a)",
      "Treatment interaction (a)",
      "Age is 4 days (a)",
      "Age is 5 days (a)",
      "Family RE correlation (a)",
      "1st gen. family RE sd. (a)",
      "2nd gen. family RE sd. (a)",
      "1st gen. em. sd. (a)",
      "2nd gen. em. sd. (a)",
      "Intercept (b)",
      "Within-gen. treatment (b)",
      "Transgen. treatment (b)",
      "Treatment interaction (b)",
      "Age is 4 days (b)",
      "Age is 5 days (b)",
      "Emergence weight (b)",
      "Family RE correlation (b)",
      "1st gen. family RE sd. (b)",
      "2nd gen. family RE sd. (b)",
      "1st gen 10-day bm. sd. (b)",
      "2nd gen 10-day bm. sd. (b)"
    )) |>
    relocate(Description)
}

repro_table_rate <- function(x, idx) {
  ex_pars <- paste0("beta[", 2:9, ", 1]")
  gennorm_pars <- c("kappa[2]", "lambda[2]")
  effects_repro <- c(ex_pars, gennorm_pars)
  transform_repro <- list(
    list(cols = gennorm_pars, fun = exp)
  )
  mcmc_summary(x,
    filter = c("phi", "sigma", "theta"),
    transforms = transform_repro, effects = effects_repro
  ) |>
    select(!Median) |>
    filter(!grepl("^beta\\[., 2\\]", Parameter)) |>
    slice(str_order(Parameter, numeric = TRUE)) |>
    mutate(Parameter = recode(
      Parameter,
      "beta[1, 1]" = "\\beta_{01}^{(idx)}",
      "beta[2, 1]" = "\\beta_{11}^{(idx)}",
      "beta[3, 1]" = "\\beta_{21}^{(idx)}",
      "beta[4, 1]" = "\\beta_{31}^{(idx)}",
      "beta[5, 1]" = "\\beta_{41}^{(idx)}",
      "beta[6, 1]" = "\\beta_{51}^{(idx)}",
      "beta[7, 1]" = "\\beta_{61}^{(idx)}",
      "beta[8, 1]" = "\\beta_{71}^{(idx)}",
      "beta[9, 1]" = "\\beta_{81}^{(idx)}",
      "kappa[1]"   = "\\kappa_{0}^{(idx)}",
      "kappa[2]"   = "\\exp(\\kappa_{1}^{(idx)})",
      "lambda[1]"  = "\\lambda_{0}^{(idx)}",
      "lambda[2]"  = "\\exp(\\lambda_{1}^{(idx)})"
    )) |>
    mutate(Parameter = gsub("idx", idx, Parameter)) |>
    mutate(Description = c(
      "Intercept (loc.)",
      "Within-gen. treatment (loc.)",
      "Transgen. treatment (loc.)",
      "Treatment interaction (loc.)",
      "Age is 4 days (loc.)",
      "Age is 5 days (loc.)",
      "2nd generation (loc.)",
      "10-day body mass (loc.)",
      "Number of batches (loc.)",
      "Intercept (scale)",
      "2nd generation (scale)",
      "Intercept (shape)",
      "2nd generation (shape)"
    )) |>
    relocate(Description)
}

repro_table_hatching <- function(x, idx) {
  ex_pars <- paste0("beta[", 2:9, ", 2]")
  pseudo_pars <- paste0("theta[", 2:9, "]")
  effects_repro <- c(ex_pars, pseudo_pars)
  transform_repro <- list(
    list(cols = effects_repro, fun = exp)
  )
  mcmc_summary(x,
    filter = c("phi", "lambda", "kappa"),
    transforms = transform_repro, effects = effects_repro
  ) |>
    select(!Median) |>
    filter(!grepl("^beta\\[., 1\\]", Parameter)) |>
    slice(str_order(Parameter, numeric = TRUE)) |>
    mutate(Parameter = recode(
      Parameter,
      "beta[1, 2]" = "\\beta_{02}^{(4)}",
      "beta[2, 2]" = "\\exp(\\beta_{12}^{(idx)})",
      "beta[3, 2]" = "\\exp(\\beta_{22}^{(idx)})",
      "beta[4, 2]" = "\\exp(\\beta_{32}^{(idx)})",
      "beta[5, 2]" = "\\exp(\\beta_{42}^{(idx)})",
      "beta[6, 2]" = "\\exp(\\beta_{52}^{(idx)})",
      "beta[7, 2]" = "\\exp(\\beta_{62}^{(idx)})",
      "beta[8, 2]" = "\\exp(\\beta_{72}^{(idx)})",
      "beta[9, 2]" = "\\exp(\\beta_{82}^{(idx)})",
      "theta[1]"   = "\\theta_{0}^{(idx)}",
      "theta[2]"   = "\\exp(\\theta_{1}^{(idx)})",
      "theta[3]"   = "\\exp(\\theta_{2}^{(idx)})",
      "theta[4]"   = "\\exp(\\theta_{3}^{(idx)})",
      "theta[5]"   = "\\exp(\\theta_{4}^{(idx)})",
      "theta[6]"   = "\\exp(\\theta_{5}^{(idx)})",
      "theta[7]"   = "\\exp(\\theta_{6}^{(idx)})",
      "theta[8]"   = "\\exp(\\theta_{7}^{(idx)})",
      "theta[9]"   = "\\exp(\\theta_{8}^{(idx)})",
    )) |>
    mutate(Parameter = gsub("idx", idx, Parameter)) |>
    mutate(Description = c(
      "Intercept (mean)",
      "Within-gen. treatment (mean)",
      "Transgen. treatment (mean)",
      "Treatment interaction (mean)",
      "Age is 4 days (mean)",
      "Age is 5 days (mean)",
      "2nd generation (mean)",
      "10-day body mass (mean)",
      "log(Number of eggs) (mean)",
      "Intercept (sample size)",
      "Within-gen. treatment (sample size)",
      "Transgen. treatment (sample size)",
      "Treatment interaction (sample size)",
      "Age is 4 days (sample size)",
      "Age is 5 days (sample size)",
      "2nd generation (sample size)",
      "10-day body mass (sample size)",
      "log(Number of eggs) (sample size)"
    )) |>
    relocate(Description)
}
