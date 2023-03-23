####################
# Modeling results #
####################

# 24h survival
surv_pred_24h <- survival_predictive(surv_pp_samples_24h, gens)
surv_gens_24h <- generational_differences(surv_samples_24h, exp, `/`)
surv_family_24h <- family_effects(surv_samples_24h)
surv_curve_24h <- survival_curve(surv_samples_24h, surv_const_24h)
surv_benefit_24h <- survival_curve(surv_samples_24h, surv_const_24h, fun = exp, prob = TRUE)

# Larva-to-adult survival
surv_pred_la <- survival_predictive(surv_pp_samples_la, gens_la)
surv_gens_la <- generational_differences(surv_samples_la, exp, `/`)
surv_family_la <- family_effects(surv_samples_la, slice_vec = surv_family_24h$family)
surv_curve_la <- survival_curve(surv_samples_la, surv_const_la)
surv_benefit_la <- survival_curve(surv_samples_la, surv_const_la, fun = exp, prob = TRUE)

# Body mass (females)
bm_pp_ad_f <- bm_pp_samples_f |>
  as.data.frame() |>
  select(matches("bm\\[.+, 1\\]"))
bm_pvals_ad_f <- bayes_p_values(bm_data_f$bm$AD_weight_mg, bm_pp_ad_f)
bm_gens_ad_f <- generational_differences(bm_samples_f, idx = 1)
bm_curve_ad_f <- body_mass_curve(bm_samples_f, idx = 1, bm_data_f)
bm_family_ad_f <- family_effects(bm_samples_f, idx = 1, slice_vec = surv_family_24h$family)

bm_pp_10d_f <- bm_pp_samples_f |>
  as.data.frame() |>
  select(matches("bm\\[.+, 2\\]"))
bm_pvals_10d_f <- bayes_p_values(bm_data_f$bm$weight_10d_mg, bm_pp_10d_f)
bm_gens_10d_f <- generational_differences(bm_samples_f, idx = 2)
bm_curve_10d_f <- body_mass_curve(bm_samples_f, idx = 2, bm_data_f)
bm_family_10d_f <- family_effects(bm_samples_f, idx = 2, slice_vec = surv_family_24h$family)

# Body mass (males)
bm_pp_ad_m <- bm_pp_samples_m |>
  as.data.frame() |>
  select(matches("bm\\[.+, 1\\]"))
bm_pvals_ad_m <- bayes_p_values(bm_data_m$bm$AD_weight_mg, bm_pp_ad_m)
bm_gens_ad_m <- generational_differences(bm_samples_m, idx = 1)
bm_curve_ad_m <- body_mass_curve(bm_samples_m, idx = 1, bm_data_m)
bm_family_ad_m <- family_effects(bm_samples_m, idx = 1, slice_vec = surv_family_24h$family)

bm_pp_10d_m <- bm_pp_samples_m |>
  as.data.frame() |>
  select(matches("bm\\[.+, 2\\]"))
bm_pvals_10d_m <- bayes_p_values(bm_data_m$bm$weight_10d_mg, bm_pp_10d_m)
bm_gens_10d_m <- generational_differences(bm_samples_m, idx = 2)
bm_curve_10d_m <- body_mass_curve(bm_samples_m, idx = 2, bm_data_m)
bm_family_10d_m <- family_effects(bm_samples_m, idx = 2, slice_vec = surv_family_24h$family)

# Reproduction
# Get indices matching complete cases in the data, because values predicted
# for missing observations are not simulated for some reason
repro_pp_df <- repro_pp_samples |> as.data.frame()
repro_pp_batches <- repro_pp_df |> select(starts_with("batches"))
repro_pp_rate <- repro_pp_df |> select(starts_with("rate"))
repro_pp_hatched <- repro_pp_df |> select(starts_with("hatched"))
repro_pp_common_ix <- get_matches("hatched\\[(.+)\\]", colnames(repro_pp_hatched)) |> as.integer()
repro_pp_batches_ix <- paste0("batches[", repro_pp_common_ix, "]")
repro_pp_rate_ix <- paste0("rate[", repro_pp_common_ix, "]")
repro_pp_total <- trunc(repro_pp_rate[, repro_pp_rate_ix] * repro_pp_batches[, repro_pp_batches_ix])
repro_pp_low <- which(as.matrix(repro_pp_total) < as.matrix(repro_pp_batches[, repro_pp_batches_ix]), arr.ind = TRUE)
repro_pp_total[repro_pp_low] <- repro_pp_batches[repro_pp_low]
repro_pp_hatched_obs <- get_matches("hatched\\[(.+)\\]", colnames(repro_pp_hatched)) |> as.integer()

repro_pvals_total <- bayes_p_values(repros$total, repro_pp_total)
repro_pvals_hatched <- bayes_p_values(repros$hatched, repro_pp_hatched)

repro_pp_cor <- sapply(1:nrow(repro_samples), function(i) {
  cor(unlist(repro_pp_total[i, ]), unlist(repro_pp_hatched[i, ]))
})
repro_pvals_cor <- mean(repro_pp_cor <= cor(repros$total, repros$hatched, use = "complete.obs"))

repro_curve_rate <- repro_curve(repro_samples, idx = 1, transfun = identity)
repro_curve_hatching <- repro_curve(repro_samples, idx = 2, transfun = expit)

save(
  surv_pred_24h, surv_gens_24h, surv_family_24h, surv_curve_24h, surv_benefit_24h,
  surv_pred_la, surv_gens_la, surv_family_la, surv_curve_la, surv_benefit_la,
  bm_pvals_ad_f, bm_gens_ad_f, bm_curve_ad_f, bm_family_ad_f,
  bm_pvals_10d_f, bm_gens_10d_f, bm_curve_10d_f, bm_family_10d_f,
  bm_pvals_ad_m, bm_gens_ad_m, bm_family_ad_m,
  bm_pvals_10d_m, bm_gens_10d_m, bm_family_10d_m,
  repro_pvals_total, repro_pvals_hatched, repro_pvals_cor,
  repro_curve_rate, repro_curve_hatching,
  file = "data/results.RData"
)
