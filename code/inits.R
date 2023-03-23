surv_df_24h <- surv_samples_24h |> as.data.frame()
surv_cov_init_c_24h <- surv_df_24h |>
  select(matches("sp_c", perl = TRUE)) |>
  cov()
surv_cov_init_t_24h <- surv_df_24h |>
  select(matches("sp_t", perl = TRUE)) |>
  cov()
surv_sp_init_c_24h <- surv_df_24h |>
  select(matches("^sp_c", perl = TRUE)) |>
  colMeans()
surv_sp_init_t_24h <- surv_df_24h |>
  select(matches("^sp_t", perl = TRUE)) |>
  colMeans()
surv_sigma_sp_init_c_24h <- surv_df_24h[, "log_sigma_sp_c"] |> mean()
surv_sigma_sp_init_t_24h <- surv_df_24h[, "log_sigma_sp_t"] |> mean()

surv_df_la <- surv_samples_la |> as.data.frame()
surv_cov_init_c_la <- surv_samples_la |>
  as.data.frame() |>
  select(matches("sp_c", perl = TRUE)) |>
  cov()
surv_cov_init_t_la <- surv_samples_la |>
  as.data.frame() |>
  select(matches("sp_t", perl = TRUE)) |>
  cov()
surv_sp_init_c_la <- surv_df_la |>
  select(matches("^sp_c", perl = TRUE)) |>
  colMeans()
surv_sp_init_t_la <- surv_df_la |>
  select(matches("^sp_t", perl = TRUE)) |>
  colMeans()
surv_sigma_sp_init_c_la <- surv_df_la[, "log_sigma_sp_c"] |> mean()
surv_sigma_sp_init_t_la <- surv_df_la[, "log_sigma_sp_t"] |> mean()

bm_df_f <- bm_samples_f |> as.data.frame()
bm_cov_init_ad_f <- bm_samples_f |>
  as.data.frame() |>
  select(matches("beta\\[[1-6]", perl = TRUE)) |>
  select(matches(", 1\\]")) |>
  cov()
bm_cov_init_10d_f <- bm_samples_f |>
  as.data.frame() |>
  select(matches("beta\\[[2-6]", perl = TRUE)) |>
  select(matches(", 2\\]")) |>
  cov()

bm_df_m <- bm_samples_m |> as.data.frame()
bm_cov_init_ad_m <- bm_samples_m |>
  as.data.frame() |>
  select(matches("beta\\[[1-6]", perl = TRUE)) |>
  select(matches(", 1\\]")) |>
  cov()
bm_cov_init_10d_m <- bm_samples_m |>
  as.data.frame() |>
  select(matches("beta\\[[2-6]", perl = TRUE)) |>
  select(matches(", 2\\]")) |>
  cov()

repro_df <- repro_samples |> as.data.frame()
repro_cov_init_rate <- repro_samples |>
  as.data.frame() |>
  select(matches("beta", perl = TRUE)) |>
  select(matches(", 1\\]")) |>
  cov()
repro_cov_init_kl <- repro_samples |>
  as.data.frame() |>
  select(matches("kappa|lambda", perl = TRUE)) |>
  cov()
repro_cov_init_hatched <- repro_samples |>
  as.data.frame() |>
  select(matches("beta", perl = TRUE)) |>
  select(matches(", 2\\]")) |>
  cov()
repro_cov_init_t <- repro_samples |>
  as.data.frame() |>
  select(matches("theta", perl = TRUE)) |>
  cov()

save(surv_cov_init_c_24h, surv_cov_init_t_24h, surv_sp_init_c_24h, surv_sp_init_t_24h, surv_sigma_sp_init_c_24h, surv_sigma_sp_init_t_24h,
  surv_cov_init_c_la, surv_cov_init_t_la, surv_sp_init_c_la, surv_sp_init_t_la, surv_sigma_sp_init_c_la, surv_sigma_sp_init_t_la,
  bm_cov_init_ad_f, bm_cov_init_10d_f, bm_cov_init_ad_m, bm_cov_init_10d_m,
  repro_cov_init_rate, repro_cov_init_kl, repro_cov_init_hatched, repro_cov_init_t,
  file = par_inits
)
