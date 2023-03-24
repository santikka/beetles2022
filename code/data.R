##############
# Data input #
##############

# Initial values for various parameters for MCMC
par_inits <- paste0(data_path, "data/inits.RData")
if (file.exists(par_inits)) {
  load(par_inits)
} else {
  surv_cov_init_c_24h <- NULL
  surv_cov_init_t_24h <- NULL
  surv_sp_init_c_24h <- NULL
  surv_sp_init_t_24h <- NULL
  surv_sigma_sp_init_c_24h <- NULL
  surv_sigma_sp_init_t_24h <- NULL
  surv_cov_init_c_la <- NULL
  surv_cov_init_t_la <- NULL
  surv_sp_init_c_la <- NULL
  surv_sp_init_t_la <- NULL
  surv_sigma_sp_init_c_la <- NULL
  surv_sigma_sp_init_t_la <- NULL
  bm_cov_init_ad_f <- NULL
  bm_cov_init_10d_f <- NULL
  bm_cov_init_ad_m <- NULL
  bm_cov_init_10d_m <- NULL
  repro_cov_init_rate <- NULL
  repro_cov_init_kl <- NULL
  repro_cov_init_hatched <- NULL
  repro_cov_init_t <- NULL
}

# Survival and body mass
gen1 <- read_spss2(data_path, "OD_1gen_2017_A.sav") |>
  rename(family1 = family) |>
  mutate(
    generation = 1L,
    experiment = 1L,
    parental_treatment = 0.0,
    Mother_ID = NA_character_,
    Father_ID = NA_character_
  )

gen2_24h <- read_spss2(data_path, "OD_2gen_2018_24h_survival_A.sav") |>
  select(!family) |>
  mutate(experiment = 3L) |>
  rename(treatment = trea) |>
  rename(parental_treatment = parental_trea)

gen2 <- read_spss2(data_path, "OD_2gen_2018_A.sav") |>
  select(!family) |>
  mutate(
    experiment = 2L,
    treatment = recode(as.numeric(treatment),
      "2" = 1.59, "3" = 2.6, "5" = 4.5
    ),
    parental_treatment = recode(as.numeric(parental_treatment),
      "2" = 1.59, "3" = 2.6, "4" = 6.0
    )
  ) |>
  bind_rows(gen2_24h) |>
  rename(larva_to_adult_survival = larvae_to_adult_survival) |>
  mutate(
    generation = 2L,
    family2 = as.integer(fct_inorder(as.factor(Mother_ID)))
  )

gens <- bind_rows(gen1, gen2) |>
  mutate(treatment_interaction = (treatment > 0 & parental_treatment > 0)) |>
  arrange(treatment_interaction)

# 1st generation families with no offspring in the 2nd generation
mis_families <- gens |>
  filter(ID %in% names(table(gens$Mother_ID))) |>
  select(family1) |>
  unique() |>
  pull(family1) |>
  setdiff(x = names(table(gens$family1))) |>
  as.integer()

gens_la <- gens |>
  filter(!is.na(larva_to_adult_survival))

gens_f <- gens |>
  filter(experiment < 3) |>
  arrange(sex, treatment_interaction)

gens_m <- gens |>
  filter(experiment < 3) |>
  arrange(desc(sex), treatment_interaction)

surv_const_24h <- generate_consts(
  x = gens,
  cov_inits = list(surv_cov_init_c_24h, surv_cov_init_t_24h),
  sp_init_c = surv_sp_init_c_24h,
  sp_init_t = surv_sp_init_t_24h,
  sigma_sp_init_c = surv_sigma_sp_init_c_24h,
  sigma_sp_init_t = surv_sigma_sp_init_t_24h
)

surv_const_la <- generate_consts(
  x = gens_la,
  cov_inits = list(surv_cov_init_c_24h, surv_cov_init_t_24h),
  sp_init_c = surv_sp_init_c_la,
  sp_init_t = surv_sp_init_t_la,
  sigma_sp_init_c = surv_sigma_sp_init_c_la,
  sigma_sp_init_t = surv_sigma_sp_init_t_la
)

bm_const_f <- generate_consts(
  x = gens_f,
  s = 1,
  cov_inits = list(bm_cov_init_ad_f, bm_cov_init_10d_f)
)

bm_const_m <- generate_consts(
  x = gens_m,
  s = 2,
  cov_inits = list(bm_cov_init_ad_m, bm_cov_init_10d_m)
)

surv_data_24h <- gens |>
  pull(survival_24h) |>
  list() |>
  setNames("surv")

surv_data_la <- gens_la |>
  pull(larva_to_adult_survival) |>
  list() |>
  setNames("surv")

bm_data_f <- gens_f |>
  filter(sex == 1) |>
  select(AD_weight_mg, weight_10d_mg) |>
  sqrt() |>
  list() |>
  setNames("bm")

bm_data_m <- gens_m |>
  filter(sex == 2) |>
  select(AD_weight_mg, weight_10d_mg) |>
  sqrt() |>
  list() |>
  setNames("bm")

# Reproduction
repro1 <- read_spss2(data_path, "OD_1gen_reproduction_A.sav") |>
  mutate(
    generation = 1L,
    parental_treatment = 0.0,
    treatment = recode(
      as.numeric(treatment),
      "1" = 0.3, "2" = 0.5, "3" = 0.8, "4" = 1.0, "5" = 1.59, "6" = 2.6,
      "7" = 3.5, "8" = 4.0, "9" = 4.5, "10" = 5.0, "11" = 5.5, "12" = 6.0,
      "13" = 6.5, "14" = 8.0
    )
  )

repro2 <- read_spss2(data_path, "OD_2gen_reproduction.sav") |>
  mutate(
    generation = 2L,
    treatment = recode(
      treatment,
      "2" = 1.59, "3" = 2.6, "4" = 4.5, "5" = 6
    )
  ) |>
  rename(parental_treatment = parental_trea) |>
  mutate(parental_treatment = recode(
    parental_treatment,
    "2" = 1.59, "3" = 2.6, "5" = 6.0
  ))

repros <- bind_rows(repro1, repro2) |>
  filter(nr_eggs_10d > 0 | is.na(nr_eggs_10d)) |>
  rename(
    total = nr_eggs_10d,
    batches = nr_batches_10d,
    hatched = nr_hatching_10d
  ) |>
  mutate(
    mother = match(Mother_ID, gens$ID, incomparables = NA),
    rate = total / batches,
    treat_date = gens[mother, ]$treat_date,
    hatching_date = gens[mother, ]$hatching_date,
    age = as.integer(treat_date - hatching_date),
    treatment_interaction = (treatment > 0 & parental_treatment > 0)
  ) |>
  left_join(gens |> select(ID, weight_10d_mg), by = c("Mother_ID" = "ID")) |>
  filter(!is.na(weight_10d_mg))

repro_const <- generate_repro_consts(
  x = repros,
  cov_inits = list(
    repro_cov_init_rate,
    repro_cov_init_kl,
    repro_cov_init_hatched,
    repro_cov_init_t
  )
)

repro_data <- repros |>
  select(batches, rate, hatched) |>
  as.list()
