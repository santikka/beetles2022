###########
# Figures #
###########

theme_set(theme_bw(base_size = 20))
update_geom_defaults("line", list(size = 0.75))

# 24h survival
gg_surv_family_24h <- family_effects_plot(
    family = surv_family_24h,
    group = experiment,
    fun = "log",
    breaks = c(64, 32, 16, 8, 4, 2, 1, 1/2, 1/4, 1/8, 1/16, 1/32),
    labels = c("1st", "2nd", "2nd (addtional)"),
    y_label = "Relative odds of 24-hour survival"
)

gg_surv_curve_24h <- curve_plot(
    curve = surv_curve_24h,
    fun = "logit",
    breaks = c(0.999, 0.99, 0.95, 0.8, 0.6, 0.2),
    y_label = "24-hour survival probability"
)

gg_surv_benefit_24h <- benefit_plot(surv_benefit_24h)

gg_surv_pred_24h <- survival_predictive_plot(
    pred = surv_pred_24h,
    surv_var = survival_24h,
    data = gens,
    y_label = "24-hour survival probability"
)

gg_surv_family_24h_pairs <- surv_family_24h |>
    arrange(family) |>
    filter(experiment != 1) |>
    mutate(family_ix = family - 25,
           pa_family = surv_const_24h$mo[family_ix],
           pa_rank = surv_family_24h |>
               filter(experiment == 1) |>
               arrange(family) |>
               slice(pa_family) |>
               pull(rank) |>
               as.integer(),
           pa_mean = surv_family_24h$logmean[pa_family]) |>
    ggplot(aes(x = pa_mean, y = logmean, color = factor(pa_rank))) +
    geom_line() +
    geom_point() +
    geom_text_repel(aes(label = rank), key_glyph = "path", max.overlaps = 16) +
    scale_colour_discrete(name = "1st gen.\nfamily index",
                          guide = guide_legend(order = 1)) +
    scale_x_continuous(limits = c(-3.0, 3), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-3.0, 3), expand = c(0, 0)) +
    labs(x = "Log-mean of 1st gen. family effect", y = "Log-mean of 2nd gen. family effect") +
    new_scale_color() +
    geom_point(data = data.frame(extinct = c("a", "b", "c", "d"), pa_rank = 1),
               mapping = aes(x = 0, y = 0, colour = extinct, alpha = I(0))) +
    scale_color_manual(name = "Extinct\n1st gen.\nfamilies",
                       labels = surv_family_24h |>
                           arrange(family) |>
                           slice(mis_families) |>
                           arrange(rank) |>
                           pull(rank),
                       values = rep("gray", 4),
                       drop = FALSE,
                       guide = guide_legend(order = 2))

# larva-to-adult survival
gg_surv_family_la <- family_effects_plot(
    family = surv_family_la,
    group = generation,
    fun = "log",
    breaks = c(8, 4, 2, 1, 1/2, 1/4, 1/8),
    labels = c("1st", "2nd"),
    y_label = "Relative odds of larva-to-adult survival"
)

gg_surv_curve_la <- curve_plot(
    curve = surv_curve_la,
    fun = "logit",
    breaks = c(0.95, 0.9, 0.8, 0.6, 0.4),
    y_label = "Larva-to-adult survival probability"
)

gg_surv_benefit_la <- benefit_plot(surv_benefit_la)

gg_surv_pred_la <- survival_predictive_plot(
    pred = surv_pred_la,
    surv_var = larva_to_adult_survival,
    data = gens_la,
    y_label = "Larva-to-adult survival probability"
)

# Body mass (females)
gg_bm_family_ad_f <- family_effects_plot(
    family = bm_family_ad_f,
    group = generation,
    labels = c("1st", "2nd"),
    y_label = expression("Emergence body mass difference ( " * sqrt(mg) * ")")
)

#gg_bm_curve_ad_f <- curve_plot(
#    curve = bm_curve_ad_f,
#    y_label = "Emergence body mass (mg)"
#)

gg_bm_family_10d_f <- family_effects_plot(
    family = bm_family_10d_f,
    group = generation,
    labels = c("1st", "2nd"),
    y_label = expression("10-day body mass difference ( " * sqrt(mg) * ")")
)

#gg_bm_curve_10d_f <- curve_plot(
#    curve = bm_curve_10d_f,
#    y_label = "10-day body mass (mg)"
#)

# Body mass (males)
gg_bm_family_ad_m <- family_effects_plot(
    family = bm_family_ad_m,
    group = generation,
    labels = c("1st", "2nd"),
    y_label = expression("Emergence body mass difference ( " * sqrt(mg) * ")")
)

#gg_bm_curve_ad_m <- curve_plot(
#    curve = bm_curve_ad_m,
#    y_label = "Emergence body mass (mg)"
#)

gg_bm_family_10d_m <- family_effects_plot(
    family = bm_family_10d_m,
    group = generation,
    labels = c("1st", "2nd"),
    y_label = expression("10-day body mass difference ( " * sqrt(mg) * ")")
)

#gg_bm_curve_10d_m <- curve_plot(
#    curve = bm_curve_10d_m,
#    y_label = "10-day body mass (mg)"
#)

##########
# Tables #
##########

tab_surv_24h <- surv_table(surv_samples_24h, idx = 1, suffix = "(a)")
tab_surv_la <- surv_table(surv_samples_la, idx = 2, suffix = "(b)")
tab_bm_f <- bm_table(bm_samples_f, idx = 3)
tab_bm_m <- bm_table(bm_samples_m, idx = 3)
tab_repro_rate <- repro_table_rate(repro_samples, idx = 4)
tab_repro_hatching <- repro_table_hatching(repro_samples, idx = 4)

save(gg_surv_family_24h, gg_surv_family_24h_pairs, gg_surv_curve_24h, gg_surv_benefit_24h,
     gg_surv_family_la, gg_surv_curve_la, gg_surv_benefit_la,
     gg_bm_family_ad_f, gg_bm_family_10d_f,
     gg_bm_family_ad_m, gg_bm_family_10d_m,
     tab_surv_24h, tab_surv_la,
     tab_bm_f, tab_bm_m,
     tab_repro_rate, tab_repro_hatching,
     file = "data/figtab.RData")
