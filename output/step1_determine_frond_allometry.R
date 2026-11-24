################################################################################
# 
# Giant kelp frond length-biomass allometry
#
# Goal:
# Estimate whole-frond wet biomass (kg) as a function of total frond length (m)
# measured underwater.
#
# Model:
# biomass = a * length^b
#
################################################################################


rm(list = ls())

librarian::shelf(
  tidyverse,
  here,
  janitor,
  googlesheets4
)

gs4_auth()


################################################################################
# set paths

figdir <- here::here("frond_allometry", "figures")

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)


################################################################################
# read kelp data

read_kelp_data <- function(sheet_id, quarter) {
  read_sheet(
    sheet_id,
    sheet = "data_entry"
  ) %>%
    clean_names() %>%
    filter(
      !is.na(total_length_cm),
      !is.na(wet_biomass_kg),
      total_length_cm > 0,
      wet_biomass_kg > 0
    ) %>%
    mutate(
      total_length_m = total_length_cm / 100,
      quarter = quarter
    )
}

kelp_dat_orig <- bind_rows(
  read_kelp_data(
    "1cguAWQF4F-1imx25-EvRWCyR3heTAANonbGufxH9y6E",
    "Q3"
  ),
  read_kelp_data(
    "1DnBrfHtwV6snrxeDRvqKWlDVHeEFc1CM6o5wgzgs3Qo",
    "Q4"
  )
) %>%
  mutate(
    quarter = factor(quarter, levels = c("Q3", "Q4"))
  )


################################################################################
# inspect data

glimpse(kelp_dat_orig)

summary(kelp_dat_orig$total_length_m)
summary(kelp_dat_orig$wet_biomass_kg)

nrow(kelp_dat_orig)


################################################################################
# raw plot

plot(
  kelp_dat_orig$total_length_m,
  kelp_dat_orig$wet_biomass_kg,
  xlab = "Total frond length (m)",
  ylab = "Wet biomass (kg)",
  pch = 16
)


################################################################################
# fit power-function models directly on the original biomass scale

fit_power <- function(dat) {
  start_fit <- lm(
    log(wet_biomass_kg) ~ log(total_length_m),
    data = dat
  )

  nls(
    wet_biomass_kg ~ a * total_length_m^b,
    data = dat,
    start = list(
      a = exp(coef(start_fit)[1]),
      b = coef(start_fit)[2]
    )
  )
}

kelp_q3 <- filter(kelp_dat_orig, quarter == "Q3")
kelp_q4 <- filter(kelp_dat_orig, quarter == "Q4")

fit_q3 <- fit_power(kelp_q3)
fit_q4 <- fit_power(kelp_q4)
fit <- fit_power(kelp_dat_orig)

start_q3 <- lm(
  log(wet_biomass_kg) ~ log(total_length_m),
  data = kelp_q3
)
start_q4 <- lm(
  log(wet_biomass_kg) ~ log(total_length_m),
  data = kelp_q4
)

fit_by_quarter <- nls(
  wet_biomass_kg ~ ifelse(
    quarter == "Q3",
    a_q3 * total_length_m^b_q3,
    a_q4 * total_length_m^b_q4
  ),
  data = kelp_dat_orig,
  start = list(
    a_q3 = exp(coef(start_q3)[1]),
    b_q3 = coef(start_q3)[2],
    a_q4 = exp(coef(start_q4)[1]),
    b_q4 = coef(start_q4)[2]
  )
)

comparison <- anova(fit, fit_by_quarter)

summary(fit_q3)
summary(fit_q4)
summary(fit)
comparison


################################################################################
# extract estimated parameters

a_est <- coef(fit)["a"]
b_est <- coef(fit)["b"]

cat("\nEstimated Parameters:\n")
cat("a:", a_est, "\n")
cat("b:", b_est, "\n")

cat(
  "\nEstimated relationship:\n",
  "Wet biomass (kg) = ",
  round(a_est, 4),
  " * frond length (m)^",
  round(b_est, 3),
  "\n",
  sep = ""
)


################################################################################
# determine fit
#
# Calculate R-squared on the original biomass scale

predicted_values <- fitted(fit)

model_residuals <-
  kelp_dat_orig$wet_biomass_kg - predicted_values

# Residual Sum of Squares
rss <- sum(model_residuals^2)

# Total Sum of Squares
mean_biomass <-
  mean(kelp_dat_orig$wet_biomass_kg)

tss <- sum(
  (kelp_dat_orig$wet_biomass_kg - mean_biomass)^2
)

# R-squared
r_squared <- 1 - (rss / tss)

cat("\nR-squared:", r_squared, "\n")


################################################################################
# prediction error

rmse <- sqrt(
  mean(model_residuals^2)
)

mae <- mean(
  abs(model_residuals)
)

cat("RMSE:", rmse, "kg\n")
cat("MAE:", mae, "kg\n")


################################################################################
# add predicted biomass and residuals to dataset

kelp_dat <- kelp_dat_orig %>%
  mutate(
    predicted_biomass_kg = predict(fit),
    residual_kg =
      wet_biomass_kg - predicted_biomass_kg
  )


################################################################################
# residual plot

g_resid <- ggplot(
  kelp_dat,
  aes(
    x = predicted_biomass_kg,
    y = residual_kg
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Predicted wet biomass (kg)",
    y = "Residual (kg)"
  ) +
  theme_classic()

g_resid


################################################################################
# plot theme

base_theme <- theme(
  axis.text = element_text(
    size = 12,
    color = "black"
  ),
  
  axis.title = element_text(
    size = 12,
    color = "black"
  ),
  
  plot.tag = element_text(
    size = 9,
    color = "black"
  ),
  
  plot.title = element_text(
    size = 12,
    color = "black",
    face = "bold"
  ),
  
  # Gridlines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  axis.line = element_line(
    colour = "black"
  ),
  
  # Legend
  legend.key.size = unit(
    0.3,
    "cm"
  ),
  
  legend.spacing.y = unit(
    0.1,
    "cm"
  ),
  
  legend.text = element_text(
    size = 8,
    color = "black"
  ),
  
  legend.title = element_blank(),
  
  # Facets
  strip.text = element_text(
    size = 10,
    face = "bold",
    color = "black",
    hjust = 0
  ),
  
  strip.background = element_blank()
)


################################################################################
# create an allometry plot for each dataset

make_allometry_plot <- function(dat, model, plot_title) {
  a_plot <- unname(coef(model)["a"])
  b_plot <- unname(coef(model)["b"])
  predicted <- predict(model, newdata = dat)
  rss_plot <- sum((dat$wet_biomass_kg - predicted)^2)
  tss_plot <- sum((dat$wet_biomass_kg - mean(dat$wet_biomass_kg))^2)
  r_squared_plot <- 1 - rss_plot / tss_plot

  pred_dat <- tibble(
    total_length_m = seq(
      min(dat$total_length_m),
      max(dat$total_length_m),
      length.out = 500
    )
  ) %>%
    mutate(
      predicted_biomass_kg = a_plot * total_length_m^b_plot
    )

  equation_text <- paste0(
    "italic(y) == ",
    round(a_plot, 3),
    " %*% italic(x)^",
    round(b_plot, 2),
    " ~~ italic(R)^2 == ",
    round(r_squared_plot, 2)
  )

  ggplot(
    dat,
    aes(
      x = total_length_m,
      y = wet_biomass_kg
    )
  ) +
    geom_point(size = 2.5) +
    geom_line(
      data = pred_dat,
      aes(
        x = total_length_m,
        y = predicted_biomass_kg
      ),
      linewidth = 1
    ) +
    labs(
      title = plot_title,
      x = "Total frond length (m)",
      y = "Whole-frond wet biomass (kg)"
    ) +
    theme_bw() +
    base_theme +
    annotate(
      "text",
      x = min(dat$total_length_m),
      y = max(dat$wet_biomass_kg),
      label = equation_text,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 5,
      color = "black"
    ) +
    annotate(
      "text",
      x = min(dat$total_length_m),
      y = max(dat$wet_biomass_kg) * 0.90,
      label = paste("n =", nrow(dat)),
      hjust = 0,
      vjust = 1,
      size = 5,
      color = "black"
    )
}


g_q3 <- make_allometry_plot(kelp_q3, fit_q3, "Q3 2026")
g_q4 <- make_allometry_plot(kelp_q4, fit_q4, "Q4 2026")
g_all <- make_allometry_plot(kelp_dat_orig, fit, "All samples: Q3 and Q4 2026")

# Build one faceted figure so the axes are labeled once and scales are shared.
make_prediction_data <- function(dat, model) {
  a_plot <- unname(coef(model)["a"])
  b_plot <- unname(coef(model)["b"])

  tibble(
    quarter = unique(dat$quarter),
    total_length_m = seq(
      min(dat$total_length_m),
      max(dat$total_length_m),
      length.out = 500
    )
  ) %>%
    mutate(
      predicted_biomass_kg = a_plot * total_length_m^b_plot
    )
}

make_panel_annotation <- function(dat, model) {
  a_plot <- unname(coef(model)["a"])
  b_plot <- unname(coef(model)["b"])
  predicted <- predict(model, newdata = dat)
  r_squared_plot <- 1 - sum((dat$wet_biomass_kg - predicted)^2) /
    sum((dat$wet_biomass_kg - mean(dat$wet_biomass_kg))^2)

  tibble(
    quarter = unique(dat$quarter),
    x = min(dat$total_length_m),
    y = max(dat$wet_biomass_kg),
    equation_text = paste0(
      "italic(y) == ",
      round(a_plot, 3),
      " %*% italic(x)^",
      round(b_plot, 2),
      " ~~ italic(R)^2 == ",
      round(r_squared_plot, 2)
    ),
    sample_size_text = paste("n =", nrow(dat))
  )
}

quarter_predictions <- bind_rows(
  make_prediction_data(kelp_q3, fit_q3),
  make_prediction_data(kelp_q4, fit_q4)
)

quarter_annotations <- bind_rows(
  make_panel_annotation(kelp_q3, fit_q3),
  make_panel_annotation(kelp_q4, fit_q4)
) %>%
  mutate(
    x = min(kelp_dat_orig$total_length_m),
    y = max(kelp_dat_orig$wet_biomass_kg) * 1.20
)

g_quarters <- ggplot(
  kelp_dat_orig,
  aes(x = total_length_m, y = wet_biomass_kg)
) +
  geom_point(size = 2.5) +
  geom_line(
    data = quarter_predictions,
    aes(y = predicted_biomass_kg, group = quarter),
    linewidth = 1
  ) +
  facet_wrap(
    ~quarter,
    nrow = 1,
    labeller = as_labeller(c(Q3 = "Q3 2026", Q4 = "Q4 2026"))
  ) +
  labs(
    x = "Total frond length (m)",
    y = "Whole-frond wet biomass (kg)"
  ) +
  theme_bw() +
  base_theme +
  geom_text(
    data = quarter_annotations,
    aes(x = x, y = y, label = equation_text),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = quarter_annotations,
    aes(x = x, y = y * 0.90, label = sample_size_text),
    hjust = 0,
    vjust = 1,
    size = 4,
    inherit.aes = FALSE
  )

g_q3
g_q4
g_quarters
g_all


################################################################################
# save the allometry figures

ggsave(
  g_q3,
  file = file.path(figdir, "giant_kelp_frond_length_biomass_Q3_2026.png"),
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 300
)

ggsave(
  g_q4,
  file = file.path(figdir, "giant_kelp_frond_length_biomass_Q4_2026.png"),
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 300
)

ggsave(
  g_quarters,
  file = file.path(figdir, "giant_kelp_frond_length_biomass_Q3_Q4_2026.png"),
  width = 12,
  height = 6.5,
  units = "in",
  dpi = 300
)

ggsave(
  g_all,
  file = file.path(figdir, "giant_kelp_frond_length_biomass_all_quarters.png"),
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 300
)


################################################################################
# print final model summary

cat(
  "\n============================================================\n",
  "GIANT KELP FROND LENGTH-BIOMASS MODELS\n",
  "============================================================\n",
  "\nQ3 sample size: n = ", nrow(kelp_q3),
  "\nQ4 sample size: n = ", nrow(kelp_q4),
  "\nCombined sample size: n = ", nrow(kelp_dat_orig),
  "\n\nModel:",
  "\nBiomass (kg) = a * Length (m)^b",
  "\n\nQ3 estimated equation:",
  "\nBiomass (kg) = ",
  round(coef(fit_q3)["a"], 4),
  " * Length (m)^",
  round(coef(fit_q3)["b"], 3),
  "\n\nQ4 estimated equation:",
  "\nBiomass (kg) = ",
  round(coef(fit_q4)["a"], 4),
  " * Length (m)^",
  round(coef(fit_q4)["b"], 3),
  "\n\nCombined estimated equation:",
  "\nBiomass (kg) = ",
  round(a_est, 4),
  " * Length (m)^",
  round(b_est, 3),
  "\n\nR-squared = ",
  round(r_squared, 3),
  "\nRMSE = ",
  round(rmse, 3),
  " kg",
  "\nMAE = ",
  round(mae, 3),
  " kg",
  "\n\nQuarter comparison (shared vs. quarter-specific parameters):\n",
  "F = ",
  round(comparison$`F`[2], 3),
  ", p = ",
  format.pval(comparison$`Pr(>F)`[2], digits = 3),
  "\n============================================================\n",
  sep = ""
)
