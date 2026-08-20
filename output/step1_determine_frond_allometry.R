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

kelp_dat_orig <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1cguAWQF4F-1imx25-EvRWCyR3heTAANonbGufxH9y6E/edit?gid=0#gid=0",
  sheet = "data_entry"
) %>%
  clean_names() %>%
  filter(
    !(is.na(total_length_cm) | is.na(wet_biomass_kg))
  ) %>%
  filter(
    total_length_cm > 0,
    wet_biomass_kg > 0
  ) %>%
  mutate(
    total_length_m = total_length_cm / 100
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
# derive starting parameters
#
# Power model:
#
# biomass = a * length^b
#
# Taking logs:
#
# log(biomass) = log(a) + b * log(length)
#
# We use this log-log relationship ONLY to generate reasonable
# starting values for the nonlinear model below.

start_fit <- lm(
  log(wet_biomass_kg) ~ log(total_length_m),
  data = kelp_dat_orig
)

summary(start_fit)

a_init <- exp(coef(start_fit)[1])
b_init <- coef(start_fit)[2]

cat("\nInitial parameter estimates:\n")
cat("a:", a_init, "\n")
cat("b:", b_init, "\n")


################################################################################
# fit power-function model directly on original biomass scale
#
# biomass = a * length^b

fit <- nls(
  wet_biomass_kg ~ a * total_length_m^b,
  data = kelp_dat_orig,
  start = list(
    a = a_init,
    b = b_init
  )
)

summary(fit)


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
# create smooth prediction curve

pred_dat <- tibble(
  total_length_m = seq(
    min(
      kelp_dat_orig$total_length_m,
      na.rm = TRUE
    ),
    max(
      kelp_dat_orig$total_length_m,
      na.rm = TRUE
    ),
    length.out = 500
  )
) %>%
  mutate(
    predicted_biomass_kg =
      a_est * total_length_m^b_est
  )


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
# generate equation text
#
# Character string + parse = TRUE avoids the annotate() error
# caused by passing an expression object directly.

equation_text <- paste0(
  "italic(y) == ",
  round(a_est, 3),
  " %*% italic(x)^",
  round(b_est, 2),
  " ~~ italic(R)^2 == ",
  round(r_squared, 2)
)


################################################################################
# sample size

n <- nrow(kelp_dat_orig)

sample_size_text <- paste(
  "n =",
  n
)


################################################################################
# create main allometry plot

g <- ggplot(
  kelp_dat_orig,
  aes(
    x = total_length_m,
    y = wet_biomass_kg
  )
) +
  
  geom_point(
    size = 2.5
  ) +
  
  geom_line(
    data = pred_dat,
    aes(
      x = total_length_m,
      y = predicted_biomass_kg
    ),
    linewidth = 1
  ) +
  
  labs(
    x = "Total frond length (m)",
    y = "Whole-frond wet biomass (kg)"
  ) +
  
  theme_bw() +
  base_theme +
  
  annotate(
    "text",
    x = min(
      kelp_dat_orig$total_length_m,
      na.rm = TRUE
    ),
    y = max(
      kelp_dat_orig$wet_biomass_kg,
      na.rm = TRUE
    ),
    label = equation_text,
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 5,
    color = "black"
  ) +
  
  annotate(
    "text",
    x = min(
      kelp_dat_orig$total_length_m,
      na.rm = TRUE
    ),
    y = max(
      kelp_dat_orig$wet_biomass_kg,
      na.rm = TRUE
    ) * 0.90,
    label = sample_size_text,
    hjust = 0,
    vjust = 1,
    size = 5,
    color = "black"
  )

g


################################################################################
# save main figure

ggsave(
  g,
  file = file.path(
    figdir,
    "giant_kelp_frond_length_biomass.png"
  ),
  width = 6.5,
  height = 6.5,
  units = "in",
  dpi = 300
)


################################################################################
# print final model summary

cat(
  "\n============================================================\n",
  "GIANT KELP FROND LENGTH-BIOMASS MODEL\n",
  "============================================================\n",
  "\nSample size: n = ", n,
  "\n\nModel:",
  "\nBiomass (kg) = a * Length (m)^b",
  "\n\nEstimated equation:",
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
  "\n============================================================\n",
  sep = ""
)
