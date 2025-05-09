# Title: Power Analysis and Simulations for PGS–Cognition Studies
# Author: Cameron Watson
# Purpose: Estimate minimum detectable effects, simulate power, and evaluate models across multiple datasets and scenarios

# =======================
# 1. SETUP
# =======================
library(pwr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(lavaan)
library(semPower)

# =======================
# 2. SAMPLE SIZES
# =======================
n_ukb <- 450000        # UK Biobank (PGS–Cognition)
n_ukb_neuro <- 40000   # UKB neuroimaging subset
n_clinical <- 1000     # GAP/EUGEI clinical cohort
n_gap <- 600           # GAP cohort only
n_meg_alspac <- 203    # ALSPAC RbG MEG sample
n_meg_hub <- 600       # Brain & Genomics Hub MEG sample

# =======================
# 3. MINIMUM DETECTABLE R²
# =======================
calculate_effect_size <- function(n, power_level) {
  pwr.f2.test(u = 1, v = n - 2, f2 = NULL, sig.level = 0.05, power = power_level)$f2
}

effect_sizes <- data.frame(
  Dataset = c("UK Biobank (PGS-Cognition)", "UKB Neuroimaging", "Clinical Cohorts", "GAP",
              "MEG ALSPAC", "MEG Hub"),
  Sample_Size = c(n_ukb, n_ukb_neuro, n_clinical, n_gap, n_meg_alspac, n_meg_hub),
  Effect_Size_80 = sapply(c(n_ukb, n_ukb_neuro, n_clinical, n_meg_alspac, n_meg_hub, n_gap), calculate_effect_size, power_level = 0.8),
  Effect_Size_90 = sapply(c(n_ukb, n_ukb_neuro, n_clinical, n_meg_alspac, n_meg_hub, n_gap), calculate_effect_size, power_level = 0.9)
)

effect_sizes$R2_80 <- effect_sizes$Effect_Size_80 / (1 + effect_sizes$Effect_Size_80)
effect_sizes$R2_90 <- effect_sizes$Effect_Size_90 / (1 + effect_sizes$Effect_Size_90)

print(effect_sizes)

# =======================
# 4. UKB R² vs Sample Size and Predictors
# =======================
calculate_minimum_effect_size <- function(sample_size, alpha = 0.001, power = 0.80, num_predictors) {
  pwr.f2.test(u = num_predictors, v = sample_size - num_predictors - 1, sig.level = alpha, power = power)$f2
}

calculate_minimum_R2 <- function(num_predictors, sample_size) {
  f2 <- calculate_minimum_effect_size(sample_size, num_predictors)
  return(f2 / (1 + f2))
}

sample_sizes <- seq(1000, 100000, by = 2000)

r2_1 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 1) * 100
r2_5 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 5) * 100
r2_10 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 10) * 100

ukb_plot_data <- data.frame(
  SampleSize = rep(sample_sizes, 3),
  MinR2 = c(r2_1, r2_5, r2_10),
  Predictors = factor(rep(c("1 Predictor", "5 Predictors", "10 Predictors"), each = length(sample_sizes)))
)

ukb_plot <- ggplot(ukb_plot_data, aes(x = MinR2, y = SampleSize, color = Predictors)) +
  geom_line(size = 1) +
  geom_point(shape = 16) +
  labs(title = "Power = 80%, Alpha = 0.001",
       x = expression("Minimum Detectable R"^2 ~ "(%)"),
       y = "Sample Size") +
  theme_minimal(base_size = 15) +
  scale_color_manual(values = c("#DD2F92", "#7A40E5", "#1E90FF")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_y_continuous(labels = scales::comma)

# =======================
# 5. GAP/EUGEI R² vs Sample Size and Predictors
# =======================
sample_sizes <- seq(100, 1500, by = 100)

r2_1 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 1) * 100
r2_5 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 5) * 100
r2_10 <- sapply(sample_sizes, calculate_minimum_R2, num_predictors = 10) * 100

gap_plot_data <- data.frame(
  SampleSize = rep(sample_sizes, 3),
  MinR2 = c(r2_1, r2_5, r2_10),
  Predictors = factor(rep(c("1 Predictor", "5 Predictors", "10 Predictors"), each = length(sample_sizes)))
)

gap_plot <- ggplot(gap_plot_data, aes(x = MinR2, y = SampleSize, color = Predictors)) +
  geom_line(size = 1) +
  geom_point(shape = 16) +
  labs(title = "Power = 80%, Alpha = 0.01",
       x = expression("Minimum Detectable R"^2 ~ "(%)"),
       y = "Sample Size") +
  theme_minimal(base_size = 15) +
  scale_color_manual(values = c("#DD2F92", "#7A40E5", "#1E90FF")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_y_continuous(labels = scales::comma)

# =======================
# 6. Combine UKB and GAP Plots
# =======================
combined_plot <- plot_grid(ukb_plot, gap_plot, labels = c("Study 1 - UK Biobank", "Study 2 - GAP/EUGEI"), ncol = 2)
print(combined_plot)

# =======================
# 7. Connectivity Regression: Minimum R² for 10 predictors at N = 40,000
# =======================
n <- 40000
u <- 10
v <- n - u - 1
result <- pwr.f2.test(u = u, v = v, f2 = NULL, sig.level = 0.001, power = 0.80)
f2 <- result$f2
R2 <- f2 / (1 + f2)
cat("Connectivity model minimum detectable R²:", round(R2, 6), "\n")

# =======================
# 8. Mediation Power (semPower)
# =======================
powerMed <- semPower.powerMediation(
  type = "a-priori",
  alpha = 0.001,
  power = 0.80,
  bYX = 0.02,
  bMX = 0.05,
  bYM = 0.05,
  nullEffect = "ind = 0",
  nIndicator = c(1, 1, 1),
  loadM = c(1, 1, 1),
  standardized = TRUE
)
summary(powerMed)

# =======================
# 9. Post-hoc Mediation Power
# =======================
powerMed_posthoc <- semPower.powerMediation(
  type = "post-hoc",
  alpha = 0.001,
  N = 40000,
  bYX = 0.01,
  bMX = 0.01,
  bYM = 0.01,
  nullEffect = "ind = 0",
  nIndicator = c(1, 1, 1),
  loadM = c(1, 1, 1),
  standardized = TRUE
)
summary(powerMed_posthoc)

# =======================
# 10. RMSEA Sensitivity (Model Fit)
# =======================
comp <- semPower.compromise(
  effect = NULL,
  effect.measure = "RMSEA",
  alpha = NULL,
  power = 0.80,
  N = 40000,
  df = 1
)
summary(comp)

# End of Script