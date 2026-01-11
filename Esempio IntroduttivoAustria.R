# =========================================================
# Austria — Plots from precomputed models (RDS)
# - 1) House effects by polling institute (points + CI) — BINOMIAL BIAS model
# - 2) MSE decomposition (Bias² + Variance) by institute — GAUSSIAN (free variance)
# - 3) Party trends over time — BINOMIAL model (no bias)
# =========================================================

library(tidyverse)
library(lubridate)
library(ggrepel)
library(stringr)

# --- Load precomputed object ---
modelli <- readRDS("modelli_precotti_AT.rds")
partiti <- names(modelli)

# Party colors
party_colors <- c(
  "ÖVP"   = "black",
  "SPÖ"   = "red3",
  "FPÖ"   = "cornflowerblue",
  "GRÜNE" = "forestgreen",
  "NEOS"  = "deeppink"
)

# =========================================================
# 1) HOUSE EFFECTS (bias in percentage points) – BINOMIAL BIAS model
# =========================================================
bias_all <- map_dfr(partiti, function(p) {
  modelli[[p]]$binbias$bias_stimati %>%
    transmute(
      Institute = Istituto,
      Party = p,
      Bias = Bias_pp,
      Low = Low_95_pp,
      High = High_95_pp
    )
})

ggplot(bias_all, aes(x = Bias, y = Party, colour = Party)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = .4) +
  geom_point() +
  geom_errorbarh(aes(xmin = Low, xmax = High), height = .15) +
  scale_colour_manual(values = party_colors, guide = "none") +
  facet_wrap(~ Institute, ncol = 3) +
  labs(
    title = "Estimated institute bias (house effect)",
    subtitle = "Points = estimate; bars = 95% credible interval; vertical line = zero bias",
    x = "Bias (percentage points)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot")

# =========================================================
# 2) MSE decomposition per institute = Bias² + E[σ² / n]
#    (σ² from Gaussian free-variance model; n from raw dataset)
# =========================================================

# Load sample sizes to compute E[1/n]
dt_n <- read.csv("DatiAustria.csv") %>%
  mutate(
    Istituto = str_squish(as.character(Istituto)),
    Partito  = as.character(Partito),
    Data     = as.Date(Data)
  )

inv_n <- dt_n %>%
  filter(!is.na(Campione), Campione > 0) %>%
  group_by(Partito, Istituto) %>%
  summarise(inv_n_mean = mean(1 / Campione, na.rm = TRUE), .groups = "drop")

mse_long <- map_dfr(partiti, function(p) {
  tabV <- modelli[[p]]$bfv$tabella_varianze %>%
    transmute(
      Institute = str_squish(as.character(Istituto)),
      sigma2   = as.numeric(Varianza)
    ) %>%
    distinct(Institute, .keep_all = TRUE)
  
  tabB <- modelli[[p]]$bfv$tabella_bias %>%
    transmute(
      Institute = str_squish(as.character(Istituto)),
      Bias_pp  = as.numeric(Bias)
    ) %>%
    distinct(Institute, .keep_all = TRUE)
  
  tabN <- inv_n %>%
    filter(Partito == p) %>%
    select(Institute = Istituto, inv_n_mean)
  
  tabB %>%
    inner_join(tabV, by = "Institute") %>%
    inner_join(tabN, by = "Institute") %>%
    mutate(
      Party = p,
      Bias2   = (Bias_pp)^2,                 # δ_j² (in percentage points²)
      VarComp = sigma2 * inv_n_mean          # E[σ_j² / n_{j,t}]
    ) %>%
    pivot_longer(c(Bias2, VarComp),
                 names_to = "Component", values_to = "MSE")
})

ggplot(mse_long, aes(x = MSE, y = Party, fill = Component)) +
  geom_col(width = .7) +
  facet_wrap(~ Institute, ncol = 3, scales = "free_x") +
  scale_fill_manual(
    values = c("Bias2" = "#F8766D", "VarComp" = "#00BFC4"),
    labels = c("Bias2" = "Bias²", "VarComp" = "Variance (σ² / n)")
  ) +
  labs(
    title = "Mean squared error decomposition by institute",
    subtitle = "MSE ≈ Bias² + E[σ² / n]",
    x = "MSE (percentage points²)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot")

# =========================================================
# 3) Trends over time — BINOMIAL model (no bias)
# =========================================================
trend_all <- map_dfr(partiti, function(p) {
  modelli[[p]]$bin$stime_credibilita %>%
    mutate(Party = p)
})

ggplot(trend_all, aes(x = Data, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = Party), alpha = 0.20) +
  geom_line(aes(color = Party), linewidth = 0.9) +
  scale_color_manual(values = party_colors, guide = "none") +
  scale_fill_manual(values = party_colors, guide = "none") +
  facet_wrap(~ Party, ncol = 2, scales = "free_y") +
  labs(
    title = "Voting intentions in Austria",
    subtitle = "Line = estimated mean; shaded area = 95% credible interval",
    x = NULL,
    y = "Share (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot")
