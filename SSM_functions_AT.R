

# ---- 1) free_var_mod ----
free_var_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_party <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), by = "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_samples <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  
  # imputing NA solo su colonne numeriche (escludo Data)
  num_cols <- names(dt_samples)[-1]
  for (cc in num_cols) {
    dt_samples[[cc]][is.na(dt_samples[[cc]])] <- min(dt_samples[[cc]], na.rm = TRUE)
  }
  
  Y <- as.matrix(dt_party[,-1])
  X <- as.matrix(dt_samples[,-1])
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  mZ <- matrix(1, p, 1)
  mH <- array(0, c(p, p, n))
  mT <- matrix(1)
  mQ <- matrix(NA_real_, 1, 1)
  mR <- matrix(1)
  va1 <- matrix(20)
  mP1 <- matrix(56.25)
  mP1inf <- matrix(0)
  
  mod <- SSModel(Y ~ 0 + SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf, state_names = "mu"),
                 H = mH)
  
  updt <- function(pars, model, samples) {
    model$Q[1, 1, 1] <- exp(pars[1])
    var_eps <- exp(pars[-1])
    for (t in 1:nrow(samples)) {
      s <- pmax(samples[t, ], 1)  # guardia contro zeri
      diag(model$H[,,t]) <- var_eps / s
    }
    model
  }
  
  fit <- fitSSM(model = mod,
                inits = log(c(var_eta = 0.1, var_eps = rep(2000, p))),
                updatefn = updt,
                update_args = list(samples = X),
                hessian = TRUE, method = 'BFGS')
  
  kfs <- KFS(fit$model, filtering = "state", smoothing = "state")
  
  dt_party_trend <- tibble(Data = dt_party$Data,
                           Trend = kfs$alphahat[, 1])
  
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
    ggtitle(paste("Voting intentions for", par)) +
    theme_minimal(base_size = 13)
  
  # Intervalli per varianze istituto-specifiche
  param_estimates <- fit$optim.out$par
  H <- fit$optim.out$hessian
  H <- 0.5 * (H + t(H))
  stopifnot(all(is.finite(H)))
  
  requireNamespace("MASS", quietly = TRUE)
  vcov_matrix <- MASS::ginv(H)
  std_errors  <- sqrt(diag(vcov_matrix))
  z_value     <- qnorm(0.975)
  
  ci_lower <- exp(param_estimates - z_value * std_errors)[-1]
  ci_upper <- exp(param_estimates + z_value * std_errors)[-1]
  
  confidence_intervals <- data.frame(
    Istituto = colnames(Y),
    Varianza = exp(param_estimates[-1]),
    Lower95 = ci_lower,
    Upper95 = ci_upper
  )
  
  g2 <- confidence_intervals %>%
    mutate(Istituto = fct_reorder(Istituto, Varianza)) %>%
    ggplot(aes(x = Istituto, y = Varianza)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .25) +
    coord_flip() +
    labs(title = paste("Institute-specific variances –", par),
         x = "Institute",
         y = expression(hat(sigma)^2)) +
    theme_minimal(base_size = 13)
  
  return(list(
    grafico_trend = g1,
    grafico_varianze = g2,
    intervalli_varianza = confidence_intervals,
    trend_temporale = dt_party_trend
  ))
}


# ---- 2) bias_mod (senza confronto elezioni) ----
bias_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(Matrix)
  library(lubridate)
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_party <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_samples <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  
  # imputing NA solo su colonne numeriche (escludo Data)
  num_cols <- names(dt_samples)[-1]
  for (cc in num_cols) {
    dt_samples[[cc]][is.na(dt_samples[[cc]])] <- min(dt_samples[[cc]], na.rm = TRUE)
  }
  
  Y <- as.matrix(dt_party[, -1])
  X <- as.matrix(dt_samples[, -1])
  p <- ncol(Y)
  n <- nrow(Y)
  
  # Struttura con bias (delta) con vincolo somma-zero (ultima riga)
  mZ <- cbind(matrix(1, p, 1), diag(1, p))
  mZ <- mZ[, -ncol(mZ)]
  mZ[p, ] <- c(1, rep(-1, p - 1))
  mH <- array(0, c(p, p, n))
  
  mT <- diag(1, p)
  mQ <- matrix(NA_real_, 1, 1)
  mR <- rbind(matrix(1, 1, 1), matrix(0, p - 1, 1))
  
  va1 <- matrix(c(20, rep(0, p - 1)), p, 1)
  mP1 <- as.matrix(bdiag(matrix(56.25, 1, 1), diag(10, p - 1)))
  mP1inf <- matrix(0, p, p)
  
  s_names <- c("mu", paste0("delta_", 1:(p - 1)))
  
  mod <- SSModel(Y ~ 0 + SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf,
                                   state_names = s_names),
                 H = mH)
  
  updt <- function(pars, model, samples) {
    model$Q[1, 1, 1] <- exp(pars[1])
    var_eps <- exp(pars[2])
    for (t in 1:nrow(samples)) {
      s <- pmax(samples[t, ], 1)
      diag(model$H[,,t]) <- var_eps / s
    }
    model
  }
  
  set.seed(123)
  fit <- fitSSM(model = mod,
                inits = log(c(var_eta = 0.1, var_eps = 2000)),
                updatefn = updt,
                update_args = list(samples = X))
  
  kfs <- KFS(fit$model, filtering = "state", smoothing = "state")
  
  dt_party_trend <- tibble(Data = dt_party$Data,
                           Trend = kfs$alphahat[, 1])
  
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
    ggtitle(paste("Voting intentions for", par)) +
    theme_minimal(base_size = 13)
  
  # Stima bias per istituto (al tempo finale)
  delta_idx <- 2:p
  last_t <- nrow(kfs$alphahat)
  
  delta_hat <- kfs$alphahat[last_t, delta_idx]
  var_delta <- diag(kfs$V[delta_idx, delta_idx, last_t])
  se_delta <- sqrt(var_delta)
  
  z <- qnorm(0.975)
  ci_low  <- delta_hat - z * se_delta
  ci_high <- delta_hat + z * se_delta
  
  bias_tbl <- tibble(
    Istituto = colnames(Y)[1:(p - 1)],
    Bias = delta_hat,
    Lower95 = ci_low,
    Upper95 = ci_high
  ) %>% arrange(desc(Bias))
  
  bias_ultimo <- -sum(delta_hat)
  se_ultimo   <- sqrt(sum(var_delta))
  
  bias_tbl <- bind_rows(
    bias_tbl,
    tibble(
      Istituto = colnames(Y)[p],
      Bias = bias_ultimo,
      Lower95 = bias_ultimo - z * se_ultimo,
      Upper95 = bias_ultimo + z * se_ultimo
    )
  ) %>% arrange(desc(Bias))
  
  g2 <- bias_tbl %>%
    ggplot(aes(x = reorder(Istituto, Bias), y = Bias)) +
    geom_point() +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .2) +
    coord_flip() +
    labs(x = "Institute", y = "Bias (%)",
         title = "Estimated institute bias (house effect) with 95% CI") +
    theme_minimal(base_size = 13)
  
  return(list(
    grafico_trend = g1,
    grafico_bias_stimato = g2,
    bias_stimato = bias_tbl,
    trend_temporale = dt_party_trend
  ))
}


# ---- 3) bias_free_var_mod (senza confronto elezioni) ----
bias_free_var_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(Matrix)
  library(lubridate)
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_party <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), by = "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_samples <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  
  # imputing NA solo su colonne numeriche (escludo Data)
  num_cols <- names(dt_samples)[-1]
  for (cc in num_cols) {
    dt_samples[[cc]][is.na(dt_samples[[cc]])] <- min(dt_samples[[cc]], na.rm = TRUE)
  }
  
  Y <- as.matrix(dt_party[, -1])
  X <- as.matrix(dt_samples[, -1])
  p <- ncol(Y)
  n <- nrow(Y)
  
  mZ <- cbind(matrix(1, p, 1), diag(1, p))
  mZ <- mZ[, -ncol(mZ)]
  mZ[p, ] <- c(1, rep(-1, p - 1))
  mH <- array(0, c(p, p, n))
  mT <- diag(1, p)
  mQ <- matrix(NA_real_, 1, 1)
  mR <- rbind(matrix(1, 1, 1), matrix(0, p - 1, 1))
  
  va1 <- matrix(c(20, rep(0, p - 1)), p, 1)
  mP1 <- as.matrix(bdiag(matrix(56.25, 1, 1), diag(10, p - 1)))
  mP1inf <- matrix(0, p, p)
  
  s_names <- c("mu", paste0("delta_", 1:(p - 1)))
  mod <- SSModel(Y ~ 0 + SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf,
                                   state_names = s_names),
                 H = mH)
  
  updt <- function(pars, model, samples) {
    model$Q[1, 1, 1] <- exp(pars[1])
    var_eps <- exp(pars[-1])
    for (t in 1:nrow(samples)) {
      s <- pmax(samples[t, ], 1)
      diag(model$H[,,t]) <- var_eps / s
    }
    model
  }
  
  fit <- fitSSM(model       = mod,
                inits       = log(c(var_eta = 0.1, var_eps = rep(2000, p))),
                updatefn    = updt,
                update_args = list(samples = X),
                hessian     = TRUE,
                method      = "BFGS")
  
  kfs   <- KFS(fit$model, filtering = "state", smoothing = "state")
  trend <- tibble(Data = dt_party$Data, Trend = kfs$alphahat[, 1])
  
  # TREND
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(data = trend, aes(x = Data, y = Trend), linewidth = 1) +
    ggtitle(paste("Voting intentions for", par)) +
    theme_minimal(base_size = 13)
  
  # VARIANZE per istituto (+ CI)
  param_estimates <- fit$optim.out$par
  H <- fit$optim.out$hessian
  H <- 0.5 * (H + t(H))
  stopifnot(all(is.finite(H)))
  
  requireNamespace("MASS", quietly = TRUE)
  vcov_matrix <- MASS::ginv(H)
  
  std_errors <- sqrt(diag(vcov_matrix))
  z          <- qnorm(0.975)
  
  ci_lower <- exp(param_estimates - z * std_errors)[-1]
  ci_upper <- exp(param_estimates + z * std_errors)[-1]
  
  var_tbl <- data.frame(
    Istituto = colnames(Y),
    Varianza = exp(param_estimates[-1]),
    Lower95  = ci_lower,
    Upper95  = ci_upper
  )
  
  g2 <- var_tbl %>%
    mutate(Istituto = fct_reorder(Istituto, Varianza)) %>%
    ggplot(aes(x = Istituto, y = Varianza)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .25) +
    coord_flip() +
    labs(title = "Institute-specific variances – 95% CI",
         x = "Institute", y = expression(hat(sigma)^2)) +
    theme_minimal(base_size = 13)
  
  # BIAS stimati (ultimo tempo)
  delta_idx <- 2:p
  last_t    <- nrow(kfs$alphahat)
  delta_hat <- kfs$alphahat[last_t, delta_idx]
  var_delta <- diag(kfs$V[delta_idx, delta_idx, last_t])
  se_delta  <- sqrt(var_delta)
  
  ci_low  <- delta_hat - z * se_delta
  ci_high <- delta_hat + z * se_delta
  
  bias_tbl <- tibble(
    Istituto = colnames(Y)[1:(p - 1)],
    Bias     = delta_hat,
    Lower95  = ci_low,
    Upper95  = ci_high
  ) %>% arrange(desc(Bias))
  
  bias_ultimo <- -sum(delta_hat)
  se_ultimo   <- sqrt(sum(var_delta))
  
  bias_tbl <- bind_rows(
    bias_tbl,
    tibble(
      Istituto = colnames(Y)[p],
      Bias     = bias_ultimo,
      Lower95  = bias_ultimo - z * se_ultimo,
      Upper95  = bias_ultimo + z * se_ultimo
    )
  ) %>% arrange(desc(Bias))
  
  g3 <- bias_tbl %>%
    ggplot(aes(x = reorder(Istituto, Bias), y = Bias)) +
    geom_point() +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .2) +
    coord_flip() +
    labs(x = "Institute", y = "Bias (%)",
         title = "Estimated institute bias (house effect) with 95% CI") +
    theme_minimal()
  
  list(
    trend            = g1,
    varianze         = g2,
    bias_stimato     = g3,
    tabella_varianze = var_tbl,
    tabella_bias     = bias_tbl,
    trend_dati       = trend
  )
}


# ---- 4) bin_mod ----
bin_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_party <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_samples <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data)
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  
  # imputing NA solo su colonne numeriche (escludo Data)
  num_cols <- names(dt_samples)[-1]
  for (cc in num_cols) {
    dt_samples[[cc]][is.na(dt_samples[[cc]])] <- min(dt_samples[[cc]], na.rm = TRUE)
  }
  
  # Modello binomiale
  Y_pct    <- as.matrix(dt_party[, -1]) / 100
  U_trials <- as.matrix(round(dt_samples[, -1]))
  Y        <- round(Y_pct * U_trials)
  
  p <- ncol(Y)
  n <- nrow(dt_party)
  
  mZ  <- matrix(1, p, 1)
  mT  <- matrix(1)
  mQ  <- matrix(NA_real_, 1, 1)
  mR  <- matrix(1)
  va1 <- matrix(-1.4)
  mP1 <- matrix(0.3)
  mP1inf <- matrix(0)
  
  mod <- SSModel(Y ~ 0 + SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ,
                                   a1 = va1, P1 = mP1, P1inf = mP1inf,
                                   state_names = "mu"),
                 distribution = "binomial", u = U_trials)
  
  updt <- function(pars, model) {
    model$Q[1, 1, 1] <- exp(pars[1])
    model
  }
  
  fit <- fitSSM(model = mod,
                inits = log(0.001),
                updatefn = updt,
                method = 'BFGS')
  
  kfs <- KFS(fit$model, smoothing = "mean")
  
  dt_party_trend <- tibble(Data = dt_party$Data,
                           Trend = kfs$muhat[, 1] * 100)
  
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(data = dt_party_trend, aes(x = Data, y = Trend), linewidth = 1) +
    ggtitle(paste("Voting intentions for", par)) +
    theme_minimal(base_size = 13)
  
  # Importance sampling
  imp <- importanceSSM(fit$model, type = "state", antithetics = TRUE)
  w <- imp$weights / sum(imp$weights)
  
  mean_imp <- low_imp <- up_imp <- numeric(nrow(Y_pct))
  
  for (i in 1:nrow(Y_pct)) {
    ilogit_imp <- plogis(imp$samples[i, 1, ])
    mean_imp[i] <- sum(ilogit_imp * w)
    oo <- order(ilogit_imp)
    low_imp[i] <- ilogit_imp[oo][which.min(abs(cumsum(w[oo]) - 0.025))]
    up_imp[i]  <- ilogit_imp[oo][which.min(abs(cumsum(w[oo]) - 0.975))]
  }
  
  res <- tibble(
    Data = dt_party$Data,
    mean = mean_imp * 100,
    low  = low_imp * 100,
    high = up_imp * 100
  )
  
  g2 <- ggplot(res, aes(x = Data, y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = co, alpha = 0.25) +
    geom_line(size = 1, colour = co) +
    labs(title    = paste("Estimated share of", par, "over time"),
         subtitle = "Binomial state-space model – smoothed estimate via importance sampling",
         y        = "(%)",
         x        = NULL,
         caption  = "Band = 95% credible interval") +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot",
          plot.caption.position = "plot")
  
  list(
    grafico_trend_lineare = g1,
    grafico_importance_ci = g2,
    stime_credibilita = res
  )
}


# ---- 5) bin_bias_mod (fix 2: niente read.csv2; fix 3: niente confronto elezioni) ----
bin_bias_mod <- function(party, co) {
  library(tidyverse)
  library(Matrix)
  library(KFAS)
  library(lubridate)
  
  # Usa dt e nomi_istituti globali
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  # Wide Percentuale
  dt_party <- dt_filtro %>%
    filter(Partito == party) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto, values_from = Percentuale,
                values_fill = NA, values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- left_join(dates, dt_party, by = "Data")
  
  # Wide Campione
  dt_samples <- dt_filtro %>%
    filter(Partito == party) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto, values_from = Campione,
                values_fill = NA, values_fn = mean) %>%
    arrange(Data)
  
  dt_samples <- left_join(dates, dt_samples, by = "Data")
  
  # imputing NA solo su colonne numeriche (escludo Data)
  num_cols <- names(dt_samples)[-1]
  for (cc in num_cols) {
    dt_samples[[cc]][is.na(dt_samples[[cc]])] <- min(dt_samples[[cc]], na.rm = TRUE)
  }
  
  # Modello binomiale con delta (house effects) e vincolo somma-zero
  Y_pct    <- as.matrix(dt_party[, -1]) / 100
  U_trials <- as.matrix(round(dt_samples[, -1]))
  Y        <- round(Y_pct * U_trials)
  p <- ncol(Y)
  
  mZ <- rbind(cbind(1, diag(p - 1)), c(1, rep(-1, p - 1)))
  mT <- diag(p)
  mR <- rbind(1, matrix(0, p - 1, 1))
  mQ <- matrix(NA_real_, 1, 1)
  va1 <- matrix(c(-1.3, rep(0, p - 1)), p, 1)
  mP1 <- as.matrix(bdiag(1, diag(10, p - 1)))
  mP1inf <- matrix(0, p, p)
  s_names <- c("mu", paste0("delta_", 1:(p - 1)))
  
  mod <- SSModel(Y ~ -1 + SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ,
                                    a1 = va1, P1 = mP1, P1inf = mP1inf,
                                    state_names = s_names),
                 distribution = "binomial", u = U_trials)
  
  fit <- fitSSM(mod, inits = log(0.05), updatefn = \(pars, model) {
    model$Q[1, 1, 1] <- exp(pars[1]); model
  }, method = "BFGS")
  
  kfs <- KFS(fit$model, smoothing = c("state", "signal", "mean"))
  
  # Importance sampling per CI sullo share
  imp <- importanceSSM(fit$model, type = "state", antithetics = TRUE)
  w <- imp$weights / sum(imp$weights)
  
  extract_ci <- function(i) {
    ilogit <- plogis(imp$samples[i, 1, ])
    oo <- order(ilogit)
    c(
      mean = sum(ilogit * w),
      low  = ilogit[oo][which.min(abs(cumsum(w[oo]) - 0.025))],
      high = ilogit[oo][which.min(abs(cumsum(w[oo]) - 0.975))]
    )
  }
  
  ci_list <- t(sapply(1:nrow(Y_pct), extract_ci))
  trend_df <- data.frame(
    Data = dt_party$Data,
    mean = ci_list[, "mean"] * 100,
    low  = ci_list[, "low"]  * 100,
    high = ci_list[, "high"] * 100
  )
  
  # Grafico trend
  g_trend <- ggplot(trend_df, aes(x = Data, y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = co, alpha = 0.25) +
    geom_line(size = 1, colour = co) +
    labs(
      title = paste("Estimated share of", party, "over time"),
      subtitle = "Binomial state-space model – smoothed estimate via importance sampling",
      y = "(%)", x = NULL,
      caption = "Band = 95% credible interval"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot", plot.caption.position = "plot")
  
  # Bias (delta) su scala p.p. attorno alla media logit
  wquant <- function(x, w, probs = c(0.025, 0.975)) {
    oo <- order(x); x <- x[oo]; w <- w[oo] / sum(w); cw <- cumsum(w)
    sapply(probs, function(p) x[which.min(abs(cw - p))])
  }
  
  # prendi l'ultimo tempo disponibile
  idx_last <- nrow(trend_df)
  logit_mu_last <- kfs$alphahat[idx_last, 1]
  pi_ref <- plogis(logit_mu_last)
  
  delta_samples <- imp$samples[idx_last, 2:p, ]
  delta_samples <- rbind(delta_samples, -colSums(delta_samples))
  nomi_cols <- colnames(Y)
  
  bias_df <- map_dfr(1:p, function(j) {
    d <- delta_samples[j, ]
    m <- sum(d * w)
    ci <- wquant(d, w)
    
    bias_pp <- 100 * (plogis(qlogis(pi_ref) + m) - pi_ref)
    low_pp  <- 100 * (plogis(qlogis(pi_ref) + ci[1]) - pi_ref)
    high_pp <- 100 * (plogis(qlogis(pi_ref) + ci[2]) - pi_ref)
    
    tibble(
      Istituto   = nomi_cols[j],
      Bias_pp    = bias_pp,
      Low_95_pp  = low_pp,
      High_95_pp = high_pp
    )
  }) %>% arrange(desc(Bias_pp))
  
  g_bias <- bias_df %>%
    mutate(Istituto = fct_reorder(Istituto, Bias_pp)) %>%
    ggplot(aes(x = Istituto, y = Bias_pp)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Low_95_pp, ymax = High_95_pp),
                  width = .25, linewidth = .6) +
    coord_flip() +
    labs(
      title = paste("Estimated institute bias –", party),
      subtitle = "House effect (δ) with 95% credible interval",
      x = NULL, y = "Bias (p.p.)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot")
  
  return(list(
    grafico_trend     = g_trend,
    grafico_bias      = g_bias,
    trend_stimato     = trend_df,
    bias_stimati      = bias_df
  ))
}

