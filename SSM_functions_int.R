free_var_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(plotly)
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_party
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_samples
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)
  
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
      diag(model$H[,,t]) <- var_eps / samples[t, ]
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
  
  # GRAFICO TREND INTERATTIVO
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
    ggtitle(paste("Intenzioni di voto per", par)) +
    theme_minimal(base_size = 13)
  
  g1_plotly <- ggplotly(g1)
  
  param_estimates <- fit$optim.out$par
  vcov_matrix <- solve(fit$optim.out$hessian)
  std_errors <- sqrt(diag(vcov_matrix))
  z_value <- qnorm(0.975)
  
  ci_lower <- exp(param_estimates - z_value * std_errors)[-1]
  ci_upper <- exp(param_estimates + z_value * std_errors)[-1]
  
  confidence_intervals <- data.frame(
    Istituto = colnames(Y),
    Varianza = exp(param_estimates[-1]),
    Lower95 = ci_lower,
    Upper95 = ci_upper
  )
  
  # GRAFICO VARIANZE INTERATTIVO
  g2 <- confidence_intervals %>%
    mutate(Istituto = fct_reorder(Istituto, Varianza)) %>%
    ggplot(aes(x = Istituto, y = Varianza)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .25) +
    coord_flip() +
    labs(title = paste("Varianze specifiche degli istituti –", par),
         x = "Istituto",
         y = expression(hat(sigma)^2)) +
    theme_minimal(base_size = 13)
  
  g2_plotly <- ggplotly(g2)
  
  return(list(
    grafico_trend = g1_plotly,
    grafico_varianze = g2_plotly,
    intervalli_varianza = confidence_intervals,
    trend_temporale = dt_party_trend
  ))
}


bias_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(Matrix)
  library(readxl)
  library(lubridate)
  library(ggrepel)
  library(plotly)  # Aggiunto plotly
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_party
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_samples
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)
  
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
    var_eps <- exp(pars[2])
    for (t in 1:nrow(samples)) {
      diag(model$H[,,t]) <- var_eps / samples[t, ]
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
  
  # GRAFICO TREND (ggplot + plotly)
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
    ggtitle(paste("Intenzioni di voto per", par)) +
    theme_minimal(base_size = 13)
  
  g1_plotly <- ggplotly(g1)
  
  # Bias stimato
  delta_idx <- 2:p
  last_t <- nrow(kfs$alphahat)
  
  delta_hat <- kfs$alphahat[last_t, delta_idx]
  var_delta <- diag(kfs$V[delta_idx, delta_idx, last_t])
  se_delta <- sqrt(var_delta)
  
  z <- qnorm(0.975)
  ci_low <- delta_hat - z * se_delta
  ci_high <- delta_hat + z * se_delta
  
  bias_tbl <- tibble(
    Istituto = colnames(Y)[1:(p - 1)],
    Bias = delta_hat,
    Lower95 = ci_low,
    Upper95 = ci_high
  ) %>%
    arrange(desc(Bias))
  
  bias_ultimo <- -sum(delta_hat)
  var_ultimo <- sum(var_delta)
  se_ultimo <- sqrt(var_ultimo)
  
  bias_tbl <- bind_rows(
    bias_tbl,
    tibble(
      Istituto = colnames(Y)[p],
      Bias = bias_ultimo,
      Lower95 = bias_ultimo - z * se_ultimo,
      Upper95 = bias_ultimo + z * se_ultimo
    )
  ) %>%
    arrange(desc(Bias))
  
  # GRAFICO BIAS STIMATO (ggplot + plotly)
  g2 <- bias_tbl %>%
    ggplot(aes(x = reorder(Istituto, Bias), y = Bias)) +
    geom_point() +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .2) +
    coord_flip() +
    labs(x = "Istituto", y = "Bias (%)",
         title = "Bias stimato degli istituti (effetto 'house') con IC 95%") +
    theme_minimal(base_size = 13)
  
  g2_plotly <- ggplotly(g2)
  
  # Bias reale
  elez <- as.data.frame(read_xlsx("Elezioni2022.xlsx"))
  elez[,1] <- c(unique(dt$Partito)[1], unique(dt$Partito)[5],
                unique(dt$Partito)[4], unique(dt$Partito)[3], unique(dt$Partito)[2])
  
  data_riferimento <- as.Date("2022-09-25")
  partiti_esclusi <- c("Az", "IV")
  
  df_preelez <- dt %>%
    filter(Data <= data_riferimento, !(Partito %in% partiti_esclusi)) %>%
    group_by(Partito, Istituto) %>%
    filter(Data == max(Data)) %>%
    ungroup() %>%
    select(Partito, Istituto, Data, Percentuale)
  
  bias_reale <- df_preelez %>%
    filter(Partito == par) %>%
    mutate(Bias = Percentuale - elez[elez$Partito == par, 2]) %>%
    select(Istituto, Bias)
  
  df_estrea <- bias_reale %>%
    rename(Bias_reale = Bias) %>%
    left_join(bias_tbl %>% select(Istituto, Bias_modello = Bias),
              by = "Istituto")
  
  # GRAFICO CONFRONTO BIAS (ggplot + plotly)
  g3 <- ggplot(df_estrea, aes(x = Bias_reale, y = Bias_modello, label = Istituto)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = .4) +
    geom_smooth(method = "lm", se = FALSE, color = co, 
                linetype = "dashed", linewidth = 0.7) +
    geom_point(size = 2, colour = co) +
    geom_label_repel(size = 3, label.padding = unit(0.15, "lines"),
                     color = co, max.overlaps = Inf) +
    labs(
      x = "Bias reale (sondaggio - risultato elettorale)",
      y = "Bias stimato dal modello",
      title = paste("Confronto Bias Reale vs Stimato -", par)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot")
  
  g3_plotly <- ggplotly(g3)
  
  return(list(
    grafico_trend = g1_plotly,
    grafico_bias_stimato = g2_plotly,
    grafico_confronto_reale_stimato = g3_plotly,
    bias_stimato = bias_tbl,
    trend_temporale = dt_party_trend
  ))
}



bias_free_var_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(Matrix)
  library(readxl)
  library(lubridate)
  library(ggrepel)
  library(plotly)  # Per grafici interattivi
  
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
  dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)
  
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
      diag(model$H[,,t]) <- var_eps / samples[t, ]
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
  
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(data = trend, aes(x = Data, y = Trend), linewidth = 1) +
    ggtitle(paste("Intenzioni di voto per", par)) +
    theme_minimal(base_size = 13)
  g1_plotly <- ggplotly(g1)
  
  param_estimates <- fit$optim.out$par
  vcov_matrix     <- solve(fit$optim.out$hessian)
  std_errors      <- sqrt(diag(vcov_matrix))
  z               <- qnorm(0.975)
  
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
    labs(title = "Varianze specifiche degli istituti – IC95%",
         x = "Istituto", y = expression(hat(sigma)^2)) +
    theme_minimal(base_size = 13)
  g2_plotly <- ggplotly(g2)
  
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
    labs(x = "Istituto", y = "Bias (%)",
         title = "Bias stimato degli istituti (effetto 'house') con IC 95%") +
    theme_minimal()
  g3_plotly <- ggplotly(g3)
  
  elez <- as.data.frame(read_xlsx("Elezioni2022.xlsx"))
  elez[, 1] <- c(unique(dt$Partito)[1], unique(dt$Partito)[5],
                 unique(dt$Partito)[4], unique(dt$Partito)[3],
                 unique(dt$Partito)[2])
  
  data_rif <- as.Date("2022-09-25")
  partiti_esclusi <- c("Az", "IV")
  
  df_preelez <- dt %>%
    filter(Data <= data_rif, !(Partito %in% partiti_esclusi)) %>%
    group_by(Partito, Istituto) %>%
    filter(Data == max(Data)) %>%
    ungroup() %>%
    select(Partito, Istituto, Percentuale)
  
  bias_reale <- df_preelez %>%
    filter(Partito == par) %>%
    mutate(Bias = Percentuale - elez[elez$Partito == par, 2]) %>%
    select(Istituto, Bias)
  
  df_estrea <- bias_reale %>%
    rename(Bias_reale   = Bias) %>%
    left_join(bias_tbl %>% select(Istituto, Bias_modello = Bias),
              by = "Istituto")
  
  g4 <- ggplot(df_estrea, aes(x = Bias_reale, y = Bias_modello, label = Istituto)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = .4) +
    geom_smooth(method = "lm", se = FALSE, color = co, linetype = "dashed", linewidth = 0.7) +
    geom_point(size = 2, colour = co) +
    geom_label_repel(size = 3, label.padding = unit(0.15, "lines"), color = co, max.overlaps = Inf) +
    labs(x = "Bias reale (sondaggio - risultato elettorale)",
         y = "Bias stimato dal modello",
         title = paste("Confronto Bias Reale vs Stimato -", par)) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot")
  g4_plotly <- ggplotly(g4)
  
  var_eps_vec <- exp(fit$optim.out$par[-1])
  names(var_eps_vec) <- colnames(Y)
  
  dt_cmse <- dt[dt$Partito == par, ]
  dt_cmse$Campione[is.na(dt_cmse$Campione)] <- min(dt_cmse$Campione, na.rm = TRUE)
  avg_camp <- tapply(dt_cmse$Campione, dt_cmse$Istituto, mean)
  avg_camp <- avg_camp[colnames(Y)]
  
  var_media <- tibble(Istituto = colnames(Y), Var_media = var_eps_vec / avg_camp)
  
  mse_tbl <- bias_tbl %>%
    select(Istituto, Bias) %>%
    left_join(var_media, by = "Istituto") %>%
    mutate(
      MSE      = Bias^2 + Var_media,
      PropBias = (Bias^2) / MSE,
      PropVar  = Var_media / MSE
    )
  
  mse_tbl <- mse_tbl %>%
    mutate(
      Bias2    = Bias^2,
      Varianza = Var_media
    ) %>%
    select(Istituto, Bias2, Varianza) %>%
    pivot_longer(-Istituto,
                 names_to  = "Componente",
                 values_to = "Valore")
  
  g5 <- ggplot(mse_tbl, aes(x = reorder(Istituto, Valore), 
                            y = Valore, 
                            fill = Componente)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Composizione dell'MSE per istituto –", par),
      x     = "Istituto",
      y     = "MSE",
      fill  = ""
    )
  g5_plotly <- ggplotly(g5)
  
  list(
    trend             = g1_plotly,
    varianze          = g2_plotly,
    bias_stimato      = g3_plotly,
    bias_vs_reale     = g4_plotly,
    mse_decomposition = g5_plotly,
    tabella_varianze  = var_tbl,
    tabella_bias      = bias_tbl,
    trend_dati        = trend
  )
}



bin_mod <- function(par, co) {
  library(tidyverse)
  library(KFAS)
  library(plotly)  # Aggiunto per grafici interattivi
  
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto,
                values_from = Percentuale,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_party
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto,
                values_from = Campione,
                values_fill = NA,
                values_fn = mean) %>%
    arrange(Data) -> dt_samples
  
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)
  
  Y_pct <- as.matrix(dt_party[, -1]) / 100
  U_trials <- as.matrix(round(dt_samples[, -1]))
  Y <- round(Y_pct * U_trials)
  
  p <- length(unique(dt_filtro$Istituto))
  n <- nrow(dt_party)
  
  mZ <- matrix(1, p, 1)
  mT <- matrix(1)
  mQ <- matrix(NA_real_, 1, 1)
  mR <- matrix(1)
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
  
  # GRAFICO LINEARE INTERATTIVO
  g1 <- dt_filtro %>%
    filter(Partito == par) %>%
    ggplot(aes(x = Data, y = Percentuale)) +
    geom_point(color = co, alpha = 0.2) +
    geom_line(data = dt_party_trend, aes(x = Data, y = Trend), linewidth = 1) +
    ggtitle(paste("Intenzioni di voto per", par)) +
    theme_minimal(base_size = 13)
  g1_plotly <- ggplotly(g1)
  
  # IMPORTANCE SAMPLING
  imp <- importanceSSM(fit$model, type = "state", antithetics = TRUE)
  w <- imp$weights / sum(imp$weights)
  
  mean_imp <- low_imp <- up_imp <- numeric()
  
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
  
  # GRAFICO CREDIBILITÀ INTERATTIVO
  g2 <- ggplot(res, aes(x = Data, y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = co, alpha = 0.25) +
    geom_line(size = 1, colour = co) +
    labs(title    = paste("Quota stimata della", par, "nel tempo"),
         subtitle = "Modello a stati binomiale – stima smoothata con importance sampling",
         y        = "(%)",
         x        = NULL,
         caption  = "Banda = intervallo di credibilità al 95%") +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot",
          plot.caption.position = "plot")
  g2_plotly <- ggplotly(g2)
  
  list(
    grafico_trend_lineare = g1_plotly,
    grafico_importance_ci = g2_plotly,
    stime_credibilita = res
  )
}



bin_bias_mod <- function(party, co) {
  library(tidyverse)
  library(Matrix)
  library(KFAS)
  library(readxl)
  library(ggrepel)
  library(lubridate)
  library(plotly)  # <-- per grafici interattivi
  
  # --- Caricamento e preprocessing
  dt <- read.csv2("Sondaggi_2024-04-30.csv")
  dt$Percentuale <- as.numeric(dt$Percentuale)
  dt$Data <- as.Date(dt$Data, format = "%d/%m/%Y")
  nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:16])
  dt <- subset(dt, Istituto %in% nomi_istituti)
  
  # --- Dati wide per partito
  dt_party <- dt %>%
    filter(Partito == party) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(names_from = Istituto, values_from = Percentuale,
                values_fill = NA, values_fn = mean) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), "day"))
  dt_party <- left_join(dates, dt_party, by = "Data")
  
  dt_samples <- dt %>%
    filter(Partito == party) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(names_from = Istituto, values_from = Campione,
                values_fill = NA, values_fn = mean) %>%
    arrange(Data)
  dt_samples <- left_join(dates, dt_samples, by = "Data")
  dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)
  
  # --- Modello binomiale a stati
  Y_pct <- as.matrix(dt_party[, -1]) / 100
  U_trials <- as.matrix(round(dt_samples[, -1]))
  Y <- round(Y_pct * U_trials)
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
  
  # --- Importance sampling
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
  
  # --- GRAFICO TREND INTERATTIVO
  g_trend <- ggplot(trend_df, aes(x = Data, y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = co, alpha = 0.25) +
    geom_line(size = 1, colour = co) +
    labs(
      title = paste("Quota stimata di", party, "nel tempo"),
      subtitle = "Modello a stati binomiale – stima smoothata con importance sampling",
      y = "(%)", x = NULL,
      caption = "Banda = intervallo di credibilità al 95%"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot", plot.caption.position = "plot")
  g_trend_plotly <- ggplotly(g_trend)
  
  # --- Bias stimati
  wquant <- function(x, w, probs = c(0.025, 0.975)) {
    oo <- order(x); x <- x[oo]; w <- w[oo] / sum(w); cw <- cumsum(w)
    sapply(probs, function(p) x[which.min(abs(cw - p))])
  }
  
  data_elez <- as.Date("2022-09-25")
  idx_elez <- which(dt_party$Data == data_elez)
  
  logit_mu_elez <- kfs$alphahat[idx_elez, 1]
  pi_ref <- plogis(logit_mu_elez)
  
  delta_samples <- imp$samples[idx_elez, 2:p, ]
  delta_samples <- rbind(delta_samples, -colSums(delta_samples))
  nomi_istituti <- colnames(Y)
  
  bias_df <- map_dfr(1:p, function(j) {
    d <- delta_samples[j, ]
    m <- sum(d * w)
    ci <- wquant(d, w)
    
    bias_pp <- 100 * (plogis(qlogis(pi_ref) + m) - pi_ref)
    low_pp  <- 100 * (plogis(qlogis(pi_ref) + ci[1]) - pi_ref)
    high_pp <- 100 * (plogis(qlogis(pi_ref) + ci[2]) - pi_ref)
    
    tibble(
      Istituto   = nomi_istituti[j],
      Bias_pp    = bias_pp,
      Low_95_pp  = low_pp,
      High_95_pp = high_pp
    )
  }) %>% arrange(desc(Bias_pp))
  
  # --- GRAFICO BIAS INTERATTIVO
  g_bias <- bias_df %>%
    mutate(Istituto = fct_reorder(Istituto, Bias_pp)) %>%
    ggplot(aes(x = Istituto, y = Bias_pp)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Low_95_pp, ymax = High_95_pp),
                  width = .25, linewidth = .6) +
    coord_flip() +
    labs(
      title = paste("Bias stimato degli istituti –", party),
      subtitle = "Effetto-casa (δ) con intervallo di credibilità al 95 %",
      x = NULL, y = "Bias (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot")
  g_bias_plotly <- ggplotly(g_bias)
  
  # --- Confronto bias reale vs stimato
  elez <- as.data.frame(read_xlsx("Elezioni2022.xlsx"))
  elez[, 1] <- c(unique(dt$Partito)[1], unique(dt$Partito)[5],
                 unique(dt$Partito)[4], unique(dt$Partito)[3],
                 unique(dt$Partito)[2])
  
  data_rif <- as.Date("2022-09-25")
  partiti_esclusi <- c("Az", "IV")
  
  df_preelez <- dt %>%
    filter(Data <= data_rif, !(Partito %in% partiti_esclusi)) %>%
    group_by(Partito, Istituto) %>%
    filter(Data == max(Data)) %>%
    ungroup() %>%
    select(Partito, Istituto, Data, Percentuale)
  
  bias_reale <- df_preelez %>%
    filter(Partito == party) %>%
    mutate(Bias = Percentuale - elez[elez$Partito == party, 2]) %>%
    select(Istituto, Bias)
  
  df_confronto <- bias_reale %>%
    rename(Bias_reale = Bias) %>%
    left_join(bias_df %>% select(Istituto, Bias_modello = Bias_pp),
              by = "Istituto")
  
  # --- GRAFICO CONFRONTO INTERATTIVO
  g_confronto <- ggplot(df_confronto, aes(Bias_reale, Bias_modello, label = Istituto)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = .4) +
    geom_smooth(method = "lm", se = FALSE, colour = co,
                linetype = "dashed", linewidth = .7) +
    geom_point(size = 2, colour = co) +
    geom_label_repel(size = 3, label.padding = unit(0.15, "lines"),
                     colour = co, max.overlaps = Inf) +
    labs(
      title = paste("Confronto bias reale vs stimato –", party),
      subtitle = "Ultimo sondaggio prima del voto (25 set 2022) vs house-effect",
      x = "Bias reale (p.p.)", y = "Bias modello (p.p.)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title.position = "plot")
  g_confronto_plotly <- ggplotly(g_confronto)
  
  # --- OUTPUT finale con grafici interattivi
  return(list(
    grafico_trend     = g_trend_plotly,
    grafico_bias      = g_bias_plotly,
    grafico_confronto = g_confronto_plotly,
    trend_stimato     = trend_df,
    bias_stimati      = bias_df
  ))
}
