free_var_mod_influence <- function(par, co) {
  # co non utilizzato (lasciato per compatibilità firma)
  library(tidyverse)
  library(KFAS)
  
  # --- Prep dati ---
  dt_filtro <- subset(dt, Istituto %in% nomi_istituti)
  
  dt_party <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Campione) %>%
    pivot_wider(
      names_from = Istituto,
      values_from = Percentuale,
      values_fill = NA,
      values_fn = mean
    ) %>%
    arrange(Data)
  
  dates <- tibble(Data = seq(min(dt$Data), max(dt$Data), by = "day"))
  dt_party <- dates %>% left_join(dt_party, by = "Data")
  
  dt_samples <- dt_filtro %>%
    filter(Partito == par) %>%
    select(-Partito, -Percentuale) %>%
    pivot_wider(
      names_from = Istituto,
      values_from = Campione,
      values_fill = NA,
      values_fn = mean
    ) %>%
    arrange(Data)
  dt_samples <- dates %>% left_join(dt_samples, by = "Data")
  
  Y <- as.matrix(dt_party[, -1])
  X <- as.matrix(dt_samples[, -1])
  
  # Sanifica UNA VOLTA i campioni: NA/NaN/Inf/<=0 -> minimo positivo
  pos_min <- suppressWarnings(min(X[X > 0], na.rm = TRUE))
  if (!is.finite(pos_min)) pos_min <- 1
  X[!is.finite(X) | X <= 0] <- pos_min
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  # --- Modello SSM: local level comune ---
  mZ <- matrix(1, p, 1)
  mT <- matrix(1)
  mR <- matrix(1)
  
  # Q e H time-varying (terza dimensione = n)
  mQ <- array(NA_real_, c(1, 1, n))  # verrà riempita in updt()
  mH <- array(0, c(p, p, n))         # verrà riempita in updt()
  
  # Stato diffuso: più stabile/rapido
  va1 <- matrix(20)
  mP1 <- matrix(56.25)
  mP1inf <- matrix(0)
  
  mod <- SSModel(
    Y ~ 0 + SSMcustom(
      Z = mZ, T = mT, R = mR, Q = mQ,
      a1 = va1, P1 = mP1, P1inf = mP1inf,
      state_names = "mu"
    ),
    H = mH
  )
  
  # ------- helper: Sigma_bias = L L' (Cholesky parametrizzato) -------
  # Vettore parametri totale: c( log(var_eta), log_diag_camp (p), L_params (p*(p+1)/2) )
  build_Sigma_bias <- function(theta, p) {
    L <- matrix(0, p, p)
    k <- 1
    for (i in 1:p) {
      for (j in 1:i) {
        if (i == j) L[i, j] <- exp(theta[k]) else L[i, j] <- theta[k]
        k <- k + 1
      }
    }
    L %*% t(L)
  }
  
  # ------- updatefn: veloce e senza loop su t -------
  updt <- function(pars, model, samples) {
    p <- dim(model$H)[1]
    n <- dim(model$H)[3]
    
    # 1) varianza dello stato (una per giorno, gap=1)
    model$Q[1, 1, ] <- exp(pars[1])
    
    # 2) var campionarie per istituto (vettore)
    sigma2 <- exp(pars[2:(1 + p)])
    
    # 3) Sigma_bias costante nel tempo
    L_pars <- pars[(2 + p):length(pars)]
    Sigma_bias <- build_Sigma_bias(L_pars, p)
    
    # 4) H_t = Sigma_bias + diag(sigma2 / n_t)  (replica off-diagonali in blocco)
    model$H <- array(rep(Sigma_bias, n), dim = c(p, p, n))
    
    inv_nt <- 1 / samples  # n x p
    for (i in 1:p) {
      model$H[i, i, ] <- model$H[i, i, ] + sigma2[i] * inv_nt[, i]
    }
    model
  }
  
  # ------- inizializzazioni -------
  nL <- p * (p + 1) / 2
  init_log_var_eta   <- log(0.1)
  init_log_diag_camp <- rep(log(2000), p)
  init_L_diag_log    <- rep(log(1e-4), p)      # Sigma_bias ~ quasi 0
  init_L_offdiag     <- rep(0, nL - p)
  
  make_L_inits <- function(p, diag_log, offdiag_vals) {
    Lpars <- numeric(p * (p + 1) / 2)
    k <- 1; k_off <- 1
    for (i in 1:p) {
      for (j in 1:i) {
        if (i == j) {
          Lpars[k] <- diag_log[i]
        } else {
          Lpars[k] <- offdiag_vals[k_off]; k_off <- k_off + 1
        }
        k <- k + 1
      }
    }
    Lpars
  }
  L_inits <- make_L_inits(p, init_L_diag_log, init_L_offdiag)
  inits <- c(init_log_var_eta, init_log_diag_camp, L_inits)
  
  # ------- fit -------
  fit <- fitSSM(
    model = mod,
    inits = inits,
    updatefn = updt,
    update_args = list(samples = X),
    method = "BFGS",
    control = list(
      maxit  = 200,        # aumenta il numero massimo di iterazioni
      reltol = 1e-4,       # richiedi un miglioramento più piccolo prima di fermarti
      trace  = 1,           # stampa progresso
      REPORT = 5
    )
  )
  
  
  kfs <- KFS(fit$model, filtering = "state", smoothing = "state")
  
  # ----- output utili -----
  p_opt <- fit$optim.out$par
  L_pars_hat <- p_opt[(2 + p):length(p_opt)]
  Sigma_bias_hat <- build_Sigma_bias(L_pars_hat, p)
  colnames(Sigma_bias_hat) <- colnames(Y)
  rownames(Sigma_bias_hat) <- colnames(Y)
  Corr_bias_hat <- cov2cor(Sigma_bias_hat)
  
  list(
    Sigma_bias = Sigma_bias_hat,   # covarianze strutturali (influenze)
    Corr_bias  = Corr_bias_hat,    # correlazioni strutturali
    fit = fit
  )
}


dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data, format = "%d/%m/%Y")
nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:16])
partiti <- sort(setdiff(unique(dt$Partito), c("Az", "IV")))


library(corrplot)

# 1) Stima per tutti i partiti (con gestione errori)
results <- lapply(partiti, function(par_nome) {
  tryCatch({
    fit <- free_var_mod_influence(par_nome, NA)
    list(partito = par_nome, corr = fit$Corr_bias, ok = TRUE)
  }, error = function(e) {
    message(sprintf("➜ Errore per %s: %s", par_nome, conditionMessage(e)))
    list(partito = par_nome, corr = NULL, ok = FALSE, err = conditionMessage(e))
  })
})


library(corrplot)

# Tieni solo i fit andati a buon fine
res_ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$corr), results)

# Prendi massimo 5 partiti
res_ok5 <- head(res_ok, 5)

stopifnot(length(res_ok5) > 0)

# Ordine fisso degli istituti preso dalla prima matrice
ordine <- colnames(res_ok5[[1]]$corr)

# Reorder helper per uniformare l'ordine su TUTTE le matrici
reorder_cm <- function(cm, ord) {
  rn <- rownames(cm); cn <- colnames(cm)
  if (!identical(rn, cn)) stop("La matrice di correlazione non è quadrata/named consistent.")
  # Se mancano istituti, aggiungi NA per mantenere la griglia
  missing <- setdiff(ord, rn)
  if (length(missing) > 0) {
    # estendi con NA
    add <- matrix(NA_real_, nrow = length(missing), ncol = ncol(cm),
                  dimnames = list(missing, cn))
    cm <- rbind(cm, add)
    add <- matrix(NA_real_, nrow = nrow(cm), ncol = length(missing),
                  dimnames = list(rownames(cm), missing))
    cm <- cbind(cm, add)
  }
  cm[ord, ord, drop = FALSE]
}

# Corrplot helper (niente hclust, ordine originale!)
plot_corr_fixed <- function(cm, title = "") {
  corrplot(
    cm,
    method      = "color",
    type        = "upper",
    order       = "original",
    tl.col      = "black",
    tl.cex      = 0.7,
    addgrid.col = NA,
    diag        = FALSE,
    mar         = c(0, 0, 2, 0),
    title       = title
  )
}

# Crea PDF multipagina: 5 partiti + pagina con media
pdf("Corr_bias_top5_e_media.pdf", width = 9, height = 9)

# Pagine 1..5: partiti (in ordine)
stack <- list()
for (r in res_ok5) {
  cm <- reorder_cm(r$corr, ordine)
  plot_corr_fixed(cm, title = r$partito)
  stack[[length(stack) + 1]] <- cm
}

# Pagina finale: media (z di Fisher) sulle 5 matrici
to_z   <- function(r) 0.5 * log((1 + r) / (1 - r))
from_z <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)

# Stack come array coerente
p <- length(ordine)
arr <- array(NA_real_, dim = c(p, p, length(stack)),
             dimnames = list(ordine, ordine, seq_along(stack)))
for (i in seq_along(stack)) arr[, , i] <- stack[[i]]

z_mean    <- apply(arr, c(1, 2), function(v) mean(to_z(v), na.rm = TRUE))
corr_mean <- from_z(z_mean)

plot_corr_fixed(corr_mean, title = "Average Structural Correlations")

dev.off()
message("Creato: Corr_bias_top5_e_media.pdf")

saveRDS(results, file = "results_bias_corr.rds")

# Per ricaricarlo:
results <- readRDS("results_bias_corr.rds")
