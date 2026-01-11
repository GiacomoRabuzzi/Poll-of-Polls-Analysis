rm(list=ls())

dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data,
                   format = "%d/%m/%Y")
DT::datatable(dt)

sort(table(dt$Istituto), decreasing = TRUE)[1:18]

nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:18])

dt <- subset(dt, Istituto %in% nomi_istituti)


library(tidyverse)

par <- "Lega"

dt %>%
  filter(Partito == par) %>%
  select(-Partito, -Campione) %>%
  pivot_wider(names_from = Istituto,
              values_from = Percentuale,
              values_fill = NA,
              values_fn = mean) %>%
  arrange(Data) ->
  dt_party

# not all days are present in the table, let us enlarge the database
# so that all days are present, possibly with missing observations

dates <- tibble(Data = seq(dt$Data[1],
                           dt$Data[nrow(dt)],
                           "day"))

dates %>% left_join(dt_party, by = "Data") -> dt_party

dt %>%
  filter(Partito == par) %>%
  select(-Partito, -Percentuale) %>%
  pivot_wider(names_from = Istituto,
              values_from = Campione,
              values_fill = NA,
              values_fn = mean) %>%
  arrange(Data) ->
  dt_samples

dates %>% left_join(dt_samples, by = "Data") -> dt_samples

dt_samples[is.na(dt_samples)] <- min(dt_samples[, -1], na.rm = TRUE)


library(KFAS)

Y_pct <- as.matrix(dt_party[, -1]) / 100
U_trials <- as.matrix(round(dt_samples[, -1]))
Y <- round(Y_pct * U_trials)

p <- length(unique(dt$Istituto))
n <- nrow(dt_party)

# observation equation matrices
mZ <- cbind(matrix(1, p, 1),diag(NA_real_,p))

# transition equation matrices
mT <- diag(1,p+1)
mQ <- matrix(NA_real_, 1, 1)
mR <- rbind(matrix(1, 1, 1),matrix(0, p, 1))

# initial conditions
library(Matrix)
va1 <- matrix(c(-1.4,rep(1,p)), p+1, 1)
mP1 <- as.matrix(bdiag(matrix(0.6, 1, 1), diag(0, p)))
mP1inf <- matrix(0, p+1, p+1)

# let's make the SSModel
s_names <- "mu"
d_names <- paste0("delta_", 1:p)
s_names <- c(s_names, d_names)

mod <- SSModel(Y~0+SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ, a1 = va1, P1 = mP1, P1inf = mP1inf,
                             state_names = s_names), distribution = "binomial", u = U_trials)

updt <- function(pars, model) {
  model$Q[1, 1, 1] <- exp(pars[1])
  diag(model$Z[,-1,])[-p] <- pars[-1]
  diag(model$Z[,-1,])[p] <- -sum(pars[-1])
  model
}

set.seed(123)
init_Z <- rnorm(p-1,0,0.001)
fit <- fitSSM(model = mod,
              inits = c(log(c(0.001)),init_Z),
              updatefn = updt)

fit$optim.out$convergence

round(exp(fit$optim.out$par[1]), 4)
round(fit$optim.out$par[-1], 4)

kfs <- KFS(fit$model,
           filtering = "state",
           smoothing = "state")

# let us make a dataset with the trend
dt_party_trend <- tibble(Data = dt_party$Data,
                         Trend = exp(kfs$alphahat[,1])/(1+exp(kfs$alphahat[,1]))*100)

DT::datatable(dt_party_trend)

dt %>% filter(Partito == par) %>%
  ggplot(aes(x = Data, y = Percentuale)) +
  geom_point(color = "green", alpha = 0.2) +
  geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
  ggtitle(paste("Intenzioni di voto per", par))


# Calcolo bias
bias_logit <- c(fit$optim.out$par[-1], -sum(fit$optim.out$par[-1]))
names(bias_logit) <- colnames(Y)

bias_pct <- sapply(names(bias_logit), function(istituto) {
  mean(
    (exp(kfs$alphahat[,1] + bias_logit[istituto]) / (1 + exp(kfs$alphahat[,1] + bias_logit[istituto]))) -
      (exp(kfs$alphahat[,1]) / (1 + exp(kfs$alphahat[,1])))
  ) * 100
})


df_bias <- data.frame(
  Istituto = names(bias_pct),
  BiasPercentuale = bias_pct)

df_bias


