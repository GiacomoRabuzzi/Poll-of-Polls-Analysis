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

p <- length(unique(dt$Istituto))
n <- nrow(dt_party)

# data have to turned into matrix
Y <- as.matrix(dt_party[,-1])
X <- as.matrix(dt_samples[, -1])

# observation equation matrices
mZ <- cbind(matrix(1, p, 1),diag(NA_real_,p))
mH <- array(0, c(p, p, n))

# transition equation matrices
mT <- diag(1,p+1)
mQ <- matrix(NA_real_, 1, 1)
mR <- rbind(matrix(1, 1, 1),matrix(0, p, 1))

# initial conditions
library(Matrix)
va1 <- matrix(c(20,rep(1,p)), p+1, 1)
mP1 <- as.matrix(bdiag(matrix(56.25, 1, 1), diag(0, p)))
mP1inf <- matrix(0, p+1, p+1)

# let's make the SSModel
s_names <- "mu"
d_names <- paste0("delta_", 1:p)
s_names <- c(s_names, d_names)
mod <- SSModel(Y~0+SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf,
                             state_names = s_names),
               H = mH)

# let's make the update function
updt <- function(pars, model, samples) {
  model$Q[1, 1, 1] <- exp(pars[1])
  var_eps <- exp(pars[2])
  for (t in 1:nrow(samples)) {
    diag(model$H[,,t]) <- var_eps/samples[t, ]
  }
  diag(model$Z[,-1,])[-p] <- pars[-c(1,2)]
  diag(model$Z[,-1,])[p] <- -sum(pars[-c(1,2)])
  model
}

set.seed(123)
init_Z <- rnorm(p-1,0,1)
fit <- fitSSM(model = mod,
              inits = c(log(c(var_eta = 0.1, var_eps = 2000)),init_Z),
              updatefn = updt,
              update_args = list(samples = X), hessian=T)

fit$optim.out$convergence #PROBLEMA: non va a convergenza
fit$optim.out$value

round(c(exp(fit$optim.out$par[c(1,2)]),fit$optim.out$par[-c(1,2)]), 4)
sum(fit$optim.out$par[-c(1,2)])
as.data.frame(fit$model$Z)[p,p+1]

# Bias
df_bias <- as.data.frame(fit$model$Z)[,-1]
data.frame(Istituto=rownames(df_bias),Bias=diag(as.matrix(df_bias)))


kfs <- KFS(fit$model,
           filtering = "state",
           smoothing = "state")

dt_party_trend <- tibble(Data = dt_party$Data,
                         Trend = kfs$alphahat[, 1])

DT::datatable(dt_party_trend)

dt %>% filter(Partito == par) %>%
  ggplot(aes(x = Data, y = Percentuale)) +
  geom_point(color = "green", alpha = 0.2) +
  geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
  ggtitle(paste("Intenzioni di voto per", par))


param_estimates <- fit$optim.out$par

vcov_matrix <- solve(fit$optim.out$hessian)

std_errors <- sqrt(diag(vcov_matrix))

confidence_level <- 0.95
alpha <- 1 - confidence_level

z_value <- qnorm(1 - alpha / 2)

ci_lower <- (param_estimates - z_value * std_errors)[-c(1,2)]
ci_upper <- (param_estimates + z_value * std_errors)[-c(1,2)]

confidence_intervals <- data.frame(
  Parameter = colnames(Y)[-18],
  Estimate = param_estimates[-c(1,2)],
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
)

rownames(confidence_intervals) <- NULL
print(confidence_intervals)


# L'errore standard di una combinazione lineare si calcola con la varianza delle somme
bias_cov <- vcov_matrix[-c(1,2), -c(1,2)]

var_bias_last <- sum(bias_cov)

se_bias_last <- sqrt(var_bias_last)

estimate_bias_last <- -sum(param_estimates[-c(1,2)])

ci_lower_last <- estimate_bias_last - z_value * se_bias_last
ci_upper_last <- estimate_bias_last + z_value * se_bias_last

conf_int_last_bias <- data.frame(
  Parameter = colnames(Y)[18],
  Estimate = estimate_bias_last,
  CI_Lower = ci_lower_last,
  CI_Upper = ci_upper_last
)

rbind(confidence_intervals,conf_int_last_bias)

