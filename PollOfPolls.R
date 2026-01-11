rm(list=ls())

dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data,
                   format = "%d/%m/%Y")
DT::datatable(dt)

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

p <- length(unique(dt$Istituto)) # number of institutes (40)
n <- nrow(dt_party)

# data have to turned into matrix
Y <- as.matrix(dt_party[,-1]) # just eliminate the date
X <- as.matrix(dt_samples[, -1])

# observation equation matrices
mZ <- matrix(1, p, 1)
mH <- array(0, c(p, p, n)) # this changes every day => 3-dim array

# transition equation matrices
mT <- matrix(1, 1, 1)
mQ <- matrix(NA_real_, 1, 1)
mR <- matrix(1, 1, 1)

# initial conditions
va1 <- matrix(20, 1, 1)
mP1 <- matrix(56.25, 1, 1)
mP1inf <- matrix(0, 1, 1)

# let's make the SSModel
mod <- SSModel(Y~0+SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf,
                             state_names = "mu"),
               H = mH)

# let's make the update function
updt <- function(pars, model, samples) {
  model$Q[1, 1, 1] <- exp(pars[1])
  var_eps <- exp(pars[2])
  for (t in 1:nrow(samples)) {
    diag(model$H[,,t]) <- var_eps/samples[t, ]
  }
  model
}

fit <- fitSSM(model = mod,
              inits = log(c(var_eta = 0.1, var_eps = 2000)),
              updatefn = updt,
              update_args = list(samples = X))

fit$optim.out$convergence

round(exp(fit$optim.out$par), 4)

kfs <- KFS(fit$model,
           filtering = "state",
           smoothing = "state")

# let us make a dataset with the trend
dt_party_trend <- tibble(Data = dt_party$Data,
                         Trend = kfs$alphahat[, 1])

DT::datatable(dt_party_trend)

dt %>% filter(Partito == par) %>%
  ggplot(aes(x = Data, y = Percentuale)) +
  geom_point(color = "green", alpha = 0.2) +
  geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
  ggtitle(paste("Intenzioni di voto per", par))






