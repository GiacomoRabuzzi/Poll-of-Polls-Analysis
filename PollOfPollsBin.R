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
mZ <- matrix(1, p, 1)

# transition equation matrices
mT <- matrix(1, 1, 1)
mQ <- matrix(NA_real_, 1, 1)
mR <- matrix(1, 1, 1)

# initial conditions
va1 <- matrix(-1.4, 1, 1)
mP1 <- matrix(0.3, 1, 1)
mP1inf <- matrix(0, 1, 1)

# let's make the SSModel
mod <- SSModel(Y~0+SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ, a1 = va1, P1 = mP1, P1inf = mP1inf,
                             state_names = "mu"), distribution = "binomial", u = U_trials)

updt <- function(pars, model) {
  model$Q[1, 1, 1] <- exp(pars[1])
  model
}

fit <- fitSSM(model = mod,
              inits = log(c(var_eta = 0.001)),
              updatefn = updt,
              method='BFGS')

fit$optim.out$convergence

round(exp(fit$optim.out$par), 4)

kfs <- KFS(fit$model,smoothing = "mean")

dt_party_trend <- tibble(Data = dt_party$Data,
                         Trend = kfs$muhat[,1]*100)

DT::datatable(dt_party_trend)

dt %>% filter(Partito == par) %>%
  ggplot(aes(x = Data, y = Percentuale)) +
  geom_point(color = "green", alpha = 0.2) +
  geom_line(aes(x = Data, y = Trend), data = dt_party_trend, linewidth = 1) +
  ggtitle(paste("Intenzioni di voto per", par))


# Importance sampling
imp <- importanceSSM(fit$model,type= "state",antithetics = TRUE)

w <- imp$weights/sum(imp$weights)

mean_imp <- numeric()
low_imp <- numeric()
up_imp <- numeric()

for(i in 1:nrow(Y_pct)){
  ilogit_imp <- plogis(imp$samples[i,1,])
  mean_imp[i]<-sum(ilogit_imp*w)
  oo <- order(ilogit_imp)
  low_imp[i] <- ilogit_imp[oo][which.min(abs(cumsum(w[oo]) - 0.025))]
  up_imp[i] <- ilogit_imp[oo][which.min(abs(cumsum(w[oo]) - 0.975))]
}

res <- data.frame(
  Data = dt_party$Data,
  mean = mean_imp*100,
  low  = low_imp*100,
  high = up_imp*100
)


ggplot(res, aes(x = Data, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "steelblue",  alpha = 0.25) +
  geom_line(size = 1, colour = "steelblue") +
  labs(title   = "Quota stimata della Lega nel tempo",
       subtitle = "Modello a stati binomiale – stima smoothata con importance sampling",
       y       = "(%)",
       x       = NULL,
       caption = "Banda = intervallo di credibilità al 95%") +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot",
        plot.caption.position = "plot")




