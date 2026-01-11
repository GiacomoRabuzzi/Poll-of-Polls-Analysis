rm(list=ls())

dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data,
                   format = "%d/%m/%Y")
DT::datatable(dt)

sort(table(dt$Istituto), decreasing = TRUE)[1:16]

nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:16])

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

p <- length(unique(dt$Istituto))
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
  var_eps <- exp(pars[-1])
  for (t in 1:nrow(samples)) {
    diag(model$H[,,t]) <- var_eps/samples[t, ]
  }
  model
}

fit <- fitSSM(model = mod,
              inits = log(c(var_eta = 0.1, var_eps = rep(2000,p))),
              updatefn = updt,
              update_args = list(samples = X), hessian=T, method='BFGS')

fit$optim.out$convergence
fit$optim.out$value

round(exp(fit$optim.out$par), 4)

vr <- data.frame(Istituti=colnames(Y),Vars=round(exp(fit$optim.out$par), 4)[-1])
rownames(vr) <- NULL
vr

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


param_estimates <- fit$optim.out$par

vcov_matrix <- solve(fit$optim.out$hessian)

std_errors <- sqrt(diag(vcov_matrix))

confidence_level <- 0.95
alpha <- 1 - confidence_level

z_value <- qnorm(1 - alpha / 2)

ci_lower <- exp(param_estimates - z_value * std_errors)[-1]
ci_upper <- exp(param_estimates + z_value * std_errors)[-1]

confidence_intervals <- data.frame(
  Parameter = colnames(Y),
  Estimate = exp(param_estimates[-1]),
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
)

rownames(confidence_intervals) <- NULL
print(confidence_intervals)

confidence_intervals %>% 
  rename(Istituto = Parameter,
         Varianza = Estimate,
         Lower95  = CI_Lower,
         Upper95  = CI_Upper) %>% 
  mutate(Istituto = fct_reorder(Istituto, Varianza)) %>% 
  
  ggplot(aes(x = Istituto, y = Varianza)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .25) +
  coord_flip() +
  labs(title = "Varianze specifiche degli istituti – IC95%",
       x      = "Istituto",
       y      = expression(hat(sigma)^2)) +
  theme_minimal(base_size = 13)




