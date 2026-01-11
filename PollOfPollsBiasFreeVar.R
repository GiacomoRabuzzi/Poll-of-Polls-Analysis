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
mZ <- cbind(matrix(1, p, 1),diag(1,p))
mZ <- mZ[,-(p+1)]
mZ[p,] <- c(1,rep(-1,p-1))
mH <- array(0, c(p, p, n))

# transition equation matrices
mT <- diag(1,p)
mQ <- matrix(NA_real_, 1, 1)
mR <- rbind(matrix(1, 1, 1),matrix(0, p-1, 1))

# initial conditions
library(Matrix)
va1 <- matrix(c(20,rep(0,p-1)), p, 1)
mP1 <- as.matrix(bdiag(matrix(56.25, 1, 1), diag(10, p-1)))
mP1inf <- matrix(0, p, p)

# let's make the SSModel
s_names <- "mu"
d_names <- paste0("delta_", 1:(p-1))
s_names <- c(s_names, d_names)
mod <- SSModel(Y~0+SSMcustom(mZ, mT, mR, mQ, va1, mP1, mP1inf,
                             state_names = s_names),
               H = mH)

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
              update_args = list(samples = X), hessian=T, method="BFGS")

fit$optim.out$convergence
fit$optim.out$value

round(exp(fit$optim.out$par), 4)

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


# Variance analysis

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

# Bias Analysis
delta_idx <- 2:p

last_t <- nrow(kfs$alphahat)

delta_hat <- kfs$alphahat[last_t, delta_idx] # ultimo ma tutti uguali

var_delta <- diag(kfs$V[delta_idx, delta_idx, last_t])
se_delta <- sqrt(var_delta)

z <- qnorm(0.975)
ci_low <- delta_hat - z * se_delta
ci_high <- delta_hat + z * se_delta

bias_tbl <- tibble(
  Istituto = colnames(Y)[1:(p - 1)],  # solo i primi p-1 istituti
  Bias = delta_hat,
  Lower95 = ci_low,
  Upper95 = ci_high
) %>%
  arrange(desc(Bias))

bias_ultimo <- -sum(delta_hat)
var_ultimo <- sum(kfs$V[delta_idx, delta_idx, last_t])
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

DT::datatable(
  bias_tbl,
  caption = "House‑effect stimati con intervalli di confidenza al 95%"
)

bias_tbl %>%
  ggplot(aes(x = reorder(Istituto, Bias), y = Bias)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .2) +
  coord_flip() +
  labs(x = "Istituto", y = "Bias (%)",
       title = "Bias stimato degli istituti (effetto 'house') con IC 95%")


# Elezioni
library(readxl)
elez <- as.data.frame(read_xlsx("Elezioni2022.xlsx"))
elez[,1] <- c(unique(dt$Partito)[1],unique(dt$Partito)[5],unique(dt$Partito)[4],
              unique(dt$Partito)[3],unique(dt$Partito)[2])

library(lubridate)
data_riferimento <- as.Date("2022-09-25")

partiti_esclusi <- c("Az", "IV")

df_preelez <- dt %>%
  filter(Data <= data_riferimento,
         !(Partito %in% partiti_esclusi)) %>%
  group_by(Partito, Istituto) %>%
  filter(Data == max(Data)) %>%
  ungroup() %>%
  select(Partito, Istituto, Data, Percentuale) %>%
  arrange(Partito, Istituto)

df_preelez%>%
  filter(Partito==par)%>%
  mutate(Bias=Percentuale-elez[elez$Partito==par,2])%>%
  select(Istituto,Bias) %>%
  arrange(Istituto)-> bias_reale

df_estrea <- bias_reale %>%
  rename(Bias_reale = Bias) %>%
  left_join(bias_tbl %>% select(Istituto, Bias_modello = Bias),
            by = "Istituto")

library(ggrepel)
ggplot(df_estrea, aes(x = Bias_reale, y = Bias_modello, label = Istituto)) +
  geom_smooth(method = "lm", se = FALSE, color = "green4", 
              linetype = "dashed",linewidth=0.6) +
  geom_label_repel(size = 3, label.padding = unit(0.15, "lines"), color = "green3") +
  
  labs(
    x = "Bias reale (sondaggio - risultato elettorale)",
    y = "Bias stimato dal modello",
    title = paste("Confronto Bias Reale vs Stimato -", par)
  ) +
  theme_minimal()


# Scomposizione MSE
ci <- confidence_intervals %>%
  rename(Istituto = Parameter,
         Varianza_mod = Estimate)

var_eps_vec <- exp(fit$optim.out$par[-1])
names(var_eps_vec) <- colnames(Y)  

dt_cmse <- dt[dt$Partito==par,]
dt_cmse$Campione[is.na(dt_cmse$Campione)] <- min(dt_cmse$Campione,na.rm=T)
avg_camp <- tapply(dt_cmse$Campione,dt_cmse$Istituto,mean)
  
avg_camp <- avg_camp[colnames(Y)]

var_media <- tibble(Istituto=colnames(Y), Var_media=var_eps_vec/avg_camp)

mse_tbl <- bias_tbl %>%
  select(Istituto, Bias) %>%
  left_join(var_media, by = "Istituto") %>%
  mutate(
    MSE        = Bias^2 + Var_media,
    PropBias   = (Bias^2) / MSE,
    PropVar    = Var_media / MSE
  )

mse_tbl <- mse_tbl %>%
  mutate(
    Bias2     = Bias^2,
    Varianza  = Var_media
  ) %>%
  select(Istituto, Bias2, Varianza) %>%
  pivot_longer(-Istituto,
               names_to = "Componente",
               values_to = "Valore")

# plot stacked
ggplot(mse_tbl, aes(x = reorder(Istituto, Valore), 
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

