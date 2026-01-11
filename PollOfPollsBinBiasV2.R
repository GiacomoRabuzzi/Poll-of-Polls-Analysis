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

Y_pct <- as.matrix(dt_party[, -1]) / 100
U_trials <- as.matrix(round(dt_samples[, -1]))
Y <- round(Y_pct * U_trials)

p <- length(unique(dt$Istituto))
n <- nrow(dt_party)

# observation equation matrices
mZ <- cbind(matrix(1, p-1, 1),diag(1,p-1))
mZ <- rbind(mZ,c(1,rep(-1,p-1)))

# transition equation matrices
mT <- diag(1,p)
mQ <- matrix(NA_real_, 1, 1)
mR <- rbind(matrix(1, 1, 1),matrix(0, p-1, 1))

# initial conditions
library(Matrix)
va1   <- matrix(c(-1.3, rep(0, p-1)), p, 1)
mP1   <- as.matrix(bdiag(matrix(1, 1, 1), diag(10, p-1)))
mP1inf <- matrix(0, p, p)

# let's make the SSModel
s_names <- "mu"
d_names <- paste0("delta_", 1:(p-1))
s_names <- c(s_names, d_names)

mod <- SSModel(Y~0+SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ, a1 = va1, P1 = mP1, P1inf = mP1inf,
                             state_names = s_names), distribution = "binomial", u = U_trials)

updt <- function(pars, model) {
  model$Q[1, 1, 1] <- exp(pars[1])
  model
}

fit <- fitSSM(model = mod,
              inits = log(c(0.05)),
              updatefn = updt,
              method='BFGS')

fit$optim.out$convergence

round(exp(fit$optim.out$par), 4)

kfs <- KFS(fit$model,
           smoothing = c("state","signal","mean"))

kfs$alphahat
kfs$thetahat
kfs$muhat

kfs$alphahat[,1]


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



# Bias Analysis

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
  d  <- delta_samples[j, ]
  m  <- sum(d * w)
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


library(forcats)

bias_df %>% 
  mutate(Istituto = fct_reorder(Istituto, Bias_pp)) %>% 
  ggplot(aes(x = Istituto, y = Bias_pp)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Low_95_pp, ymax = High_95_pp),
                width = .25, linewidth = .6) +
  coord_flip() +
  labs(
    title = paste("Bias stimato degli istituti –", par),
    subtitle = "Effetto-casa (δ) con intervallo di credibilità al 95 %",
    x = NULL,
    y = "Bias (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot")


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
  left_join(bias_df %>% select(Istituto, Bias_modello = Bias_pp),
            by = "Istituto")

library(ggrepel)

ggplot(df_estrea, aes(Bias_reale, Bias_modello, label = Istituto)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = .4) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkgreen",
              linetype = "dashed", linewidth = .7) +
  geom_point(size = 2, colour = "forestgreen") +
  geom_label_repel(size = 3, label.padding = unit(0.15, "lines"),
                   colour = "darkgreen", max.overlaps = Inf) +
  labs(
    title = paste("Confronto bias reale vs stimato –", par),
    subtitle = "Ultimo sondaggio prima del voto (25 set 2022) vs house-effect",
    x = "Bias reale (p.p.)", y = "Bias modello (p.p.)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title.position = "plot")



