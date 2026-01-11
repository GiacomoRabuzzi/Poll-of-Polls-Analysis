# Precalcola Austria

library(tidyverse)
library(KFAS)
library(Matrix)
library(readxl)
library(lubridate)
library(ggrepel)

dt <- read.csv("DatiAustria.csv")
dt <- dt[,-1]
dt$Partito <- as.factor((dt$Partito))
dt$Istituto <- as.factor((dt$Istituto))
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data, format = "%Y-%m-%d")


nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:6])
partiti <- unique(dt$Partito)

# Colori partiti
party_colors <- c(
  "ÖVP"   = "black",        # tradizionale colore del Partito Popolare Austriaco
  "SPÖ"   = "red3",         # Socialdemocratici, rosso
  "FPÖ"   = "cornflowerblue", # Partito della Libertà, blu
  "GRÜNE" = "forestgreen",  # Verdi
  "NEOS"  = "deeppink"      # NEOS – La Nuova Austria, magenta
)



source("SSM_functions_AT.R")  # questo contiene free_var_mod, bias_mod, ecc.

modelli <- list()
for (p in partiti) {
  cat("Calcolo modelli per", p, "...\n")
  co <- party_colors[[p]] %||% "black"
  modelli[[p]] <- list(
    bfv     = bias_free_var_mod(par = p, co = co),
    bin     = bin_mod(par = p, co = co),
    binbias = bin_bias_mod(party = p, co = co)
  )
}

saveRDS(modelli, file = "modelli_precotti_AT.rds")
