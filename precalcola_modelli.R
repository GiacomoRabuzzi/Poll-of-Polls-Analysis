# precalcola_modelli.R

library(tidyverse)
library(KFAS)
library(Matrix)
library(readxl)
library(lubridate)
library(ggrepel)

dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data, format = "%d/%m/%Y")
nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:16])
partiti <- sort(setdiff(unique(dt$Partito), c("Az", "IV")))

# Colori partiti
party_colors <- c(
  "Lega"         = "forestgreen",
  "PD"           = "red3",
  "FdI"          = "black",
  "FI"           = "blue3",
  "M5S"          = "goldenrod",
  "Verdi-SI"     = "darkgreen",
  "Più Europa"   = "darkorange",
  "Noi Moderati" = "deepskyblue",
  "Italexit"     = "darkred",
  "Unione Popolare" = "purple",
  "Art.1"        = "darkolivegreen",
  "Altri"        = "grey40"
)

# Carica le funzioni dei modelli (sorgente dal tuo file R)
source("SSM_functions.R")  # questo contiene free_var_mod, bias_mod, ecc.

modelli <- list()
for (p in partiti) {
  cat("Calcolo modelli per", p, "...\n")
  co <- party_colors[[p]] %||% "black"
  modelli[[p]] <- list(
    freevar = free_var_mod(par = p, co = co),
    bias    = bias_mod(par = p, co = co),
    bfv     = bias_free_var_mod(par = p, co = co),
    bin     = bin_mod(par = p, co = co),
    binbias = bin_bias_mod(party = p, co = co)
  )
}

saveRDS(modelli, file = "modelli_precotti.rds")
