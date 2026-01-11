rm(list=ls())

library(tidyverse)

dt <- read.csv2("Sondaggi_2024-04-30.csv")
dt$Percentuale <- as.numeric(dt$Percentuale)
dt$Data <- as.Date(dt$Data, format = "%d/%m/%Y")
nomi_istituti <- names(sort(table(dt$Istituto), decreasing = TRUE)[1:16])
partiti <- sort(setdiff(unique(dt$Partito), c("Az", "IV")))

sort(table(dt$Istituto),decreasing=T)

dt_filtro <- subset(dt, Istituto %in% nomi_istituti)

tab_cam <- dt_filtro %>%
  group_by(Istituto) %>%
  summarise(
    Media = mean(Campione, na.rm = TRUE),
    Mediana = median(Campione, na.rm = TRUE),
    Minimo = min(Campione, na.rm = TRUE),
    Massimo = max(Campione, na.rm = TRUE),
    Deviazione_Std = sd(Campione, na.rm = TRUE),
    Perc_Null = sum(is.na(Campione))/n()*100,
    Conteggio_Totale = n(),
    min_date = min(Data),
    max_date = max(Data)
  )

as.data.frame(tab_cam[order(tab_cam$Conteggio_Totale,decreasing = T),])


library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

multiTimeline <- read_csv("C:/Users/giaco/Downloads/multiTimeline.csv", skip = 2)

multiTimeline <- multiTimeline %>%
  mutate(Giorno = as.Date(paste0(Mese, "-01"), format = "%Y-%m-%d"))

long <- multiTimeline %>%
  select(-matches("isPartial", ignore.case = TRUE)) %>%
  pivot_longer(cols = -c(Mese, Giorno), names_to = "termine", values_to = "valore") %>%
  mutate(valore = as.numeric(valore))

eu_dates  <- as.Date(c("2004-06-01","2009-06-01","2014-05-01","2019-05-01","2024-06-01"))
pol_dates <- as.Date(c("2006-04-01","2008-04-01","2013-02-01","2018-03-01","2022-09-01"))

vlines <- tibble(
  date = c(eu_dates, pol_dates),
  tipo = c(rep("European elections", length(eu_dates)), 
           rep("National elections", length(pol_dates)))
)


ggplot(long, aes(x = Giorno, y = valore, color = termine)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  geom_vline(data = vlines, aes(xintercept = date, linetype = tipo),
             color = "grey20", linewidth = 0.5, alpha = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = NULL, y = "Interest (Google Trends)",
    color = "Search term",
    linetype = "Elections",
    title = "Google Trends interest for the word Sondaggi (Polls) in Italy",
    subtitle = "Vertical lines: European (---) and National (–)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


library(ggplot2)
library(dplyr)

# dati dal grafico originale
polls <- tibble(
  year = c(1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004),
  number = c(12, 10, 22, 28, 105, 85, 195, 245, 165)
)

ggplot(polls, aes(x = factor(year), y = number)) +
  geom_col(fill = "black", width = 0.7) +
  geom_hline(yintercept = seq(50, 250, 50), 
             color = "grey40", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 250), expand = c(0,0)) +
  labs(
    x = NULL,
    y = "Number of Polls",
    title = "Number of polls by election year",
    caption = "Sources: Roper Center IPOLL, PollingReport.com, NationalJournal.com"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.caption = element_text(size = 9, hjust = 0),
    axis.text.x = element_text(size = 10)
  )


dt_filtro2 <- subset(dt_filtro, Partito %in% partiti)
dt2 <- select(dt_filtro2, c('Partito','Data','Percentuale'))


library(dplyr)
library(ggplot2)
library(slider)
library(lubridate)
library(purrr)

sd_roll <- dt_filtro2 %>%
  arrange(Partito, Data) %>%
  group_by(Partito) %>%
  mutate(
    sd_21d = slide_index_dbl(
      .x = Percentuale,
      .i = Data,
      .f = ~ if (length(.x) >= 2) sd(.x, na.rm = TRUE) else NA_real_,
      .before = days(21),    # finestra retrospettiva: 21 giorni
      .complete = TRUE
    )
  ) %>%
  ungroup()

elections <- vlines$date
is_in_post21 <- function(d) any(d > elections & d <= elections + days(21))

sd_roll <- sd_roll %>%
  mutate(
    post21 = map_lgl(Data, is_in_post21),
    sd_plot = if_else(post21, NA_real_, sd_21d)  # NA = buco nel grafico
  )

ggplot(sd_roll, aes(x = Data, y = sd_plot)) +
  geom_line(color = "steelblue", linewidth = 0.8, na.rm = TRUE) +
  geom_vline(data = vlines, aes(xintercept = date, linetype = tipo),
             color = "grey20", linewidth = 0.5, alpha = 0.85) +
  facet_wrap(~ Partito, scales = "free_y", ncol = 2) +
  labs(
    x = NULL, y = "Rolling standard deviation (21 days)",
    linetype = "Elections",
    title = "Convergence of pollings by party"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")



