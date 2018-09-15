pwr_A <- function(rows_sums, cols_sums,  
                  n_repeats = 1000, n_matrices = 3000, 
                  alpha = .05, dev = .6, 
                  item_pos = 2, burnIn = 300, 
                  step = 16, folder = ""){
  ###############################################################
  # INPUTS:
  # rows_sums: Zeilenrandsummen
  # cols_sums: Spaltenrandsummen
  # n_repeats: Anzahl an Power-Werte pro Szenario
  # n_matrices: Anzahl an Matrizen pro Power-Wert
  # alpha: Fehler 1. Art
  # dev: DIF-Parameter
  # item_pos: Items mit Modellabweichung
  # burnIn: Burn-In Phase fuer Rasch Sampler
  # step: Step-Parameter fuer Rasch Sampler
  # folder: Speicherort
  ##############################################################
  model <- sample(rows_sums, cols_sums, 1)
  half_length <- length(cols_sums) / 2
  groups <- c(rep(1, half_length), rep(0, half_length)) 
  dif <- rep(0, length(cols_sums))
  dif[item_pos] <- dev
  
  path <- paste0(folder, "/",
                 as.character(length(rows_sums)), "x",
                 as.character(length(cols_sums)), ".csv")
  
  mcmc <- exact <- vector("numeric", n_repeats)
  
  count(rows_sums, cols_sums)
  
  mcmc <- replicate(n_repeats, 
                    pwr_mcmc(mat = model, 
                             group = groups,
                             dif = dif,
                             repetitions = n_matrices,
                             alpha = alpha, 
                             burn = burnIn,
                             steps = step))
  exact <- replicate(n_repeats, 
                     pwr_exact(rows = rows_sums, 
                               cols = cols_sums, 
                               group = groups, 
                               dif = dif, 
                               repetitions = n_matrices, 
                               alpha = alpha))
  
  rio::export(data.frame(power = c(mcmc, exact), 
                         method = rep(c("mcmc", "exact"),
                                      each = n_repeats)), path)
  
}

pwr_BCDE <- function(itempars, n_repeats = 3000, 
                     n_matrices = 8000, alpha = .05, 
                     n_pers = 100, sd_pers = 2, 
                     dev = .6, burnIn = 300, difficulty = "moderat", 
                     step = 16, folder = ""){
  ##############################################################
  # INPUTS:
  # itempars: Itemparameter
  # n_repeats: Anzahl an Power-Werte pro Szenario
  # n_matrices: Anzahl an Matrizen pro Power-Wert
  # alpha: Fehler 1. Art
  # n_pers: Anzahl Personen
  # n_items: Anzahl Items
  # sd_pers: Standardabweichung der Personenparameter
  # dev: DIF-Parameter
  # burnIn: Burn-In Phase fuer Rasch Sampler
  # difficulty: Itemschwierigkeit ("leicht", "moderat", "schwer")
  # step: Step-Parameter fuer Rasch Sampler
  # folder: Speicherort
  #############################################################
  set.seed(123)
  personenpars <- rnorm(n = n_pers, mean = 0, sd = sd_pers)
  half_length <- length(personenpars) / 2
  groups <- c(rep(1, half_length), rep(0, half_length)) 
  model <- sim.rasch(persons = personenpars, 
                     items = itempars,
                     seed = 123)
  cols_sums <- colSums(model)
  rows_sums <- rowSums(model) 
  dif <- vector("numeric", length(cols_sums))
  mcmc <- vector("numeric", n_repeats)
  path <- paste0(folder, "/",
                 as.character(length(rows_sums)), "x",
                 as.character(length(cols_sums)), "_",
                 as.character(dev), ".csv")
  switch(
    difficulty,
    "leicht" = dif[which(cols_sums == max(cols_sums[-length(itempars)]))[1]] <- dev,
    "moderat" = dif[which(cols_sums == getMiddle(cols_sums[-length(itempars)]))[1]] <- dev,
    "schwer" =  dif[which(cols_sums == min(cols_sums[-length(itempars)]))[1]] <- dev
  )
  mcmc <- replicate(n_repeats, 
                    pwr_mcmc(mat = model, 
                             group = groups,
                             dif = dif,
                             repetitions = n_matrices,
                             alpha = alpha, 
                             burn = burnIn,
                             steps = step))
  rio::export(data.frame(power = mcmc,
                         method = rep(c("mcmc"), n_repeats)), path)
}

pwr_exact <- function(rows, cols, group, 
                      dif, repetitions, alpha) {
  ################################################
  # INPUTS:
  # rows: Zeilenrandsummen
  # cols: Spaltenrandsummen
  # group: Gruppeneinteilung
  # dif: DIF-Parameter
  # repetitions: Anzahl an Matrizen pro Power-Wert
  # alpha: Fehler 1. Art
  ################################################
  s <- sample(a = rows, b = cols, k = repetitions) 
  t <- colSums(s * group)
  e <- exp(colSums(t * dif))
  pwr <- sum(e[e >= quantile(e, 1 - alpha)]) / sum(e)
  
  return(pwr)
}

pwr_mcmc <- function(mat, group, dif, repetitions, 
                     burn, steps, alpha) {
  #################################################
  # INPUTS:
  # mat: Rasch Modell
  # group: Gruppenaufteilung
  # dif: DIF-Parameter
  # repetitions: Anzahl an Matrizen pro Power-Wert
  # burn: Burn-In Phase
  # steps: Step-Parameter
  # alpha: Fehler 1. Art
  #################################################
  s <- rsampler(mat, controls = rsctrl(n_eff = (repetitions - 1),
                                       burn_in = burn, 
                                       step = steps))
  t <- rstats(s, function(x) colSums(x * group))
  e <- exp(colSums(matrix(unlist(t), ncol = s$n_tot) * dif))
  pwr <- sum(e[e >= quantile(e, 1 - alpha)]) / sum(e)
  
  return(pwr)
}

auswertung_summary <- function(df, col1, col2){
  ##########################################################
  # Auswertung mit Summary im tidy Format
  # INPUTS:
  # df: Datensatz als data.frame
  # col1: Gruppierende Spalte als quosure
  # col2: Auszuwertende Spalte als quosure
  # z.B.: auswertung_summary(daten, quo(method), quo(power))
  ##########################################################
  df %>%
    group_by(!!col1) %>%
    summarise(min = min(!!col2, na.rm = T),
              Q.025 = quantile(!!col2, .025, na.rm = T),
              Q.25 = quantile(!!col2, .25, na.rm = T),
              median = median(!!col2, na.rm = T),
              mean = mean(!!col2, na.rm = T),
              Q.75 = quantile(!!col2, .75, na.rm = T),
              Q.975 = quantile(!!col2, .975, na.rm = T),
              max = max(!!col2, na.rm = T),
              sd = sd(!!col2, na.rm = T))
}

auswertung <- function(df, col1, col2){
  ########################################
  # Allgemeine Auswertungsfunktion
  # INPUTS:
  # df: Datensaetze als Liste
  # col1: Gruppierende Spalte
  # col2: Auszuwertende Spalte
  # z.B.: auswertung(daten, method, power)
  ########################################
  col1 <- enquo(col1)
  col2 <- enquo(col2)
  
  # Allgemeine Auswertung mit Summary im tidy Format
  a <- df %>%
    map_df(~ auswertung_summary(df = .x, col1 = col1, col2 = col2)) %>%
    mutate(szenario = rep(names(df), each = 2))
  # Vergleich der Standardabweichungen
  b <- a %>% 
    select(!!col1, sd, szenario) %>%
    spread(!!col1, sd) %>%
    mutate(mcmc_smaller = mcmc < exact)
  
  return(list(Allgemeine_Auswertung = a, 
              Vergleich_Standardabweichungen = b))
}

getMiddle <- function(x){
  ##########################
  # INPUTS:
  # x: Numerischer Vektor
  ##########################
  # Mittlere Zahl zurueckgeben
  # bei ungeraden Itemanzahlen rundet R ab
  # (z.B. bei 10.5 nimmt es den 10. Index)
  sorted <- sort(x)
  middle <- sorted[length(x) / 2]
  return(middle)
}