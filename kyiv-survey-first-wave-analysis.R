# ============================================================
# Trust Kyiv Wave 1 — Revised SEM (5 separate models)
# ============================================================
# Theoretical model:
#   Infrastructure Deprivation (additive index) → Negative Affect → DV
#   + direct path Deprivation → DV
#   Controls: SES (income), age, female
#
# 5 separate models:
#   Model 1: DV = INTTRUST        (interpersonal trust, latent, 4 items)
#   Model 2: DV = trust_president (observed indicator)
#   Model 3: DV = trust_rada      (observed indicator)
#   Model 4: DV = trust_local     (observed indicator)
#   Model 5: DV = PROSOC          (prosocial behavior, latent, 4 items)
# ============================================================

# ============================================================
# SETUP
# ============================================================

sink("output_log.txt", split = TRUE)

library(dplyr)
library(lavaan)
library(ggplot2)
library(tidyr)
library(lavaangui)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ============================================================
# STEP 1. DATA LOADING & RECODING
# ============================================================

load("Trust_Kyiv_wave1.RData")
df <- trust_wave1.0

model_df <- data.frame(id = df[["Respondent ID"]])

# --- Infrastructure deprivation components -------------------

# deprivation: 4 - number of available utilities (0-4, higher = worse)
model_df$deprivation <- 4 - (
  df[["Електропостачання"]] +
  df[["Гаряча вода з крану"]] +
  df[["Холодна вода з крану"]] +
  df[["Централізоване опалення"]]
)

# temp_cold: temperature in apartment (higher = colder = worse)
model_df$temp_cold <- case_match(
  df[["Вкажіть (хоча б приблизно), яка у Вас зараз температура у квартирі:"]],
  "22\u00b0C і більше"          ~ 0,
  "19-21\u00b0C"                ~ 1,
  "16-18\u00b0C"                ~ 2,
  "10-15\u00b0C"                ~ 3,
  "0-9\u00b0C"                  ~ 4,
  "Менше, ніж 0\u00b0C"         ~ 5,
  "Не знаю навіть приблизно" ~ NA_real_
)

# --- Negative affect (binary indicators) ---------------------

model_df$mood_sad        <- df[["Засмученим(-ою)"]]
model_df$mood_nervous    <- df[["Знервованим(-ою)"]]
model_df$mood_scared     <- df[["Наляканим(-ою)"]]
model_df$mood_aggressive <- df[["Агресивним(-ою)"]]

model_df$neg_affect_sum <- model_df$mood_sad +
                           model_df$mood_nervous +
                           model_df$mood_scared +
                           model_df$mood_aggressive

# --- Trust variables (0-10 scale) ----------------------------

model_df$trust_neighbors <- df[["Наскільки Ви довіряєте Вашим сусідам?"]]
model_df$trust_district  <- df[["Наскільки Ви довіряєте мешканцям району, в якому Ви живете?"]]
model_df$trust_city      <- df[["Наскільки Ви довіряєте мешканцям вашого міста?"]]
model_df$trust_ukraine   <- df[["Наскільки Ви довіряєте тим, хто проживає в Україні?"]]
model_df$trust_diaspora  <- df[["Наскільки Ви довіряєте українцям, які живуть за кордоном?"]]
model_df$trust_president <- df[["Наскільки Ви довіряєте Президенту України?"]]
model_df$trust_rada      <- df[["Наскільки Ви довіряєте Верховній Раді України?"]]
model_df$trust_local     <- df[["Наскільки Ви довіряєте місцевій владі?"]]

# --- Prosocial behavior (1-7 scale) --------------------------

help_cols <- names(df)[grepl("допомогти іншій людині", names(df))]
model_df$help_friend   <- as.numeric(df[[help_cols[grepl("підтримали", help_cols)]]])
model_df$help_lost     <- as.numeric(df[[help_cols[grepl("загублену", help_cols)]]])
model_df$help_sick     <- as.numeric(df[[help_cols[grepl("хворим", help_cols)]]])
model_df$help_stranger <- as.numeric(df[[help_cols[grepl("донести", help_cols)]]])

# --- SES: income only (higher = poorer) ----------------------

model_df$income_low <- case_match(
  df[["Оцінка матеріального стану"]],
  "Можемо дозволити собі придбати майже все, що захочемо"                                                                                     ~ 0,
  "Можемо легко купуємо товари тривалого вжитку, але деякі покупки (квартира, автомобіль тощо) ми поки що не в змозі"                          ~ 1,
  "Грошей вистачає на їжу та одяг, але купівля товарів тривалого вжитку (новий холодильник або телевізор) викликає труднощі"                    ~ 2,
  "Вистачає на їжу, але купівля одягу проблематична"                                                                                           ~ 3,
  "Бракує грошей навіть на їжу"                                                                                                                ~ 4,
  "Важко відповісти"                                                                                                                           ~ NA_real_
)

# --- Controls ------------------------------------------------

model_df$age <- as.numeric(df[["Вік респондента"]])

model_df$female <- case_match(
  df[["Стать"]],
  "Жінка"   ~ 1,
  "Чоловік" ~ 0
)

cat("Data loaded and recoded:", nrow(model_df), "obs,", ncol(model_df), "vars\n\n")

# ============================================================
# STEP 2. ADDITIVE INDICES
# ============================================================

# Deprivation index (formative: sum of 2 components)
model_df$depriv_index <- model_df$deprivation + model_df$temp_cold

cat("========================================\n")
cat("ADDITIVE INDICES — DESCRIPTIVES\n")
cat("========================================\n\n")

idx_vars <- c("depriv_index", "income_low", "neg_affect_sum")
for (v in idx_vars) {
  vals <- model_df[[v]]
  cat(sprintf("  %-18s N=%d  M=%.2f  SD=%.2f  [%.1f, %.1f]\n",
              v,
              sum(!is.na(vals)),
              mean(vals, na.rm = TRUE),
              sd(vals, na.rm = TRUE),
              min(vals, na.rm = TRUE),
              max(vals, na.rm = TRUE)))
}

# ============================================================
# STEP 3. CFA — MEASUREMENT MODEL
# ============================================================

# --- Individual CFA models per factor -----------------------

cfa_specs <- list(
  list(
    name  = "Interpersonal Trust",
    label = "INTTRUST",
    model = 'INTTRUST =~ trust_district + trust_city + trust_ukraine + trust_diaspora
              trust_district ~~ trust_city'
  ),
  list(
    name  = "Prosocial Behavior",
    label = "PROSOC",
    model = 'PROSOC =~ help_friend + help_lost + help_sick + help_stranger'
  )
)

# Full CFA (all latent factors together)
cfa_model <- '
  INTTRUST  =~ trust_district + trust_city + trust_ukraine + trust_diaspora
  PROSOC    =~ help_friend + help_lost + help_sick + help_stranger
  trust_district ~~ trust_city
'

cat("\n========================================\n")
cat("STEP 3: CFA — MEASUREMENT MODEL\n")
cat("========================================\n\n")

compute_omega <- function(fit, factor_name) {
  sl <- standardizedSolution(fit)
  lam <- sl[sl$op == "=~" & sl$lhs == factor_name, "est.std"]
  error_var <- 1 - lam^2
  omega <- sum(lam)^2 / (sum(lam)^2 + sum(error_var))
  return(round(omega, 3))
}

# --- 3a. Per-factor CFA fit ----------------------------------

cat("--- Per-Factor CFA Fit ---\n\n")

cfa_factor_fits <- list()

for (spec in cfa_specs) {
  fit_i <- cfa(spec$model, data = model_df, estimator = "MLR")
  fi <- fitMeasures(fit_i, c("chisq", "df", "pvalue",
                              "cfi", "tli", "rmsea",
                              "rmsea.ci.lower", "rmsea.ci.upper",
                              "srmr"))
  om <- compute_omega(fit_i, spec$label)
  cfa_factor_fits[[spec$label]] <- fi

  cat(sprintf("  %s:\n", spec$name))
  cat(sprintf("    Chi-sq = %.2f, df = %.0f, p = %.4f\n",
              fi["chisq"], fi["df"], fi["pvalue"]))
  cat(sprintf("    CFI = %.3f  TLI = %.3f  RMSEA = %.3f [%.3f, %.3f]  SRMR = %.3f\n",
              fi["cfi"], fi["tli"], fi["rmsea"],
              fi["rmsea.ci.lower"], fi["rmsea.ci.upper"], fi["srmr"]))
  cat(sprintf("    Omega = %.3f\n\n", om))
}

# --- 3b. Full CFA (all factors) ------------------------------

cat("--- Full CFA (all factors) ---\n\n")

cfa_fit <- cfa(cfa_model, data = model_df, estimator = "MLR")

summary(cfa_fit, fit.measures = TRUE, standardized = TRUE)

fit_cfa <- fitMeasures(cfa_fit, c("chisq", "df", "pvalue",
                                   "cfi", "tli", "rmsea",
                                   "rmsea.ci.lower", "rmsea.ci.upper",
                                   "srmr"))

cat("\n--- Full CFA Fit Indices ---\n")
cat(sprintf("  Chi-sq = %.2f, df = %.0f, p = %.4f\n",
            fit_cfa["chisq"], fit_cfa["df"], fit_cfa["pvalue"]))
cat(sprintf("  CFI    = %.3f\n", fit_cfa["cfi"]))
cat(sprintf("  TLI    = %.3f\n", fit_cfa["tli"]))
cat(sprintf("  RMSEA  = %.3f [%.3f, %.3f]\n",
            fit_cfa["rmsea"], fit_cfa["rmsea.ci.lower"], fit_cfa["rmsea.ci.upper"]))
cat(sprintf("  SRMR   = %.3f\n", fit_cfa["srmr"]))

# --- 3c. Standardized loadings ------------------------------

cat("\n--- Standardized Factor Loadings ---\n\n")
std_loadings <- standardizedSolution(cfa_fit)
loadings_only <- std_loadings[std_loadings$op == "=~", ]
print(loadings_only[, c("lhs", "rhs", "est.std", "se", "pvalue")], row.names = FALSE)

# --- 3d. Reliability (omega) for full model -----------------

cat("\n--- Factor Reliability (Omega from Full CFA) ---\n\n")

factors       <- c("INTTRUST", "PROSOC")
factor_labels <- c("Interpersonal Trust", "Prosocial Behavior")

for (i in seq_along(factors)) {
  om <- compute_omega(cfa_fit, factors[i])
  cat(sprintf("  %-30s omega = %.3f\n", factor_labels[i], om))
}

# ============================================================
# STEP 4. FIVE SEPARATE SEM MODELS
# ============================================================
# Each model: Deprivation → NegAff → DV
#             + direct Deprivation → DV
#             + controls: income_low, age, female
# Models 1 & 5: latent DV; Models 2-4: observed DV

cat("\n\n========================================\n")
cat("STEP 4: FIVE SEPARATE SEM MODELS\n")
cat("========================================\n\n")

# --- Model syntax definitions --------------------------------

sem_inttrust <- '
  INTTRUST =~ trust_district + trust_city + trust_ukraine + trust_diaspora
  trust_district ~~ trust_city

  neg_affect_sum ~ a*depriv_index + income_low
  INTTRUST       ~ b*neg_affect_sum + c*depriv_index + income_low + age + female

  indirect := a * b
'

sem_president <- '
  neg_affect_sum  ~ a*depriv_index + income_low
  trust_president ~ b*neg_affect_sum + c*depriv_index + income_low + age + female

  indirect := a * b
'

sem_rada <- '
  neg_affect_sum ~ a*depriv_index + income_low
  trust_rada     ~ b*neg_affect_sum + c*depriv_index + income_low + age + female

  indirect := a * b
'

sem_local <- '
  neg_affect_sum ~ a*depriv_index + income_low
  trust_local    ~ b*neg_affect_sum + c*depriv_index + income_low + age + female

  indirect := a * b
'

sem_prosoc <- '
  PROSOC =~ help_friend + help_lost + help_sick + help_stranger

  neg_affect_sum ~ a*depriv_index + income_low
  PROSOC         ~ b*neg_affect_sum + c*depriv_index + income_low + age + female

  indirect := a * b
'

model_specs <- list(
  list(name = "Interpersonal Trust",  label = "INTTRUST",        syntax = sem_inttrust),
  list(name = "Trust: President",     label = "trust_president",  syntax = sem_president),
  list(name = "Trust: Rada",          label = "trust_rada",       syntax = sem_rada),
  list(name = "Trust: Local Gov",     label = "trust_local",      syntax = sem_local),
  list(name = "Prosocial Behavior",   label = "PROSOC",           syntax = sem_prosoc)
)

# --- Fit all models ------------------------------------------

sem_fits     <- list()
sem_fit_idx  <- list()
all_struct   <- data.frame()
all_indirect <- data.frame()
all_r2       <- data.frame()

for (spec in model_specs) {

  cat("--------------------------------------------\n")
  cat("Model:", spec$name, "\n")
  cat("--------------------------------------------\n\n")

  fit <- sem(spec$syntax, data = model_df, estimator = "MLR")

  sem_fits[[spec$label]] <- fit

  cat("--- Summary ---\n\n")
  summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

  # Fit indices
  fi <- fitMeasures(fit, c("chisq", "df", "pvalue",
                            "cfi", "tli", "rmsea",
                            "rmsea.ci.lower", "rmsea.ci.upper",
                            "srmr"))
  sem_fit_idx[[spec$label]] <- fi

  cat(sprintf("\n  CFI = %.3f  TLI = %.3f  RMSEA = %.3f [%.3f, %.3f]  SRMR = %.3f\n\n",
              fi["cfi"], fi["tli"], fi["rmsea"],
              fi["rmsea.ci.lower"], fi["rmsea.ci.upper"], fi["srmr"]))

  # Structural paths
  std_sol <- standardizedSolution(fit)
  struct <- std_sol[std_sol$op == "~", ]
  struct$model <- spec$name
  struct$model_label <- spec$label
  all_struct <- rbind(all_struct, struct)

  # Indirect effects
  ind <- std_sol[std_sol$op == ":=", ]
  ind$model <- spec$name
  ind$model_label <- spec$label
  all_indirect <- rbind(all_indirect, ind)

  # R-squared
  r2 <- lavInspect(fit, "rsquare")
  r2_df <- data.frame(variable = names(r2), r2 = unname(r2), model = spec$name)
  all_r2 <- rbind(all_r2, r2_df)
}

# --- Modification indices for Interpersonal Trust model ------

cat("\n========================================\n")
cat("MODIFICATION INDICES: Interpersonal Trust\n")
cat("========================================\n\n")

mi_inttrust <- modindices(sem_fits[["INTTRUST"]], sort. = TRUE, minimum.value = 4)
print(mi_inttrust[, c("lhs", "op", "rhs", "mi", "epc", "sepc.all")], row.names = FALSE)

# ============================================================
# STEP 5. MODEL FIT SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("STEP 5: MODEL FIT SUMMARY\n")
cat("========================================\n\n")

fit_indices <- c("chisq", "df", "pvalue", "cfi", "tli",
                 "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")

# CFA per-factor
cat("--- CFA (Per-Factor) ---\n\n")
for (spec in cfa_specs) {
  fi <- cfa_factor_fits[[spec$label]]
  cat(sprintf("  %-25s  Chi2=%7.2f  df=%2.0f  p=%.4f  CFI=%.3f  TLI=%.3f  RMSEA=%.3f [%.3f, %.3f]  SRMR=%.3f\n",
              spec$name,
              fi["chisq"], fi["df"], fi["pvalue"],
              fi["cfi"], fi["tli"],
              fi["rmsea"], fi["rmsea.ci.lower"], fi["rmsea.ci.upper"],
              fi["srmr"]))
}

# CFA full model
cat(sprintf("\n  %-25s  Chi2=%7.2f  df=%2.0f  p=%.4f  CFI=%.3f  TLI=%.3f  RMSEA=%.3f [%.3f, %.3f]  SRMR=%.3f\n",
            "Full CFA (2 factors)",
            fit_cfa["chisq"], fit_cfa["df"], fit_cfa["pvalue"],
            fit_cfa["cfi"], fit_cfa["tli"],
            fit_cfa["rmsea"], fit_cfa["rmsea.ci.lower"], fit_cfa["rmsea.ci.upper"],
            fit_cfa["srmr"]))

# SEM models
cat("\n--- SEM (Structural) ---\n\n")
for (spec in model_specs) {
  fi <- sem_fit_idx[[spec$label]]
  cat(sprintf("  %-25s  Chi2=%7.2f  df=%2.0f  p=%.4f  CFI=%.3f  TLI=%.3f  RMSEA=%.3f [%.3f, %.3f]  SRMR=%.3f\n",
              spec$name,
              fi["chisq"], fi["df"], fi["pvalue"],
              fi["cfi"], fi["tli"],
              fi["rmsea"], fi["rmsea.ci.lower"], fi["rmsea.ci.upper"],
              fi["srmr"]))
}

# Fit comparison table
cat("\n--- Fit Comparison Table ---\n\n")
comparison <- data.frame(
  Model = sapply(model_specs, `[[`, "name"),
  ChiSq = sapply(sem_fit_idx, function(x) x["chisq"]),
  df    = sapply(sem_fit_idx, function(x) x["df"]),
  CFI   = sapply(sem_fit_idx, function(x) x["cfi"]),
  TLI   = sapply(sem_fit_idx, function(x) x["tli"]),
  RMSEA = sapply(sem_fit_idx, function(x) x["rmsea"]),
  SRMR  = sapply(sem_fit_idx, function(x) x["srmr"]),
  row.names = NULL
)
print(comparison, digits = 3)

# Key structural paths
cat("\n--- Key Structural Paths ---\n\n")

key_paths <- all_struct[all_struct$rhs %in%
  c("depriv_index", "neg_affect_sum", "income_low") &
  all_struct$lhs %in%
  c("neg_affect_sum", "INTTRUST", "trust_president", "trust_rada", "trust_local", "PROSOC"), ]

print(key_paths[, c("model", "lhs", "rhs", "est.std", "se", "pvalue")], row.names = FALSE)

# Indirect effects
cat("\n--- Indirect Effects (Depriv -> NegAff -> DV) ---\n\n")
print(all_indirect[, c("model", "lhs", "est.std", "se", "pvalue")], row.names = FALSE)

# R-squared for DVs
cat("\n--- R-squared for DVs ---\n\n")
dv_r2 <- all_r2[all_r2$variable %in%
  c("neg_affect_sum", "INTTRUST", "trust_president", "trust_rada", "trust_local", "PROSOC"), ]
print(dv_r2, row.names = FALSE)

# ============================================================
# STEP 6. VISUALIZATIONS — SUMMARY CHARTS
# ============================================================

cat("\n\n========================================\n")
cat("STEP 6: VISUALIZATIONS — SUMMARY CHARTS\n")
cat("========================================\n\n")

theme_trust <- theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 14, color = "gray40"),
    plot.margin   = margin(10, 15, 10, 15)
  )

# --- 6a. Theoretical DAG diagram (ggplot2) -------------------

cat("Generating theoretical model DAG...\n")

nodes <- data.frame(
  name = c("Infrastructure\nDeprivation", "Negative\nAffect",
           "Interpersonal\nTrust", "Trust:\nPresident", "Trust:\nRada",
           "Trust:\nLocal Gov", "Prosocial\nBehavior",
           "SES", "Age", "Sex"),
  x    = c(1,    3.5,   7,    7,    7,    7,    7,      0.2,  0.2,  0.2),
  y    = c(3,    3,     6,    4.5,  3,    1.5,  0,      6,    0,   -1.2),
  type = c("iv", "mediator", "dv", "dv", "dv", "dv", "dv",
           "control", "control", "control")
)

dv_ys <- c(6, 4.5, 3, 1.5, 0)
edges <- data.frame(
  from_x = c(1,     rep(3.5, 5),          rep(1, 5),
             rep(0.2, 5),
             rep(0.2, 5),
             rep(0.2, 5)),
  from_y = c(3,     rep(3, 5),            rep(3, 5),
             rep(6, 5),
             rep(0, 5),
             rep(-1.2, 5)),
  to_x   = c(3.5,   rep(7, 5),           rep(7, 5),
             c(3.5, rep(7, 4)),
             rep(7, 5),
             rep(7, 5)),
  to_y   = c(3,     dv_ys,               dv_ys,
             c(3, dv_ys[-1]),
             dv_ys,
             dv_ys),
  label  = c("a",   rep("b", 5),         rep("c'", 5),
             rep("", 5),
             rep("", 5),
             rep("", 5)),
  type   = c("mediation", rep("mediation", 5), rep("direct", 5),
             rep("control", 5),
             rep("control", 5),
             rep("control", 5))
)
# Fix SES control: first goes to mediator (y=3), rest to DVs
edges$to_y[edges$type == "control" & edges$from_y == 6] <- c(3, 6, 4.5, 1.5, 0)

# Shorten arrows so they don't overlap with node labels
shorten <- 0.55
edges$dx  <- edges$to_x - edges$from_x
edges$dy  <- edges$to_y - edges$from_y
edges$len <- sqrt(edges$dx^2 + edges$dy^2)
edges$ux  <- edges$dx / edges$len
edges$uy  <- edges$dy / edges$len
edges$x_start <- edges$from_x + shorten * edges$ux
edges$y_start <- edges$from_y + shorten * edges$uy
edges$x_end   <- edges$to_x   - shorten * edges$ux
edges$y_end   <- edges$to_y   - shorten * edges$uy

# Label positions — offset to avoid overlap
edges$lab_x <- (edges$from_x + edges$to_x) / 2
edges$lab_y <- (edges$from_y + edges$to_y) / 2
edges$lab_y[edges$label == "a"]  <- edges$lab_y[edges$label == "a"]  + 0.3
edges$lab_x[edges$label == "b"]  <- edges$lab_x[edges$label == "b"]  + 0.15
edges$lab_x[edges$label == "c'"] <- edges$lab_x[edges$label == "c'"] - 0.15

p_dag <- ggplot() +
  geom_segment(data = edges[edges$type == "control", ],
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "gray75", linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = edges[edges$type == "direct", ],
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "#d95f0e", linewidth = 0.9, linetype = "longdash") +
  geom_segment(data = edges[edges$type == "mediation", ],
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "#2c7fb8", linewidth = 1.1) +
  geom_label(data = edges[edges$label != "", ],
             aes(x = lab_x, y = lab_y, label = label),
             size = 4, fill = "white", label.size = 0,
             fontface = "bold.italic", color = "gray30") +
  geom_label(data = nodes[nodes$type == "iv", ],
             aes(x = x, y = y, label = name),
             size = 4.2, fontface = "bold",
             fill = "#e0f3f8", color = "#2c7fb8",
             label.padding = unit(0.5, "lines"), linewidth = 0.8) +
  geom_label(data = nodes[nodes$type == "mediator", ],
             aes(x = x, y = y, label = name),
             size = 4.2, fontface = "bold",
             fill = "#deebf7", color = "#3182bd",
             label.padding = unit(0.5, "lines"), linewidth = 0.8) +
  geom_label(data = nodes[nodes$type == "dv", ],
             aes(x = x, y = y, label = name),
             size = 3.8, fontface = "bold",
             fill = "#fee8c8", color = "#d95f0e",
             label.padding = unit(0.5, "lines"), linewidth = 0.8) +
  geom_label(data = nodes[nodes$type == "control", ],
             aes(x = x, y = y, label = name),
             size = 3, fontface = "plain",
             fill = "gray93", color = "gray45",
             label.padding = unit(0.3, "lines"), linewidth = 0.3) +
  annotate("segment", x = 0.5, xend = 1.3, y = -1.3, yend = -1.3,
           color = "#2c7fb8", linewidth = 1,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = 1.5, y = -1.3, label = "Mediation path (a \u00d7 b)",
           size = 3, hjust = 0, color = "gray30") +
  annotate("segment", x = 3.8, xend = 4.6, y = -1.3, yend = -1.3,
           color = "#d95f0e", linewidth = 0.9, linetype = "longdash",
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = 4.8, y = -1.3, label = "Direct path (c')",
           size = 3, hjust = 0, color = "gray30") +
  coord_cartesian(xlim = c(-0.8, 8.5), ylim = c(-2.5, 7.5)) +
  labs(title = "Theoretical Model",
       subtitle = "Infrastructure Deprivation \u2192 Negative Affect \u2192 Social Outcomes") +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5),
    plot.margin   = margin(15, 15, 15, 15)
  )

cairo_pdf("Trust_Kyiv_DAG.pdf", width = 12, height = 8)
print(p_dag)
dev.off()
cat("  Saved: Trust_Kyiv_DAG.pdf\n")

# --- 6b. Factor loadings bar chart ---------------------------

cat("Generating factor loadings chart...\n")

loadings_df <- loadings_only[, c("lhs", "rhs", "est.std")]
names(loadings_df) <- c("Factor", "Indicator", "Loading")

loadings_df$Factor_label <- factor(loadings_df$Factor,
  levels = factors,
  labels = factor_labels
)

loadings_df$Indicator <- factor(loadings_df$Indicator,
  levels = rev(loadings_df$Indicator))

p_loadings <- ggplot(loadings_df,
                     aes(x = Indicator, y = Loading, fill = Factor_label)) +
  geom_col(width = 0.7, show.legend = TRUE) +
  geom_text(aes(label = sprintf("%.2f", Loading)),
            hjust = -0.1, size = 3.2) +
  coord_flip() +
  facet_wrap(~ Factor_label, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#7fcdbb", "#e6550d")) +
  labs(title = "CFA Factor Loadings (Standardized)",
       x = NULL, y = "Standardized Loading") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text  = element_text(face = "bold", size = 11),
    plot.title  = element_text(face = "bold", size = 13)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     limits = c(0, 1.1))

cairo_pdf("Trust_Kyiv_CFA_loadings.pdf", width = 10, height = 8)
print(p_loadings)
dev.off()
cat("  Saved: Trust_Kyiv_CFA_loadings.pdf\n")

# --- 6c. Forest plots: one per model -------------------------

cat("Generating structural coefficients forest plots...\n")

forest_all <- key_paths[, c("model", "lhs", "rhs", "est.std", "se", "pvalue")]
names(forest_all) <- c("Model", "DV", "IV", "Beta", "SE", "p")

iv_labels <- c(
  "depriv_index"   = "Deprivation",
  "neg_affect_sum" = "Neg. Affect",
  "income_low"     = "SES (income)"
)

forest_all$path  <- paste(iv_labels[forest_all$IV], "\u2192",
                          ifelse(forest_all$DV == "neg_affect_sum", "Neg. Affect", "DV"))
forest_all$ci_lo <- forest_all$Beta - 1.96 * forest_all$SE
forest_all$ci_hi <- forest_all$Beta + 1.96 * forest_all$SE
forest_all$sig   <- ifelse(forest_all$p < 0.001, "***",
                   ifelse(forest_all$p < 0.01,  "**",
                   ifelse(forest_all$p < 0.05,  "*", "ns")))

model_titles <- c(
  "Interpersonal Trust" = "Model 1: Interpersonal Trust",
  "Trust: President"    = "Model 2: Trust in President",
  "Trust: Rada"         = "Model 3: Trust in Rada",
  "Trust: Local Gov"    = "Model 4: Trust in Local Government",
  "Prosocial Behavior"  = "Model 5: Prosocial Behavior"
)

model_files <- c(
  "Interpersonal Trust" = "Trust_Kyiv_SEM_forest_intrust.pdf",
  "Trust: President"    = "Trust_Kyiv_SEM_forest_president.pdf",
  "Trust: Rada"         = "Trust_Kyiv_SEM_forest_rada.pdf",
  "Trust: Local Gov"    = "Trust_Kyiv_SEM_forest_local.pdf",
  "Prosocial Behavior"  = "Trust_Kyiv_SEM_forest_prosoc.pdf"
)

for (mod_name in names(model_titles)) {
  fdf <- forest_all[forest_all$Model == mod_name, ]
  fdf$path <- factor(fdf$path, levels = rev(unique(fdf$path)))

  p_f <- ggplot(fdf, aes(x = path, y = Beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi, color = Beta > 0),
                    size = 0.9, linewidth = 0.9) +
    geom_text(aes(label = sprintf("\u03b2 = %.2f%s", Beta, sig)),
              hjust = -0.2, size = 3.5) +
    coord_flip() +
    scale_color_manual(values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
                       guide = "none") +
    labs(title = model_titles[mod_name],
         subtitle = "Standardized estimates with 95% CI | * p<.05, ** p<.01, *** p<.001",
         x = NULL, y = "Standardized Beta") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.35)))

  cairo_pdf(model_files[mod_name], width = 10, height = 5)
  print(p_f)
  dev.off()
  cat("  Saved:", model_files[mod_name], "\n")
}

# --- 6d. Indirect effects plot -------------------------------

cat("Generating indirect effects plot...\n")

indirect_df <- all_indirect[, c("model", "lhs", "est.std", "se", "pvalue")]
names(indirect_df) <- c("Model", "Effect", "Beta", "SE", "p")

indirect_df$Label <- paste("Depriv \u2192 NegAff \u2192", indirect_df$Model)
indirect_df$ci_lo <- indirect_df$Beta - 1.96 * indirect_df$SE
indirect_df$ci_hi <- indirect_df$Beta + 1.96 * indirect_df$SE
indirect_df$sig   <- ifelse(indirect_df$p < 0.001, "***",
                    ifelse(indirect_df$p < 0.01,  "**",
                    ifelse(indirect_df$p < 0.05,  "*", "ns")))

indirect_df$Label <- factor(indirect_df$Label, levels = rev(indirect_df$Label))

p_indirect <- ggplot(indirect_df, aes(x = Label, y = Beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi, color = Beta > 0),
                  size = 0.7, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.3f%s", Beta, sig)),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_color_manual(values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
                     guide = "none") +
  labs(title = "Indirect (Mediation) Effects Across Models",
       subtitle = "Standardized estimates with 95% CI | Depriv \u2192 NegAff \u2192 DV",
       x = NULL, y = "Standardized Indirect Effect") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "gray40")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.35)))

cairo_pdf("Trust_Kyiv_SEM_indirect.pdf", width = 10, height = 5)
print(p_indirect)
dev.off()
cat("  Saved: Trust_Kyiv_SEM_indirect.pdf\n")

# ============================================================
# STEP 7. DEPRIVATION x TRUST PLOTS
# ============================================================

cat("\n\n========================================\n")
cat("STEP 7: DEPRIVATION x TRUST PLOTS\n")
cat("========================================\n\n")

dv_items <- c("trust_district", "trust_city",
              "trust_ukraine", "trust_diaspora",
              "trust_president", "trust_rada", "trust_local",
              "help_friend", "help_lost", "help_sick", "help_stranger")

dv_labels <- c(
  "trust_district"  = "Мешк. району (0-10)",
  "trust_city"      = "Мешк. міста (0-10)",
  "trust_ukraine"   = "Мешк. України (0-10)",
  "trust_diaspora"  = "Діаспора (0-10)",
  "trust_president" = "Президент (0-10)",
  "trust_rada"      = "Верх. Рада (0-10)",
  "trust_local"     = "Місц. влада (0-10)",
  "help_friend"     = "Допомога знайомому (1-7)",
  "help_lost"       = "Допомога знайти річ (1-7)",
  "help_sick"       = "Догляд за хворим (1-7)",
  "help_stranger"   = "Допомога незнайомцю (1-7)"
)

dv_ymin <- c(rep(0, 7), rep(1, 4))
dv_ymax <- c(rep(10, 7), rep(7, 4))
names(dv_ymin) <- dv_items
names(dv_ymax) <- dv_items

# Helper: Spearman rho + p-value annotation per facet
compute_spearman <- function(data, x_var, dv_var) {
  data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[dv_var]])) %>%
    summarise(
      rho = cor(.data[[x_var]], .data[[dv_var]], method = "spearman"),
      p   = cor.test(.data[[x_var]], .data[[dv_var]], method = "spearman",
                     exact = FALSE)$p.value
    ) %>%
    mutate(
      sig   = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
              ifelse(p < 0.05, "*", "ns"))),
      label = sprintf("rho = %.2f%s", rho, sig)
    )
}

# Helper: build Spearman annotation data frame for a given x variable
build_spearman_df <- function(data, x_var) {
  do.call(rbind, lapply(dv_items, function(dv) {
    res <- compute_spearman(data, x_var, dv)
    res$dv_var   <- dv
    res$dv_label <- dv_labels[dv]
    res$y_pos    <- dv_ymax[dv] * 0.98
    res
  }))
}

# Helper: reshape to long format for a given x variable
to_long <- function(data, x_var) {
  data %>%
    select(all_of(x_var), all_of(dv_items)) %>%
    pivot_longer(cols = all_of(dv_items),
                 names_to = "dv_var",
                 values_to = "dv_value") %>%
    filter(!is.na(.data[[x_var]]), !is.na(dv_value)) %>%
    mutate(dv_label = dv_labels[dv_var])
}

# Helper: compute group means
group_means <- function(long_df, x_var) {
  long_df %>%
    group_by(comp_val = .data[[x_var]], dv_var, dv_label) %>%
    summarise(
      mean_dv = mean(dv_value, na.rm = TRUE),
      se_dv   = sd(dv_value, na.rm = TRUE) / sqrt(n()),
      n       = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 5)
}

# Helper: save plot to both PDF and PNG
save_plot <- function(p, file_base, width, height) {
  cairo_pdf(paste0(file_base, ".pdf"), width = width, height = height)
  print(p)
  dev.off()
  png(paste0(file_base, ".png"), width = width, height = height, units = "in", res = 500)
  print(p)
  dev.off()
  cat("  Saved:", paste0(file_base, ".pdf / .png\n"))
}

# --- 7a. Deprivation INDEX vs each DV ------------------------

cat("Generating deprivation index vs DV plots...\n")

spearman_idx <- build_spearman_df(model_df, "depriv_index")

dv_long <- to_long(model_df, "depriv_index")

depriv_means <- group_means(dv_long, "depriv_index")

depriv_means$dv_label  <- factor(depriv_means$dv_label,  levels = dv_labels)
spearman_idx$dv_label  <- factor(spearman_idx$dv_label,  levels = dv_labels)
spearman_idx$face      <- ifelse(spearman_idx$p < 0.05, "bold.italic", "italic")

scale_bounds_7a <- rbind(
  data.frame(depriv_index = 0, mean_dv = dv_ymin[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels)),
  data.frame(depriv_index = 0, mean_dv = dv_ymax[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels))
)

p_depriv_dv <- ggplot(depriv_means, aes(x = depriv_index, y = mean_dv)) +
  geom_point(size = 2.5, color = "#2c7fb8") +
  geom_line(color = "#2c7fb8", linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean_dv - se_dv, ymax = mean_dv + se_dv),
                width = 0.2, color = "#2c7fb8", alpha = 0.6) +
  geom_blank(data = scale_bounds_7a, aes(x = depriv_index, y = mean_dv)) +
  geom_text(aes(label = n), vjust = -1.5, size = 3.5, color = "gray50") +
  geom_text(data = spearman_idx,
            aes(x = Inf, y = y_pos, label = label, fontface = face),
            hjust = 1.1, vjust = 1, size = 5.5,
            color = "gray30", inherit.aes = FALSE) +
  facet_wrap(~ dv_label, ncol = 4, scales = "free") +
  labs(
    title    = "Індекс інфраструктурної депривації та залежні змінні",
    subtitle = "Середнє за рівнем депривації | SE | N у групі | Spearman rho | * p<.05, ** p<.01, *** p<.001",
    x        = "Індекс депривації (вищий = гірше)",
    y        = NULL
  ) +
  theme_trust +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text  = element_text(face = "bold", size = 13)
  )

save_plot(p_depriv_dv, "Trust_Kyiv_deprivation_trust", width = 16, height = 12)

# --- 7b. Components of deprivation vs DVs -------------------

cat("Generating deprivation components vs DV plots...\n")

depriv_components <- c("deprivation", "temp_cold")
depriv_comp_labels <- c(
  "deprivation" = "Втрата комунальних послуг (0=всі є, 4=нічого немає)",
  "temp_cold"   = "Холодна температура (0=22\u00b0C+, 5=<0\u00b0C)"
)

scale_bounds_comp <- rbind(
  data.frame(comp_val = 0, mean_dv = dv_ymin[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels)),
  data.frame(comp_val = 0, mean_dv = dv_ymax[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels))
)

# Build all component plots (reused for both PDF and PNG)
component_plots <- lapply(depriv_components, function(comp) {

  spearman_comp <- build_spearman_df(model_df, comp)
  comp_long     <- to_long(model_df, comp)
  comp_means    <- group_means(comp_long, comp)

  comp_means$dv_label    <- factor(comp_means$dv_label,    levels = dv_labels)
  spearman_comp$dv_label <- factor(spearman_comp$dv_label, levels = dv_labels)
  spearman_comp$face     <- ifelse(spearman_comp$p < 0.05, "bold.italic", "italic")

  ggplot(comp_means, aes(x = comp_val, y = mean_dv)) +
    geom_point(size = 2.5, color = "#d95f0e") +
    geom_line(color = "#d95f0e", linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean_dv - se_dv, ymax = mean_dv + se_dv),
                  width = 0.15, color = "#d95f0e", alpha = 0.6) +
    geom_blank(data = scale_bounds_comp, aes(x = comp_val, y = mean_dv)) +
    geom_text(aes(label = n), vjust = -1.5, size = 3.5, color = "gray50") +
    geom_text(data = spearman_comp,
              aes(x = Inf, y = y_pos, label = label, fontface = face),
              hjust = 1.1, vjust = 1, size = 5.5,
              color = "gray30", inherit.aes = FALSE) +
    facet_wrap(~ dv_label, ncol = 4, scales = "free") +
    labs(
      title    = depriv_comp_labels[comp],
      subtitle = paste0("Середнє за рівнем | Spearman rho | N = ",
                        sum(!is.na(model_df[[comp]]))),
      x = NULL, y = NULL
    ) +
    theme_trust +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text  = element_text(face = "bold", size = 13)
    )
})

cairo_pdf("Trust_Kyiv_deprivation_components_trust.pdf", width = 16, height = 12)
for (p in component_plots) print(p)
dev.off()

png("Trust_Kyiv_deprivation_components_trust.png", width = 16, height = 12, units = "in", res = 500)
for (p in component_plots) print(p)
dev.off()

cat("  Saved: Trust_Kyiv_deprivation_components_trust.pdf / .png\n")

# --- 7c. Neg Affect INDEX vs each DV -------------------------

cat("Generating neg affect index vs DV plots...\n")

spearman_neg <- build_spearman_df(model_df, "neg_affect_sum")

neg_long  <- to_long(model_df, "neg_affect_sum")
neg_means <- group_means(neg_long, "neg_affect_sum")

neg_means$dv_label    <- factor(neg_means$dv_label,    levels = dv_labels)
spearman_neg$dv_label <- factor(spearman_neg$dv_label, levels = dv_labels)
spearman_neg$face     <- ifelse(spearman_neg$p < 0.05, "bold.italic", "italic")

scale_bounds_7c <- rbind(
  data.frame(neg_affect_sum = 0, mean_dv = dv_ymin[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels)),
  data.frame(neg_affect_sum = 0, mean_dv = dv_ymax[names(dv_labels)],
             dv_label = factor(dv_labels, levels = dv_labels))
)

p_neg_dv <- ggplot(neg_means, aes(x = neg_affect_sum, y = mean_dv)) +
  geom_point(size = 2.5, color = "#7570b3") +
  geom_line(color = "#7570b3", linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean_dv - se_dv, ymax = mean_dv + se_dv),
                width = 0.2, color = "#7570b3", alpha = 0.6) +
  geom_blank(data = scale_bounds_7c, aes(x = neg_affect_sum, y = mean_dv)) +
  geom_text(aes(label = n), vjust = -1.5, size = 3.5, color = "gray50") +
  geom_text(data = spearman_neg,
            aes(x = Inf, y = y_pos, label = label, fontface = face),
            hjust = 1.1, vjust = 1, size = 5.5,
            color = "gray30", inherit.aes = FALSE) +
  facet_wrap(~ dv_label, ncol = 4, scales = "free") +
  labs(
    title    = "\u0406\u043d\u0434\u0435\u043a\u0441 \u043d\u0435\u0433\u0430\u0442\u0438\u0432\u043d\u043e\u0433\u043e \u0430\u0444\u0435\u043a\u0442\u0443 \u0442\u0430 \u0437\u0430\u043b\u0435\u0436\u043d\u0456 \u0437\u043c\u0456\u043d\u043d\u0456",
    subtitle = "\u0421\u0435\u0440\u0435\u0434\u043d\u0454 \u0437\u0430 \u0440\u0456\u0432\u043d\u0435\u043c \u043d\u0435\u0433. \u0430\u0444\u0435\u043a\u0442\u0443 | SE | N \u0443 \u0433\u0440\u0443\u043f\u0456 | Spearman rho | * p<.05, ** p<.01, *** p<.001",
    x        = "\u0406\u043d\u0434\u0435\u043a\u0441 \u043d\u0435\u0433\u0430\u0442\u0438\u0432\u043d\u043e\u0433\u043e \u0430\u0444\u0435\u043a\u0442\u0443 (0\u20134, \u0432\u0438\u0449\u0438\u0439 = \u0433\u0456\u0440\u0448\u0435)",
    y        = NULL
  ) +
  theme_trust +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text  = element_text(face = "bold", size = 13)
  )

save_plot(p_neg_dv, "Trust_Kyiv_negaffect_trust", width = 16, height = 12)

# ============================================================
# STEP 8. SEM PATH DIAGRAMS (ggplot2)
# ============================================================

cat("\n\n========================================\n")
cat("STEP 8: SEM PATH DIAGRAMS (ggplot2)\n")
cat("========================================\n\n")

# --- Shared helpers ------------------------------------------

fmt_coef <- function(est, pval) {
  star <- ifelse(pval < 0.001, "***",
          ifelse(pval < 0.01,  "**",
          ifelse(pval < 0.05,  "*", "")))
  sprintf("%.2f%s", est, star)
}

# Shorten an arrow segment so it doesn't overlap node labels
make_edge <- function(fx, fy, tx, ty, lbl, etype,
                      lbl_nudge_x = 0, lbl_nudge_y = 0.3,
                      shorten_amt = 0.6) {
  dx <- tx - fx; dy <- ty - fy
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len; uy <- dy / len
  data.frame(
    x_start = fx + shorten_amt * ux,
    y_start = fy + shorten_amt * uy,
    x_end   = tx - shorten_amt * ux,
    y_end   = ty - shorten_amt * uy,
    lab_x   = (fx + tx) / 2 + lbl_nudge_x,
    lab_y   = (fy + ty) / 2 + lbl_nudge_y,
    label   = lbl,
    etype   = etype,
    stringsAsFactors = FALSE
  )
}

# --- Draw SEM path diagram for a LATENT DV model -------------

draw_sem_latent <- function(fit, model_title, dv_label, indicator_names,
                            indicator_labels, file_prefix) {

  std    <- standardizedSolution(fit)
  struct <- std[std$op == "~", ]

  # Factor loadings
  loadings    <- std[std$op == "=~" & std$lhs == dv_label, ]
  load_labels <- sapply(seq_len(nrow(loadings)), function(i)
    fmt_coef(loadings$est.std[i], loadings$pvalue[i]))

  # Key structural paths
  a_row <- struct[struct$lhs == "neg_affect_sum" & struct$rhs == "depriv_index", ]
  b_row <- struct[struct$lhs == dv_label         & struct$rhs == "neg_affect_sum", ]
  c_row <- struct[struct$lhs == dv_label         & struct$rhs == "depriv_index", ]

  a_lab <- paste0("a = ", fmt_coef(a_row$est.std, a_row$pvalue))
  b_lab <- paste0("b = ", fmt_coef(b_row$est.std, b_row$pvalue))
  c_lab <- paste0("c = ", fmt_coef(c_row$est.std, c_row$pvalue))

  # Control paths
  ctrl_vars      <- c("income_low", "age", "female")
  ctrl_labels_ua <- c("низький\nдохід", "вік", "жіноча\nстать")
  ctrl_coefs <- sapply(ctrl_vars, function(cv) {
    r <- struct[struct$lhs == dv_label & struct$rhs == cv, ]
    if (nrow(r) > 0) fmt_coef(r$est.std, r$pvalue) else ""
  })

  inc_neg_row <- struct[struct$lhs == "neg_affect_sum" & struct$rhs == "income_low", ]
  inc_neg_lab <- if (nrow(inc_neg_row) > 0) fmt_coef(inc_neg_row$est.std, inc_neg_row$pvalue) else ""

  ind     <- std[std$op == ":=", ]
  ind_lab <- if (nrow(ind) > 0) fmt_coef(ind$est.std, ind$pvalue) else ""

  # Node positions
  n_ind      <- length(indicator_names)
  ind_spacing <- 2.2
  ind_total  <- (n_ind - 1) * ind_spacing

  iv_x  <- 1;   iv_y  <- 1.5
  med_x <- 4;   med_y <- 3.5
  dv_x  <- 8;   dv_y  <- 3.5
  ctrl_x <- 13
  ctrl_ys <- c(5, 3.5, 2)
  ind_xs <- seq(dv_x - ind_total / 2, dv_x + ind_total / 2, length.out = n_ind)
  ind_y  <- 7

  nodes_df <- data.frame(
    x     = c(iv_x, med_x, dv_x, ind_xs, ctrl_x, ctrl_x, ctrl_x),
    y     = c(iv_y, med_y, dv_y, rep(ind_y, n_ind), ctrl_ys),
    label = c("Індекс\nдепривації", "Сума\nнегативних\nемоцій", dv_label,
              indicator_labels, ctrl_labels_ua),
    type  = c("iv", "mediator", "dv", rep("indicator", n_ind),
              "control", "control", "control"),
    stringsAsFactors = FALSE
  )

  edges_list <- list(
    make_edge(iv_x, iv_y, med_x, med_y, a_lab, "main",   0,   0.4),
    make_edge(med_x, med_y, dv_x, dv_y, b_lab, "main",   0,   0.35),
    make_edge(iv_x, iv_y,  dv_x, dv_y,  c_lab, "direct", 0,  -0.4),
    make_edge(ctrl_x, ctrl_ys[1], med_x, med_y, inc_neg_lab, "control", 0, 0.3)
  )
  for (i in seq_len(n_ind))
    edges_list[[length(edges_list) + 1]] <- make_edge(
      dv_x, dv_y, ind_xs[i], ind_y, load_labels[i], "loading", 0, 0.25)
  for (i in seq_along(ctrl_vars))
    edges_list[[length(edges_list) + 1]] <- make_edge(
      ctrl_x, ctrl_ys[i], dv_x, dv_y, ctrl_coefs[i], "control", 0.3, 0.2)

  edges_df <- do.call(rbind, edges_list)

  fi         <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
  r2         <- lavInspect(fit, "rsquare")
  dv_r2_val  <- if (dv_label %in% names(r2)) r2[dv_label] else NA
  fit_text   <- sprintf("CFI = %.3f | TLI = %.3f | RMSEA = %.3f | SRMR = %.3f | R\u00b2(%s) = %.3f",
                        fi["cfi"], fi["tli"], fi["rmsea"], fi["srmr"], dv_label, dv_r2_val)
  ind_text   <- paste0("Indirect (a\u00d7b) = ", ind_lab)

  p <- ggplot() +
    geom_segment(data = edges_df[edges_df$etype == "control", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                 color = "gray65", linewidth = 0.5, linetype = "dashed") +
    geom_segment(data = edges_df[edges_df$etype == "direct", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "#d95f0e", linewidth = 1, linetype = "longdash") +
    geom_segment(data = edges_df[edges_df$etype == "main", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "#2c7fb8", linewidth = 1.2) +
    geom_segment(data = edges_df[edges_df$etype == "loading", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 color = "#31a354", linewidth = 0.8) +
    geom_label(data = edges_df[edges_df$label != "", ],
               aes(x = lab_x, y = lab_y, label = label),
               size = 3.5, fill = "white", label.size = 0,
               fontface = "bold", color = "gray20") +
    geom_label(data = nodes_df[nodes_df$type == "iv", ],
               aes(x = x, y = y, label = label),
               size = 4, fontface = "bold",
               fill = "#fc8d59", color = "white",
               label.padding = unit(0.5, "lines"), linewidth = 0.8) +
    geom_label(data = nodes_df[nodes_df$type == "mediator", ],
               aes(x = x, y = y, label = label),
               size = 3.8, fontface = "bold",
               fill = "#6baed6", color = "white",
               label.padding = unit(0.5, "lines"), linewidth = 0.8) +
    annotate("point", x = dv_x, y = dv_y, size = 22,
             shape = 21, fill = "#31a354", color = "#31a354") +
    annotate("text", x = dv_x, y = dv_y, label = gsub("_", "\n", dv_label),
             size = 3.5, fontface = "bold", color = "white") +
    geom_label(data = nodes_df[nodes_df$type == "indicator", ],
               aes(x = x, y = y, label = label),
               size = 3.2, fontface = "bold",
               fill = "#74c476", color = "white",
               label.padding = unit(0.45, "lines"), linewidth = 0.6) +
    geom_label(data = nodes_df[nodes_df$type == "control", ],
               aes(x = x, y = y, label = label),
               size = 3, fontface = "plain",
               fill = "gray88", color = "gray35",
               label.padding = unit(0.35, "lines"), linewidth = 0.3) +
    annotate("text", x = 7, y =  0.2, label = fit_text,
             size = 3.2, color = "gray40", hjust = 0.5) +
    annotate("text", x = 7, y = -0.3, label = ind_text,
             size = 3.5, color = "#2c7fb8", fontface = "bold", hjust = 0.5) +
    coord_cartesian(xlim = c(-0.5, 15), ylim = c(-1, 8.5)) +
    labs(title    = model_title,
         subtitle = "Standardized estimates | * p<.05, ** p<.01, *** p<.001") +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
      plot.margin   = margin(10, 10, 10, 10)
    )

  save_plot(p, file_prefix, width = 14, height = 9)
}

# --- Draw SEM path diagram for an OBSERVED DV model ----------

draw_sem_observed <- function(fit, model_title, dv_var, dv_label_ua,
                              file_prefix) {

  std    <- standardizedSolution(fit)
  struct <- std[std$op == "~", ]

  # Key structural paths
  a_row <- struct[struct$lhs == "neg_affect_sum" & struct$rhs == "depriv_index", ]
  b_row <- struct[struct$lhs == dv_var           & struct$rhs == "neg_affect_sum", ]
  c_row <- struct[struct$lhs == dv_var           & struct$rhs == "depriv_index", ]

  a_lab <- paste0("a = ", fmt_coef(a_row$est.std, a_row$pvalue))
  b_lab <- paste0("b = ", fmt_coef(b_row$est.std, b_row$pvalue))
  c_lab <- paste0("c = ", fmt_coef(c_row$est.std, c_row$pvalue))

  # Control paths
  ctrl_vars      <- c("income_low", "age", "female")
  ctrl_labels_ua <- c("низький\nдохід", "вік", "жіноча\nстать")
  ctrl_coefs <- sapply(ctrl_vars, function(cv) {
    r <- struct[struct$lhs == dv_var & struct$rhs == cv, ]
    if (nrow(r) > 0) fmt_coef(r$est.std, r$pvalue) else ""
  })

  inc_neg_row <- struct[struct$lhs == "neg_affect_sum" & struct$rhs == "income_low", ]
  inc_neg_lab <- if (nrow(inc_neg_row) > 0) fmt_coef(inc_neg_row$est.std, inc_neg_row$pvalue) else ""

  ind     <- std[std$op == ":=", ]
  ind_lab <- if (nrow(ind) > 0) fmt_coef(ind$est.std, ind$pvalue) else ""

  # Node positions
  iv_x  <- 1;   iv_y  <- 1.5
  med_x <- 4.5; med_y <- 3.5
  dv_x  <- 9;   dv_y  <- 3.5
  ctrl_x <- 13
  ctrl_ys <- c(5, 3.5, 2)

  nodes_df <- data.frame(
    x     = c(iv_x, med_x, dv_x, ctrl_x, ctrl_x, ctrl_x),
    y     = c(iv_y, med_y, dv_y, ctrl_ys),
    label = c("Індекс\nдепривації", "Сума\nнегативних\nемоцій", dv_label_ua,
              ctrl_labels_ua),
    type  = c("iv", "mediator", "dv", "control", "control", "control"),
    stringsAsFactors = FALSE
  )

  edges_list <- list(
    make_edge(iv_x, iv_y,  med_x, med_y, a_lab, "main",   0,   0.4),
    make_edge(med_x, med_y, dv_x, dv_y,  b_lab, "main",   0,   0.35),
    make_edge(iv_x, iv_y,   dv_x, dv_y,  c_lab, "direct", 0,  -0.4),
    make_edge(ctrl_x, ctrl_ys[1], med_x, med_y, inc_neg_lab, "control", 0, 0.3)
  )
  for (i in seq_along(ctrl_vars))
    edges_list[[length(edges_list) + 1]] <- make_edge(
      ctrl_x, ctrl_ys[i], dv_x, dv_y, ctrl_coefs[i], "control", 0.3, 0.2)

  edges_df <- do.call(rbind, edges_list)

  fi        <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
  r2        <- lavInspect(fit, "rsquare")
  dv_r2_val <- if (dv_var %in% names(r2)) r2[dv_var] else NA
  fit_text  <- sprintf("CFI = %.3f | TLI = %.3f | RMSEA = %.3f | SRMR = %.3f | R\u00b2(%s) = %.3f",
                       fi["cfi"], fi["tli"], fi["rmsea"], fi["srmr"], dv_var, dv_r2_val)
  ind_text  <- paste0("Indirect (a\u00d7b) = ", ind_lab)

  p <- ggplot() +
    geom_segment(data = edges_df[edges_df$etype == "control", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                 color = "gray65", linewidth = 0.5, linetype = "dashed") +
    geom_segment(data = edges_df[edges_df$etype == "direct", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "#d95f0e", linewidth = 1, linetype = "longdash") +
    geom_segment(data = edges_df[edges_df$etype == "main", ],
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "#2c7fb8", linewidth = 1.2) +
    geom_label(data = edges_df[edges_df$label != "", ],
               aes(x = lab_x, y = lab_y, label = label),
               size = 3.8, fill = "white", label.size = 0,
               fontface = "bold", color = "gray20") +
    geom_label(data = nodes_df[nodes_df$type == "iv", ],
               aes(x = x, y = y, label = label),
               size = 4.2, fontface = "bold",
               fill = "#fc8d59", color = "white",
               label.padding = unit(0.55, "lines"), linewidth = 0.8) +
    geom_label(data = nodes_df[nodes_df$type == "mediator", ],
               aes(x = x, y = y, label = label),
               size = 4, fontface = "bold",
               fill = "#6baed6", color = "white",
               label.padding = unit(0.55, "lines"), linewidth = 0.8) +
    geom_label(data = nodes_df[nodes_df$type == "dv", ],
               aes(x = x, y = y, label = label),
               size = 4.2, fontface = "bold",
               fill = "#31a354", color = "white",
               label.padding = unit(0.6, "lines"), linewidth = 0.8) +
    geom_label(data = nodes_df[nodes_df$type == "control", ],
               aes(x = x, y = y, label = label),
               size = 3.2, fontface = "plain",
               fill = "gray88", color = "gray35",
               label.padding = unit(0.35, "lines"), linewidth = 0.3) +
    annotate("text", x = 7, y =  0.3, label = fit_text,
             size = 3.2, color = "gray40", hjust = 0.5) +
    annotate("text", x = 7, y = -0.2, label = ind_text,
             size = 3.8, color = "#2c7fb8", fontface = "bold", hjust = 0.5) +
    coord_cartesian(xlim = c(-0.5, 15), ylim = c(-1, 6)) +
    labs(title    = model_title,
         subtitle = "Standardized estimates | * p<.05, ** p<.01, *** p<.001") +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
      plot.margin   = margin(10, 10, 10, 10)
    )

  save_plot(p, file_prefix, width = 14, height = 7)
}

# --- Draw all 5 models ---------------------------------------

cat("Generating SEM path diagrams...\n")

draw_sem_latent(
  fit              = sem_fits[["INTTRUST"]],
  model_title      = "Model 1: \u041c\u0456\u0436\u043e\u0441\u043e\u0431\u0438\u0441\u0442\u0456\u0441\u043d\u0430 \u0434\u043e\u0432\u0456\u0440\u0430 (Interpersonal Trust)",
  dv_label         = "INTTRUST",
  indicator_names  = c("trust_district", "trust_city", "trust_ukraine", "trust_diaspora"),
  indicator_labels = c("Жителі\nрайону", "Жителі\nміста", "Люди в\nУкраїні", "Українці\nзакордоном"),
  file_prefix      = "Trust_Kyiv_SEM_path_inttrust"
)

draw_sem_observed(
  fit          = sem_fits[["trust_president"]],
  model_title  = "Model 2: \u0414\u043e\u0432\u0456\u0440\u0430 \u041f\u0440\u0435\u0437\u0438\u0434\u0435\u043d\u0442\u0443",
  dv_var       = "trust_president",
  dv_label_ua  = "\u0414\u043e\u0432\u0456\u0440\u0430\n\u041f\u0440\u0435\u0437\u0438\u0434\u0435\u043d\u0442\u0443",
  file_prefix  = "Trust_Kyiv_SEM_path_president"
)

draw_sem_observed(
  fit          = sem_fits[["trust_rada"]],
  model_title  = "Model 3: \u0414\u043e\u0432\u0456\u0440\u0430 \u0412\u0435\u0440\u0445\u043e\u0432\u043d\u0456\u0439 \u0420\u0430\u0434\u0456",
  dv_var       = "trust_rada",
  dv_label_ua  = "\u0414\u043e\u0432\u0456\u0440\u0430\n\u0412\u0435\u0440\u0445\u043e\u0432\u043d\u0456\u0439 \u0420\u0430\u0434\u0456",
  file_prefix  = "Trust_Kyiv_SEM_path_rada"
)

draw_sem_observed(
  fit          = sem_fits[["trust_local"]],
  model_title  = "Model 4: \u0414\u043e\u0432\u0456\u0440\u0430 \u043c\u0456\u0441\u0446\u0435\u0432\u0456\u0439 \u0432\u043b\u0430\u0434\u0456",
  dv_var       = "trust_local",
  dv_label_ua  = "\u0414\u043e\u0432\u0456\u0440\u0430\n\u043c\u0456\u0441\u0446\u0435\u0432\u0456\u0439 \u0432\u043b\u0430\u0434\u0456",
  file_prefix  = "Trust_Kyiv_SEM_path_local"
)

draw_sem_latent(
  fit              = sem_fits[["PROSOC"]],
  model_title      = "Model 5: \u041f\u0440\u043e\u0441\u043e\u0446\u0456\u0430\u043b\u044c\u043d\u0430 \u043f\u043e\u0432\u0435\u0434\u0456\u043d\u043a\u0430 (Prosocial Behavior)",
  dv_label         = "PROSOC",
  indicator_names  = c("help_friend", "help_lost", "help_sick", "help_stranger"),
  indicator_labels = c("Допомога\nзнайомому", "Допомога\nзнайти річ", "Догляд за\nхворим", "Допомога\nнезнайомцю"),
  file_prefix      = "Trust_Kyiv_SEM_path_prosoc"
)

cat("\nAll SEM path diagrams done.\n")

# ============================================================
# STEP 9. SAVE RESULTS
# ============================================================

sem_results <- list(
  cfa_fit      = cfa_fit,
  sem_fits     = sem_fits,
  fit_cfa      = fit_cfa,
  sem_fit_idx  = sem_fit_idx,
  std_loadings = loadings_only,
  all_struct   = all_struct,
  all_indirect = all_indirect,
  all_r2       = all_r2,
  model_df     = model_df
)

save(sem_results, file = "Trust_Kyiv_SEM_results.RData")
cat("\nResults saved: Trust_Kyiv_SEM_results.RData\n")

sink()
