# -------------------------------------------------------------------------
# Necessary libraries
# -------------------------------------------------------------------------
library(plm)
library(ggstatsplot)
library(ggplot2)
library(ggcharts)
library(patchwork)

# -------------------------------------------------------------------------
# Data input
# -------------------------------------------------------------------------
# Load the prepared dataset
# All analyses below use the panel dataset `trust_both_waves2`
load("Data_trust_Kyiv.RData")

# Rename the wave labels to English so plots and model summaries are in English
levels(trust_both_waves2$wave) <- c("First wave", "Second wave")

# Save the two survey wave labels once and reuse them later when the data are split into the first and second wave
wave_levels <- levels(trust_both_waves2$wave)

# -------------------------------------------------------------------------
# 0. Helper functions
# -------------------------------------------------------------------------
# Convert selected variables from long panel format to wide format
# In the output, each respondent has one row with separate columns for wave 1 and wave 2. This format is needed for paired comparisons
make_wide_data <- function(vars) {
  wave1_data <- as.data.frame(
    trust_both_waves2[trust_both_waves2$wave == wave_levels[1], c("id", vars)]
  )
  wave2_data <- as.data.frame(
    trust_both_waves2[trust_both_waves2$wave == wave_levels[2], c("id", vars)]
  )

  merge(wave1_data, wave2_data, by = "id", suffixes = c("_w1", "_w2"))
}

# Format p-values so they can be inserted directly into figure labels
format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }

  if (p_value < 0.001) {
    return("p < 0.001")
  }

  sprintf("p = %.3f", p_value)
}

# Extract the coefficients needed for plotting from a `plm` model
# The table is reshaped to match the format expected by `ggcoefstats` when it receives a precomputed coefficient table instead of a model object
extract_model_terms <- function(model, pattern) {
  coefficient_table <- as.data.frame(summary(model)$coefficients)
  confidence_intervals <- as.data.frame(confint(model))

  coefficient_table$term <- rownames(coefficient_table)
  confidence_intervals$term <- rownames(confidence_intervals)

  model_data <- merge(coefficient_table, confidence_intervals, by = "term")
  model_data <- model_data[grep(pattern, model_data$term), , drop = FALSE]

  if (nrow(model_data) == 0) {
    stop("No model terms matched the requested pattern: ", pattern)
  }

  ci_names <- names(confidence_intervals)[names(confidence_intervals) != "term"]

  data.frame(
    term = factor(model_data$term, levels = model_data$term),
    estimate = model_data$Estimate,
    std.error = model_data$`Std. Error`,
    statistic = model_data$`t-value`,
    p.value = model_data$`Pr(>|t|)`,
    conf.low = model_data[[ci_names[1]]],
    conf.high = model_data[[ci_names[2]]],
    conf.level = 0.95,
    conf.method = "Wald",
    stringsAsFactors = FALSE
  )
}

# Calculate mean values for the first and second wave separately
# The result is a two-column data frame that can be used directly in the descriptive dumbbell plots
calculate_wave_means <- function(var_range) {
  means_matrix <- sapply(
    trust_both_waves2[, var_range],
    function(x) tapply(x, trust_both_waves2$wave, mean, na.rm = TRUE)
  ) |>
    round(2) |>
    t() |>
    as.data.frame()

  colnames(means_matrix) <- c("wave1", "wave2")
  means_matrix
}

# Use one shared theme for all coefficient plots so they have the same appearance and can be compared more easily
coef_theme <- theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(size = 11, color = "black")
)

# Build a dumbbell chart that compares mean values between the two waves
# The same function is used for trust indicators and prosocial indicators
build_dumbbell_plot <- function(data, title, y_limits, y_breaks, point_size, wave2_hjust, wave1_hjust) {
  dumbbell_chart(
    data = data,
    x = category,
    y1 = wave1,
    y2 = wave2,
    sort = FALSE,
    legend = FALSE,
    line_color = data$colors,
    point_colors = c("transparent", "navy"),
    line_size = 5.5,
    point_size = point_size
  ) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    xlab("") +
    labs(title = title) +
    ylab("") +
    theme(
      axis.text.x = element_text(size = 10.5),
      plot.title = element_text(size = 13.5, color = "black", face = "bold"),
      axis.text.y = element_text(size = 11.5)
    ) +
    geom_text(
      data = data,
      aes(x = category, y = wave2, label = sprintf("%.2f", wave2)),
      hjust = wave2_hjust,
      size = 3.8,
      fontface = "bold",
      color = "navy"
    ) +
    geom_text(
      data = data,
      aes(x = category, y = wave1, label = sprintf("%.2f", wave1)),
      hjust = wave1_hjust,
      size = 3.8,
      fontface = "bold",
      color = "grey45"
    )
}

# Build a coefficient plot for models that include survey wave and `sumgoods`
# Only the coefficient for wave 2 and the `sumgoods` category coefficients are kept in the final plot
build_sumgoods_coef_plot <- function(model, title, point_colors) {
  model_terms <- extract_model_terms(model, "^wave|^as.factor\\(sumgoods\\)")
  desired_order <- c(
    "as.factor(sumgoods)4",
    "as.factor(sumgoods)3",
    "as.factor(sumgoods)2",
    "as.factor(sumgoods)1",
    "waveSecond wave"
  )
  term_labels <- c(
    "as.factor(sumgoods)4" = "Total number of services: 4",
    "as.factor(sumgoods)3" = "Total number of services: 3",
    "as.factor(sumgoods)2" = "Total number of services: 2",
    "as.factor(sumgoods)1" = "Total number of services: 1",
    "waveSecond wave" = "Wave: Second wave"
  )
  present_terms <- desired_order[desired_order %in% model_terms$term]
  model_terms <- model_terms[match(present_terms, model_terms$term), , drop = FALSE]
  model_terms$term <- factor(
    unname(term_labels[present_terms]),
    levels = unname(term_labels[present_terms])
  )

  ggcoefstats(
    model_terms,
    point.args = list(size = 3.2, color = point_colors, na.rm = TRUE)
  ) +
    xlab(title) +
    ylab("") +
    coef_theme
}

# -------------------------------------------------------------------------
# 1. Services: percentage distributions by wave and McNemar tests
# -------------------------------------------------------------------------
# For each of the five goods, calculate the percentage of respondents who reported access to that good in wave 1 and wave 2. Then run a separate
# McNemar test to evaluate whether the paired proportion changed over time
# In the Vox article this plot is built using Datawrapper 
goods_vars <- c("goods1", "goods4", "goods3", "goods2", "goods5")
goods_labels <- c(
  "Electricity",
  "Heating",
  "Cold water",
  "Hot water",
  "None of them"
)
goods_wide <- make_wide_data(goods_vars)

goods_distribution_data <- do.call(
  rbind,
  lapply(seq_along(goods_vars), function(i) {
    variable <- goods_vars[i]

    data.frame(
      good = goods_labels[i],
      wave = factor(wave_levels, levels = wave_levels),
      percent = c(
        mean(
          trust_both_waves2[trust_both_waves2$wave == wave_levels[1], variable][[1]],
          na.rm = TRUE
        ),
        mean(
          trust_both_waves2[trust_both_waves2$wave == wave_levels[2], variable][[1]],
          na.rm = TRUE
        )
      ) * 100
    )
  })
)
goods_distribution_data$good <- factor(goods_distribution_data$good, levels = goods_labels)

# Create one McNemar label for each good so the test result is shown directly on the descriptive chart
goods_mcnemar_labels <- data.frame(
  good = goods_labels,
  percent = rep(max(goods_distribution_data$percent) + 4, length(goods_vars)),
  label = sapply(seq_along(goods_vars), function(i) {
    test_table <- table(
      goods_wide[[paste0(goods_vars[i], "_w1")]],
      goods_wide[[paste0(goods_vars[i], "_w2")]]
    )

    paste0("McNemar: ", format_p_value(mcnemar.test(test_table)$p.value))
  })
)
goods_mcnemar_labels$good <- factor(goods_mcnemar_labels$good, levels = goods_labels)

ggplot(goods_distribution_data, aes(x = good, y = percent, fill = wave)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.5
  ) +
  geom_text(
    data = goods_mcnemar_labels,
    aes(x = good, y = percent, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = -2,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("grey70", "navy")) +
  scale_y_continuous(
    limits = c(0, max(goods_mcnemar_labels$percent) + 10),
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Availability of goods in the first and second waves",
    x = "",
    y = "Percent",
    fill = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# -------------------------------------------------------------------------
# 2. Total number of services: percentage distribution by wave and McNemar test
# -------------------------------------------------------------------------
# Compare the full distribution of the total number of available goods between wave 1 and wave 2
# The paired McNemar test is reported in the subtitle of the figure.
# In the Vox article this plots is built using Datawrapper
sumgoods_wide <- make_wide_data("sumgoods")
sumgoods_mcnemar <- mcnemar.test(table(sumgoods_wide$sumgoods_w1, sumgoods_wide$sumgoods_w2))

sumgoods_distribution_data <- as.data.frame(
  prop.table(table(trust_both_waves2$sumgoods, trust_both_waves2$wave), margin = 2)
)

colnames(sumgoods_distribution_data) <- c("sumgoods", "wave", "percent")
sumgoods_distribution_data$percent <- sumgoods_distribution_data$percent * 100
sumgoods_distribution_data$sumgoods <- factor(
  sumgoods_distribution_data$sumgoods,
  levels = sort(unique(sumgoods_distribution_data$sumgoods))
)

ggplot(sumgoods_distribution_data, aes(x = sumgoods, y = percent, fill = wave)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = c("grey70", "navy")) +
  scale_y_continuous(
    limits = c(0, max(sumgoods_distribution_data$percent) + 8),
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Distribution of total goods in the first and second waves",
    subtitle = paste0("McNemar chi-square test: ", format_p_value(sumgoods_mcnemar$p.value)),
    x = "Total number of services",
    y = "Percent",
    fill = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# -------------------------------------------------------------------------
# 3. Mean comparisons between waves
# -------------------------------------------------------------------------
# Compute mean levels of trust in wave 1 and wave 2 for the nine trust items
# Green color indicates increase in trust, coral - decrease, transparent - no changes
trust_means <- calculate_wave_means(12:20)
trust_means <- data.frame(
  category = c(
    "Family members",
    "Neighbors",
    "Residents of one's district",
    "Residents of the city",
    "Residents of Ukraine",
    "Ukrainians abroad",
    "President of Ukraine",
    "Verkhovna Rada",
    "Local authorities"
  ),
  colors = c(
    "coral",
    "transparent",
    "transparent",
    "seagreen3",
    "transparent",
    "transparent",
    "coral",
    "transparent",
    "coral"
  ),
  trust_means
)

# Set the display order of trust categories manually so the figure follows the intended substantive sequence instead of the default alphabetic order
trust_means$category <- factor(
  trust_means$category,
  levels = names(table(trust_means$category))[c(2, 9, 4, 8:5, 3, 1)]
)

# Compute mean levels of prosocial behavior in wave 1 and wave 2
# Green color indicates increase in prosocial behavior, coral - decrease, transparent - no changes
prosocial_means <- calculate_wave_means(21:24)
prosocial_means <- data.frame(
  category = c(
    "Supporting an acquaintance",
    "Helping a stranger look \n for a lost item",
    "Helping care for a sick person",
    "Helping a stranger with \n a small task"
  ),
  colors = c("transparent", "coral", "transparent", "coral"),
  prosocial_means
)

prosocial_means$category <- factor(
  prosocial_means$category,
  levels = names(table(prosocial_means$category))[c(2, 3, 1, 4)]
)

# Combine the two descriptive figures into one patchwork figure
# Grey labels represent mean values from the first wave, while blue labels represent those from the second
# Lower values are positioned to the left, and higher values to the right 
# A green rectangle indicates an increase in the indicator during the second wave, whereas a coral rectangle indicates a decrease. The absence of a rectangle, with only a blue dot (representing the second-wave mean), identifies no significant changes
build_dumbbell_plot(
  data = trust_means,
  title = "Trust",
  y_limits = c(0, 10),
  y_breaks = seq(0, 10, 1),
  point_size = 2.1,
  wave2_hjust = c(rep(1.5, 2), -0.5, -0.5, -0.5, -0.5, rep(1.5, 3)),
  wave1_hjust = c(-0.3, -0.3, 1.4, 1.3, 1.4, 1.5, -0.4, -0.3, -0.3)
) /
  build_dumbbell_plot(
    data = prosocial_means,
    title = "Prosocial behavior",
    y_limits = c(1, 7),
    y_breaks = seq(1, 7, 1),
    point_size = 3.1,
    wave2_hjust = 1.5,
    wave1_hjust = -0.4
  )

# -------------------------------------------------------------------------
# 4. Respondent fixed-effects models: wave only
# -------------------------------------------------------------------------
# Estimate within-respondent change using `plm` with individual fixed effects
wave_only_models <- list(
  plm(trfam ~ wave, data = trust_both_waves2, index = c("id", "wave"), model = "within", effect = "individual"),
  plm(trpres ~ wave, data = trust_both_waves2, index = c("id", "wave"), model = "within", effect = "individual"),
  plm(trloc ~ wave, data = trust_both_waves2, index = c("id", "wave"), model = "within", effect = "individual"),
  plm(soc2 ~ wave, data = trust_both_waves2, index = c("id", "wave"), model = "within", effect = "individual"),
  plm(soc4 ~ wave, data = trust_both_waves2, index = c("id", "wave"), model = "within", effect = "individual")
)

wave_effects_table <- do.call(
  rbind,
  lapply(wave_only_models, extract_model_terms, pattern = "^wave")
)

# Replace technical model labels with readable labels for the final figure
wave_effects_table[[1]] <- c(
  "Trust in family members",
  "Trust in the President of Ukraine",
  "Trust in local authorities",
  "Helping a stranger look \n for a lost item",
  "Helping a stranger with \n a small task"
) |>
  as.factor()

wave_effects_table <- wave_effects_table[5:1, ]

# Red color indicates decrease in trust/prosocial behavior, dark gray - no changes
ggcoefstats(
  wave_effects_table,
  point.args = list(
    size = 3.2,
    color = c(rep("firebrick2", 3), "grey20", "firebrick2"),
    na.rm = TRUE
  )
) +
  xlab("") +
  ylab("") +
  coef_theme +
  scale_x_continuous(limits = c(-0.8, 0.2), breaks = seq(-1, 0.2, 0.2))

# -------------------------------------------------------------------------
# 5. Respondent fixed-effects models: wave + total number of services
# -------------------------------------------------------------------------
# Add the total number of available goods in the apartment (`sumgoods`) to
# the fixed-effects models to test whether this variable is associated with
# trust and prosocial behavior after controlling for respondent constants
family_trust_model <- plm(
  trfam ~ wave + as.factor(sumgoods),
  data = trust_both_waves2,
  index = c("id", "wave"),
  model = "within",
  effect = "individual"
)
local_trust_model <- plm(
  trloc ~ wave + as.factor(sumgoods),
  data = trust_both_waves2,
  index = c("id", "wave"),
  model = "within",
  effect = "individual"
)
lost_item_help_model <- plm(
  soc2 ~ wave + as.factor(sumgoods),
  data = trust_both_waves2,
  index = c("id", "wave"),
  model = "within",
  effect = "individual"
)
small_task_help_model <- plm(
  soc4 ~ wave + as.factor(sumgoods),
  data = trust_both_waves2,
  index = c("id", "wave"),
  model = "within",
  effect = "individual"
)

# Keep model summaries in the script output for anyone who reruns the file and wants to inspect the full regression results in the console
summary(family_trust_model)
summary(local_trust_model)
summary(lost_item_help_model)
summary(small_task_help_model)

# Combine plots using patchwork library
# Red color indicates decrease in trust/prosocial behavior, green - increase, dark gray - no changes
(
  build_sumgoods_coef_plot(
    model = family_trust_model,
    title = "Trust in family members",
    point_colors = c("springgreen3", "grey20", "springgreen3", "springgreen3", "firebrick2")
  ) +
    build_sumgoods_coef_plot(
      model = local_trust_model,
      title = "Trust in local authorities",
      point_colors = c(rep("grey20", 4), "firebrick2")
    )
) /
  (
    build_sumgoods_coef_plot(
      model = lost_item_help_model,
      title = "Helping a stranger look \n for a lost item",
      point_colors = c(rep("grey20", 4), "firebrick2")
    ) +
      build_sumgoods_coef_plot(
        model = small_task_help_model,
        title = "Helping a stranger with \n a small task",
        point_colors = "grey20"
      )
  )
