# ============================================================
# File: 03_plots.R
# Purpose:
#   (A) IRT plots: Item Characteristic Curves (ICC) for 2PL models
#   (B) Modeling plots: ggpredict() from lme1 for key terms
# Notes:
#   - Requires objects: m2pl_e4, m2pl_rec (ltm 2PL models), lme1 (mixed model)
#   - Exports Figure_2.png and Figure_3.jpeg
# Author: Diego Rivera
# Date: 2025-10-26
# ============================================================

# --- 0) Packages and data -------------------------------------------------------------
suppressPackageStartupMessages({
  library(ltm)        # ICC extraction (coef, 2PL)
  library(ggplot2)    # plotting
  library(dplyr)      # bind, mutate
  library(ggeffects)  # ggpredict()
  library(cowplot)    # plot_grid
})

# --- 0) Load model objects ---------------------------------------------------
# Ensure required objects are available (lme1 and irt_results)
# Both were saved previously in the analysis pipeline.

# Load mixed-effects model (LMM)
lme1 <- readRDS("outputs/lmm_selected_model.rds")

# Load IRT results (Rasch + 2PL fits for each Trial)
irt_results <- readRDS("outputs/irt_results.rds")

# --- 1) Output dir -----------------------------------------------------------
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# (A) IRT PLOTS — ICC for Trial 4 and Recognition
# ============================================================

# Helper: extract tidy ICC curves from a 2PL model
extract_icc <- function(model, label, theta_range = c(-4, 4), n_points = 400) {
  cf <- coef(model)  # expected columns: "Dffclt", "Dscrmn"
  if (!all(c("Dffclt", "Dscrmn") %in% colnames(cf))) {
    stop("coef(model) must contain columns 'Dffclt' and 'Dscrmn'.")
  }
  theta <- seq(theta_range[1], theta_range[2], length.out = n_points)
  
  # Build long data.frame with ICC for each item
  out <- lapply(seq_len(nrow(cf)), function(i) {
    a <- cf[i, "Dscrmn"]
    b <- cf[i, "Dffclt"]
    p <- plogis(a * (theta - b))
    data.frame(Item = paste0("Item ", i), Theta = theta, Prob = p, stringsAsFactors = FALSE)
  })
  df <- dplyr::bind_rows(out)
  df$Condition <- label
  df
}

# Build ICC data
df_e4  <- extract_icc(irt_results$Trial_4$pl2,  "Trial 4")
df_rec <- extract_icc(irt_results$Recognition$pl2, "Recognition")
df_all <- dplyr::bind_rows(df_e4, df_rec)
df_all$Condition <- factor(df_all$Condition, levels = c("Trial 4", "Recognition"))

# Dynamic linetypes: repeat a vector to cover all items
n_items <- length(unique(df_all$Item))
linetypes <- rep_len(c(1, 2, 3, 4, 5, 6), n_items)

p_icc <- ggplot(df_all, aes(x = Theta, y = Prob, group = Item, linetype = Item)) +
  geom_line(color = "black", linewidth = 0.55) +
  facet_wrap(~ Condition, ncol = 2) +
  scale_linetype_manual(values = linetypes) +
  labs(
    title = NULL,
    x = "Ability (θ)",
    y = "Probability of Correct Response"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave(filename = "outputs/Figure_2.png", plot = p_icc,
       width = 10, height = 5, dpi = 300)

# ============================================================
# (B) MODELING PLOTS — ggpredict from lme1
# ============================================================

if (!exists("lme1")) stop("`lme1` not found.")

custom_theme <- theme(
  plot.title   = element_text(size = 14, color = "black"),
  axis.title.x = element_text(size = 10, color = "black"),
  axis.title.y = element_text(size = 12, color = "black"),
  axis.text.x  = element_text(size = 10, color = "black"),
  axis.text.y  = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 10, color = "black"),
  legend.text  = element_text(size = 9,  color = "black")
)

# Note: ggpredict() returns a 'ggeffects' object with a plot() method
# We keep CI off (ci_level = FALSE) to match your original.

p1 <- ggpredict(lme1, terms = c("MPE",  "country"), ci_level = FALSE) %>%
  plot(colors = "bw") +
  labs(
    x = "MPE (years)",
    y = "Predicted mean ability-score",
    title = NULL,
    linetype = "Country"
  ) +
  scale_x_continuous(breaks = seq(0, 25, by = 2)) +
  custom_theme

p2 <- ggpredict(lme1, terms = c("age",  "Trial"),   ci_level = FALSE) %>%
  plot(colors = "bw") +
  labs(
    x = "Age (years)",
    y = "Predicted mean ability-score",
    title = NULL,
    linetype = "Trial"
  ) +
  scale_x_continuous(breaks = seq(6, 17, by = 1)) +
  custom_theme +
  theme(legend.position = "right")

p3 <- ggpredict(lme1, terms = c("country", "sex"),  ci_level = FALSE) %>%
  plot(colors = "bw") +
  labs(
    x = "Country",
    y = "Predicted mean theta-score",
    title = NULL,
    linetype = "Sex"
  ) +
  custom_theme

p4 <- ggpredict(lme1, terms = c("Trial",   "country"), ci_level = FALSE) %>%
  plot(colors = "bw") +
  labs(
    x = "Trial",
    y = "Predicted mean theta-score",
    title = NULL,
    linetype = "Country"
  ) +
  custom_theme

# Arrange into a single figure and export (high DPI)
figure_3 <- cowplot::plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol = 2)
ggsave("outputs/Figure_3.jpeg", plot = figure_3,
       width = 16, height = 10, dpi = 600, units = "in", device = "jpeg")
