# ============================================================
# Title: TAMVI — IRT Analysis (Rasch vs 2PL) and Normative Data
# Author: Rivera et al., 2025
# Date: 2025-10-27
# Description:
#   Clean, reproducible script to fit Rasch and 2PL (ltm) models
#   across predefined item blocks, compare models (ANOVA),
#   and append EAP theta scores per case.
# ============================================================

# --- 1) Packages -------------------------------------------------------------

# Note: 'rasch', 'ltm', 'factor.scores' come from the ltm package
suppressPackageStartupMessages({
  library(ltm)
})

# --- 2) Load data --------------------------------------------------------------

load("00_tamvi.RData")
tamvi <- na.exclude(tamvi)

# Expect: an object named `tamvi` already in memory with item columns.
# If not, stop with a clear message.
if (!exists("tamvi")) {
  stop("Object `tamvi` not found. Please load your data before running this script.")
}

# ============================================================
# Title: TAMVI — IRT looped analysis (Rasch vs 2PL) without plots
# Author: Diego Rivera
# Date: 2025-10-25
# Description:
#   Loop over item blocks, fit Rasch and 2PL (ltm) models,
#   compare models (ANOVA), compute EAP theta, and append
#   theta columns to `tamvi` (Trial_1 ... Trial_4).
#   No plots are produced.
# ============================================================

# --- 1) Packages -------------------------------------------------------------
suppressPackageStartupMessages({
  library(ltm)
  library(dplyr)
  library(purrr)
  library(tibble)
})

# --- 2) Create directory if it doesn’t exist ---------------------------------
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# --- 3) Preconditions --------------------------------------------------------
if (!exists("tamvi")) {
  stop("Object `tamvi` not found. Please load your data before running this script.")
}

# --- 4) Block specification --------------------------------------------------
# NOTE: update column ranges if your item layout changes.
blocks <- list(
  Trial_1 = 7:18,
  Trial_2 = 20:31,
  Trial_3 = 33:44,
  Trial_4 = 46:57,
  Delayed_Recall = 59:70,
  Recognition = 72:83
)

# --- 5) Helper to fit a single block ----------------------------------------
fit_block <- function(data, cols) {
  # Defensive checks
  stopifnot(is.data.frame(data))
  stopifnot(all(cols %in% seq_len(ncol(data))))
  
  items <- data[, cols]
  
  # Rasch (1PL)
  rasch_fit <- ltm::rasch(
    items,
    constraint = NULL,
    IRT.param  = TRUE,
    start.val  = NULL,
    na.action  = NULL,
    control    = list(),
    Hessian    = TRUE
  )
  
  # 2PL
  pl2_fit <- ltm::ltm(
    items ~ z1,
    IRT.param = TRUE,
    na.action = stats::na.exclude
  )
  
  # Model comparison
  comp <- stats::anova(rasch_fit, pl2_fit)
  
  # EAP factor scores from 2PL
  fs <- ltm::factor.scores(
    pl2_fit,
    resp.patterns = items,
    method = "EAP"
  )
  theta <- fs$score.dat$z1
  
  list(
    rasch = rasch_fit,
    pl2   = pl2_fit,
    anova = comp,
    theta = theta
  )
}

# --- 6) Run over all blocks --------------------------------------------------
irt_results <- imap(
  blocks,
  ~ {
    message("Fitting block: ", .y, " | cols: ", paste0(range(.x), collapse = ":"))
    fit_block(tamvi, .x)
  }
)

# --- 7) Append theta to `tamvi` ---------------------------------------------
# Each element’s name (e.g., "Trial_1") becomes a new column with its theta.
for (nm in names(irt_results)) {
  tamvi[[nm]] <- irt_results[[nm]]$theta
}

# --- 8) Compact summaries to inspect in console ------------------------------
# (a) Print model objects (optional; comment out if you want a quieter log)
walk(irt_results, ~ { print(.x$rasch); print(.x$pl2) })

# (b) Build a tidy ANOVA-like summary using LRT, AIC, BIC
lr_row <- function(block_name, rasch_fit, pl2_fit) {
  ll1 <- as.numeric(logLik(rasch_fit))
  ll2 <- as.numeric(logLik(pl2_fit))
  df1 <- attr(logLik(rasch_fit), "df")
  df2 <- attr(logLik(pl2_fit),  "df")
  
  # Likelihood-ratio test
  lr_stat <- 2 * (ll2 - ll1)
  df_diff <- df2 - df1
  p_val   <- stats::pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
  
  tibble::tibble(
    Block      = block_name,
    Model_1    = "Rasch (1PL)",
    Model_2    = "2PL",
    LL_1       = ll1,
    LL_2       = ll2,
    df_1       = df1,
    df_2       = df2,
    LR_stat    = lr_stat,
    df_diff    = df_diff,
    p_value    = p_val,
    AIC_1      = AIC(rasch_fit),
    AIC_2      = AIC(pl2_fit),
    BIC_1      = BIC(rasch_fit),
    BIC_2      = BIC(pl2_fit),
    Better_AIC = c("Rasch","2PL")[which.min(c(AIC(rasch_fit), AIC(pl2_fit)))]
  )
}

anova_tbl <- purrr::imap_dfr(
  irt_results,
  ~ lr_row(.y, .x$rasch, .x$pl2)
)

print(anova_tbl)

# --- Save outputs ------------------------------------------------------------

saveRDS(irt_results, "outputs/irt_results.rds")


# ============================================================
# Add EAP theta scores from 2PL fits to `tamvi`
# ============================================================

# Define the column ranges for each trial
blocks <- list(
  Trial_1        = 7:18,
  Trial_2        = 20:31,
  Trial_3        = 33:44,
  Trial_4        = 46:57,
  Delayed_Recall = 59:70,
  Recognition    = 72:83
)

# Loop over blocks to compute and append theta (EAP)
for (nm in names(blocks)) {
  message("Computing EAP scores for: ", nm)
  
  # Extract item responses for this block
  items_block <- tamvi[, blocks[[nm]]]
  
  # Extract corresponding 2PL fit from irt_results
  model_2pl <- irt_results[[nm]]$pl2
  
  # Compute EAP factor scores
  fs <- ltm::factor.scores(
    model_2pl,
    resp.patterns = items_block,
    method = "EAP"
  )
  
  # Add theta (z1) to tamvi
  tamvi[[nm]] <- fs$score.dat$z1
  
  # Optional: also store SE of theta
  # tamvi[[paste0(nm, "_se")]] <- fs$score.dat$se.z1
}

# Verify that all new columns were added
names(tamvi)[names(tamvi) %in% names(blocks)]

# ============================================================
# Save df2 subset from tamvi
# ============================================================

# Create subset: columns 1–5 (sociodemographic / IDs, etc.) + 86–91 (Trial_1–Trial_6)
df2 <- tamvi[, c(1:5, 86:91)]

# --- Save outputs ------------------------------------------------------------

# Save as native R format
saveRDS(df2, file = "data/df2.rds")


# --- Reproducibility -----------------------------------------------------
sessionInfo()
