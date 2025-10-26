# ============================================================
# File: 02_variable_selection.R
# Purpose:
#   Variable selection (best-subset) and mixed-effects modeling
#   for theta using tidy, reproducible practices.
# Author: Rivera et al., 2025
# Date: 2025-10-26
# ============================================================

# --- 0) Packages -------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(leaps)        # regsubsets()
  library(lme4)         # lmer()
  library(lmerTest)     # p-values for lmer
  library(broom)        # tidy() / glance() for models
  library(rlang)
  # Optional: library(merTools) # REsim(), if you want RE CIs
})


# --- 1) Inputs & basic checks ------------------------------------------------
df2 <- readRDS("data/df2.rds")
colnames(df2)

if (!exists("df2")) {
  stop("`df2` not found in the environment. Please load your analysis dataset.")
}

# --- 1) Rename columns -------------------------------------------------------
# Rename Delayed_Recall -> Trial_5, Recognition -> Trial_6
df2 <- df2 %>%
  rename(
    Trial_5 = Delayed_Recall,
    Trial_6 = Recognition
  )


# Convertir a formato largo manteniendo el ID
df2 <- df2 %>%
  pivot_longer(-c(ID, country, MPE, age, sex),  # Mantener solo estas columnas
               names_to = "Trial", # Nombre de la nueva columna 
               values_to = "theta") 

df2 <- na.exclude(df2)

# --- 2) Best-subset variable selection (no random effects) -------------------
# Model space with interactions up to 2-way; cap size via nvmax to avoid explosion.
# NOTE: leaps::regsubsets works on a fixed-effects linear model space.
sel_formula <- as.formula(
  theta ~ (poly(age, 2) + log(MPE + 1) + sex + country + Trial)^2
)

# Control flags
nvmax_val    <- 20    # maximum number of terms in the subset
criterion    <- "bic" 
plot_selection <- FALSE  # set TRUE if you want simple base R plots

message("Running best-subset selection with criterion = ", toupper(criterion))

mod_sel <- regsubsets(
  sel_formula,
  data = df2,
  really.big = TRUE,
  nvmax = nvmax_val,
  method = "seqrep"
)

sum_sel <- summary(mod_sel)

# Choose best model according to the criterion
best_idx <- switch(
  tolower(criterion),
  "bic"   = which.min(sum_sel$bic),
  "adjr2" = which.max(sum_sel$adjr2),
  "cp"    = which.min(sum_sel$cp),
  {
    warning("Unknown criterion; defaulting to BIC.")
    which.min(sum_sel$bic)
  }
)

message("Best subset index (", toupper(criterion), "): ", best_idx)

# Optional quick-look plots (disabled by default)
if (isTRUE(plot_selection)) {
  op <- par(mfrow = c(1, 3))
  plot(sum_sel$adjr2, type = "b", main = "Adj R^2"); abline(v = which.max(sum_sel$adjr2), col = 2)
  plot(sum_sel$cp,    type = "b", main = "Mallows Cp"); abline(v = which.min(sum_sel$cp), col = 2)
  plot(sum_sel$bic,   type = "b", main = "BIC");        abline(v = which.min(sum_sel$bic), col = 2)
  par(op)
}

# Extract selected coefficients (names correspond to model.matrix columns)
coefs_best <- coef(mod_sel, id = best_idx)
selected_terms <- names(coefs_best)[-1] # drop intercept label
message("Selected fixed-effect terms (no RE):")
print(selected_terms)


# --- 3) Mixed-effects model (random intercept for ID) ------------------------
# We append (1|ID) to the selected fixed-effects formula.
# If you need cross-level interactions not captured by best-subset, add them here manually.
lme1 <- lmer(theta ~ 
               poly(age, 2)+ 
               log(MPE + 1)+
               Trial+
               sex+
               country+
               poly(age,2):Trial+
               log(MPE + 1):country +
               sex:country +
               country:Trial+
               (1|ID), 
             data=df2)

# Summaries
cat("\n--- lmer summary ---\n")
print(summary(lme1))
cat("\n--- Random effects (variance components) ---\n")
print(lme4::VarCorr(lme1))


# --- 4) Prediction and normative percentile example -------------------------
# Goal: Given a new case, predict expected theta (population-level) from lme1,
#       compute an observed theta from IRT (EAP), and convert to a percentile.
irt_results <- readRDS("outputs/irt_results.rds")

# 1) New case (EDIT values/levels as needed)
newdat <- data.frame(
  age     = 17,
  MPE     = 18,
  Trial   = "Trial_6",  # must match a valid level in df2$Trial
  sex     = "Girl",     # must match a valid level in df2$sex
  country = "Spain"     # must match a valid level in df2$country
)

# 2) Population-level prediction from the mixed model (no random effects)
#    re.form = NA -> ignores subject-specific random intercepts
pred <- predict(lme1, newdata = newdat, re.form = NA)

# 3) Residual SD of the mixed model (used for percentile conversion)
sigma_eps <- sigma(lme1)

# 4) IRT observed theta (EAP) from the 2PL model for Trial_6 (recognition)
#    NOTE: The response pattern must have the same number of items as m2pl_rec.
#    Here we give a single subject with 12 items (0/1 responses).
theta <- ltm::factor.scores(
  irt_results$Recognition$pl2,
  resp.patterns = rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
  method = "EAP"
)$score.dat$z1

# 5) Percentile using z-score under Normal(pred, sigma_eps)
z   <- as.numeric((theta - pred) / sigma_eps)
pct <- pnorm(z) * 100

# 6) (Equivalent) Percentile via CDF form P(Theta <= theta_obs | N(pred, sigma_eps))
prob <- pnorm(theta, mean = pred, sd = sigma_eps) * 100

# 7) Print a compact summary
cat("\n--- Normative conversion ---\n",
    "Observed theta (IRT):    ", round(theta, 3), "\n",
    "Predicted theta (LMM):   ", round(pred, 3), "\n",
    "Residual sigma (LMM):    ", round(sigma_eps, 3), "\n",
    "Percentile (z method):   ", round(pct, 1), "%\n",
    "Percentile (CDF method): ", round(prob, 1), "%\n", sep = "")


# --- 5) (Optional) Save artifacts -------------------------------------------
# readr::write_csv(broom::tidy(lme1, effects = "fixed", conf.int = TRUE), "outputs/lmm_fixed_effects.csv")
# saveRDS(lme1, "outputs/lmm_selected_model.rds")
# saveRDS(list(model_selection = sum_sel, selected_coefs = coefs_best), "outputs/selection_artifacts.rds")

# --- End of file -------------------------------------------------------------
