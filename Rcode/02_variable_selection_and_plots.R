############
# Packages #
############
library(tidyverse)
library(leaps)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(rlang)
library(ltm)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(cowplot)

#########################
# Inputs & basic checks #
#########################
source("Rcode/auxiliary_functions.R")
df2 <- readRDS("data/df2.rds")
colnames(df2)

if (!exists("df2")) {
  stop("`df2` not found in the environment. Please load your analysis dataset.")
}

##################
# Rename columns #
##################
df2 <- df2 %>%
  rename(Trial_5 = Delayed_Recall, Trial_6 = Recognition)

##########################
# convert in long format #
##########################
df2 <- df2 %>% pivot_longer(-c(ID, country, MPE, age, sex),
                            names_to = "Trial",
                            values_to = "theta")
df2 <- na.exclude(df2)

######################################################
# Best-subset variable selection (no random effects) #
######################################################
sel_formula <- as.formula(theta ~ (poly(age, 2) + log(MPE + 1) + sex + country + Trial)^2)

#################
# Control flags #
#################
nvmax_val <- 20
criterion <- "bic" 
plot_selection <- FALSE  

message("Running best-subset selection with criterion = ", toupper(criterion))

mod_sel <- regsubsets(
  sel_formula,
  data = df2,
  really.big = TRUE,
  nvmax = nvmax_val,
  method = "seqrep"
)

sum_sel <- summary(mod_sel)

################################################
# Choose best model according to the criterion #
################################################
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

###################################################
# Optional quick-look plots (disabled by default) #
###################################################
if (isTRUE(plot_selection)) {
  op <- par(mfrow = c(1, 3))
  plot(sum_sel$adjr2, type = "b", main = "Adj R^2"); abline(v = which.max(sum_sel$adjr2), col = 2)
  plot(sum_sel$cp,    type = "b", main = "Mallows Cp"); abline(v = which.min(sum_sel$cp), col = 2)
  plot(sum_sel$bic,   type = "b", main = "BIC");        abline(v = which.min(sum_sel$bic), col = 2)
  par(op)
}

############################################################################
# Extract selected coefficients (names correspond to model.matrix columns) #
############################################################################
coefs_best <- coef(mod_sel, id = best_idx)
selected_terms <- names(coefs_best)[-1]
message("Selected fixed-effect terms (no RE):")
print(selected_terms)

#################################################
# Mixed-effects model (random intercept for ID) #
#################################################
lme1 <- lmer(
  theta ~
    poly(age, 2) +
    log(MPE + 1) +
    Trial +
    sex +
    country +
    poly(age, 2):Trial +
    log(MPE + 1):country +
    sex:country +
    country:Trial +
    (1 | ID),
  data = df2
)

print(summary(lme1))
print(lme4::VarCorr(lme1))

###############################################
# Prediction and normative percentile example #
###############################################
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

#############
# IRT plots #
#############

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

width <- 7
height <- width / (16/9)
units <- "in"
dpi <- 300
ggsave(filename = "outputs/Figure_2.tiff", plot = p_icc, units = units,
       width = width, height = height, dpi = dpi, compression = "lzw")

########################################
# MODELING PLOTS — ggpredict from lme1 #
########################################

if (!exists("lme1")) stop("`lme1` not found.")

custom_theme <- theme(
  plot.title   = element_text(size = 8, color = "black"),
  axis.title.x = element_text(size = 8, color = "black"),
  axis.title.y = element_text(size = 8, color = "black"),
  axis.text.x  = element_text(size = 8, color = "black"),
  axis.text.y  = element_text(size = 8, color = "black"),
  legend.title = element_text(size = 8, color = "black"),
  legend.text  = element_text(size = 8,  color = "black")
)

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


# Arrange into a single figure and export #
figure_3 <- cowplot::plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol = 2)
width <- 7.5
height <- width / (5/4)
ggsave("outputs/Figure_3.tiff", plot = figure_3, units = units,
       width = width, height = height, dpi = dpi, compression = "lzw")

