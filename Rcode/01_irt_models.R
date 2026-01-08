############
# Packages #
############
library(ltm)
library(dplyr)
library(purrr)
# library(readr)
library(tibble)

#############
# Load data #
#############
tamvi <-  read_delim("https://zenodo.org/record/18176818/files/tamvi.csv?download=1", delim = ";")
source("Rcode/auxiliary_functions.R")
tamvi <- na.exclude(tamvi)
set.seed(1234)
tamvi <- tamvi %>%
  mutate(ID = replicate(n(), paste0(sample(c(0:9, LETTERS, letters), 8, replace = TRUE), collapse = "")))

if (!exists("tamvi")) {
  stop("Object `tamvi` not found. Please load your data before running this script.")
}

########################################
# Create directory if it doesn’t exist #
########################################
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

#######################
# Block specification #
#######################
blocks <- list(
  Trial_1 = grep("E1_", colnames(tamvi)),
  Trial_2 = grep("E2_", colnames(tamvi)),
  Trial_3 = grep("E3_", colnames(tamvi)),
  Trial_4 = grep("E4_", colnames(tamvi)),
  Delayed_Recall = grep("ET", colnames(tamvi)),
  Recognition = grep("R_", colnames(tamvi))
)

#######################
# Run over all blocks #
#######################
irt_results <- purrr::imap(blocks, ~ {
  message("Fitting block: ", .y, " | cols: ", paste0(range(.x), collapse = ":"))
  fit_block(tamvi, .x)
})

###########################################################
# Append theta to `tamvi`                                 #
# Note that each element’s name (e.g., "Trial_1") becomes #
# a new column with its theta.                            #
###########################################################
for (nm in names(irt_results)) {
  tamvi[[nm]] <- irt_results[[nm]]$theta
}

###########################################
# Compact summaries to inspect in console #
###########################################
walk(irt_results, ~ { print(.x$rasch); print(.x$pl2) })

#######################################################
# Build a tidy ANOVA-like summary using LRT, AIC, BIC #
#######################################################
anova_tbl <- purrr::imap_dfr(irt_results, ~ lr_row(.y, .x$rasch, .x$pl2))
print(anova_tbl)

################
# Save outputs #
################
saveRDS(irt_results, "outputs/irt_results.rds")

#################################################
# Add EAP theta scores from 2PL fits to `tamvi` #
#################################################

######################################################
# Loop over blocks to compute and append theta (EAP) #
######################################################
for (nm in names(blocks)) {
  message("Computing EAP scores for: ", nm)
  items_block <- tamvi[, blocks[[nm]]]
  model_2pl <- irt_results[[nm]]$pl2
  fs <- ltm::factor.scores(model_2pl, resp.patterns = items_block, method = "EAP")
  tamvi[[nm]] <- fs$score.dat$z1
}

##########################################
# Verify that all new columns were added #
##########################################
names(tamvi)[names(tamvi) %in% names(blocks)]

##############################
# Save df2 subset from tamvi #
##############################
df2 <- tamvi[, c(1:5, 86:91)]
saveRDS(df2, file = "data/df2.rds")

###################
# Reproducibility #
###################
sessionInfo()

