# tamvi-irt-lmm-scripts

This repository contains the reproducible code used in the manuscript: _Normative Data for Learning and Memory Test (TAMV-I) in Latin American and Spanish Children: An Item Response Theory and Linear Mixed Models Approach_

## Table of contents

-   [R code](#R-code)
-   [Reproducible workflow](#Reproducible workflow)
-   [References](#References)

# R-code

This folder contains the following files:

-   [**01_irt_models.R**](https://github.com/diegoriveraps/tamvi-irt-lmm-scripts/Rcode/01_irt_models.R)

R script Fits **Rasch (1PL)** and **2PL** models for each trial using `ltm::rasch()` and `ltm::ltm()`, compares models via log-likelihood ratio tests, AIC, and BIC, extracts individual **EAP theta scores** (`factor.scores(method = "EAP")`), and appends theta estimates (`Trial_1` to `Trial_6`) to the dataset `tamvi`.

-   [**02_variable_selection_and_plots.R**](https://github.com/diegoriveraps/tamvi-irt-lmm-scripts/Rcode/02_variable_selection_and_plots.R)

R scripts performs **best-subset variable selection** using `leaps::regsubsets()` across polynomial age, log-transformed MPE, and interaction terms; selects the optimal model based on **BIC** (default); refits the selected structure in a **Linear Mixed-Effects Model** (`lmer()`) with random intercepts for `ID`; generates predictions for new cases and computes **normative percentiles** by comparing observed IRT scores ($\theta$) to model-based expectations; and generate Figures 2 and 3 showing IRT and modeling results at 300 dpi.

-   [**auxiliary_functions.R**](https://github.com/diegoriveraps/tamvi-irt-lmm-scripts/Rcode/auxiliary_functions.R)

R script containing functions that execute procedures used in the main scripts.

# Reproducible workflow

Run the scripts sequentially:

```
source("scripts/01_irt_models.R")
source("scripts/02_variable_selection_and_plots.R")
```
# References
