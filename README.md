# tamvi-irt-lmm-scripts
R scripts used for the psychometric and normative analyses of the TAMV-I Learning and Memory Test. This repository contains the reproducible code used in the manuscript “Establishing Normative Data for the TAMV-I Learning and Memory Test: An Approach Using Item Response Theory and Linear Mixed Models"

## Script overview

### `01_irt_models.R` – IRT modeling and theta extraction
- Fits **Rasch (1PL)** and **2PL** models for each trial using `ltm::rasch()` and `ltm::ltm()`.
- Compares models via log-likelihood ratio tests, AIC, and BIC.
- Extracts individual **EAP theta scores** (`factor.scores(method = "EAP")`).
- Appends theta estimates (`Trial_1` to `Trial_6`) to the dataset `tamvi`.

---

### `02_variable_selection.R` – Variable selection and mixed modeling
- Performs **best-subset variable selection** using `leaps::regsubsets()`  
  across polynomial age, log-transformed MPE, and interaction terms.
- Selects the optimal model based on **BIC** (default).
- Refits the selected structure in a **Linear Mixed-Effects Model** (`lmer()`)  
  with random intercepts for `ID`.
- Generates predictions for new cases and computes **normative percentiles**  
  comparing observed IRT scores (θ) to model-based expectations.

---

### `03_plots.R` – Figures for IRT and modeling results
- Loads the saved models (`irt_results` and `lme1`).
- (A) **IRT plots:** Item Characteristic Curves (ICCs) for selected 2PL models.  
  - Exports `outputs/Figure_2.png`
- (B) **Modeling plots:** Marginal predictions using `ggeffects::ggpredict()`  
  across key covariates (MPE, age, Trial, sex, country).  
  - Exports `outputs/Figure_3.jpeg` (high-resolution, 600 dpi)

---

### Reproducible workflow

Run the scripts sequentially:

```r
source("scripts/01_irt_models.R")
source("scripts/02_variable_selection.R")
source("scripts/03_plots.R")
