# --- 4) Helper to fit a single block ----------------------------------------
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