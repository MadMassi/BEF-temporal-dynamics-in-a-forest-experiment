## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Design models
# Author: Massimo Martini

# Script description: This is the first script in the analysis workflow. It
# reads plot_data.csv, prepares the plot-year variables and response-specific
# datasets, and fits design GLMMs for the five focal responses. Biological
# fixed effects are specified a priori; model comparisons evaluate suitable
# repeated-measures and dispersion structures, alongside declared sensitivity
# checks. The script saves the prepared data, selected design models, and model-
# selection results used by the subsequent analyses.

# Run this script from the project root.

# Reproducibility: run with R 4.6.1 and use package versions recorded
# in renv.lock. The supplied .Rprofile activates renv;
# if needed run:

# source("renv/activate.R")
# renv::restore(retry = FALSE)
# renv::status()

# For simplicity, we leave out plotting procedures

# Approximate runtime: 7-8 minutes

# ------------------------------------------------------------------------------
# 1. Packages and helper functions ####
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(glmmTMB)
library(tibble)
library(tidyr)

source("analysis_functions.R")

# ------------------------------------------------------------------------------
# 2. Read and prepare data ####
# ------------------------------------------------------------------------------

data_file <- "plot_data.csv"
if (!file.exists(data_file)) {
  stop("Cannot find plot_data.csv. Run the script from the repository root.")
}

plot_data <- read_csv(data_file, show_col_types = FALSE)

# Fail early if the released input is incomplete or internally inconsistent.
required_columns <- c(
  "site", "plot_id", "plot_yr", "year", "tree_r", "tree_fd",
  "stand_volume", "elevation", "eastness", "northness", "slope",
  "annual_temperature", "annual_humidity", "total_cells", "overall_cells",
  "host_rich", "par_abund", "par_rich", "network_size", "n_links",
  "linkage_density", "h2", "mean_dprime_hl", "niche_overlap_hl",
  "interaction_evenness", "robustness_hl", "focal_rich"
)
missing_columns <- setdiff(required_columns, names(plot_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(plot_data) != 792L || anyDuplicated(plot_data$plot_yr)) {
  stop("plot_data.csv must contain 792 unique plot-year records.")
}
if (any(plot_data$tree_r <= 0, na.rm = TRUE) ||
    any(plot_data$par_abund < 0, na.rm = TRUE) ||
    any(plot_data$par_abund > plot_data$total_cells, na.rm = TRUE) ||
    any(plot_data$par_abund > plot_data$overall_cells, na.rm = TRUE)) {
  stop("Invalid richness or brood-cell counts in plot_data.csv.")
}

# Group four 24-species plots with the 16-species level, then express richness
# on a doubling scale across the experimental gradient.
plot_data$tree_r[plot_data$tree_r > 16] <- 16
plot_data$log_tr <- log2(plot_data$tree_r)

plot_data$par_rich10 <- plot_data$par_rich
plot_data$par_rich10[plot_data$par_rich > 10] <- 10
# Cap the four communities above 10 parasitoid species. This predictor is used
# only in the later parasitism models to avoid extrapolation at the sparse tail.

plot_data <- plot_data |>
  arrange(plot_id, year)

plot_data$for_age <- ifelse(plot_data$site == "A",
                            plot_data$year - 2009,
                            plot_data$year - 2010)
# Stand age accounts for the one-year establishment difference between sites.

year_levels <- seq.int(min(plot_data$year), max(plot_data$year))
plot_data$year <- factor(
  plot_data$year,
  levels = year_levels
) # Keep the unsampled 2017 level so observation spacing remains explicit.

plot_data <- plot_data |>
  mutate(
    cell_prate = if_else(
      total_cells > 0,
      par_abund / total_cells,
      NA_real_,
      missing = NA_real_
    )
  )

plot_data$p_failure <- (plot_data$total_cells - plot_data$par_abund)
# Primary failures use successfully developed brood cells.

plot_data$p_failure2 <- (plot_data$overall_cells - plot_data$par_abund)
# The mortality sensitivity adds confirmed premature deaths as failures.

plot_data <- plot_data |>
  mutate(
    rtsv = sqrt(stand_volume),
    rtfd = sqrt(tree_fd),
    rtcells = sqrt(total_cells),
    rtovcells = sqrt(overall_cells)
  )

pd_all <- plot_data

# Scale continuous predictors once on the full 792-row dataset. Downstream
# subsets therefore retain the same reference means and standard deviations.
name_map <- c(
  log_tr               = "sc_logtr",
  stand_volume         = "sc_sv",
  rtsv                 = "sc_rtsv",
  tree_fd              = "sc_fd",
  rtfd                 = "sc_rtfd",
  for_age              = "sc_fa",
  elevation            = "sc_elev",
  eastness             = "sc_east",
  northness            = "sc_north",
  slope                = "sc_slope",
  annual_temperature   = "sc_temp",
  annual_humidity      = "sc_humid",
  total_cells          = "sc_cells",
  rtcells              = "sc_rtcells",
  overall_cells        = "sc_overall_cells",
  rtovcells            = "sc_rtovcells",
  host_rich            = "sc_hr",
  par_rich             = "sc_pr",
  par_rich10           = "sc_pr10",
  network_size         = "sc_netsize",
  n_links              = "sc_links",
  linkage_density      = "sc_linkdense",
  h2                   = "sc_h2",
  mean_dprime_hl       = "sc_mdprime",
  niche_overlap_hl     = "sc_niche",
  interaction_evenness = "sc_intev",
  robustness_hl        = "sc_robust",
  focal_rich           = "sc_frich"
)

# Scale shared predictors once using the complete dataset.
for (old in names(name_map)) {
  new <- name_map[[old]]
  pd_all[[new]] <- as.numeric(scale(pd_all[[old]]))
}

pd_climate <- pd_all |>
  filter(!is.na(annual_temperature)) # Microclimate subset: 448 plot-years.

pd_parasitism <- pd_all |>
  filter(total_cells > 0) # Brood-cell-positive subset: 786 plot-years.

pd_net <- pd_all |>
  filter(!is.na(network_size)) # Network-eligible subset: 448 plot-years.

if (nrow(pd_parasitism) != 786L ||
    nrow(pd_climate) != 448L ||
    nrow(pd_net) != 448L) {
  stop("Unexpected analytical sample sizes in plot_data.csv.")
}

# Descriptive correlations complement the later model-specific VIF checks.
vars <- pd_all[, c("log_tr", "stand_volume", "tree_fd", "for_age",
                   "elevation", "eastness", "northness", "slope",
                   "annual_temperature", "annual_humidity", "total_cells",
                   "host_rich", "par_abund", "par_rich", "focal_rich")]

cortable <- cor(vars, use = "complete.obs", method = "pearson")

cortable |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "cor") |>
  filter(var1 < var2, abs(cor) > 0.60) |>
  arrange(desc(abs(cor)))

# Screen network metrics separately because they use a smaller subset.
vars2 <- pd_net[, c("network_size", "n_links", "linkage_density", "h2",
                       "mean_dprime_hl", "niche_overlap_hl",
                       "interaction_evenness", "robustness_hl", "focal_rich")]

cortable2 <- cor(vars2, use = "complete.obs", method = "pearson")

cortable2 |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "cor") |>
  filter(var1 < var2, abs(cor) > 0.60) |>
  arrange(desc(abs(cor))) # many network metrics correlated with one another

# ------------------------------------------------------------------------------
# 3. GENERALIZED LINEAR MIXED MODELS (GLMMs) ####
# ------------------------------------------------------------------------------

## 3.1 Model-selection workflow ##

# The primary hypothesis is a tree-richness-by-forest-age interaction. The
# biological fixed effects are unchanged across candidates. Selection concerns
# only repeated-measures and dispersion structures and considers numerical
# validity, information criteria, residual checks, and focal-effect stability.

# ----------------------------
## 3.2 Parasitism
# ----------------------------

## 3.2.1 Linear design model and candidate variance structures ##

prt_design_ri <- glmmTMB(
  cbind(par_abund, p_failure) ~
    sc_logtr + sc_fa +
      sc_logtr:sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = pd_parasitism,
  family = betabinomial(link = "logit"),
  dispformula = ~sc_cells
)

# Random-intercept reference for the hypothesized linear change in the
# experimental tree-richness effect through stand age. Cell-dependent
# dispersion accommodates changing heterogeneity with the trial denominator.

prt_design_rs <- update(
  prt_design_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)
# Uncorrelated plot-specific intercepts and stand-age slopes.

prt_design_ar <- update(
  prt_design_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)
# Conditional AR(1) deviations within plots, indexed by successive observed
# sampling occasions.

prt_design_disp_ar <- update(
  prt_design_ri,
  dispformula = ~
    sc_cells +
    ar1(year + 0 | site:plot_id)
)
# Plot-specific AR(1) structure in the beta-binomial dispersion model.
# This tests whether residual heterogeneity, beyond its relationship with the
# brood-cell denominator, is temporally correlated within plots.

prt_repeated_candidates <- list(
  ri = prt_design_ri,
  rs = prt_design_rs,
  ar = prt_design_ar
)

# Report numerical validity and information criteria for each candidate.
prt_repeated_fit <- compare_candidate_models(prt_repeated_candidates)

# Simulate residuals and summarize uniformity, dispersion, outliers, and zeros.
prt_repeated_diagnostics <- run_candidate_diagnostics(
  prt_repeated_candidates,
  n_sim = 1000L,
  seed = 123L
)

# Record the focal coefficients and Wald intervals across candidates.
prt_repeated_effects <- summarize_fixed_terms(
  prt_repeated_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

prt_repeated_fit
prt_repeated_diagnostics$summary
prt_repeated_effects

prt_dispersion_candidates <- list(
  ri = prt_design_ri,
  disp_ar = prt_design_disp_ar
)

prt_dispersion_fit <- compare_candidate_models(prt_dispersion_candidates)
prt_dispersion_diagnostics <- run_candidate_diagnostics(
  prt_dispersion_candidates,
  n_sim = 1000L,
  seed = 123L
)
prt_dispersion_effects <- summarize_fixed_terms(
  prt_dispersion_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

prt_dispersion_fit
prt_dispersion_diagnostics$summary
prt_dispersion_effects

## 3.2.2 Select primary design specification ##

prt_design <- prt_design_ri
prt_design_cbind <- prt_design

# Retain the random-intercept model. The alternatives do not give a clear
# information-criterion or diagnostic improvement, and the focal interaction
# is stable across specifications.

# Fit the equivalent rate-and-weights form used by later path analyses.
prt_design_sem <- glmmTMB(
  cell_prate ~
    sc_logtr * sc_fa +
      (1 | site/plot_id) + (1 | year),
  weights = total_cells,
  data = pd_parasitism,
  family = betabinomial(link = "logit"),
  dispformula = ~ sc_cells
)

## 3.2.3 Nonlinear temporal sensitivities of the selected primary model ##

if (!is.null(prt_design)) {
  prt_design_age2 <- update(
    prt_design,
    . ~ . + I(sc_fa^2)
  )

  prt_design_age2_interaction <- update(
    prt_design_age2,
    . ~ . + sc_logtr:I(sc_fa^2)
  )

  prt_nonlinear_sensitivities <- list(
    linear = prt_design,
    quadratic = prt_design_age2,
    quadratic_interaction = prt_design_age2_interaction
  )

  prt_nonlinear_fit <- compare_candidate_models(prt_nonlinear_sensitivities)
  prt_nonlinear_diagnostics <- run_candidate_diagnostics(
    prt_nonlinear_sensitivities,
    n_sim = 1000L,
    seed = 123L
  )
  prt_nonlinear_effects <- summarize_fixed_terms(
    prt_nonlinear_sensitivities,
    terms = c(
      "sc_logtr",
      "sc_fa",
      "I(sc_fa^2)",
      "sc_logtr:sc_fa",
      "sc_logtr:I(sc_fa^2)"
    )
  )
}

# Curvature is supported in the overall age trend, but not in the temporal
# change in the tree-richness effect. The linear interaction remains primary.

prt_design_present <- update(
  prt_design,
  data = subset(pd_parasitism, par_rich > 0)
) # Sensitivity excluding parasitoid absences.

prt_design_mortality <- update(
  prt_design,
  cbind(par_abund, p_failure2) ~ .,
  dispformula = ~ sc_overall_cells
) # Sensitivity to inclusion of non-parasitism mortality.

# ----------------------------
## 3.3 Intermediate responses
# ----------------------------

# Fit one design model for each downstream outcome.

## ----------------------------------
## 3.3.1 Brood cells (host abundance)

## 3.3.1.1 Linear design model and candidate variance structures ##

ha_design_ri <- glmmTMB(
  total_cells ~
    sc_logtr + sc_fa +
      sc_logtr:sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = pd_all,
  family = nbinom1()
)

ha_design_rs <- update(
  ha_design_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

ha_design_ar <- update(
  ha_design_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

ha_design_disp_ar <- update(
  ha_design_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

ha_repeated_candidates <- list(
  ri = ha_design_ri,
  rs = ha_design_rs,
  ar = ha_design_ar
)

ha_repeated_fit <- compare_candidate_models(ha_repeated_candidates)
ha_repeated_diagnostics <- run_candidate_diagnostics(
  ha_repeated_candidates,
  n_sim = 1000L,
  seed = 123L
)
ha_repeated_effects <- summarize_fixed_terms(
  ha_repeated_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

ha_repeated_fit
ha_repeated_diagnostics$summary
ha_repeated_effects

ha_dispersion_candidates <- list(
  ri = ha_design_ri,
  disp_ar = ha_design_disp_ar
)

ha_dispersion_fit <- compare_candidate_models(ha_dispersion_candidates)
ha_dispersion_diagnostics <- run_candidate_diagnostics(
  ha_dispersion_candidates,
  n_sim = 1000L,
  seed = 123L
)
ha_dispersion_effects <- summarize_fixed_terms(
  ha_dispersion_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

ha_dispersion_fit
ha_dispersion_diagnostics$summary
ha_dispersion_effects

## 3.3.1.2 Select primary design specification ##

ha_design <- ha_design_disp_ar
# Retain dispersion AR(1): it improves information criteria and resolves the
# dispersion and excess-zero discrepancies of the random-intercept model.

## 3.3.1.3 Nonlinear temporal sensitivities of the selected primary model ##

if (!is.null(ha_design)) {
  ha_design_age2 <- update(
    ha_design,
    . ~ . + I(sc_fa^2)
  )

  ha_design_age2_interaction <- update(
    ha_design_age2,
    . ~ . + sc_logtr:I(sc_fa^2)
  )

  ha_nonlinear_sensitivities <- list(
    primary_linear = ha_design,
    quadratic_age = ha_design_age2,
    quadratic_age_interaction = ha_design_age2_interaction
  )

  ha_nonlinear_fit <- compare_candidate_models(ha_nonlinear_sensitivities)
  ha_nonlinear_diagnostics <- run_candidate_diagnostics(
    ha_nonlinear_sensitivities,
    n_sim = 1000L,
    seed = 211L
  )
  ha_nonlinear_effects <- summarize_fixed_terms(
    ha_nonlinear_sensitivities,
    terms = c(
      "sc_logtr",
      "sc_fa",
      "I(sc_fa^2)",
      "sc_logtr:sc_fa",
      "sc_logtr:I(sc_fa^2)"
    )
  )
}

ha_nonlinear_fit
ha_nonlinear_diagnostics
ha_nonlinear_effects

## ----------------------------------
## 3.3.2 Host species richness

## 3.3.2.1 Linear design model and candidate variance structures ##

hr_design_ri <- glmmTMB(
  host_rich ~
    sc_logtr + sc_fa +
      sc_logtr : sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = pd_all,
  family = genpois()
)

hr_design_rs <- update(
  hr_design_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

hr_design_ar <- update(
  hr_design_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

hr_design_disp_ar <- update(
  hr_design_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

hr_repeated_candidates <- list(
  ri = hr_design_ri,
  rs = hr_design_rs,
  ar = hr_design_ar
)

hr_repeated_fit <- compare_candidate_models(hr_repeated_candidates)
hr_repeated_diagnostics <- run_candidate_diagnostics(
  hr_repeated_candidates,
  n_sim = 1000L,
  seed = 301L
)
hr_repeated_effects <- summarize_fixed_terms(
  hr_repeated_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

hr_dispersion_candidates <- list(
  ri = hr_design_ri,
  disp_ar = hr_design_disp_ar
)

hr_dispersion_fit <- compare_candidate_models(hr_dispersion_candidates)
hr_dispersion_diagnostics <- run_candidate_diagnostics(
  hr_dispersion_candidates,
  n_sim = 1000L,
  seed = 311L
)
hr_dispersion_effects <- summarize_fixed_terms(
  hr_dispersion_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

hr_repeated_fit
hr_repeated_diagnostics$summary
hr_repeated_effects

hr_dispersion_fit
hr_dispersion_diagnostics$summary
hr_dispersion_effects

## 3.3.2.2 Select primary design specification ##

hr_design <- hr_design_ri
# Retain the simpler random-intercept model. More complex candidates offer no
# clear improvement and conditional AR(1) produces poor residual diagnostics.

## 3.3.2.3 Nonlinear temporal sensitivities of the selected primary model ##

if (!is.null(hr_design)) {
  hr_design_age2 <- update(
    hr_design,
    . ~ . + I(sc_fa^2)
  )

  hr_design_age2_interaction <- update(
    hr_design_age2,
    . ~ . + sc_logtr:I(sc_fa^2)
  )

  hr_nonlinear_sensitivities <- list(
    primary_linear = hr_design,
    quadratic_age = hr_design_age2,
    quadratic_age_interaction = hr_design_age2_interaction
  )

  hr_nonlinear_fit <- compare_candidate_models(hr_nonlinear_sensitivities)
  hr_nonlinear_diagnostics <- run_candidate_diagnostics(
    hr_nonlinear_sensitivities,
    n_sim = 1000L,
    seed = 321L
  )
  hr_nonlinear_effects <- summarize_fixed_terms(
    hr_nonlinear_sensitivities,
    terms = c(
      "sc_logtr",
      "sc_fa",
      "I(sc_fa^2)",
      "sc_logtr:sc_fa",
      "sc_logtr:I(sc_fa^2)"
    )
  )
}

## ----------------------------------
## 3.3.3 Realized parasitism (or parasitoid abundance)

## 3.3.3.1 Linear design model and candidate variance structures ##

pa_design_ri <- glmmTMB(
  par_abund ~
    sc_logtr + sc_fa +
      sc_logtr : sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = pd_parasitism,
  family = nbinom1(),
  ziformula = ~1
)

pa_design_rs <- update(
  pa_design_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

pa_design_ar <- update(
  pa_design_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

pa_design_disp_ar <- update(
  pa_design_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

pa_repeated_candidates <- list(
  ri = pa_design_ri,
  rs = pa_design_rs,
  ar = pa_design_ar
)

pa_repeated_fit <- compare_candidate_models(pa_repeated_candidates)
pa_repeated_diagnostics <- run_candidate_diagnostics(
  pa_repeated_candidates,
  n_sim = 1000L,
  seed = 401L
)
pa_repeated_effects <- summarize_fixed_terms(
  pa_repeated_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

pa_dispersion_candidates <- list(
  ri = pa_design_ri,
  disp_ar = pa_design_disp_ar
)

pa_dispersion_fit <- compare_candidate_models(pa_dispersion_candidates)
pa_dispersion_diagnostics <- run_candidate_diagnostics(
  pa_dispersion_candidates,
  n_sim = 1000L,
  seed = 411L
)
pa_dispersion_effects <- summarize_fixed_terms(
  pa_dispersion_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

pa_repeated_fit
pa_repeated_diagnostics$summary
pa_repeated_effects

pa_dispersion_fit
pa_dispersion_diagnostics$summary
pa_dispersion_effects

pa_design_ri2 <- update(
  pa_design_ri,
  . ~ . - (1 | site/plot_id) +
    site + (1 | plot_id)
) # Treat site as fixed to avoid a singular site-level random effect.

## 3.3.3.2 Select primary design specification ##

pa_design <- pa_design_ri2
# Retain the random-intercept structure; the random slope has an invalid
# Hessian and neither AR(1) alternative improves fit or diagnostics.

## 3.3.3.3 Nonlinear temporal sensitivities of the selected primary model ##

if (!is.null(pa_design)) {
  pa_design_age2 <- update(
    pa_design,
    . ~ . + I(sc_fa^2)
  )

  pa_design_age2_interaction <- update(
    pa_design_age2,
    . ~ . + sc_logtr:I(sc_fa^2)
  )

  pa_nonlinear_sensitivities <- list(
    primary_linear = pa_design,
    quadratic_age = pa_design_age2,
    quadratic_age_interaction = pa_design_age2_interaction
  )

  pa_nonlinear_fit <- compare_candidate_models(pa_nonlinear_sensitivities)
  pa_nonlinear_diagnostics <- run_candidate_diagnostics(
    pa_nonlinear_sensitivities,
    n_sim = 1000L,
    seed = 421L
  )
  pa_nonlinear_effects <- summarize_fixed_terms(
    pa_nonlinear_sensitivities,
    terms = c(
      "sc_logtr",
      "sc_fa",
      "I(sc_fa^2)",
      "sc_logtr:sc_fa",
      "sc_logtr:I(sc_fa^2)"
    )
  )
}

pa_design_present <- update(
  pa_design,
  data = subset(pd_parasitism, par_rich > 0),
  family = truncated_nbinom1()
) # Sensitivity excluding parasitoid absences.

## ----------------------------------
## 3.3.4 Parasitoid species richness

## 3.3.4.1 Linear design model and candidate variance structures ##

pr_design_ri <- glmmTMB(
  par_rich ~
    sc_logtr + sc_fa +
      sc_logtr : sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = pd_parasitism,
  family = genpois(),
  ziformula = ~1
)

pr_design_rs <- update(
  pr_design_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

pr_design_ar <- update(
  pr_design_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

pr_design_disp_ar <- update(
  pr_design_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

pr_repeated_candidates <- list(
  ri = pr_design_ri,
  rs = pr_design_rs,
  ar = pr_design_ar
)

pr_repeated_fit <- compare_candidate_models(pr_repeated_candidates)
pr_repeated_diagnostics <- run_candidate_diagnostics(
  pr_repeated_candidates,
  n_sim = 1000L,
  seed = 501L
)
pr_repeated_effects <- summarize_fixed_terms(
  pr_repeated_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

pr_dispersion_candidates <- list(
  ri = pr_design_ri,
  disp_ar = pr_design_disp_ar
)

pr_dispersion_fit <- compare_candidate_models(pr_dispersion_candidates)
pr_dispersion_diagnostics <- run_candidate_diagnostics(
  pr_dispersion_candidates,
  n_sim = 1000L,
  seed = 511L
)
pr_dispersion_effects <- summarize_fixed_terms(
  pr_dispersion_candidates,
  terms = c("sc_logtr", "sc_fa", "sc_logtr:sc_fa")
)

pr_repeated_fit
pr_repeated_diagnostics$summary
pr_repeated_effects

pr_dispersion_fit
pr_dispersion_diagnostics$summary
pr_dispersion_effects

pr_design_disp_ar2 <- update(
  pr_design_disp_ar,
  . ~ . - (1 | site/plot_id) +
    site + (1 | plot_id)
) # Treat site as fixed to avoid a singular site-level random effect.

## 3.3.4.2 Select primary design specification ##

pr_design <- pr_design_disp_ar2
# Retain dispersion AR(1), which improves AIC and dispersion diagnostics. The
# conditional AR(1) candidate does not converge.

## 3.3.4.3 Nonlinear temporal sensitivities of the selected primary model ##

if (!is.null(pr_design)) {
  pr_design_age2 <- update(
    pr_design,
    . ~ . + I(sc_fa^2)
  )

  pr_design_age2_interaction <- update(
    pr_design_age2,
    . ~ . + sc_logtr:I(sc_fa^2)
  )

  pr_nonlinear_sensitivities <- list(
    primary_linear = pr_design,
    quadratic_age = pr_design_age2,
    quadratic_age_interaction = pr_design_age2_interaction
  )

  pr_nonlinear_fit <- compare_candidate_models(pr_nonlinear_sensitivities)
  pr_nonlinear_diagnostics <- run_candidate_diagnostics(
    pr_nonlinear_sensitivities,
    n_sim = 1000L,
    seed = 521L
  )
  pr_nonlinear_effects <- summarize_fixed_terms(
    pr_nonlinear_sensitivities,
    terms = c(
      "sc_logtr",
      "sc_fa",
      "I(sc_fa^2)",
      "sc_logtr:sc_fa",
      "sc_logtr:I(sc_fa^2)"
    )
  )
}

pr_nonlinear_fit
pr_nonlinear_diagnostics$summary
pr_nonlinear_effects

design_models <- list(
  parasitism = prt_design,
  realized_parasitism = pa_design,
  host_abundance = ha_design,
  host_richness = hr_design,
  parasitoid_richness = pr_design
) # primary model list

# Retain complete term-level VIF tables for the five selected design models.
# These records support publication reporting without recalculating diagnostics
# in a table-rendering script.
design_collinearity <- lapply(design_models, function(model) {
  tryCatch(
    as.data.frame(performance::check_collinearity(model)),
    error = function(e) data.frame(error = conditionMessage(e))
  )
})

analysis_data <- list(
  pd_all = pd_all,
  pd_parasitism = pd_parasitism,
  pd_climate = pd_climate,
  pd_net = pd_net
)

# Machine-readable record of the complete design-model selection workflow.
# Candidate model objects are retained so that additional read-only checks can
# be run without refitting the output-producing analysis.
design_selection_results <- list(
  parasitism = list(
    current_selection = prt_design,
    collinearity = design_collinearity$parasitism,
    sem_equivalent = prt_design_sem,
    repeated_candidates = prt_repeated_candidates,
    repeated_comparison = prt_repeated_fit,
    repeated_diagnostics = prt_repeated_diagnostics$summary,
    repeated_effects = prt_repeated_effects,
    dispersion_candidates = prt_dispersion_candidates,
    dispersion_comparison = prt_dispersion_fit,
    dispersion_diagnostics = prt_dispersion_diagnostics$summary,
    dispersion_effects = prt_dispersion_effects,
    nonlinear_candidates = prt_nonlinear_sensitivities,
    nonlinear_comparison = prt_nonlinear_fit,
    nonlinear_diagnostics = prt_nonlinear_diagnostics$summary,
    nonlinear_effects = prt_nonlinear_effects,
    sensitivities = list(
      presence = prt_design_present,
      mortality = prt_design_mortality
    )
  ),
  host_abundance = list(
    current_selection = ha_design,
    collinearity = design_collinearity$host_abundance,
    repeated_candidates = ha_repeated_candidates,
    repeated_comparison = ha_repeated_fit,
    repeated_diagnostics = ha_repeated_diagnostics$summary,
    repeated_effects = ha_repeated_effects,
    dispersion_candidates = ha_dispersion_candidates,
    dispersion_comparison = ha_dispersion_fit,
    dispersion_diagnostics = ha_dispersion_diagnostics$summary,
    dispersion_effects = ha_dispersion_effects,
    nonlinear_candidates = ha_nonlinear_sensitivities,
    nonlinear_comparison = ha_nonlinear_fit,
    nonlinear_diagnostics = ha_nonlinear_diagnostics$summary,
    nonlinear_effects = ha_nonlinear_effects
  ),
  host_richness = list(
    current_selection = hr_design,
    collinearity = design_collinearity$host_richness,
    repeated_candidates = hr_repeated_candidates,
    repeated_comparison = hr_repeated_fit,
    repeated_diagnostics = hr_repeated_diagnostics$summary,
    repeated_effects = hr_repeated_effects,
    dispersion_candidates = hr_dispersion_candidates,
    dispersion_comparison = hr_dispersion_fit,
    dispersion_diagnostics = hr_dispersion_diagnostics$summary,
    dispersion_effects = hr_dispersion_effects,
    nonlinear_candidates = hr_nonlinear_sensitivities,
    nonlinear_comparison = hr_nonlinear_fit,
    nonlinear_diagnostics = hr_nonlinear_diagnostics$summary,
    nonlinear_effects = hr_nonlinear_effects
  ),
  realized_parasitism = list(
    current_selection = pa_design,
    collinearity = design_collinearity$realized_parasitism,
    repeated_candidates = pa_repeated_candidates,
    repeated_comparison = pa_repeated_fit,
    repeated_diagnostics = pa_repeated_diagnostics$summary,
    repeated_effects = pa_repeated_effects,
    dispersion_candidates = pa_dispersion_candidates,
    dispersion_comparison = pa_dispersion_fit,
    dispersion_diagnostics = pa_dispersion_diagnostics$summary,
    dispersion_effects = pa_dispersion_effects,
    nonlinear_candidates = pa_nonlinear_sensitivities,
    nonlinear_comparison = pa_nonlinear_fit,
    nonlinear_diagnostics = pa_nonlinear_diagnostics$summary,
    nonlinear_effects = pa_nonlinear_effects,
    sensitivities = list(presence = pa_design_present)
  ),
  parasitoid_richness = list(
    current_selection = pr_design,
    collinearity = design_collinearity$parasitoid_richness,
    repeated_candidates = pr_repeated_candidates,
    repeated_comparison = pr_repeated_fit,
    repeated_diagnostics = pr_repeated_diagnostics$summary,
    repeated_effects = pr_repeated_effects,
    dispersion_candidates = pr_dispersion_candidates,
    dispersion_comparison = pr_dispersion_fit,
    dispersion_diagnostics = pr_dispersion_diagnostics$summary,
    dispersion_effects = pr_dispersion_effects,
    nonlinear_candidates = pr_nonlinear_sensitivities,
    nonlinear_comparison = pr_nonlinear_fit,
    nonlinear_diagnostics = pr_nonlinear_diagnostics$summary,
    nonlinear_effects = pr_nonlinear_effects
  )
)

design_model_output_dir <- file.path("model_outputs", "design_models")
dir.create(design_model_output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  design_models,
  file.path(design_model_output_dir, "design_models.rds")
)
saveRDS(
  analysis_data,
  file.path(design_model_output_dir, "analysis_data.rds")
)
saveRDS(
  design_selection_results,
  file.path(design_model_output_dir, "design_selection_results.rds")
)

message(
  "Selected design models, selection results, and prepared analysis data exported to: ",
  design_model_output_dir
)

message("Script ran successfully ✔️")
