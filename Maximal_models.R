## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Maximal models
# Author: Massimo Martini

# Script description: This script fits and validates the five maximal
# GLMMs of hypothesized conditional associations used in the pathway analyses.
# It adds the prespecified forest and community predictors to the approved
# design-model structures. The script records the nuisance-structure review,
# checks the final models, tests the declared predictor and interaction blocks,
# and runs the prespecified sensitivities. It saves the final models together
# with their validation and inference results.

# Run this script once from the project root after Design_models.R.
# It first records the prespecified structure-selection audit,
# then applies the approved structures and completes validation,
# inference, and sensitivity analyses.

# Reproducibility: run with R 4.6.1 and the package versions in renv.lock.

# Approximate runtime: 18-20 minutes

structure_selection_n_sim <- 500L
final_validation_n_sim <- 1000L
sensitivity_validation_n_sim <- 500L

if (any(c(
  structure_selection_n_sim,
  final_validation_n_sim,
  sensitivity_validation_n_sim
) > 1000L)) {
  stop("Residual simulations must not exceed 1000 per model.")
}

# ------------------------------------------------------------------------------
# 1. Packages and helper functions ####
# ------------------------------------------------------------------------------

library(dplyr)
library(glmmTMB)
library(tibble)
library(tidyr)

source("analysis_functions.R")

# ------------------------------------------------------------------------------
# 2. Read prepared analysis data and design models ####
# ------------------------------------------------------------------------------

model_output_dir <- "model_outputs"
design_model_output_dir <- file.path(model_output_dir, "design_models")
maximal_model_output_dir <- file.path(model_output_dir, "maximal_models")
dir.create(maximal_model_output_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  analysis_data = file.path(design_model_output_dir, "analysis_data.rds"),
  design_models = file.path(design_model_output_dir, "design_models.rds")
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input file(s): ",
    paste(basename(missing_inputs), collapse = ", "),
    ". Run Design_models.R first."
  )
}

analysis_data <- readRDS(input_paths[["analysis_data"]])
design_models <- readRDS(input_paths[["design_models"]])

required_datasets <- c("pd_all", "pd_parasitism", "pd_climate", "pd_net")
required_design_models <- c(
  "parasitism",
  "realized_parasitism",
  "host_abundance",
  "host_richness",
  "parasitoid_richness"
)

if (!all(required_datasets %in% names(analysis_data))) {
  stop("analysis_data.rds does not contain all required analysis datasets.")
}

if (!all(required_design_models %in% names(design_models))) {
  stop("design_models.rds does not contain all five authoritative models.")
}

pd_all <- analysis_data$pd_all
pd_parasitism <- analysis_data$pd_parasitism
pd_climate <- analysis_data$pd_climate
pd_net <- analysis_data$pd_net
rm(analysis_data)

# Every candidate for a response is fitted to one explicit complete-case
# dataset. This prevents nuisance-structure comparisons from changing samples.
prt_maximal_data <- pd_parasitism |>
  tidyr::drop_na(
    cell_prate, par_abund, p_failure, total_cells, sc_cells,
    sc_logtr, sc_fa, sc_rtfd, sc_rtsv, sc_rtcells, sc_hr, sc_pr10,
    site, plot_id, year
  )

pa_maximal_data <- pd_parasitism |>
  tidyr::drop_na(
    par_abund, sc_logtr, sc_fa, sc_rtfd, sc_rtsv, sc_rtcells, sc_hr,
    sc_pr10, site, plot_id, year
  )

pr_maximal_data <- pd_parasitism |>
  tidyr::drop_na(
    par_rich, sc_logtr, sc_fa, sc_rtfd, sc_rtsv, sc_rtcells, sc_hr,
    site, plot_id, year
  )

hr_maximal_data <- pd_all |>
  tidyr::drop_na(
    host_rich, sc_logtr, sc_fa, sc_rtfd, sc_rtsv, sc_slope, sc_elev,
    sc_east, sc_north, sc_rtcells, site, plot_id, year
  )

ha_maximal_data <- pd_all |>
  tidyr::drop_na(
    total_cells, sc_logtr, sc_fa, sc_rtfd, sc_rtsv, sc_slope, sc_elev,
    sc_east, sc_north, site, plot_id, year
  )

if (any(vapply(
  list(
    parasitism = prt_maximal_data,
    realized_parasitism = pa_maximal_data,
    parasitoid_richness = pr_maximal_data,
    host_richness = hr_maximal_data,
    host_abundance = ha_maximal_data
  ),
  nrow,
  integer(1)
) == 0L)) {
  stop("At least one maximal-model dataset has no complete observations.")
}

# ------------------------------------------------------------------------------
# 3. Maximal fixed effects and nuisance-structure candidates ####
# ------------------------------------------------------------------------------

# Candidate diagnostics use 500 simulations per model. The selected models
# receive 1,000-simulation validation later in this script.

# ----------------------------
## 3.1 Parasitism ####
# ----------------------------

prt_design_cbind <- design_models$parasitism

# The weighted-rate response is numerically equivalent to the cbind response
# and is retained because piecewiseSEM requires a scalar endogenous response.
prt_maximal_ri <- update(
  prt_design_cbind,
  cell_prate ~
    sc_logtr + sc_fa +
      sc_rtfd + sc_rtsv +
      sc_rtcells + sc_hr + sc_pr10 +
      sc_logtr:sc_fa +
      sc_rtsv:sc_fa +
      sc_rtcells:sc_fa +
      sc_hr:sc_fa +
      sc_pr10:sc_fa +
      sc_pr10:sc_rtsv +
      sc_hr:sc_rtcells +
      (1 | site/plot_id) + (1 | year),
  weights = total_cells,
  data = prt_maximal_data,
  dispformula = ~ sc_cells
)

prt_maximal_rs <- update(
  prt_maximal_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

prt_maximal_cond_ar <- update(
  prt_maximal_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

prt_maximal_disp_ar <- update(
  prt_maximal_ri,
  dispformula = ~
    sc_cells + ar1(year + 0 | site:plot_id)
)

prt_maximal_candidates <- list(
  random_intercept = prt_maximal_ri,
  random_slope = prt_maximal_rs,
  conditional_ar1 = prt_maximal_cond_ar,
  dispersion_ar1 = prt_maximal_disp_ar
)

# Compile fit, residual, temporal, singularity, and coefficient audits.
prt_maximal_structure <- evaluate_maximal_structure_candidates(
  models = prt_maximal_candidates,
  seed = 701L,
  n_sim = structure_selection_n_sim,
  temporal_data = prt_maximal_data
)

# ----------------------------------
## 3.2 Realized parasitism (parasitoid abundance) ####
# ----------------------------------

pa_design <- design_models$realized_parasitism

pa_maximal_ri <- update(
  pa_design,
  par_abund ~
    sc_logtr + sc_fa +
      sc_rtfd + sc_rtsv +
      sc_rtcells + sc_hr + sc_pr10 +
      sc_logtr:sc_fa +
      sc_rtsv:sc_fa +
      sc_rtcells:sc_fa +
      sc_hr:sc_fa +
      sc_pr10:sc_fa +
      sc_pr10:sc_rtsv +
      sc_hr:sc_rtcells +
      site + (1 | plot_id) + (1 | year),
  data = pa_maximal_data,
  ziformula = ~ 0,
  dispformula = ~ 1
)

pa_maximal_rs <- update(
  pa_maximal_ri,
  . ~ . - (1 | plot_id) + diag(1 + sc_fa | plot_id)
)

pa_maximal_cond_ar <- update(
  pa_maximal_ri,
  . ~ . + ar1(year + 0 | plot_id)
)

pa_maximal_disp_ar <- update(
  pa_maximal_ri,
  dispformula = ~ ar1(year + 0 | plot_id)
)

pa_maximal_candidates <- list(
  random_intercept = pa_maximal_ri,
  random_slope = pa_maximal_rs,
  conditional_ar1 = pa_maximal_cond_ar,
  dispersion_ar1 = pa_maximal_disp_ar
)

pa_maximal_structure <- evaluate_maximal_structure_candidates(
  models = pa_maximal_candidates,
  seed = 711L,
  n_sim = structure_selection_n_sim,
  temporal_data = pa_maximal_data
)

# ----------------------------------
## 3.3 Parasitoid species richness ####
# ----------------------------------

pr_design <- design_models$parasitoid_richness

# The maximal model omits zero inflation following the completed zero-process
# assessment. The no-plot, constant-dispersion model is the simplest candidate.
pr_maximal_no_plot <- update(
  pr_design,
  par_rich ~
    sc_logtr + sc_fa +
      sc_rtfd + sc_rtsv +
      sc_rtcells + sc_hr +
      sc_logtr:sc_fa +
      sc_rtsv:sc_fa +
      sc_rtcells:sc_fa +
      sc_hr:sc_fa +
      sc_hr:sc_rtcells +
      site + (1 | year),
  data = pr_maximal_data,
  ziformula = ~ 0,
  dispformula = ~ 1
)

pr_maximal_ri <- update(
  pr_maximal_no_plot,
  . ~ . + (1 | plot_id)
)

pr_maximal_rs <- update(
  pr_maximal_no_plot,
  . ~ . + diag(1 + sc_fa | plot_id)
)

pr_maximal_cond_ar <- update(
  pr_maximal_no_plot,
  . ~ . + ar1(year + 0 | plot_id)
)

pr_maximal_disp_ar <- update(
  pr_maximal_no_plot,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

# Retain the original maximal formulation with both the plot intercept and
# dispersion AR(1), so no previously fitted structure is lost from review.
pr_maximal_ri_disp_ar <- update(
  pr_maximal_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

pr_maximal_candidates <- list(
  no_plot = pr_maximal_no_plot,
  random_intercept = pr_maximal_ri,
  random_slope = pr_maximal_rs,
  conditional_ar1 = pr_maximal_cond_ar,
  dispersion_ar1 = pr_maximal_disp_ar,
  random_intercept_dispersion_ar1 = pr_maximal_ri_disp_ar
)

pr_maximal_structure <- evaluate_maximal_structure_candidates(
  models = pr_maximal_candidates,
  seed = 721L,
  n_sim = structure_selection_n_sim,
  temporal_data = pr_maximal_data
)

# ----------------------------------
## 3.4 Host species richness ####
# ----------------------------------

hr_design <- design_models$host_richness

hr_maximal_linear_ri <- update(
  hr_design,
  host_rich ~
    sc_logtr + sc_fa +
      sc_rtfd + sc_rtsv +
      sc_slope + sc_elev + sc_east + sc_north +
      sc_rtcells +
      sc_logtr:sc_fa +
      sc_rtsv:sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = hr_maximal_data,
  dispformula = ~ 1
)

hr_maximal_linear_rs <- update(
  hr_maximal_linear_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

hr_maximal_linear_cond_ar <- update(
  hr_maximal_linear_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

hr_maximal_linear_disp_ar <- update(
  hr_maximal_linear_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

hr_maximal_quadratic_ri <- update(
  hr_maximal_linear_ri,
  . ~ . + I(sc_rtcells^2)
)

hr_maximal_quadratic_rs <- update(
  hr_maximal_quadratic_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

hr_maximal_quadratic_cond_ar <- update(
  hr_maximal_quadratic_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

hr_maximal_quadratic_disp_ar <- update(
  hr_maximal_linear_disp_ar,
  . ~ . + I(sc_rtcells^2)
)

hr_maximal_quadratic_no_plot <- update(
  hr_maximal_quadratic_ri,
  . ~ . - (1 | site/plot_id) + (1 | site)
)

hr_maximal_quadratic_no_plot_disp_ar <- update(
  hr_maximal_quadratic_no_plot,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

# This combined set retains the previous temporal candidates and completes the
# required linear/quadratic brood-cell by constant/AR(1)-dispersion comparison.
hr_maximal_candidates <- list(
  linear_random_intercept = hr_maximal_linear_ri,
  linear_random_slope = hr_maximal_linear_rs,
  linear_conditional_ar1 = hr_maximal_linear_cond_ar,
  linear_dispersion_ar1 = hr_maximal_linear_disp_ar,
  quadratic_random_intercept = hr_maximal_quadratic_ri,
  quadratic_random_slope = hr_maximal_quadratic_rs,
  quadratic_conditional_ar1 = hr_maximal_quadratic_cond_ar,
  quadratic_dispersion_ar1 = hr_maximal_quadratic_disp_ar,
  quadratic_no_plot = hr_maximal_quadratic_no_plot,
  quadratic_no_plot_dispersion_ar1 = hr_maximal_quadratic_no_plot_disp_ar
)

hr_maximal_structure <- evaluate_maximal_structure_candidates(
  models = hr_maximal_candidates,
  seed = 731L,
  n_sim = structure_selection_n_sim,
  temporal_data = hr_maximal_data
)

# ----------------------------------
## 3.5 Host abundance ####
# ----------------------------------

ha_design <- design_models$host_abundance

ha_maximal_ri <- update(
  ha_design,
  total_cells ~
    sc_logtr + sc_fa +
      sc_rtfd + sc_rtsv +
      sc_slope + sc_elev + sc_east + sc_north +
      sc_logtr:sc_fa +
      sc_rtsv:sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = ha_maximal_data,
  dispformula = ~ 1
)

ha_maximal_rs <- update(
  ha_maximal_ri,
  . ~ . - (1 | site/plot_id) +
    (1 | site) + diag(1 + sc_fa | site:plot_id)
)

ha_maximal_cond_ar <- update(
  ha_maximal_ri,
  . ~ . + ar1(year + 0 | site:plot_id)
)

ha_maximal_disp_ar <- update(
  ha_maximal_ri,
  dispformula = ~ ar1(year + 0 | site:plot_id)
)

ha_maximal_candidates <- list(
  random_intercept = ha_maximal_ri,
  random_slope = ha_maximal_rs,
  conditional_ar1 = ha_maximal_cond_ar,
  dispersion_ar1 = ha_maximal_disp_ar
)

ha_maximal_structure <- evaluate_maximal_structure_candidates(
  models = ha_maximal_candidates,
  seed = 741L,
  n_sim = structure_selection_n_sim,
  temporal_data = ha_maximal_data
)

# ------------------------------------------------------------------------------
# 4. Export structure-selection review ####
# ------------------------------------------------------------------------------

maximal_structure_comparisons <- list(
  parasitism = prt_maximal_structure,
  realized_parasitism = pa_maximal_structure,
  parasitoid_richness = pr_maximal_structure,
  host_richness = hr_maximal_structure,
  host_abundance = ha_maximal_structure
)

structure_selection_overview <- dplyr::bind_rows(
  dplyr::left_join(
    dplyr::left_join(
      prt_maximal_structure$comparison,
      prt_maximal_structure$diagnostics$summary,
      by = "model"
    ),
    prt_maximal_structure$temporal_summary,
    by = "model"
  ) |>
    dplyr::mutate(response = "parasitism", .before = 1),
  dplyr::left_join(
    dplyr::left_join(
      pa_maximal_structure$comparison,
      pa_maximal_structure$diagnostics$summary,
      by = "model"
    ),
    pa_maximal_structure$temporal_summary,
    by = "model"
  ) |>
    dplyr::mutate(response = "realized_parasitism", .before = 1),
  dplyr::left_join(
    dplyr::left_join(
      pr_maximal_structure$comparison,
      pr_maximal_structure$diagnostics$summary,
      by = "model"
    ),
    pr_maximal_structure$temporal_summary,
    by = "model"
  ) |>
    dplyr::mutate(response = "parasitoid_richness", .before = 1),
  dplyr::left_join(
    dplyr::left_join(
      hr_maximal_structure$comparison,
      hr_maximal_structure$diagnostics$summary,
      by = "model"
    ),
    hr_maximal_structure$temporal_summary,
    by = "model"
  ) |>
    dplyr::mutate(response = "host_richness", .before = 1),
  dplyr::left_join(
    dplyr::left_join(
      ha_maximal_structure$comparison,
      ha_maximal_structure$diagnostics$summary,
      by = "model"
    ),
    ha_maximal_structure$temporal_summary,
    by = "model"
  ) |>
    dplyr::mutate(response = "host_abundance", .before = 1)
)

# Record the analyses that follow the structure-selection stage.
deferred_final_analysis <- list(
  model_stages = c(
    "design-to-forest comparisons",
    "forest-model drop-one tests",
    "nuisance-aligned forest-to-maximal comparisons",
    "strict optimizer refits",
    "final maximal-model validation",
    "single-interaction likelihood-ratio tests",
    "joint pathway-block likelihood-ratio tests",
    "conditional slopes and predictions"
  ),
  sensitivities = list(
    prt_maximal_mortality = c(
      "p_failure2 response failures",
      "sc_rtovcells abundance predictor",
      "sc_overall_cells dispersion predictor",
      "sc_hr:sc_rtovcells host-richness interaction"
    ),
    prt_maximal_present = "par_rich > 0 subset",
    pa_maximal_mortality = c(
      "replace sc_rtcells with sc_rtovcells",
      "replace sc_hr:sc_rtcells with sc_hr:sc_rtovcells"
    ),
    pa_maximal_present = c(
      "par_rich > 0 subset",
      "truncated_nbinom1 family"
    ),
    microclimate = "matched baseline, additive climate, and climate-by-age models",
    comparison = "full-data versus parasitoid-present estimates and predictions"
  ),
  exports = c(
    "selected scalar and cbind parasitism components",
    "five authoritative mechanistic models",
    "final validation records",
    "nested and pathway-block tests",
    "validated sensitivity models"
  )
)

structure_selection_results <- list(
  phase = "structure_selection",
  residual_simulations_per_model = structure_selection_n_sim,
  sample_accounting = tibble::tibble(
    response = c(
      "parasitism",
      "realized_parasitism",
      "parasitoid_richness",
      "host_richness",
      "host_abundance"
    ),
    eligible_rows = c(
      nrow(pd_parasitism),
      nrow(pd_parasitism),
      nrow(pd_parasitism),
      nrow(pd_all),
      nrow(pd_all)
    ),
    complete_rows = c(
      nrow(prt_maximal_data),
      nrow(pa_maximal_data),
      nrow(pr_maximal_data),
      nrow(hr_maximal_data),
      nrow(ha_maximal_data)
    )
  ) |>
    dplyr::mutate(excluded_incomplete = eligible_rows - complete_rows),
  overview = structure_selection_overview,
  comparisons = maximal_structure_comparisons,
  deferred_final_analysis = deferred_final_analysis
)

structure_selection_path <- file.path(
  maximal_model_output_dir,
  "maximal_structure_selection.rds"
)

saveRDS(structure_selection_results, structure_selection_path)

print(structure_selection_overview)
message("Structure-selection review exported to: ", structure_selection_path)
message("Structure-selection stage complete; continuing to final analysis.")

# ------------------------------------------------------------------------------
# 3. Read approved structures and strictly refit final models ####
# ------------------------------------------------------------------------------

structure_selection_path <- file.path(
  maximal_model_output_dir,
  "maximal_structure_selection.rds"
)

if (!file.exists(structure_selection_path)) {
  stop(
    "The structure-selection artifact was not created successfully: ",
    structure_selection_path
  )
}

structure_selection_results <- readRDS(structure_selection_path)
if (!identical(structure_selection_results$phase, "structure_selection")) {
  stop("The structure-selection artifact has an unexpected analysis phase.")
}

approved_structure_names <- c(
  parasitism = "random_intercept",
  realized_parasitism = "random_intercept",
  parasitoid_richness = "no_plot",
  host_richness = "quadratic_no_plot_dispersion_ar1",
  host_abundance = "dispersion_ar1"
)

selection_comparisons <- structure_selection_results$comparisons
if (!all(names(approved_structure_names) %in% names(selection_comparisons))) {
  stop("The structure-selection artifact does not contain all five responses.")
}

for (response_name in names(approved_structure_names)) {
  candidate_name <- approved_structure_names[[response_name]]
  if (!candidate_name %in%
      names(selection_comparisons[[response_name]]$candidates)) {
    stop(
      "Approved candidate `", candidate_name, "` is missing for `",
      response_name, "`."
    )
  }
}

selected_candidates <- list(
  parasitism = selection_comparisons$parasitism$candidates[[
    approved_structure_names[["parasitism"]]
  ]],
  realized_parasitism =
    selection_comparisons$realized_parasitism$candidates[[
      approved_structure_names[["realized_parasitism"]]
    ]],
  parasitoid_richness =
    selection_comparisons$parasitoid_richness$candidates[[
      approved_structure_names[["parasitoid_richness"]]
    ]],
  host_richness = selection_comparisons$host_richness$candidates[[
    approved_structure_names[["host_richness"]]
  ]],
  host_abundance = selection_comparisons$host_abundance$candidates[[
    approved_structure_names[["host_abundance"]]
  ]]
)

strict_control <- glmmTMB::glmmTMBControl(
  optCtrl = list(iter.max = 10000L, eval.max = 10000L)
)

prt_mechanistic_sem <- update(
  selected_candidates$parasitism,
  control = strict_control
)
pa_mechanistic <- update(
  selected_candidates$realized_parasitism,
  control = strict_control
)
pr_mechanistic <- update(
  selected_candidates$parasitoid_richness,
  control = strict_control
)
hr_mechanistic <- update(
  selected_candidates$host_richness,
  control = strict_control
)
ha_mechanistic <- update(
  selected_candidates$host_abundance,
  control = strict_control
)

# Retain the transparent success/failure encoding alongside the verified-
# equivalent scalar weighted-rate model needed by piecewiseSEM.
prt_design_cbind <- design_models$parasitism
prt_cbind_formula <- stats::update.formula(
  stats::formula(prt_mechanistic_sem),
  cbind(par_abund, p_failure) ~ .
)
prt_mechanistic_cbind <- update(
  prt_design_cbind,
  prt_cbind_formula,
  data = prt_maximal_data,
  dispformula = stats::formula(prt_mechanistic_sem, component = "disp"),
  control = strict_control
)

prt_encoding_equivalence <- tibble::tibble(
  metric = c(
    "logLik",
    "AIC",
    "maximum conditional fixed-effect difference",
    "maximum dispersion fixed-effect difference",
    "maximum fitted-value difference"
  ),
  difference = c(
    as.numeric(logLik(prt_mechanistic_sem) - logLik(prt_mechanistic_cbind)),
    AIC(prt_mechanistic_sem) - AIC(prt_mechanistic_cbind),
    max(abs(
      glmmTMB::fixef(prt_mechanistic_sem)$cond -
        glmmTMB::fixef(prt_mechanistic_cbind)$cond
    )),
    max(abs(
      glmmTMB::fixef(prt_mechanistic_sem)$disp -
        glmmTMB::fixef(prt_mechanistic_cbind)$disp
    )),
    max(abs(fitted(prt_mechanistic_sem) - fitted(prt_mechanistic_cbind)))
  ),
  tolerance = c(1e-6, 1e-6, 1e-4, 1e-4, 1e-4)
) |>
  dplyr::mutate(
    within_tolerance = abs(difference) <= tolerance
  )

if (any(!prt_encoding_equivalence$within_tolerance)) {
  print(prt_encoding_equivalence)
  stop(
    "The scalar and cbind parasitism encodings differ beyond the declared ",
    "numerical tolerances."
  )
}

# Compare each approved candidate with its strict optimizer refit.
optimizer_stability <- list(
  parasitism = compare_candidate_models(list(
    selection_fit = selected_candidates$parasitism,
    strict_refit = prt_mechanistic_sem
  )),
  realized_parasitism = compare_candidate_models(list(
    selection_fit = selected_candidates$realized_parasitism,
    strict_refit = pa_mechanistic
  )),
  parasitoid_richness = compare_candidate_models(list(
    selection_fit = selected_candidates$parasitoid_richness,
    strict_refit = pr_mechanistic
  )),
  host_richness = compare_candidate_models(list(
    selection_fit = selected_candidates$host_richness,
    strict_refit = hr_mechanistic
  )),
  host_abundance = compare_candidate_models(list(
    selection_fit = selected_candidates$host_abundance,
    strict_refit = ha_mechanistic
  ))
)

# Check that fixed-effect estimates remain stable after the strict refit.
optimizer_fixed_effects <- list(
  parasitism = summarize_fixed_terms(
    list(selection_fit = selected_candidates$parasitism,
         strict_refit = prt_mechanistic_sem),
    names(glmmTMB::fixef(prt_mechanistic_sem)$cond)
  ),
  realized_parasitism = summarize_fixed_terms(
    list(selection_fit = selected_candidates$realized_parasitism,
         strict_refit = pa_mechanistic),
    names(glmmTMB::fixef(pa_mechanistic)$cond)
  ),
  parasitoid_richness = summarize_fixed_terms(
    list(selection_fit = selected_candidates$parasitoid_richness,
         strict_refit = pr_mechanistic),
    names(glmmTMB::fixef(pr_mechanistic)$cond)
  ),
  host_richness = summarize_fixed_terms(
    list(selection_fit = selected_candidates$host_richness,
         strict_refit = hr_mechanistic),
    names(glmmTMB::fixef(hr_mechanistic)$cond)
  ),
  host_abundance = summarize_fixed_terms(
    list(selection_fit = selected_candidates$host_abundance,
         strict_refit = ha_mechanistic),
    names(glmmTMB::fixef(ha_mechanistic)$cond)
  )
)

# ------------------------------------------------------------------------------
# 4. Final validation of selected maximal models ####
# ------------------------------------------------------------------------------

# Run the complete numerical and residual validation for each selected model.
maximal_validations <- list(
  parasitism = validate_maximal_model(
    prt_mechanistic_sem,
    model_name = "prt_mechanistic_sem",
    n_sim = final_validation_n_sim,
    seed = 801L,
    temporal_data = prt_maximal_data
  ),
  realized_parasitism = validate_maximal_model(
    pa_mechanistic,
    model_name = "pa_mechanistic",
    n_sim = final_validation_n_sim,
    seed = 811L,
    temporal_data = pa_maximal_data
  ),
  parasitoid_richness = validate_maximal_model(
    pr_mechanistic,
    model_name = "pr_mechanistic",
    n_sim = final_validation_n_sim,
    seed = 821L,
    temporal_data = pr_maximal_data
  ),
  host_richness = validate_maximal_model(
    hr_mechanistic,
    model_name = "hr_mechanistic",
    n_sim = final_validation_n_sim,
    seed = 831L,
    temporal_data = hr_maximal_data
  ),
  host_abundance = validate_maximal_model(
    ha_mechanistic,
    model_name = "ha_mechanistic",
    n_sim = final_validation_n_sim,
    seed = 841L,
    temporal_data = ha_maximal_data
  )
)

# ------------------------------------------------------------------------------
# 5. Forest-structure and community-block comparisons ####
# ------------------------------------------------------------------------------

# Design-aligned forest models test the forest-structure block without changing
# the nuisance structure selected for the authoritative design model.
prt_design_reference <- update(
  prt_design_cbind,
  cell_prate ~ .,
  weights = total_cells,
  data = prt_maximal_data,
  control = strict_control
)
prt_forest_structure <- update(
  prt_design_reference,
  . ~ . + sc_rtfd + sc_rtsv + sc_rtsv:sc_fa
)

pa_design_reference <- update(
  design_models$realized_parasitism,
  data = pa_maximal_data,
  control = strict_control
)
pa_forest_structure <- update(
  pa_design_reference,
  . ~ . + sc_rtfd + sc_rtsv + sc_rtsv:sc_fa
)

pr_design_reference <- update(
  design_models$parasitoid_richness,
  data = pr_maximal_data,
  control = strict_control
)
pr_forest_structure <- update(
  pr_design_reference,
  . ~ . + sc_rtfd + sc_rtsv + sc_rtsv:sc_fa
)

hr_design_reference <- update(
  design_models$host_richness,
  data = hr_maximal_data,
  control = strict_control
)
hr_forest_structure <- update(
  hr_design_reference,
  . ~ . +
    sc_rtfd + sc_rtsv +
    sc_slope + sc_elev + sc_east + sc_north +
    sc_rtsv:sc_fa
)

ha_design_reference <- update(
  design_models$host_abundance,
  data = ha_maximal_data,
  control = strict_control
)
ha_forest_structure <- update(
  ha_design_reference,
  . ~ . +
    sc_rtfd + sc_rtsv +
    sc_slope + sc_elev + sc_east + sc_north +
    sc_rtsv:sc_fa
)

stopifnot(
  nobs(prt_design_reference) == nobs(prt_forest_structure),
  nobs(pa_design_reference) == nobs(pa_forest_structure),
  nobs(pr_design_reference) == nobs(pr_forest_structure),
  nobs(hr_design_reference) == nobs(hr_forest_structure),
  nobs(ha_design_reference) == nobs(ha_forest_structure)
)

forest_block_comparisons <- list(
  parasitism = anova(
    prt_design_reference, prt_forest_structure, test = "Chisq"
  ),
  realized_parasitism = anova(
    pa_design_reference, pa_forest_structure, test = "Chisq"
  ),
  parasitoid_richness = anova(
    pr_design_reference, pr_forest_structure, test = "Chisq"
  ),
  host_richness = anova(
    hr_design_reference, hr_forest_structure, test = "Chisq"
  ),
  host_abundance = anova(
    ha_design_reference, ha_forest_structure, test = "Chisq"
  )
)

forest_drop1 <- list(
  parasitism = drop1(prt_forest_structure, test = "Chisq"),
  realized_parasitism = drop1(pa_forest_structure, test = "Chisq"),
  parasitoid_richness = drop1(pr_forest_structure, test = "Chisq"),
  host_richness = drop1(hr_forest_structure, test = "Chisq"),
  host_abundance = drop1(ha_forest_structure, test = "Chisq")
)

# These reduced references inherit the selected maximal nuisance structures, so
# each comparison isolates only the added community block.
prt_mechanistic_forest_reference <- update(
  prt_mechanistic_sem,
  . ~ . -
    sc_rtcells - sc_hr - sc_pr10 -
    sc_rtcells:sc_fa - sc_hr:sc_fa - sc_pr10:sc_fa -
    sc_pr10:sc_rtsv - sc_hr:sc_rtcells
)

pa_mechanistic_forest_reference <- update(
  pa_mechanistic,
  . ~ . -
    sc_rtcells - sc_hr - sc_pr10 -
    sc_rtcells:sc_fa - sc_hr:sc_fa - sc_pr10:sc_fa -
    sc_pr10:sc_rtsv - sc_hr:sc_rtcells
)

pr_mechanistic_forest_reference <- update(
  pr_mechanistic,
  . ~ . -
    sc_rtcells - sc_hr -
    sc_rtcells:sc_fa - sc_hr:sc_fa - sc_hr:sc_rtcells
)

hr_mechanistic_forest_reference <- update(
  hr_mechanistic,
  . ~ . - sc_rtcells - I(sc_rtcells^2)
)

stopifnot(
  nobs(prt_mechanistic_forest_reference) == nobs(prt_mechanistic_sem),
  nobs(pa_mechanistic_forest_reference) == nobs(pa_mechanistic),
  nobs(pr_mechanistic_forest_reference) == nobs(pr_mechanistic),
  nobs(hr_mechanistic_forest_reference) == nobs(hr_mechanistic)
)

community_block_comparisons <- list(
  parasitism = anova(
    prt_mechanistic_forest_reference,
    prt_mechanistic_sem,
    test = "Chisq"
  ),
  realized_parasitism = anova(
    pa_mechanistic_forest_reference,
    pa_mechanistic,
    test = "Chisq"
  ),
  parasitoid_richness = anova(
    pr_mechanistic_forest_reference,
    pr_mechanistic,
    test = "Chisq"
  ),
  host_richness = anova(
    hr_mechanistic_forest_reference,
    hr_mechanistic,
    test = "Chisq"
  )
)

# ------------------------------------------------------------------------------
# 6. Final interaction and pathway-block tests ####
# ------------------------------------------------------------------------------

maximal_drop1 <- list(
  parasitism = drop1(prt_mechanistic_sem, test = "Chisq"),
  realized_parasitism = drop1(pa_mechanistic, test = "Chisq"),
  parasitoid_richness = drop1(pr_mechanistic, test = "Chisq"),
  host_richness = drop1(hr_mechanistic, test = "Chisq"),
  host_abundance = drop1(ha_mechanistic, test = "Chisq")
)

# Remove each predictor with all related interactions and compare by LRT.
interaction_blocks <- list(
  parasitism = test_interaction_blocks(
    prt_mechanistic_sem,
    exclude = c("sc_logtr", "sc_fa")
  ),
  realized_parasitism = test_interaction_blocks(
    pa_mechanistic,
    exclude = c("sc_logtr", "sc_fa")
  ),
  parasitoid_richness = test_interaction_blocks(
    pr_mechanistic,
    exclude = c("sc_logtr", "sc_fa")
  ),
  host_richness = test_interaction_blocks(
    hr_mechanistic,
    exclude = c("sc_logtr", "sc_fa")
  ),
  host_abundance = test_interaction_blocks(
    ha_mechanistic,
    exclude = c("sc_logtr", "sc_fa")
  )
)

# ------------------------------------------------------------------------------
# 7. Declared sensitivities ####
# ------------------------------------------------------------------------------

prt_mortality_data <- prt_maximal_data |>
  tidyr::drop_na(
    p_failure2, overall_cells, sc_overall_cells, sc_rtovcells
  ) |>
  dplyr::mutate(
    cell_prate_mortality = dplyr::if_else(
      overall_cells > 0,
      par_abund / overall_cells,
      NA_real_
    )
  )

pa_mortality_data <- pa_maximal_data |>
  tidyr::drop_na(sc_rtovcells)

prt_present_data <- prt_maximal_data |>
  dplyr::filter(par_rich > 0)
pa_present_data <- pa_maximal_data |>
  dplyr::filter(par_rich > 0)

prt_mechanistic_mortality_cbind <- update(
  prt_mechanistic_cbind,
  cbind(par_abund, p_failure2) ~
    . -
    sc_rtcells - sc_rtcells:sc_fa - sc_hr:sc_rtcells +
    sc_rtovcells + sc_rtovcells:sc_fa + sc_hr:sc_rtovcells,
  data = prt_mortality_data,
  dispformula = ~ sc_overall_cells,
  control = strict_control
)

prt_mechanistic_mortality_sem <- update(
  prt_mechanistic_mortality_cbind,
  cell_prate_mortality ~ .,
  weights = overall_cells,
  data = prt_mortality_data,
  control = strict_control
)

prt_mechanistic_present <- update(
  prt_mechanistic_sem,
  data = prt_present_data,
  control = strict_control
)

pa_mechanistic_mortality <- update(
  pa_mechanistic,
  . ~ . -
    sc_rtcells - sc_rtcells:sc_fa - sc_hr:sc_rtcells +
    sc_rtovcells + sc_rtovcells:sc_fa + sc_hr:sc_rtovcells,
  data = pa_mortality_data,
  control = strict_control
)

pa_mechanistic_present <- update(
  pa_mechanistic,
  data = pa_present_data,
  family = truncated_nbinom1(),
  control = strict_control
)

sensitivity_validations <- list(
  parasitism_mortality = validate_maximal_model(
    prt_mechanistic_mortality_sem,
    model_name = "prt_mechanistic_mortality_sem",
    n_sim = sensitivity_validation_n_sim,
    seed = 901L,
    temporal_data = prt_mortality_data
  ),
  parasitism_present = validate_maximal_model(
    prt_mechanistic_present,
    model_name = "prt_mechanistic_present",
    n_sim = sensitivity_validation_n_sim,
    seed = 911L,
    temporal_data = prt_present_data
  ),
  realized_parasitism_mortality = validate_maximal_model(
    pa_mechanistic_mortality,
    model_name = "pa_mechanistic_mortality",
    n_sim = sensitivity_validation_n_sim,
    seed = 921L,
    temporal_data = pa_mortality_data
  ),
  realized_parasitism_present = validate_maximal_model(
    pa_mechanistic_present,
    model_name = "pa_mechanistic_present",
    n_sim = sensitivity_validation_n_sim,
    seed = 931L,
    temporal_data = pa_present_data
  )
)

presence_drop1 <- list(
  parasitism = drop1(prt_mechanistic_present, test = "Chisq"),
  realized_parasitism = drop1(pa_mechanistic_present, test = "Chisq")
)

presence_interaction_blocks <- list(
  parasitism = test_interaction_blocks(
    prt_mechanistic_present,
    exclude = c("sc_logtr", "sc_fa")
  ),
  realized_parasitism = test_interaction_blocks(
    pa_mechanistic_present,
    exclude = c("sc_logtr", "sc_fa")
  )
)

sensitivity_effects <- list(
  parasitism_mortality = summarize_fixed_terms(
    list(
      original = prt_mechanistic_cbind,
      confirmed_mortality = prt_mechanistic_mortality_cbind
    ),
    unique(c(
      names(glmmTMB::fixef(prt_mechanistic_cbind)$cond),
      names(glmmTMB::fixef(prt_mechanistic_mortality_cbind)$cond)
    ))
  ),
  parasitism_presence = summarize_fixed_terms(
    list(
      all_data = prt_mechanistic_sem,
      parasitoid_present = prt_mechanistic_present
    ),
    names(glmmTMB::fixef(prt_mechanistic_sem)$cond)
  ),
  realized_parasitism_mortality = summarize_fixed_terms(
    list(
      original = pa_mechanistic,
      confirmed_mortality = pa_mechanistic_mortality
    ),
    unique(c(
      names(glmmTMB::fixef(pa_mechanistic)$cond),
      names(glmmTMB::fixef(pa_mechanistic_mortality)$cond)
    ))
  ),
  realized_parasitism_presence = summarize_fixed_terms(
    list(
      all_data = pa_mechanistic,
      parasitoid_present = pa_mechanistic_present
    ),
    names(glmmTMB::fixef(pa_mechanistic)$cond)
  )
)

# ------------------------------------------------------------------------------
# 7.1 Microclimate sensitivity ####
# ------------------------------------------------------------------------------

# Refit baselines on identical rows before adding climate;
# retain all original scaling and model structures.
# sc_temp and sc_humid were standardized once with as.numeric(scale(...))
# in Design_models.R; retain those values across all response subsets.
microclimate_data <- pd_all |>
  filter(is.finite(annual_temperature), is.finite(annual_humidity))
microclimate_centres <- c(
  temperature = mean(pd_all$annual_temperature, na.rm = TRUE),
  humidity = mean(pd_all$annual_humidity, na.rm = TRUE)
)
microclimate_scales <- c(
  temperature = sd(pd_all$annual_temperature, na.rm = TRUE),
  humidity = sd(pd_all$annual_humidity, na.rm = TRUE)
)
microclimate_parasitism <- microclimate_data |>
  filter(plot_yr %in% pd_parasitism$plot_yr)

prt_climate_baseline <- update(
  prt_mechanistic_sem, data = microclimate_parasitism, weights = total_cells,
  na.action = na.fail
)
pa_climate_baseline <- update(
  pa_mechanistic, data = microclimate_parasitism, na.action = na.fail
)
pr_climate_baseline <- update(
  pr_mechanistic, data = microclimate_parasitism, na.action = na.fail
)
hr_climate_baseline <- update(
  hr_mechanistic, data = microclimate_data, na.action = na.fail
)
ha_climate_baseline <- update(
  ha_mechanistic, data = microclimate_data, na.action = na.fail
)

prt_climate_additive <- update(prt_climate_baseline, . ~ . + sc_temp + sc_humid)
pa_climate_additive <- update(pa_climate_baseline, . ~ . + sc_temp + sc_humid)
pr_climate_additive <- update(pr_climate_baseline, . ~ . + sc_temp + sc_humid)
hr_climate_additive <- update(hr_climate_baseline, . ~ . + sc_temp + sc_humid)
ha_climate_additive <- update(ha_climate_baseline, . ~ . + sc_temp + sc_humid)

prt_climate_age <- update(prt_climate_additive, . ~ . + sc_fa:sc_temp + sc_fa:sc_humid)
pa_climate_age <- update(pa_climate_additive, . ~ . + sc_fa:sc_temp + sc_fa:sc_humid)
pr_climate_age <- update(pr_climate_additive, . ~ . + sc_fa:sc_temp + sc_fa:sc_humid)
hr_climate_age <- update(hr_climate_additive, . ~ . + sc_fa:sc_temp + sc_fa:sc_humid)
ha_climate_age <- update(ha_climate_additive, . ~ . + sc_fa:sc_temp + sc_fa:sc_humid)

microclimate_models <- list(
  parasitism = list(baseline = prt_climate_baseline, additive = prt_climate_additive,
                   age_dependent = prt_climate_age),
  realized_parasitism = list(baseline = pa_climate_baseline, additive = pa_climate_additive,
                            age_dependent = pa_climate_age),
  parasitoid_richness = list(baseline = pr_climate_baseline, additive = pr_climate_additive,
                           age_dependent = pr_climate_age),
  host_richness = list(baseline = hr_climate_baseline, additive = hr_climate_additive,
                      age_dependent = hr_climate_age),
  host_abundance = list(baseline = ha_climate_baseline, additive = ha_climate_additive,
                       age_dependent = ha_climate_age)
)

# Primary within this sensitivity: all four climate terms versus baseline.
# Secondary: both age interactions versus additive climate; additive versus
# baseline is descriptive. Adjust each comparison across the five responses.
# Coefficient CIs remain unadjusted 95% Wald intervals.
microclimate_fit_rows <- microclimate_effect_rows <- microclimate_test_rows <- list()
for (response_name in names(microclimate_models)) {
  climate_fits <- microclimate_models[[response_name]]
  fit_summary <- compare_candidate_models(climate_fits)
  if (!all(fit_summary$valid) || n_distinct(fit_summary$n) != 1L) {
    stop("Microclimate fits must be numerically valid and use identical samples: ", response_name)
  }
  microclimate_fit_rows[[response_name]] <- fit_summary
  microclimate_effect_rows[[response_name]] <- summarize_fixed_terms(
    climate_fits, names(fixef(climate_fits$age_dependent)$cond)
  )
  climate_tests <- list(
    joint_climate = anova(climate_fits$baseline, climate_fits$age_dependent),
    climate_by_age = anova(climate_fits$additive, climate_fits$age_dependent),
    additive_climate = anova(climate_fits$baseline, climate_fits$additive)
  )
  for (comparison_name in names(climate_tests)) {
    test <- climate_tests[[comparison_name]]
    microclimate_test_rows[[paste(response_name, comparison_name)]] <- tibble(
      response = response_name, comparison = comparison_name,
      likelihood_ratio = test$Chisq[2], df = test$`Chi Df`[2],
      p_value = test$`Pr(>Chisq)`[2], delta_AIC = test$AIC[2] - test$AIC[1]
    )
  }
}
microclimate_fit_summary <- bind_rows(microclimate_fit_rows, .id = "response")
microclimate_effects <- bind_rows(microclimate_effect_rows, .id = "response")
microclimate_comparisons <- bind_rows(microclimate_test_rows) |>
  group_by(comparison) |>
  mutate(p_holm = p.adjust(p_value, method = "holm", n = 5L)) |>
  ungroup()

# ------------------------------------------------------------------------------
# 8. Export final model, validation, and inference objects ####
# ------------------------------------------------------------------------------

selected_maximal_models <- list(
  parasitism = prt_mechanistic_sem,
  realized_parasitism = pa_mechanistic,
  host_abundance = ha_mechanistic,
  host_richness = hr_mechanistic,
  parasitoid_richness = pr_mechanistic
)

maximal_models <- list(
  maximal_models = selected_maximal_models,
  response_encodings = list(
    parasitism_cbind = prt_mechanistic_cbind
  ),
  forest_models = list(
    design_aligned = list(
      parasitism = prt_forest_structure,
      realized_parasitism = pa_forest_structure,
      parasitoid_richness = pr_forest_structure,
      host_richness = hr_forest_structure,
      host_abundance = ha_forest_structure
    ),
    maximal_aligned_references = list(
      parasitism = prt_mechanistic_forest_reference,
      realized_parasitism = pa_mechanistic_forest_reference,
      parasitoid_richness = pr_mechanistic_forest_reference,
      host_richness = hr_mechanistic_forest_reference
    )
  ),
  sensitivities = list(
    parasitism_mortality = prt_mechanistic_mortality_sem,
    parasitism_mortality_cbind = prt_mechanistic_mortality_cbind,
    parasitism_presence = prt_mechanistic_present,
    realized_parasitism_mortality = pa_mechanistic_mortality,
    realized_parasitism_presence = pa_mechanistic_present
  ),
  approved_structure_names = approved_structure_names,
  microclimate = list(
    models = microclimate_models,
    data = list(all = microclimate_data, parasitism = microclimate_parasitism),
    centres = microclimate_centres,
    scales = microclimate_scales
  ),
  sample_accounting = structure_selection_results$sample_accounting
)

maximal_validation_results <- list(
  settings = list(
    final_residual_simulations = final_validation_n_sim,
    sensitivity_residual_simulations = sensitivity_validation_n_sim,
    optimizer_control = list(iter.max = 10000L, eval.max = 10000L)
  ),
  maximal = maximal_validations,
  sensitivities = sensitivity_validations,
  optimizer_stability = optimizer_stability,
  optimizer_fixed_effects = optimizer_fixed_effects,
  parasitism_encoding_equivalence = prt_encoding_equivalence
)

maximal_inference_results <- list(
  forest_block_comparisons = forest_block_comparisons,
  community_block_comparisons = community_block_comparisons,
  forest_drop1 = forest_drop1,
  maximal_drop1 = maximal_drop1,
  interaction_blocks = interaction_blocks,
  presence_drop1 = presence_drop1,
  presence_interaction_blocks = presence_interaction_blocks,
  sensitivity_effects = sensitivity_effects,
  microclimate = list(
    fit_summary = microclimate_fit_summary,
    coefficients = microclimate_effects,
    comparisons = microclimate_comparisons
  )
)

final_output_paths <- c(
  models = file.path(maximal_model_output_dir, "maximal_models.rds"),
  validation = file.path(
    maximal_model_output_dir,
    "maximal_validation_results.rds"
  ),
  inference = file.path(
    maximal_model_output_dir,
    "maximal_inference_results.rds"
  )
)

saveRDS(maximal_models, final_output_paths[["models"]])
saveRDS(maximal_validation_results, final_output_paths[["validation"]])
saveRDS(maximal_inference_results, final_output_paths[["inference"]])

message(
  "Final maximal models, validation, and inference exported to: ",
  paste(final_output_paths, collapse = ", ")
)
message("Final maximal-model analysis phase ran successfully ✔️")
