## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Parsimonious structural equation models
# Author: Massimo Martini

# Script description: This script builds the primary and realized-parasitism
# parsimonious structural equation models from the completed maximal SEMs. It
# applies the recorded component reductions while retaining the selected
# response families, nuisance structures, random effects, and analysis samples.
# It then reconstructs each pathway graph, tests the removed paths, summarizes
# overall graph fit, and recalculates standardized pathway effects and component
# R2. The script saves complete RDS handoffs and readable CSV result tables.

# Run this script once from the project root after Maximal_SEM.R.

# The maximal models remain the confirmatory reference. The parsimonious SEMs
# are exploratory pathway summaries and do not independently establish causal
# direction.

# Reproducibility: run with R 4.6.1 and the package versions in renv.lock.

# Approximate runtime: 5 minutes

# ------------------------------------------------------------------------------
# 1. Load libraries and shared functions ####
# ------------------------------------------------------------------------------

library(dplyr)
library(glmmTMB)

source("analysis_functions.R")

# Use the same strict optimizer settings as the maximal component models.
strict_control <- glmmTMB::glmmTMBControl(
  optCtrl = list(iter.max = 10000L, eval.max = 10000L)
)

# ------------------------------------------------------------------------------
# 2. Import the maximal-SEM handoff ####
# ------------------------------------------------------------------------------

# The handoffs contain the fitted components, aligned datasets, node scales,
# relevant ranges, and graph metadata required below.
model_output_dir <- "model_outputs"
maximal_sem_output_dir <- file.path(model_output_dir, "maximal_sem")
parsimonious_sem_output_dir <- file.path(
  model_output_dir,
  "parsimonious_sem"
)
handoff_path <- file.path(
  maximal_sem_output_dir,
  "maximal_sem_handoff.rds"
)
realized_handoff_path <- file.path(
  maximal_sem_output_dir,
  "maximal_realized_sem_handoff.rds"
)

missing_handoffs <- c(handoff_path, realized_handoff_path)[
  !file.exists(c(handoff_path, realized_handoff_path))
]
if (length(missing_handoffs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(basename(missing_handoffs), collapse = ", "),
    ". ",
    "Run Maximal_SEM.R first."
  )
}

maximal_sem_handoff <- readRDS(handoff_path)
realized_maximal_sem_handoff <- readRDS(realized_handoff_path)

required_handoff_objects <- c(
  "sem_variant",
  "primary_sem_models",
  "node_reference_datasets",
  "node_reference_table",
  "reference_values",
  "relevant_ranges_ready",
  "primary_component_metadata",
  "primary_sem_term_results",
  "primary_edge_spec",
  "primary_standardization_spec",
  "primary_standardized_effects",
  "primary_fd_sv_association"
)
required_realized_handoff_objects <- c(
  "sem_variant",
  "primary_sem_models",
  "primary_component_metadata",
  "primary_edge_spec"
)

missing_handoff_objects <- setdiff(
  required_handoff_objects,
  names(maximal_sem_handoff)
)
missing_realized_handoff_objects <- setdiff(
  required_realized_handoff_objects,
  names(realized_maximal_sem_handoff)
)

if (length(missing_handoff_objects) > 0L ||
    length(missing_realized_handoff_objects) > 0L) {
  stop(
    "A maximal-SEM handoff is incomplete. Primary missing: ",
    paste(missing_handoff_objects, collapse = ", "),
    "; realized missing: ",
    paste(missing_realized_handoff_objects, collapse = ", "),
    ". Rerun Maximal_SEM.R first."
  )
}

# ------------------------------------------------------------------------------
# 3. Unpack and validate the maximal-SEM inputs ####
# ------------------------------------------------------------------------------

maximal_sem_models <- maximal_sem_handoff$primary_sem_models
sem_data <- maximal_sem_handoff$node_reference_datasets
node_reference_table <- maximal_sem_handoff$node_reference_table
reference_values <- maximal_sem_handoff$reference_values
relevant_ranges_ready <- maximal_sem_handoff$relevant_ranges_ready
maximal_component_metadata <- maximal_sem_handoff$primary_component_metadata
maximal_sem_term_results <- maximal_sem_handoff$primary_sem_term_results
maximal_edge_spec <- maximal_sem_handoff$primary_edge_spec
maximal_standardization_spec <-
  maximal_sem_handoff$primary_standardization_spec
maximal_standardized_effects <-
  maximal_sem_handoff$primary_standardized_effects
maximal_fd_sv_association <-
  maximal_sem_handoff$primary_fd_sv_association

realized_maximal_sem_models <-
  realized_maximal_sem_handoff$primary_sem_models
realized_component_metadata <-
  realized_maximal_sem_handoff$primary_component_metadata
realized_maximal_edge_spec <-
  realized_maximal_sem_handoff$primary_edge_spec

required_models <- c(
  "prt_sem",
  "pr_sem",
  "hr_sem",
  "ha_sem",
  "sv_sem",
  "fd_sem"
)
required_data <- c(
  "sem_data_master",
  "sem_data_parasitism",
  "sem_data_fd"
)

model_data_ids <- stats::setNames(
  maximal_component_metadata$data_id,
  maximal_component_metadata$source_model
)

# Registry checks prevent silent changes in component order or observational
# grain; every reduction refit must use the same native rows as its parent.
if (!identical(names(maximal_sem_models), required_models) ||
    !all(required_data %in% names(sem_data)) ||
    !identical(names(model_data_ids), required_models) ||
    !all(model_data_ids %in% names(sem_data))) {
  stop("The maximal-SEM model, data, or component registries are inconsistent.")
}

valid_maximal_inputs <- vapply(required_models, function(model_id) {
  model <- maximal_sem_models[[model_id]]
  inherits(model, "glmmTMB") &&
    isTRUE(model$fit$convergence == 0L) &&
    isTRUE(model$sdr$pdHess) &&
    all(is.finite(glmmTMB::fixef(model)$cond)) &&
    stats::nobs(model) == nrow(sem_data[[model_data_ids[[model_id]]]])
}, logical(1))

if (any(!valid_maximal_inputs)) {
  stop(
    "Invalid maximal model or native sample for: ",
    paste(names(valid_maximal_inputs)[!valid_maximal_inputs], collapse = ", "),
    "."
  )
}

realized_model_data_ids <- stats::setNames(
  realized_component_metadata$data_id,
  realized_component_metadata$source_model
)
pa_sem_maximal <- realized_maximal_sem_models[["pa_sem"]]
pa_data_id <- realized_model_data_ids[["pa_sem"]]

if (!identical(maximal_sem_handoff$sem_variant, "parasitism_rate") ||
    !identical(realized_maximal_sem_handoff$sem_variant,
               "realized_parasitism") ||
    !inherits(pa_sem_maximal, "glmmTMB") ||
    !isTRUE(pa_sem_maximal$fit$convergence == 0L) ||
    !isTRUE(pa_sem_maximal$sdr$pdHess) ||
    any(!is.finite(glmmTMB::fixef(pa_sem_maximal)$cond)) ||
    length(pa_data_id) != 1L ||
    !pa_data_id %in% names(sem_data) ||
    stats::nobs(pa_sem_maximal) != nrow(sem_data[[pa_data_id]])) {
  stop("The realized-parasitism maximal component is invalid.")
}

# ------------------------------------------------------------------------------
# 4. Apply the recorded component reductions ####
# ------------------------------------------------------------------------------

# Each plan records the reviewed reduction sequence and the reason for every
# decision. At each step, one term is removed from the last accepted model while
# the family, nuisance adjustment, random effects, dispersion, and hierarchy
# remain unchanged. Selected steps also compare grouped predictive performance.

## 4.1 Parasitism ##
parasitism_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_logtr:sc_fa", "accept", FALSE,
  paste(
    "No LRT support (P = 0.654); candidate AIC was 1.80 lower;",
    "all remaining estimates changed by less than 0.09 parent-model SE."
  ),
  2L, "sc_fa:rtcells", "accept", FALSE,
  paste(
    "No LRT support (P = 0.489); candidate AIC was 1.52 lower;",
    "all remaining estimates changed by less than 0.35 parent-model SE."
  ),
  3L, "sc_fa:host_rich", "accept", FALSE,
  paste(
    "No LRT support (P = 0.583); candidate AIC was 1.70 lower;",
    "all remaining estimates changed by less than 0.19 parent-model SE."
  ),
  4L, "rtfd", "accept", TRUE,
  paste(
    "No LRT support (P = 0.414) and candidate AIC was 1.33 lower;",
    "grouped CV R2 decreased by only 0.0016;",
    "supported retained terms changed by no more than 0.11 parent-model SE.",
    "The unsupported logtr slope changed by 0.61 SE but remained unsupported."
  ),
  5L, "sc_logtr", "accept", TRUE,
  paste(
    "No LRT support (P = 0.704) and candidate AIC was 1.86 lower;",
    "all remaining estimates changed by less than 0.09 parent-model SE;",
    "grouped CV R2 increased by 0.0009."
  )
)

# Replay the recorded decisions and compare each candidate with its parent.
parasitism_reduction <- run_component_reduction_plan(
  component = "parasitism",
  maximal_model = maximal_sem_models[["prt_sem"]],
  data = sem_data[[model_data_ids[["prt_sem"]]]],
  plan = parasitism_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

print(parasitism_reduction$removal_log)

parasitism_reduction$comparison_history$step_5$model_comparison
parasitism_reduction$comparison_history$step_5$likelihood_ratio_test
parasitism_reduction$comparison_history$step_5$coefficient_stability

parasitism_reduction$current_term_results

prt_sem_parsim <- parasitism_reduction$current_model
prt_sem_parsim_terms <- parasitism_reduction$current_term_results

# Component-level add-backs are provisional checks only; the complete basis set
# is derived and retested after all component equations have been finalized.
parasitism_missing_path_spec <- tibble::tibble(
  path_id = c(
    "tree_fd__parasitism",
    "tree_richness__parasitism"
  ),
  focal_term = c("rtfd", "sc_logtr"),
  added_terms = list(
    "rtfd",
    "sc_logtr"
  ),
  conditioning_terms = list(
    # Tree richness is a parent of tree FD and must condition the FD claim.
    "sc_logtr",
    character()
  )
)

provisional_parasitism_missing_paths <- audit_component_missing_paths(
  component = "parasitism",
  final_model = prt_sem_parsim,
  data = sem_data[[model_data_ids[["prt_sem"]]]],
  path_spec = parasitism_missing_path_spec,
  reference_values = reference_values,
  control = strict_control
)

print(
  provisional_parasitism_missing_paths$summary,
  width = Inf
) # provisional d-separation tests support pathway removals

## 4.2 Realized parasitism ##
real_parasitism_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_logtr:sc_fa", "accept", FALSE,
  paste(
    "No LRT support (P = 0.864); candidate AIC was 1.97 lower;",
    "all remaining estimates changed by no more than 0.03 parent-model SE."
  ),
  2L, "sc_fa:host_rich", "accept", TRUE,
  paste(
    "No LRT support (P = 0.200); candidate AIC was 0.36 lower;",
    "grouped CV R2 decreased by only 0.00430;",
    "retained terms changed by no more than 0.12 SE.",
    "Unsupported sc_fa:rtcells interaction changed by 0.44 SE",
    "but retained its sign and remained unsupported."
  ),
  3L, "sc_fa:rtcells", "accept", TRUE,
  paste(
    "No LRT support (P = 0.133); candidate AIC was only 0.26 higher;",
    "grouped CV R2 increased by 0.0146;",
    "all retained terms preserved their signs and inferential conclusions.",
    "The supported sc_fa:par_rich10 estimate decreased by 0.83",
    "parent-model SE (~21%) but remained strongly supported",
    "(P changed from 0.0000647 to 0.000145)."
  ),
  4L, "rtfd", "accept", TRUE,
  paste(
    "No LRT support (P = 0.798); candidate AIC was 1.93 lower;",
    "grouped CV R2 increased by 0.00517;",
    "remaining estimates changed by no more than 0.19 SE;",
    "no signs or inferential conclusions changed."
  ),
  5L, "sc_logtr", "accept", TRUE,
  paste(
    "No LRT support (P = 0.807); candidate AIC was 1.94 lower;",
    "grouped CV R2 increased by 0.00192;",
    "all remaining estimates changed by no more than 0.06 parent-model SE;",
    "no signs or inferential conclusions changed."
  )
)

real_parasitism_reduction <- run_component_reduction_plan(
  component = "realized_parasitism",
  maximal_model = pa_sem_maximal,
  data = sem_data[[pa_data_id]],
  plan = real_parasitism_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

print(real_parasitism_reduction$removal_log)

real_parasitism_reduction$comparison_history$step_5$model_comparison
real_parasitism_reduction$comparison_history$step_5$likelihood_ratio_test
real_parasitism_reduction$comparison_history$step_5$coefficient_stability

real_parasitism_reduction$current_term_results

pa_sem_parsim <- real_parasitism_reduction$current_model
pa_sem_parsim_terms <- real_parasitism_reduction$current_term_results

# Component-level add-backs are provisional checks only; the complete basis set
# is derived and retested after all component equations have been finalized.
real_parasitism_missing_path_spec <- tibble::tibble(
  path_id = c(
    "tree_fd__realized_parasitism",
    "tree_richness__realized_parasitism"
  ),
  focal_term = c("rtfd", "sc_logtr"),
  added_terms = list(
    "rtfd",
    "sc_logtr"
  ),
  conditioning_terms = list(
    "sc_logtr",
    character()
  )
)

provisional_real_parasitism_missing_paths <- audit_component_missing_paths(
  component = "realized_parasitism",
  final_model = pa_sem_parsim,
  data = sem_data[[pa_data_id]],
  path_spec = real_parasitism_missing_path_spec,
  reference_values = reference_values,
  control = strict_control
)

print(
  provisional_real_parasitism_missing_paths$summary,
  width = Inf
) # provisional d-separation tests support pathway removals

## 4.3 Parasitoid richness ##
par_rich_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_fa:rtsv", "accept", FALSE,
  paste(
    "No LRT support (P = 0.709); candidate AIC was 1.86 lower;",
    "all remaining estimates changed by no more than 0.12 parent-model SE"
  ),
  2L, "sc_fa:host_rich", "accept", TRUE,
  paste(
    "No LRT support (P = 0.342); candidate AIC was 1.10 lower;",
    "grouped CV R2 increased by 0.00047;",
    "the supported sc_fa:rtcells estimate changed by 0.45 parent-model SE",
    "but retained its sign and became more strongly supported",
    "(P changed from 0.00343 to 0.000124)."
  ),
  3L, "sc_logtr:sc_fa", "accept", FALSE,
  paste(
    "No LRT support (P = 0.655); candidate AIC was 1.80 lower;",
    "all remaining estimates changed by no more than 0.11 parent-model SE."
  ),
  4L, "rtfd", "accept", TRUE,
  paste(
    "No LRT support (P = 0.715); candidate AIC was 1.87 lower;",
    "grouped CV R2 increased by 0.00141;",
    "the tree-richness estimate retained its positive sign and changed by",
    "0.27 parent-model SE, while its P value changed from 0.0687 to 0.0208"
  )
)

par_rich_reduction <- run_component_reduction_plan(
  component = "parasitoid_richness",
  maximal_model = maximal_sem_models[["pr_sem"]],
  data = sem_data[[model_data_ids[["pr_sem"]]]],
  plan = par_rich_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

par_rich_reduction$comparison_history$step_4$model_comparison
par_rich_reduction$comparison_history$step_4$likelihood_ratio_test
par_rich_reduction$comparison_history$step_4$coefficient_stability

par_rich_reduction$current_term_results

par_rich_sem_parsim <- par_rich_reduction$current_model
par_rich_sem_parsim_terms <- par_rich_reduction$current_term_results


# Provisional missing-path tests for finalized components
par_rich_missing_path_spec <- tibble::tibble(
  path_id = "tree_fd__parasitoid_richness",
  focal_term = "rtfd",
  added_terms = list(
    "rtfd"
  )
)

provisional_par_rich_missing_paths <- audit_component_missing_paths(
  component = "parasitoid_richness",
  final_model = par_rich_sem_parsim,
  data = sem_data[[model_data_ids[["pr_sem"]]]],
  path_spec = par_rich_missing_path_spec,
  reference_values = reference_values,
  control = strict_control
)

print(
  provisional_par_rich_missing_paths$summary,
  width = Inf
) # provisional d-separation tests support pathway removals

## 4.4 Host richness ##
host_rich_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_fa:rtsv", "accept", TRUE,
  paste(
    "No LRT support (P = 0.339); candidate AIC was 1.09 lower;",
    "grouped CV R2 increased by 0.00177;",
    "the stand-volume estimate changed by 0.29 parent-model SE but retained",
    "its positive sign and strong support;",
    "other retained conclusions were unchanged."
  ),
  2L, "sc_logtr:sc_fa", "accept", FALSE,
  paste(
    "No LRT support (P = 0.284); candidate AIC was 0.85 lower;",
    "all remaining estimates changed by less than 0.24 parent-model SE;",
    "no signs or inferential conclusions changed."
  ),
  3L, "rtfd", "accept", TRUE,
  paste(
    "No LRT support (P = 0.449); candidate AIC was 1.43 lower;",
    "grouped CV R2 increased by 0.00507;",
    "the tree-richness estimate decreased by 0.56 parent-model SE",
    "(approximately 24%) but retained its positive sign and became more",
    "strongly supported (P changed from 0.0168 to 0.00659)."
  ),
  4L, "sc_fa", "reject", TRUE,
  paste(
    "Although the fixed age trend lacked LRT support (P = 0.196) and its",
    "removal lowered AIC by 0.33, grouped CV R2 decreased by 0.0251;",
    "the fixed age term was therefore retained for temporal adjustment and",
    "out-of-sample predictive performance."
  )
)

host_rich_reduction <- run_component_reduction_plan(
  component = "host_richness",
  maximal_model = maximal_sem_models[["hr_sem"]],
  data = sem_data[[model_data_ids[["hr_sem"]]]],
  plan = host_rich_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

host_rich_reduction$comparison_history$step_4$model_comparison
host_rich_reduction$comparison_history$step_4$likelihood_ratio_test
host_rich_reduction$comparison_history$step_4$coefficient_stability

host_rich_reduction$current_term_results

host_rich_sem_parsim <- host_rich_reduction$current_model
host_rich_sem_parsim_terms <- host_rich_reduction$current_term_results


# Provisional missing-path tests for finalized components
host_rich_missing_path_spec <- tibble::tibble(
  path_id = "tree_fd__host_richness",
  focal_term = "rtfd",
  added_terms = list(
    "rtfd"
  )
)

provisional_host_rich_missing_paths <- audit_component_missing_paths(
  component = "host_richness",
  final_model = host_rich_sem_parsim,
  data = sem_data[[model_data_ids[["hr_sem"]]]],
  path_spec = host_rich_missing_path_spec,
  reference_values = reference_values,
  control = strict_control
)

print(
  provisional_host_rich_missing_paths$summary,
  width = Inf
) # provisional d-separation tests support pathway removals

## 4.5 Host abundance ##
host_abund_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_logtr:sc_fa", "accept", FALSE,
  paste(
    "No LRT support (P = 0.593); candidate AIC was 1.71 lower;",
    "all remaining estimates changed by no more than 0.15 parent-model SE;",
    "the supported stand-volume-by-age interaction retained its sign and",
    "magnitude."
  ),
  2L, "sc_logtr", "accept", FALSE,
  paste(
    "No LRT support (P = 0.878); candidate AIC was 1.98 lower;",
    "all remaining estimates changed by no more than 0.12 parent-model SE"
  )
)

host_abund_reduction <- run_component_reduction_plan(
  component = "host_abundance",
  maximal_model = maximal_sem_models[["ha_sem"]],
  data = sem_data[[model_data_ids[["ha_sem"]]]],
  plan = host_abund_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

host_abund_reduction$comparison_history$step_2$model_comparison
host_abund_reduction$comparison_history$step_2$likelihood_ratio_test
host_abund_reduction$comparison_history$step_2$coefficient_stability

host_abund_reduction$current_term_results

host_abund_sem_parsim <- host_abund_reduction$current_model
host_abund_sem_parsim_terms <- host_abund_reduction$current_term_results

# Provisional missing-path tests for finalized components
host_abund_missing_path_spec <- tibble::tibble(
  path_id = "tree_richness__host_abundance",
  focal_term = "sc_logtr",
  added_terms = list(
    "sc_logtr"
  )
)

provisional_host_abund_missing_paths <- audit_component_missing_paths(
  component = "host_abundance",
  final_model = host_abund_sem_parsim,
  data = sem_data[[model_data_ids[["ha_sem"]]]],
  path_spec = host_abund_missing_path_spec,
  reference_values = reference_values,
  control = strict_control
)

print(
  provisional_host_abund_missing_paths$summary,
  width = Inf
) # provisional d-separation tests support pathway removals

## 4.6 Stand volume ##
# No stand-volume term was removed, so its maximal component is carried forward.
sv_sem_parsim <- maximal_sem_models[["sv_sem"]]

sv_sem_parsim_terms <- maximal_sem_term_results |>
  dplyr::filter(model_id == "sv_sem")

## 4.7 Tree FD ##
# No tree-FD term was removed, so its 88-plot component is carried forward.
fd_sem_parsim <- maximal_sem_models[["fd_sem"]]

fd_sem_parsim_terms <- maximal_sem_term_results |>
  dplyr::filter(model_id == "fd_sem")

# This registry is the finalized set of component equations used below.
parsimonious_sem_models <- list(
  prt_sem = prt_sem_parsim,
  pr_sem = par_rich_sem_parsim,
  hr_sem = host_rich_sem_parsim,
  ha_sem = host_abund_sem_parsim,
  sv_sem = sv_sem_parsim,
  fd_sem = fd_sem_parsim
)

realized_parsimonious_sem_models <- list(
  pa_sem = pa_sem_parsim,
  pr_sem = par_rich_sem_parsim,
  hr_sem = host_rich_sem_parsim,
  ha_sem = host_abund_sem_parsim,
  sv_sem = sv_sem_parsim,
  fd_sem = fd_sem_parsim
)

# Calculate final-fit marginal/conditional Nakagawa R2 for the seven unique
# components once, then split the summaries by SEM variant.
final_component_models <- list(
  parasitism = prt_sem_parsim,
  realized_parasitism = pa_sem_parsim,
  parasitoid_richness = par_rich_sem_parsim,
  host_richness = host_rich_sem_parsim,
  host_abundance = host_abund_sem_parsim,
  stand_volume = sv_sem_parsim,
  tree_fd = fd_sem_parsim
)
final_component_model_ids <- c(
  parasitism = "prt_sem",
  realized_parasitism = "pa_sem",
  parasitoid_richness = "pr_sem",
  host_richness = "hr_sem",
  host_abundance = "ha_sem",
  stand_volume = "sv_sem",
  tree_fd = "fd_sem"
)
final_component_data_ids <- c(
  parasitism = "sem_data_parasitism",
  realized_parasitism = "sem_data_parasitism",
  parasitoid_richness = "sem_data_parasitism",
  host_richness = "sem_data_master",
  host_abundance = "sem_data_master",
  stand_volume = "sem_data_master",
  tree_fd = "sem_data_fd"
)
final_component_data <- stats::setNames(
  lapply(final_component_data_ids, function(data_id) sem_data[[data_id]]),
  names(final_component_data_ids)
)

final_component_r2 <- dplyr::bind_rows(lapply(
  names(final_component_models),
  function(component) {
    tibble::as_tibble(sem_component_nakagawa_r2(
      final_component_models[[component]],
      final_component_data[[component]]
    )) |>
      dplyr::mutate(
        component = component,
        model_id = final_component_model_ids[[component]],
        data_id = final_component_data_ids[[component]],
        .before = 1
      )
  }
))

stopifnot(
  all(is.finite(final_component_r2$r2_marginal)),
  all(is.finite(final_component_r2$r2_conditional))
)

parsimonious_component_r2 <- final_component_r2 |>
  dplyr::filter(component != "realized_parasitism")
realized_parsimonious_component_r2 <- final_component_r2 |>
  dplyr::filter(component != "parasitism")

# ------------------------------------------------------------------------------
# 5. Assemble the graph and test the final basis set ####
# ------------------------------------------------------------------------------

# The adjusted plot-level FD--stand-volume association is descriptive rather
# than a directed path; it is unchanged because both components were retained.
parsimonious_fd_sv_association <- maximal_fd_sv_association

parsimonious_nondirectional_associations <- tibble::tibble(
  association_id = "tree_fd__stand_volume",
  node_1 = "tree_fd",
  node_2 = "stand_volume",
  association_type = "residual_correlation",
  native_n = parsimonious_fd_sv_association$n_plots,
  estimate = parsimonious_fd_sv_association$correlation,
  conditioning = parsimonious_fd_sv_association$conditioning
)

expected_removed_edges <- c(
  "tree_fd__parasitism",
  "tree_richness__parasitism",
  "tree_fd__parasitoid_richness",
  "tree_fd__host_richness",
  "tree_richness__host_abundance"
)

# Rebuild edge presence from the finalized model coefficients and assert that
# the graph differs from the maximal graph by exactly the five declared paths.
parsimonious_graph <- assemble_reduced_sem_graph(
  models = parsimonious_sem_models,
  maximal_edge_spec = maximal_edge_spec,
  expected_removed_edges = expected_removed_edges,
  nondirectional_associations = parsimonious_nondirectional_associations
)

parsimonious_edge_spec <- parsimonious_graph$directed_edges
removed_edges <- parsimonious_graph$removed_edges
parsimonious_model_terms <- parsimonious_graph$model_terms

# Each basis claim conditions on the relevant parents of both variables.
# Component models already contain most conditioning variables plus all age,
# nuisance, random, and dispersion terms; `conditioning_terms` lists only terms
# that must be added to construct the claim-specific null model.
basis_conditioning <- tibble::tibble(
  edge_id = expected_removed_edges,
  conditioning_nodes = list(
    c(
      "tree_richness", "stand_volume", "host_abundance",
      "host_richness", "parasitoid_richness"
    ),
    c(
      "stand_volume", "host_abundance",
      "host_richness", "parasitoid_richness"
    ),
    c(
      "tree_richness", "stand_volume",
      "host_abundance", "host_richness"
    ),
    c("tree_richness", "stand_volume", "host_abundance"),
    c("tree_fd", "stand_volume")
  ),
  conditioning_terms = list(
    "sc_logtr",
    character(),
    character(),
    character(),
    character()
  )
)

parsimonious_basis_set <- derive_sem_basis_set(
  maximal_edge_spec = maximal_edge_spec,
  removed_edges = removed_edges,
  conditioning_spec = basis_conditioning,
  models = parsimonious_sem_models
)

# Test every missing path in the correct response direction and native dataset.
# The conditioned null and its focal add-back share the final family and random
# structure; each LRT therefore tests one missing continuous pathway (1 df).
final_basis_audit <- audit_sem_basis_set(
  basis_set = parsimonious_basis_set,
  models = parsimonious_sem_models,
  data_sets = sem_data,
  model_data_ids = model_data_ids,
  component_metadata = maximal_component_metadata,
  reference_values = reference_values,
  control = strict_control
)

final_missing_path_audits <- final_basis_audit$audits
final_missing_path_results <- final_basis_audit$results
# Fisher's C combines the five valid claim P values; it is an exploratory
# post-selection goodness-of-fit summary, not a cross-grain global likelihood.
parsimonious_global_fit <- sem_fisher_c(final_missing_path_results)

# One publication-facing row per retained directed path. Standardizations are
# recomputed from the parsimonious fits: response-scale average partial effect
# (primary), Menard OE (sensitivity), and relevant-range contrast (sensitivity).
# Fixed-effect predictions exclude realized random effects; age-moderated paths
# are evaluated at centered age, and standardized effects are point estimates.
parsimonious_path_table <- build_sem_path_table(
  edge_spec = parsimonious_edge_spec,
  models = parsimonious_sem_models,
  data_sets = sem_data,
  node_reference = node_reference_table,
  reference_values = reference_values,
  include_relevant_range = relevant_ranges_ready
)

# Fail before export if an edge is missing, a standardization failed, or a
# declared relevant range falls outside the observed predictor support.
stopifnot(
  nrow(parsimonious_path_table) == nrow(parsimonious_edge_spec),
  identical(parsimonious_path_table$edge_id, parsimonious_edge_spec$edge_id),
  all(is.na(parsimonious_path_table$standardization_error)),
  all(is.finite(parsimonious_path_table$standardized_response_effect)),
  all(is.finite(parsimonious_path_table$standardized_oe_effect)),
  !relevant_ranges_ready || all(
    is.finite(parsimonious_path_table$standardized_relevant_range_effect) &
      parsimonious_path_table$rr_x_range_supported
  )
)

# The RDS handoff preserves the fitted models and nested graph objects. The CSV
# files report retained paths, the finalized graph, removed-path tests, overall
# graph fit, the complete reduction log, and component R2 values.
parsimonious_export <- export_parsimonious_sem_handoff(
  models = parsimonious_sem_models,
  graph = parsimonious_graph,
  basis_set = parsimonious_basis_set,
  missing_path_results = final_missing_path_results,
  global_fit = parsimonious_global_fit,
  path_table = parsimonious_path_table,
  reduction_logs = list(
    parasitism = parasitism_reduction$removal_log,
    parasitoid_richness = par_rich_reduction$removal_log,
    host_richness = host_rich_reduction$removal_log,
    host_abundance = host_abund_reduction$removal_log
  ),
  component_r2 = parsimonious_component_r2,
  output_dir = parsimonious_sem_output_dir,
  file_stem = "parsimonious_sem"
)

parsimonious_handoff <- parsimonious_export$handoff

# Console summaries provide a compact end-of-run audit of graph, coefficients,
# missing-path tests, and global directed-separation fit.
print(
  parsimonious_edge_spec |>
    dplyr::select(edge_id, model_id, predictor, age_moderated),
  n = Inf
)
print(parsimonious_graph$nondirectional_associations, width = Inf)
print(parsimonious_path_table, width = Inf)
print(final_missing_path_results, width = Inf)
print(parsimonious_global_fit, width = Inf)
# Fisher’s C = 5, df = 10, P = 0.891
# Fail to reject the parsimonious path diagram

# ------------------------------------------------------------------------------
# 6. Assemble and export the realized-parasitism SEM ####
# ------------------------------------------------------------------------------

# Only the terminal response differs; all finalized upstream components,
# conditioning sets, and the nondirectional FD--stand-volume association remain
# unchanged.
realized_expected_removed_edges <- sub(
  "__parasitism$",
  "__realized_parasitism",
  expected_removed_edges
)

realized_parsimonious_graph <- assemble_reduced_sem_graph(
  models = realized_parsimonious_sem_models,
  maximal_edge_spec = realized_maximal_edge_spec,
  expected_removed_edges = realized_expected_removed_edges,
  nondirectional_associations = parsimonious_nondirectional_associations
)

realized_parsimonious_edge_spec <-
  realized_parsimonious_graph$directed_edges

realized_basis_conditioning <- basis_conditioning
realized_basis_conditioning$edge_id <- realized_expected_removed_edges

realized_parsimonious_basis_set <- derive_sem_basis_set(
  maximal_edge_spec = realized_maximal_edge_spec,
  removed_edges = realized_parsimonious_graph$removed_edges,
  conditioning_spec = realized_basis_conditioning,
  models = realized_parsimonious_sem_models
)

realized_final_basis_audit <- audit_sem_basis_set(
  basis_set = realized_parsimonious_basis_set,
  models = realized_parsimonious_sem_models,
  data_sets = sem_data,
  model_data_ids = realized_model_data_ids,
  component_metadata = realized_component_metadata,
  reference_values = reference_values,
  control = strict_control
)

realized_final_missing_path_results <- realized_final_basis_audit$results
realized_parsimonious_global_fit <- sem_fisher_c(
  realized_final_missing_path_results
)

realized_parsimonious_path_table <- build_sem_path_table(
  edge_spec = realized_parsimonious_edge_spec,
  models = realized_parsimonious_sem_models,
  data_sets = sem_data,
  node_reference = node_reference_table,
  reference_values = reference_values,
  include_relevant_range = relevant_ranges_ready
)

stopifnot(
  nrow(realized_parsimonious_path_table) ==
    nrow(realized_parsimonious_edge_spec),
  identical(
    realized_parsimonious_path_table$edge_id,
    realized_parsimonious_edge_spec$edge_id
  ),
  all(is.na(realized_parsimonious_path_table$standardization_error)),
  all(is.finite(
    realized_parsimonious_path_table$standardized_response_effect
  )),
  all(is.finite(realized_parsimonious_path_table$standardized_oe_effect)),
  !relevant_ranges_ready || all(
    is.finite(
      realized_parsimonious_path_table$standardized_relevant_range_effect
    ) & realized_parsimonious_path_table$rr_x_range_supported
  )
)

realized_parsimonious_export <- export_parsimonious_sem_handoff(
  models = realized_parsimonious_sem_models,
  graph = realized_parsimonious_graph,
  basis_set = realized_parsimonious_basis_set,
  missing_path_results = realized_final_missing_path_results,
  global_fit = realized_parsimonious_global_fit,
  path_table = realized_parsimonious_path_table,
  reduction_logs = list(
    realized_parasitism = real_parasitism_reduction$removal_log,
    parasitoid_richness = par_rich_reduction$removal_log,
    host_richness = host_rich_reduction$removal_log,
    host_abundance = host_abund_reduction$removal_log
  ),
  component_r2 = realized_parsimonious_component_r2,
  output_dir = parsimonious_sem_output_dir,
  file_stem = "parsimonious_realized_sem"
)

realized_parsimonious_handoff <- realized_parsimonious_export$handoff

print(
  realized_parsimonious_edge_spec |>
    dplyr::select(edge_id, model_id, predictor, age_moderated),
  n = Inf
)
print(realized_parsimonious_path_table, width = Inf)
print(realized_final_missing_path_results, width = Inf)
print(realized_parsimonious_global_fit, width = Inf)

message("Primary and realized parsimonious SEMs ran successfully ✔️")
