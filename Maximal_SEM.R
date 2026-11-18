## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Maximal structural equation models
# Author: Massimo Martini

# Script description: This script builds the maximal parasitism and realized-
# parasitism structural equation models from the validated maximal GLMMs. It
# aligns the response-specific datasets and assembles the component models into
# directed pathways. The script summarizes pathway estimates, standardized
# effects, component-model R2, and the plot-level association between
# tree functional diversity and stand volume. It saves readable result tables
# and complete RDS handoffs for downstream figure and reporting scripts.

# Run this script once from the project root after Maximal_models.R.

# Reproducibility: run with R 4.6.1 and the package versions in renv.lock.

# Approximate runtime: 1 minute

# ------------------------------------------------------------------------------
# 1. Packages and helper functions ####
# ------------------------------------------------------------------------------

library(dplyr)
library(glmmTMB)

source("analysis_functions.R")

# Shared optimizer control inherited by the SEM-compatible component refits.
strict_control <- glmmTMB::glmmTMBControl(
  optCtrl = list(iter.max = 10000L, eval.max = 10000L)
)

# ------------------------------------------------------------------------------
# 2. Read prepared data and final maximal models ####
# ------------------------------------------------------------------------------

model_output_dir <- "model_outputs"
design_model_output_dir <- file.path(model_output_dir, "design_models")
maximal_model_output_dir <- file.path(model_output_dir, "maximal_models")
maximal_sem_output_dir <- file.path(model_output_dir, "maximal_sem")

input_paths <- c(
  analysis_data = file.path(design_model_output_dir, "analysis_data.rds"),
  maximal_models = file.path(maximal_model_output_dir, "maximal_models.rds")
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input file(s): ",
    paste(basename(missing_inputs), collapse = ", "),
    ". Run Design_models.R and Maximal_models.R first."
  )
}

analysis_data <- readRDS(input_paths[["analysis_data"]])
maximal_results <- readRDS(input_paths[["maximal_models"]])

required_datasets <- c(
  "pd_all",
  "pd_parasitism"
)

if (!all(required_datasets %in% names(analysis_data))) {
  stop(
    "analysis_data.rds does not contain the required datasets: ",
    paste(required_datasets, collapse = ", "),
    "."
  )
}

if (!"maximal_models" %in% names(maximal_results)) {
  stop("maximal_models.rds does not contain `maximal_models`.")
}

required_maximal_models <- c(
  "parasitism",
  "realized_parasitism",
  "host_abundance",
  "host_richness",
  "parasitoid_richness"
)

maximal_models <- maximal_results$maximal_models

if (!all(required_maximal_models %in% names(maximal_models))) {
  stop(
    "`maximal_models` does not contain all five authoritative mechanistic ",
    "models."
  )
}

pd_all <- analysis_data$pd_all
pd_parasitism <- analysis_data$pd_parasitism

# Authoritative fitted mechanistic components. These objects provide the
# selected response families and nuisance structures used as templates for the
# SEM-compatible refits below.
# `maximal_models$parasitism` is the scalar weighted parasitism-rate model used
# for the SEM component refit. The equivalent cbind response is stored
# separately under `response_encodings` and is intentionally not used here.
prt_mechanistic_sem <- maximal_models$parasitism
pa_mechanistic <- maximal_models$realized_parasitism
ha_mechanistic <- maximal_models$host_abundance
hr_mechanistic <- maximal_models$host_richness
pr_mechanistic <- maximal_models$parasitoid_richness

if (!identical(all.vars(formula(prt_mechanistic_sem)[[2L]]), "cell_prate")) {
  stop(
    "The imported parasitism component is not the scalar weighted ",
    "`cell_prate` model required for this SEM."
  )
}
if (!identical(all.vars(formula(pa_mechanistic)[[2L]]), "par_abund")) {
  stop(
    "The imported realized-parasitism component does not use the required ",
    "scalar `par_abund` response."
  )
}

# ------------------------------------------------------------------------------
# 3. Prepare the aligned SEM datasets ####
# ------------------------------------------------------------------------------

if (!all(c("plot_yr", "total_cells", "par_abund") %in% names(pd_all))) {
  stop("`pd_all` lacks one or more variables required for SEM alignment.")
}

if (!"plot_yr" %in% names(pd_parasitism)) {
  stop("`pd_parasitism` must contain `plot_yr`.")
}

if (anyNA(pd_all$plot_yr) || anyDuplicated(pd_all$plot_yr)) {
  stop("`plot_yr` must uniquely identify every row in `pd_all`.")
}

if (anyNA(pd_parasitism$plot_yr) || anyDuplicated(pd_parasitism$plot_yr)) {
  stop("`plot_yr` must uniquely identify every row in `pd_parasitism`.")
}

if (!all(pd_parasitism$plot_yr %in% pd_all$plot_yr)) {
  stop("Every `pd_parasitism` row must also occur in `pd_all`.")
}

# This master table retains every plot-year. Parasitism is undefined for rows
# with no brood cells, but those rows remain present for component models whose
# responses are defined there.
sem_data_master <- pd_all

# Explicit response-specific subset corresponding to plot-years on which
# parasitism is defined. Other component models retain the model frames stored
# in their fitted objects rather than being forced onto this subset.
sem_data_parasitism <- sem_data_master |>
  filter(plot_yr %in% pd_parasitism$plot_yr)

# Match FD data with 2016 tree volume
sem_data_fd <- sem_data_master |>
  filter(as.character(year) == "2016") |>
  arrange(plot_id)

stopifnot(
  nrow(sem_data_master) == nrow(pd_all),
  identical(sem_data_master$plot_yr, pd_all$plot_yr),
  identical(sem_data_parasitism$plot_yr, pd_parasitism$plot_yr),
  nrow(sem_data_fd) == dplyr::n_distinct(sem_data_master$plot_id),
  !anyDuplicated(sem_data_fd$plot_id),
  all(!is.na(sem_data_fd$rtfd)),
  all(!is.na(sem_data_fd$rtsv)),
  all(sem_data_parasitism$total_cells > 0),
  all(!is.na(sem_data_parasitism$cell_prate)),
  all(sem_data_parasitism$cell_prate >= 0),
  all(sem_data_parasitism$cell_prate <= 1),
  all(!is.na(sem_data_parasitism$par_abund)),
  all(sem_data_parasitism$par_abund >= 0)
)

# Define one reference population and fitted scale per conceptual SEM node.

node_reference_datasets <- list(
  sem_data_fd = sem_data_fd,
  sem_data_master = sem_data_master,
  sem_data_parasitism = sem_data_parasitism
)

node_reference_spec <- tibble::tribble(
  ~node, ~model_variable, ~biological_variable, ~transformation,
  ~reference_dataset, ~observational_grain,
  "tree_richness", "sc_logtr", "tree_r", "cap at 16; log2; z-scale",
  "sem_data_fd", "unique plot",
  "tree_fd", "rtfd", "tree_fd", "square root",
  "sem_data_fd", "unique plot",
  "stand_volume", "rtsv", "stand_volume", "square root",
  "sem_data_master", "plot-year",
  "host_abundance", "rtcells", "total_cells", "square root",
  "sem_data_master", "plot-year",
  "host_richness", "host_rich", "host_rich", "identity",
  "sem_data_master", "plot-year",
  "parasitoid_richness", "par_rich10", "par_rich", "cap at 10",
  "sem_data_parasitism", "brood-cell-positive plot-year",
  "parasitism", "cell_prate", "cell_prate", "identity (proportion)",
  "sem_data_parasitism", "brood-cell-positive plot-year",
  "realized_parasitism", "par_abund", "par_abund", "identity",
  "sem_data_parasitism", "brood-cell-positive plot-year",
  "forest_age", "sc_fa", "for_age", "z-scale",
  "sem_data_master", "plot-year"
)

node_reference_stats <- Map(
  function(dataset_name, variable_name, biological_name) {
    reference_data <- node_reference_datasets[[dataset_name]]
    if (!variable_name %in% names(reference_data)) {
      stop("Missing node-reference variable: `", variable_name, "`.")
    }
    if (!biological_name %in% names(reference_data)) {
      stop("Missing biological node variable: `", biological_name, "`.")
    }
    values <- reference_data[[variable_name]]
    biological_values <- reference_data[[biological_name]]
    tibble::tibble(
      n_reference = sum(!is.na(values)),
      mean_model = mean(values, na.rm = TRUE),
      sd_model = stats::sd(values, na.rm = TRUE),
      observed_low_model = min(values, na.rm = TRUE),
      observed_high_model = max(values, na.rm = TRUE),
      observed_low_biological = min(biological_values, na.rm = TRUE),
      observed_high_biological = max(biological_values, na.rm = TRUE)
    )
  },
  node_reference_spec$reference_dataset,
  node_reference_spec$model_variable,
  node_reference_spec$biological_variable
)

scaling_metadata <- tibble::tibble(
  model_variable = c("sc_logtr", "sc_fa"),
  scaling_base_variable = c("log_tr", "for_age")
)
scaling_metadata$scaling_center <- vapply(
  scaling_metadata$scaling_base_variable,
  function(variable_name) mean(sem_data_master[[variable_name]], na.rm = TRUE),
  numeric(1)
)
scaling_metadata$scaling_divisor <- vapply(
  scaling_metadata$scaling_base_variable,
  function(variable_name) stats::sd(sem_data_master[[variable_name]], na.rm = TRUE),
  numeric(1)
)

# ----------------------------------------------------------------------------
# Declared relevant ranges
# ----------------------------------------------------------------------------
# Endpoints are expressed in the original biological units; transformations to
# fitted scales are calculated reproducibly below.
relevant_range_spec <- tibble::tribble(
  ~node, ~rr_low_biological, ~rr_high_biological,
  ~rr_rationale,
  "tree_richness", 1, 16, "experimental limits",
  "tree_fd", 0, 5.593435, "10th-90th percentiles",
  "stand_volume", 0.6813006, 183.6210250, "10th-90th percentiles",
  "host_abundance", 16, 165.4333, "10th-90th percentiles",
  "host_richness", 3, 9, "10th-90th percentiles",
  "parasitoid_richness", 1, 6, "10th-90th percentiles",
  "parasitism", 0.02061154, 0.29416221, "10th-90th percentiles",
  "realized_parasitism", 1, 27, "10th-90th percentiles",
  "forest_age", 4, 14, "experimental limits"
)

node_reference_table <- node_reference_spec |>
  dplyr::bind_cols(dplyr::bind_rows(node_reference_stats)) |>
  dplyr::left_join(scaling_metadata, by = "model_variable")

# Validate the biological ranges and convert them to each fitted node scale.
range_setup <- prepare_sem_relevant_ranges(
  node_reference = node_reference_table,
  range_spec = relevant_range_spec,
  required_nodes = node_reference_spec$node
)
node_reference_table <- range_setup$node_reference
relevant_ranges_ready <- range_setup$ready

stopifnot(
  !anyDuplicated(node_reference_table$node),
  !anyDuplicated(node_reference_table$model_variable),
  all(node_reference_table$n_reference >= 2L),
  all(is.finite(node_reference_table$mean_model)),
  all(is.finite(node_reference_table$sd_model)),
  all(node_reference_table$sd_model > 0)
)

# ------------------------------------------------------------------------------
# 4. Component model assembly ####
# ------------------------------------------------------------------------------

prt_sem <- update(
  prt_mechanistic_sem,
  cell_prate ~
    sc_logtr + sc_fa +
      rtfd + rtsv +
      rtcells + host_rich + par_rich10 +
      sc_logtr:sc_fa +
      rtsv:sc_fa +
      rtcells:sc_fa +
      host_rich:sc_fa +
      par_rich10:sc_fa +
      par_rich10:rtsv +
      host_rich:rtcells +
      (1 | site/plot_id) + (1 | year),
  weights = total_cells,
  data = sem_data_parasitism,
  dispformula = ~ sc_cells,
  control = strict_control
)

pa_sem <- update(
  pa_mechanistic,
  par_abund ~
    sc_logtr + sc_fa +
      rtfd + rtsv +
      rtcells + host_rich + par_rich10 +
      sc_logtr:sc_fa +
      rtsv:sc_fa +
      rtcells:sc_fa +
      host_rich:sc_fa +
      par_rich10:sc_fa +
      par_rich10:rtsv +
      host_rich:rtcells +
      site + (1 | plot_id) + (1 | year),
  data = sem_data_parasitism,
  dispformula = ~ 1,
  control = strict_control
)

pr_sem <- update(
  pr_mechanistic,
  par_rich10 ~
    sc_logtr + sc_fa +
      rtfd + rtsv +
      rtcells + host_rich +
      sc_logtr:sc_fa +
      rtsv:sc_fa +
      rtcells:sc_fa +
      host_rich:sc_fa +
      rtcells:host_rich +
      site + (1 | year),
  data = sem_data_parasitism,
  dispformula = ~ 1,
  control = strict_control
)

hr_sem <- update(
  hr_mechanistic,
  host_rich ~
    sc_logtr + sc_fa +
      rtfd + rtsv +
      sc_slope + sc_elev + sc_east + sc_north +
      rtcells +
      sc_logtr:sc_fa +
      rtsv:sc_fa +
      (1 | site) + (1 | year),
  data = sem_data_master,
  dispformula = ~ ar1(year + 0 | site:plot_id),
  control = strict_control
)

ha_sem <- update(
  ha_mechanistic,
  rtcells ~
    sc_logtr + sc_fa +
      rtfd + rtsv +
      sc_slope + sc_elev + sc_east + sc_north +
      sc_logtr:sc_fa +
      rtsv:sc_fa +
      (1 | site/plot_id) + (1 | year),
  data = sem_data_master,
  family = gaussian(link = "identity"),
  dispformula = ~ ar1(year + 0 | site:plot_id),
  control = strict_control
)

sv_sem <- glmmTMB(
  rtsv ~
    sc_logtr * sc_fa +
      site +
      diag(1 + sc_fa | plot_id) +
      (1 | year),
  data = sem_data_master,
  family = t_family(link = "identity"),
  control = strict_control
)

fd_sem <- glmmTMB(
  rtfd ~
    sc_logtr + site,
  data = sem_data_fd,
  family = tweedie(link = "log"),
  control = strict_control
)

secondary_sem_models <- list(
  realized_parasitism = pa_sem,
  parasitoid_richness = pr_sem,
  host_richness = hr_sem,
  host_abundance = ha_sem,
  stand_volume = sv_sem,
  tree_fd = fd_sem
)

primary_edge_models <- list(
  parasitism = prt_sem,
  parasitoid_richness = pr_sem,
  host_richness = hr_sem,
  host_abundance = ha_sem,
  stand_volume = sv_sem,
  tree_fd = fd_sem
)

valid_primary_fits <- vapply(primary_edge_models, function(model) {
  isTRUE(model$fit$convergence == 0L) &&
    isTRUE(model$sdr$pdHess) &&
    all(is.finite(glmmTMB::fixef(model)$cond))
}, logical(1))
if (any(!valid_primary_fits)) {
  stop(
    "Invalid primary SEM component fit(s): ",
    paste(names(valid_primary_fits)[!valid_primary_fits], collapse = ", "),
    "."
  )
}

primary_component_metadata <- tibble::tibble(
  component = names(primary_edge_models),
  Response = vapply(
    primary_edge_models,
    function(model) all.vars(stats::formula(model))[1],
    character(1)
  ),
  source_model = c("prt_sem", "pr_sem", "hr_sem", "ha_sem", "sv_sem", "fd_sem"),
  data_id = c(
    "sem_data_parasitism", "sem_data_parasitism",
    "sem_data_master", "sem_data_master", "sem_data_master", "sem_data_fd"
  ),
  n = vapply(primary_edge_models, stats::nobs, numeric(1)),
  family = vapply(primary_edge_models, function(model) family(model)$family, character(1)),
  link = vapply(primary_edge_models, function(model) family(model)$link, character(1))
)

primary_edge_table_raw <- lapply(primary_edge_models, function(model) {
  coefficients <- summary(model)$coefficients$cond
  tibble::tibble(
    Response = all.vars(stats::formula(model))[[1]],
    Predictor = rownames(coefficients),
    Estimate = unname(coefficients[, "Estimate"]),
    Std.Error = unname(coefficients[, "Std. Error"]),
    Crit.Value = unname(coefficients[, "z value"]),
    P.Value = unname(coefficients[, "Pr(>|z|)"])
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::left_join(
  primary_component_metadata,
  by = "Response"
  )

nuisance_predictors <- c(
  "site", "sc_slope", "sc_elev", "sc_east", "sc_north"
)

primary_edge_table <- primary_edge_table_raw |>
  filter(
    Predictor != "(Intercept)",
    !Predictor %in% nuisance_predictors,
    !startsWith(Predictor, "site")
  ) |>
  mutate(
    term_type = if_else(
      grepl(":", Predictor, fixed = TRUE),
      "interaction",
      "direct_path"
    ),
    interaction_type = case_when(
      term_type == "direct_path" ~ NA_character_,
      grepl("sc_fa", Predictor, fixed = TRUE) ~ "age_moderation",
      TRUE ~ "ecological_moderation"
    ),
    conf_low_link = Estimate - stats::qnorm(0.975) * Std.Error,
    conf_high_link = Estimate + stats::qnorm(0.975) * Std.Error
  )

reference_values <- stats::setNames(
  node_reference_table$mean_model,
  node_reference_table$model_variable
)
# Reference values must be on the fitted model scale. Explicitly centered and
# scaled predictors use zero; raw and transformed predictors retain their
# declared model-scale node means.
centered_reference_variables <- scaling_metadata$model_variable
reference_values[intersect(
  names(reference_values), centered_reference_variables
)] <- 0

primary_sem_models <- stats::setNames(
  unname(primary_edge_models),
  primary_component_metadata$source_model
)
# Test each pathway at its declared moderator reference values.
primary_edge_table <- primary_edge_table |>
  dplyr::bind_cols(dplyr::bind_rows(Map(
    function(model_id, term, term_type) {
      sem_term_wald_test(
        primary_sem_models[[model_id]], term, term_type, reference_values
      )
    },
    primary_edge_table$source_model,
    primary_edge_table$Predictor,
    primary_edge_table$term_type
  )))

stopifnot(
  all(is.finite(primary_edge_table$inference_estimate)),
  all(is.finite(primary_edge_table$inference_standard_error)),
  all(primary_edge_table$inference_standard_error > 0),
  all(is.finite(primary_edge_table$inference_p_value))
)

# Computational adapter for the standardization functions.
# Stand age is retained in the component models and represented through pathway
# moderation markers, not as a focal ecological arrow.
age_moderated_edges <- primary_edge_table |>
  filter(interaction_type == "age_moderation") |>
  transmute(
    Response,
    Predictor = sub("(^sc_fa:|:sc_fa$)", "", Predictor),
    age_moderated = TRUE
  )

source_node_lookup <- node_reference_table |>
  select(source_node = node, Predictor = model_variable)

destination_node_lookup <- node_reference_table |>
  select(destination_node = node, Response = model_variable)

forest_age_variable <- node_reference_table |>
  filter(node == "forest_age") |>
  pull(model_variable)

primary_edge_spec <- primary_edge_table |>
  filter(
    term_type == "direct_path",
    Predictor != forest_age_variable
  ) |>
  left_join(source_node_lookup, by = "Predictor") |>
  left_join(destination_node_lookup, by = "Response") |>
  left_join(age_moderated_edges, by = c("Response", "Predictor")) |>
  mutate(
    age_moderated = dplyr::coalesce(age_moderated, FALSE),
    edge_id = paste(source_node, destination_node, sep = "__"),
    model_id = source_model,
    predictor = Predictor,
    response = Response,
    at = lapply(
      age_moderated,
      function(is_moderated) {
        if (is_moderated) list(sc_fa = 0) else list()
      }
    ),
    average_weight_variable = NA_character_
  ) |>
  select(
    edge_id, source_node, destination_node,
    model_id, data_id, predictor, response,
    age_moderated, at, average_weight_variable
  )

forest_age_effect_spec <- primary_edge_table |>
  filter(
    term_type == "direct_path",
    Predictor == forest_age_variable
  ) |>
  left_join(destination_node_lookup, by = "Response") |>
  mutate(
    edge_id = paste("forest_age", destination_node, sep = "__"),
    source_node = "forest_age",
    model_id = source_model,
    predictor = Predictor,
    response = Response,
    age_moderated = NA,
    at = rep(list(list()), dplyr::n()),
    average_weight_variable = NA_character_
  ) |>
  select(
    edge_id, source_node, destination_node,
    model_id, data_id, predictor, response,
    age_moderated, at, average_weight_variable
  )

primary_standardization_spec <- bind_rows(
  primary_edge_spec,
  forest_age_effect_spec
)

stopifnot(
  length(forest_age_variable) == 1L,
  !anyNA(primary_standardization_spec$source_node),
  !anyNA(primary_standardization_spec$destination_node),
  !anyDuplicated(primary_standardization_spec$edge_id),
  all(primary_standardization_spec$model_id %in%
        primary_component_metadata$source_model),
  all(primary_standardization_spec$data_id %in% names(node_reference_datasets))
)

# Calculate standardized point estimates for the focal direct
# paths and forest-age main effects. Interacting covariates retain their
# observed values when forest age is the focal predictor. These standardized
# quantities summarize effect magnitude; the Wald results above remain the
# declared source of primary SEM inference.
primary_standardized_effects <- standardize_sem_edges(
  edge_spec = primary_standardization_spec,
  models = primary_sem_models,
  data_sets = node_reference_datasets,
  node_reference = node_reference_table,
  include_relevant_range = relevant_ranges_ready,
  include_random = FALSE,
  on_error = "stop"
)

stopifnot(
  nrow(primary_standardized_effects) == nrow(primary_standardization_spec),
  identical(
    primary_standardized_effects$edge_id,
    primary_standardization_spec$edge_id
  ),
  all(is.na(primary_standardized_effects$standardization_error)),
  all(is.finite(primary_standardized_effects$standardized_response_effect)),
  all(is.finite(primary_standardized_effects$standardized_oe_effect)),
  all(is.finite(primary_standardized_effects$oe_r2)),
  all(primary_standardized_effects$random_effects == "excluded")
)

if (relevant_ranges_ready) {
  stopifnot(
    all(is.finite(
      primary_standardized_effects$standardized_relevant_range_effect
    )),
    all(primary_standardized_effects$rr_x_range_supported)
  )
}

primary_term_standardization_results <- primary_standardized_effects |>
  left_join(
    primary_standardization_spec |>
      select(edge_id, age_moderated, average_weight_variable),
    by = "edge_id"
  )

# One publication-facing table for all focal component-model terms. For main
# pathways, the primary inference is the exact link-scale slope at the declared
# moderator reference conditions. For interactions, it is the exact fitted
# interaction coefficient test. Raw component-model coefficients are retained
# separately, and main effects also receive the manual standardizations.
primary_sem_term_results <- primary_edge_table |>
  mutate(
    term_role = case_when(
      term_type == "interaction" & interaction_type == "age_moderation" ~
        "age_moderation",
      term_type == "interaction" ~ "ecological_interaction",
      Predictor == forest_age_variable ~ "forest_age_main_effect",
      TRUE ~ "main_ecological_path"
    )
  ) |>
  transmute(
    component,
    model_id = source_model,
    data_id,
    response = Response,
    predictor = Predictor,
    term_role,
    interaction_type,
    native_n = n,
    family,
    link,
    coefficient_estimate = Estimate,
    coefficient_std_error = Std.Error,
    coefficient_conf_low = conf_low_link,
    coefficient_conf_high = conf_high_link,
    coefficient_test_statistic = Crit.Value,
    coefficient_p_value = P.Value,
    inference_estimand,
    inference_method,
    reference_conditioning,
    estimate = inference_estimate,
    std_error = inference_standard_error,
    conf_low = inference_conf_low,
    conf_high = inference_conf_high,
    df = NA_real_,
    test_statistic = inference_test_statistic,
    p_value = inference_p_value
  ) |>
  left_join(
    primary_term_standardization_results |>
      select(
        model_id, data_id, response, predictor,
        edge_id, source_node, destination_node, age_moderated,
        average_weight_variable, random_effects, averaging, conditioning,
        n_evaluation, derivative_step, average_partial_effect,
        sd_x_reference, sd_y_reference, standardized_response_effect,
        oe_link_slope, oe_r2, oe_variance_eta, oe_sd_y,
        standardized_oe_effect,
        rr_x_low_model, rr_x_high_model, rr_y_low_model, rr_y_high_model,
        rr_response_contrast, standardized_relevant_range_effect,
        rr_x_low_supported, rr_x_high_supported, rr_x_range_supported,
        standardization_error
      ),
    by = c("model_id", "data_id", "response", "predictor")
  ) |>
  mutate(
    standardization_status = case_when(
      term_role %in% c("main_ecological_path", "forest_age_main_effect") &
        is.na(standardization_error) ~
        "available",
      term_role %in% c("main_ecological_path", "forest_age_main_effect") ~
        paste0("failed: ", standardization_error),
      TRUE ~ "not calculated: interaction"
    ),
    .after = age_moderated
  )

stopifnot(
  nrow(primary_sem_term_results) == nrow(primary_edge_table),
  all(is.finite(primary_sem_term_results$coefficient_estimate)),
  all(is.finite(primary_sem_term_results$coefficient_std_error)),
  all(is.finite(primary_sem_term_results$estimate)),
  all(is.finite(primary_sem_term_results$std_error)),
  all(is.na(primary_sem_term_results$df)),
  sum(primary_sem_term_results$term_role %in%
        c("main_ecological_path", "forest_age_main_effect")) ==
    nrow(primary_standardization_spec),
  all(is.finite(primary_sem_term_results$standardized_response_effect[
    primary_sem_term_results$term_role %in%
      c("main_ecological_path", "forest_age_main_effect")
  ])),
  all(is.na(primary_sem_term_results$standardized_response_effect[
    primary_sem_term_results$term_role %in%
      c("age_moderation", "ecological_interaction")
  ])),
  all(is.finite(primary_sem_term_results$standardized_relevant_range_effect[
    primary_sem_term_results$term_role %in%
      c("main_ecological_path", "forest_age_main_effect")
  ]))
)

interaction_rows <- primary_sem_term_results$term_role %in%
  c("age_moderation", "ecological_interaction")
stopifnot(
  isTRUE(all.equal(
    primary_sem_term_results$estimate[interaction_rows],
    primary_sem_term_results$coefficient_estimate[interaction_rows],
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    primary_sem_term_results$std_error[interaction_rows],
    primary_sem_term_results$coefficient_std_error[interaction_rows],
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    primary_sem_term_results$p_value[interaction_rows],
    primary_sem_term_results$coefficient_p_value[interaction_rows],
    tolerance = 1e-12
  ))
)

# ------------------------------------------------------------------------------
# 5. Component-model R2 summaries ####
# ------------------------------------------------------------------------------

sem_component_models <- list(
  parasitism = prt_sem,
  realized_parasitism = pa_sem,
  parasitoid_richness = pr_sem,
  host_richness = hr_sem,
  host_abundance = ha_sem,
  stand_volume = sv_sem,
  tree_fd = fd_sem
)

sem_component_data <- list(
  parasitism = sem_data_parasitism,
  realized_parasitism = sem_data_parasitism,
  parasitoid_richness = sem_data_parasitism,
  host_richness = sem_data_master,
  host_abundance = sem_data_master,
  stand_volume = sem_data_master,
  tree_fd = sem_data_fd
)

# Reuse one plot-to-fold assignment across all native component datasets.
component_r2_settings <- list(folds = 5L, seed = 1201L)
plot_fold_registry <- sem_data_master |>
  distinct(site, plot_id)
plot_fold_registry$fold <- make_grouped_sem_folds(
  plot_fold_registry,
  group = "plot_id",
  strata = "site",
  folds = component_r2_settings$folds,
  seed = component_r2_settings$seed
)
plot_fold_key <- paste(plot_fold_registry$site, plot_fold_registry$plot_id)

# Calculate marginal and conditional Nakagawa R2 for each component.
component_nakagawa_r2 <- lapply(names(sem_component_models), function(component) {
  sem_component_nakagawa_r2(
    sem_component_models[[component]],
    sem_component_data[[component]]
  )
}) |>
  stats::setNames(names(sem_component_models)) |>
  dplyr::bind_rows(.id = "component")

# Refit components by whole-plot fold for predictive response-scale R2.
component_cv_results <- lapply(names(sem_component_models), function(component) {
  component_data <- sem_component_data[[component]]
  fold_id <- plot_fold_registry$fold[
    match(paste(component_data$site, component_data$plot_id), plot_fold_key)
  ]
  sem_component_predictive_r2(
    model = sem_component_models[[component]],
    data = component_data,
    group = "plot_id",
    strata = "site",
    fold_id = fold_id,
    update_arguments = list(control = strict_control),
    on_refit_failure = "record"
  )
}) |>
  stats::setNames(names(sem_component_models))

component_predictive_r2 <- lapply(
  component_cv_results,
  `[[`,
  "summary"
) |>
  dplyr::bind_rows(.id = "component")

component_cv_fold_stability <- lapply(
  component_cv_results,
  `[[`,
  "folds"
) |>
  dplyr::bind_rows(.id = "component") |>
  dplyr::group_by(component) |>
  dplyr::summarise(
    n_folds = dplyr::n(),
    n_folds_available = sum(refit_status == "available"),
    n_folds_refit_failed = sum(refit_status == "failed"),
    r2_response_cv_fold_min = if (any(refit_status == "available")) {
      min(r2_response[refit_status == "available"])
    } else {
      NA_real_
    },
    r2_response_cv_fold_max = if (any(refit_status == "available")) {
      max(r2_response[refit_status == "available"])
    } else {
      NA_real_
    },
    n_negative_r2_folds = sum(r2_response < 0, na.rm = TRUE),
    .groups = "drop"
  )

component_r2_table <- component_nakagawa_r2 |>
  left_join(
    component_predictive_r2,
    by = c("component", "family", "link", "n")
  ) |>
  left_join(
    component_cv_fold_stability,
    by = "component"
  )

allowed_cv_unavailable <- "stand_volume"
unexpected_cv_failures <- component_r2_table |>
  dplyr::filter(
    cv_status != "available",
    !component %in% allowed_cv_unavailable
  )
if (nrow(unexpected_cv_failures) > 0L) {
  stop(
    "Unexpected grouped-CV refit failure(s): ",
    paste(unexpected_cv_failures$component, collapse = ", "),
    "."
  )
}

stopifnot(
  identical(component_r2_table$component, names(sem_component_models)),
  all(is.finite(component_r2_table$r2_marginal)),
  all(is.finite(component_r2_table$r2_conditional)),
  all(component_r2_table$r2_conditional >= component_r2_table$r2_marginal),
  all(component_r2_table$cv_status %in% c("available", "unavailable")),
  all(
    is.finite(component_r2_table$r2_response_cv[
      component_r2_table$cv_status == "available"
    ])
  ),
  all(
    is.na(component_r2_table$r2_response_cv[
      component_r2_table$cv_status == "unavailable"
    ])
  )
)

# ------------------------------------------------------------------------------
# 6. FD--stand-volume plot-level association ####
# ------------------------------------------------------------------------------

# Correlate FD residuals with the SV random intercept, which is each plot's
# adjusted SV deviation at centred age. This uses all SV years without copying
# the 88 FD observations into a 792-row correlation.
sv_plot_intercepts <- glmmTMB::ranef(sv_sem)$cond$plot_id |>
  as.data.frame() |>
  tibble::rownames_to_column("plot_id") |>
  dplyr::transmute(
    plot_id = as.character(plot_id),
    sv_plot_intercept = .data[["(Intercept)"]]
  )

fd_sv_association_data <- sem_data_fd |>
  dplyr::transmute(
    plot_id = as.character(plot_id),
    fd_residual = residuals(fd_sem, type = "response")
  ) |>
  dplyr::left_join(sv_plot_intercepts, by = "plot_id")

primary_fd_sv_association <- fd_sv_association_data |>
  dplyr::summarise(
    n_plots = dplyr::n(),
    correlation = stats::cor(fd_residual, sv_plot_intercept),
    conditioning = "FD: tree richness and site; SV: tree richness, age, site, and year"
  )

stopifnot(
  nrow(fd_sv_association_data) == nrow(sem_data_fd),
  !anyNA(fd_sv_association_data),
  nrow(primary_fd_sv_association) == 1L,
  primary_fd_sv_association$n_plots == 88L,
  is.finite(primary_fd_sv_association$correlation)
)

# Append the nondirectional association to the SEM table
primary_sem_term_results <- primary_sem_term_results |>
  tibble::add_row(
    component = "tree_fd_stand_volume_association",
    model_id = "fd_sem__sv_sem",
    data_id = "sem_data_fd__sv_sem_plot_intercepts",
    response = "rtfd",
    predictor = "rtsv",
    term_role = "nondirectional_association",
    native_n = primary_fd_sv_association$n_plots,
    inference_estimand = "plot_level_correlation",
    inference_method = paste(
      "Pearson correlation of FD response residuals and",
      "SV conditional plot random intercepts"
    ),
    reference_conditioning = primary_fd_sv_association$conditioning,
    estimate = primary_fd_sv_association$correlation,
    edge_id = "tree_fd__stand_volume",
    source_node = "tree_fd",
    destination_node = "stand_volume",
    standardization_status =
      "not applicable: nondirectional association"
  )

fd_sv_sem_term <- primary_sem_term_results |>
  dplyr::filter(term_role == "nondirectional_association")
stopifnot(
  nrow(primary_sem_term_results) == nrow(primary_edge_table) + 1L,
  nrow(fd_sv_sem_term) == 1L,
  fd_sv_sem_term$native_n == 88L,
  identical(
    fd_sv_sem_term$estimate,
    primary_fd_sv_association$correlation
  ),
  all(is.na(fd_sv_sem_term[c(
    "std_error", "conf_low", "conf_high", "df", "test_statistic", "p_value",
    "standardized_response_effect", "standardized_oe_effect",
    "standardized_relevant_range_effect"
  )]))
)

# ------------------------------------------------------------------------------
# 7. Realized-parasitism maximal SEM ####
# ------------------------------------------------------------------------------

# The realized-parasitism SEM differs from the primary SEM only in its terminal
# response. Reuse all upstream component results and calculate only the pa_sem
# terminal terms, edges, and standardizations.
realized_terminal <- build_sem_component_results(
  model = pa_sem,
  component = "realized_parasitism",
  model_id = "pa_sem",
  data_id = "sem_data_parasitism",
  data_sets = node_reference_datasets,
  node_reference = node_reference_table,
  reference_values = reference_values,
  include_relevant_range = relevant_ranges_ready
)

realized_sem_models <- c(
  realized_terminal$model,
  primary_sem_models[setdiff(names(primary_sem_models), "prt_sem")]
)
realized_component_metadata <- dplyr::bind_rows(
  realized_terminal$component_metadata,
  primary_component_metadata |>
    dplyr::filter(component != "parasitism")
)
realized_edge_spec <- dplyr::bind_rows(
  realized_terminal$edge_spec,
  primary_edge_spec |>
    dplyr::filter(destination_node != "parasitism")
)
realized_standardization_spec <- dplyr::bind_rows(
  realized_terminal$standardization_spec,
  primary_standardization_spec |>
    dplyr::filter(destination_node != "parasitism")
)
realized_standardized_effects <- dplyr::bind_rows(
  realized_terminal$standardized_effects,
  primary_standardized_effects |>
    dplyr::filter(destination_node != "parasitism")
)
realized_sem_term_results <- dplyr::bind_rows(
  realized_terminal$term_results,
  primary_sem_term_results |>
    dplyr::filter(component != "parasitism")
)
primary_component_r2_table <- component_r2_table |>
  dplyr::filter(component != "realized_parasitism")
primary_component_cv_fold_stability <- component_cv_fold_stability |>
  dplyr::filter(component != "realized_parasitism")
realized_component_r2_table <- component_r2_table |>
  dplyr::filter(component != "parasitism")
realized_component_cv_fold_stability <- component_cv_fold_stability |>
  dplyr::filter(component != "parasitism")

stopifnot(
  identical(
    names(realized_sem_models),
    realized_component_metadata$source_model
  ),
  identical(
    primary_component_r2_table$component,
    primary_component_metadata$component
  ),
  identical(
    realized_component_r2_table$component,
    realized_component_metadata$component
  ),
  !"parasitism" %in% realized_sem_term_results$component,
  "realized_parasitism" %in% realized_sem_term_results$component,
  !anyDuplicated(realized_edge_spec$edge_id),
  setequal(
    realized_standardization_spec$edge_id,
    realized_standardized_effects$edge_id
  ),
  sum(realized_sem_term_results$term_role ==
        "nondirectional_association") == 1L
)

# ------------------------------------------------------------------------------
# 8. Machine-readable path-analysis handoffs ####
# ------------------------------------------------------------------------------

# Keep most important columns
primary_sem_term_results_export <- primary_sem_term_results |>
  select(
    component,
    model_id,
    data_id,
    response,
    predictor,
    term_role,
    native_n,
    family,
    link,
    inference_estimand,
    inference_method,
    reference_conditioning,
    estimate,
    std_error,
    conf_low,
    conf_high,
    test_statistic,
    p_value,
    edge_id,
    source_node,
    destination_node,
    age_moderated,
    standardized_response_effect,
    standardized_oe_effect,
    standardized_relevant_range_effect
  )

realized_sem_term_results_export <- realized_sem_term_results |>
  select(
    component,
    model_id,
    data_id,
    response,
    predictor,
    term_role,
    native_n,
    family,
    link,
    inference_estimand,
    inference_method,
    reference_conditioning,
    estimate,
    std_error,
    conf_low,
    conf_high,
    test_statistic,
    p_value,
    edge_id,
    source_node,
    destination_node,
    age_moderated,
    standardized_response_effect,
    standardized_oe_effect,
    standardized_relevant_range_effect
  )

maximal_sem_handoff <- list(
  sem_variant = "parasitism_rate",
  primary_sem_models = primary_sem_models,
  node_reference_datasets = node_reference_datasets,
  node_reference_table = node_reference_table,
  reference_values = reference_values,
  relevant_ranges_ready = relevant_ranges_ready,
  primary_component_metadata = primary_component_metadata,
  primary_sem_term_results = primary_sem_term_results,
  primary_edge_spec = primary_edge_spec,
  primary_standardization_spec = primary_standardization_spec,
  primary_standardized_effects = primary_standardized_effects,
  component_r2_table = primary_component_r2_table,
  component_cv_fold_stability = primary_component_cv_fold_stability,
  primary_fd_sv_association = primary_fd_sv_association
)

maximal_realized_sem_handoff <- list(
  sem_variant = "realized_parasitism",
  primary_sem_models = realized_sem_models,
  node_reference_datasets = node_reference_datasets,
  node_reference_table = node_reference_table,
  reference_values = reference_values,
  relevant_ranges_ready = relevant_ranges_ready,
  primary_component_metadata = realized_component_metadata,
  primary_sem_term_results = realized_sem_term_results,
  primary_edge_spec = realized_edge_spec,
  primary_standardization_spec = realized_standardization_spec,
  primary_standardized_effects = realized_standardized_effects,
  component_r2_table = realized_component_r2_table,
  component_cv_fold_stability = realized_component_cv_fold_stability,
  primary_fd_sv_association = primary_fd_sv_association
)

if (!dir.exists(maximal_sem_output_dir)) {
  dir.create(maximal_sem_output_dir, recursive = TRUE)
}
if (!dir.exists(maximal_sem_output_dir)) {
  stop("Could not create the maximal-SEM output directory.")
}

saveRDS(
  maximal_sem_handoff,
  file.path(maximal_sem_output_dir, "maximal_sem_handoff.rds")
)
utils::write.csv(
  primary_sem_term_results_export,
  file.path(maximal_sem_output_dir, "maximal_sem_terms.csv"),
  row.names = FALSE
)
utils::write.csv(
  primary_component_r2_table,
  file.path(maximal_sem_output_dir, "maximal_sem_component_r2.csv"),
  row.names = FALSE
)
saveRDS(
  maximal_realized_sem_handoff,
  file.path(maximal_sem_output_dir, "maximal_realized_sem_handoff.rds")
)
utils::write.csv(
  realized_sem_term_results_export,
  file.path(maximal_sem_output_dir, "maximal_realized_sem_terms.csv"),
  row.names = FALSE
)
utils::write.csv(
  realized_component_r2_table,
  file.path(
    maximal_sem_output_dir,
    "maximal_realized_sem_component_r2.csv"
  ),
  row.names = FALSE
)

message("Maximal SEM outputs exported to: ", maximal_sem_output_dir)
message("Primary and realized-parasitism maximal SEMs ran successfully ✔️")
