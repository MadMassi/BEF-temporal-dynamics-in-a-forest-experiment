## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Interaction-network structural equation models
# Author: Massimo Martini

# Script description: This script builds the maximal and parsimonious
# interaction-network SEMs for parasitism and realized parasitism. It fits the
# total- and focal-richness, niche-overlap, and terminal-response components,
# calculates standardized pathway effects and component R2, applies the
# recorded reductions, and evaluates the reduced graphs with d-separation.
# Hierarchical parametric bootstraps provide model-based uncertainty for the
# bounded focal-richness component. A final sensitivity analysis repeats the
# focal-richness adjustment using the full parasitism dataset. Complete RDS
# handoffs and readable CSV tables are saved for every analysis branch.

# Run this script once from the project root after Maximal_models.R. By default,
# the first run performs two 1,000-replicate focal-richness bootstraps; later
# runs validate and reuse the saved checkpoints.

# Network metrics are available for plot-years with at least two host species,
# two parasitoid species, and two trophic links. Their existing scaling and the
# network-eligibility rule are retained throughout.

# Reproducibility: run with R 4.6.1 and the package versions in renv.lock.

# Approximate runtime: 4-5 hours when quick_review <- FALSE
#                      30 minutes when quick_review <- TRUE

# Keep TRUE to use 100 bootstrap replicates for a quicker review of the
# workflow. Use FALSE for the publication analysis, which uses 1,000.
# Review-mode results are saved separately and are not intended for inference.
quick_review <- TRUE

if (!is.logical(quick_review) || length(quick_review) != 1L ||
    is.na(quick_review)) {
  stop("`quick_review` must be either TRUE or FALSE.")
}
bootstrap_replicates <- if (quick_review) 100L else 1000L

# ------------------------------------------------------------------------------
# 1. Load libraries ####
# ------------------------------------------------------------------------------

library(dplyr)
library(glmmTMB)
library(gamlss)
library(gamlss.dist)

# ------------------------------------------------------------------------------
# 2. Read prepared analysis data and prepare the modelling dataset ####
# ------------------------------------------------------------------------------

model_output_dir <- "model_outputs"
design_model_output_dir <- file.path(model_output_dir, "design_models")
network_output_dir <- file.path(model_output_dir, "network_sem")
if (quick_review) {
  network_output_dir <- file.path(network_output_dir, "quick_review")
  message(
    "Quick-review mode: using 100 bootstrap replicates and writing to ",
    network_output_dir
  )
}

source("analysis_functions.R")

analysis_data_path <- file.path(design_model_output_dir, "analysis_data.rds")
if (!file.exists(analysis_data_path)) {
  stop("Missing analysis_data.rds. Run Design_models.R first.")
}
analysis_data <- readRDS(analysis_data_path)

required_datasets <- c("pd_all", "pd_net", "pd_parasitism")
missing_datasets <- setdiff(required_datasets, names(analysis_data))
if (length(missing_datasets) > 0L) {
  stop(
    "analysis_data.rds is missing: ",
    paste(missing_datasets, collapse = ", "),
    ". Rerun Design_models.R."
  )
}

# Use the saved network-eligible sample; do not reapply the network threshold
# using community richness, which is distinct from richness within networks.
# Keep exact zero niche-overlap values for response-family assessment.
pd_net <- analysis_data$pd_net |>
  arrange(plot_id, year) |>
  droplevels()

# focal_rich counts the four focal parasitoid species present (0--4).
# par_rich is total observed parasitoid richness; par_rich10 is the capped
# predictor used in the existing parasitism models. All remain available here.

strict_control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000L, eval.max = 10000L)
)

# ------------------------------------------------------------------------------
# 3. Component model assembly ####
# ------------------------------------------------------------------------------

# Parasitism
prt_sem <- glmmTMB(
  cell_prate ~
    sc_logtr + par_rich + focal_rich + niche_overlap_hl +
    sc_fa + sc_hr + sc_rtcells +
    sc_logtr:sc_fa +
    par_rich:sc_fa +
    focal_rich:sc_fa +
    sc_hr:sc_fa +
    sc_rtcells:sc_fa +
    (1 | site / plot_id) + (1 | year),
  weights = total_cells,
  data = pd_net,
  dispformula = ~sc_cells,
  family = betabinomial(link = "logit"),
  control = strict_control
)

# Realized parasitism (parasitoid abundance)
pa_sem <- glmmTMB(
  par_abund ~
    sc_logtr + par_rich + focal_rich + niche_overlap_hl +
    sc_fa + sc_hr + sc_rtcells +
    sc_logtr:sc_fa +
    par_rich:sc_fa +
    focal_rich:sc_fa +
    sc_hr:sc_fa +
    sc_rtcells:sc_fa +
    site + (1 | plot_id) + (1 | year),
  data = pd_net,
  family = nbinom2(link = "log"),
  ziformula = ~ 0,
  dispformula = ~ 1,
  control = strict_control
)

# Niche overlap
niche_sem <- glmmTMB(
  niche_overlap_hl ~
    sc_logtr + par_rich + focal_rich + sc_fa +
    sc_hr + sc_rtcells +
    sc_logtr:sc_fa +
    par_rich:sc_fa +
    focal_rich:sc_fa +
    sc_hr:sc_fa +
    sc_rtcells:sc_fa +
    (1 | site / plot_id) + (1 | year),
  data = pd_net,
  family = beta_family(link = "logit"),
  ziformula = ~1,
  dispformula = ~ par_rich + (1 | plot_id) + (1 | year),
  control = strict_control
) # plot and year account for residual heterogeneity in dispersion

# Focal parasitoid species richness
# The double-binomial family permits underdispersion and respects
# the observation-specific maximum of min(4, par_rich).
# Dedicated dataset containing only variables used by this component.
focal_data <- as.data.frame(pd_net)[, c(
  "focal_rich", "par_rich", "sc_logtr", "sc_fa",
  "sc_hr", "sc_rtcells", "site", "plot_id", "year"
), drop = FALSE]

focal_data$focal_max <- pmin(4L, focal_data$par_rich)
focal_data$plot_id <- factor(focal_data$plot_id)
focal_data$year <- factor(focal_data$year)

# Preserve the sample: do not silently remove observations.
if (anyNA(focal_data)) {
  print(colSums(is.na(focal_data)))
  stop("Missing values in focal-model variables; inspect before fitting.")
}

# Focal richness conditional on total parasitoid richness.
# Double binomial allows underdispersion within the biological bounds.
frich_sem <- gamlss(
  cbind(focal_rich, focal_max - focal_rich) ~
    sc_logtr + par_rich + sc_fa +
    sc_hr + sc_rtcells + site +
    sc_logtr:sc_fa +
    par_rich:sc_fa +
    sc_hr:sc_fa +
    sc_rtcells:sc_fa +
    random(plot_id) + random(year),
  sigma.formula = ~1,
  family = DBI(mu.link = "logit", sigma.link = "log"),
  sigma.start = 0.4,
  data = focal_data,
  control = gamlss.control(
    n.cyc = 300,
    c.crit = 1e-6,
    trace = FALSE
  )
)

# Hierarchical bootstrap of focal-component coefficients on the location-logit
# scale. The checkpoint preserves all simulated responses, retries, and fit
# audits. Matching reruns reuse it; a mismatch stops without overwriting it.
focal_bootstrap_settings <- list(
  times = bootstrap_replicates, seed = 123L, confidence = 0.95,
  minimum_success = 0.95, progress_every = 50L,
  cycle_limits = c(300L, 1000L, 1000L), tolerance = 1e-6
)
focal_bootstrap_reference_dir <- file.path(
  network_output_dir, "focal_richness_bootstrap",
  paste0("maximal_", focal_bootstrap_settings$times)
)
focal_bootstrap_reference_file <- file.path(
  focal_bootstrap_reference_dir, "focal_dbi_bootstrap.rds"
)
focal_bootstrap_signature <- list(
  settings = focal_bootstrap_settings,
  bootstrap_function = deparse(bootstrap_focal_dbi),
  success_function = deparse(.sem_required_successes),
  package_versions = vapply(c("gamlss", "gamlss.dist"), function(pkg) {
    as.character(utils::packageVersion(pkg))
  }, character(1))
)
frich_bootstrap <- if (file.exists(focal_bootstrap_reference_file)) {
  saved <- readRDS(focal_bootstrap_reference_file)
  matches <- tryCatch({
    stopifnot(identical(saved$cache_signature, focal_bootstrap_signature))
    # Checks formula, ordered data, location coefficients and dispersion.
    sem_focal_dbi_inference_adapter(frich_sem, focal_data, saved)
    # Also check fitted random contributions and their generating variances.
    stopifnot(
      identical(frich_sem$family, saved$reference_model$family),
      identical(frich_sem$mu.link, saved$reference_model$mu.link),
      identical(frich_sem$sigma.link, saved$reference_model$sigma.link),
      isTRUE(frich_sem$converged),
      isTRUE(all.equal(frich_sem$mu.fv, saved$reference_model$mu.fv,
        tolerance = 1e-8)),
      isTRUE(all.equal(
        vapply(1:2, function(j) gamlss::getSmo(frich_sem, "mu", which = j)$sigb,
          numeric(1)),
        unname(c(saved$generating_parameters$plot_SD,
          saved$generating_parameters$year_SD)), tolerance = 1e-8
      ))
    )
    TRUE
  }, error = function(e) FALSE)
  if (!matches) stop(
    paste0(
      "Saved focal bootstrap does not match the current model, ",
      "data, settings or implementation. "
    ),
    "Reference preserved at: ", focal_bootstrap_reference_file,
    paste0(
      ". Inspect the mismatch; choose a new reference ",
      "directory for an intentional new run."
    )
  )
  message("Reusing maximal focal bootstrap: ", focal_bootstrap_reference_file)
  saved
} else {
  result <- do.call(bootstrap_focal_dbi, c(
    list(model = frich_sem, data = focal_data), focal_bootstrap_settings
  ))
  result$cache_signature <- focal_bootstrap_signature
  export_focal_dbi_bootstrap(
    result, focal_bootstrap_reference_dir,
    source_paths = c(
      analysis_data_path,
      "Network_SEM.R",
      "analysis_functions.R"
    )
  )
  result
}

# Parasitoid richness
pr_sem <- glmmTMB(
  par_rich ~
    sc_logtr + sc_fa +
    sc_rtcells + sc_hr +
    sc_logtr:sc_fa +
    sc_hr:sc_fa +
    sc_rtcells:sc_fa +
    site + (1 | plot_id) + (1 | year),
  data = pd_net,
  family = genpois(link = "log"),
  ziformula = ~0,
  dispformula = ~1,
  control = strict_control
)

# ------------------------------------------------------------------------------
# 4. Node reference scales and relevant ranges ####
# ------------------------------------------------------------------------------

# Retain upstream predictor scaling. Network node SDs describe the native
# network sample; tree richness uses one observation per plot. Range endpoints
# are in biological units, with observed integer quantiles for richness nodes.
network_models <- list(prt_sem = prt_sem, niche_sem = niche_sem,
  frich_sem = frich_sem, pr_sem = pr_sem)
network_components <- c(prt_sem = "parasitism", niche_sem = "niche_overlap",
  frich_sem = "focal_richness",
  pr_sem = "parasitoid_richness")
network_data_ids <- c(prt_sem = "pd_net", niche_sem = "pd_net",
  frich_sem = "focal_data", pr_sem = "pd_net")
node_reference_datasets <- list(
  pd_net = pd_net, focal_data = focal_data,
  network_plots = pd_net |> distinct(plot_id, .keep_all = TRUE)
)

node_reference_spec <- tibble::tribble(
  ~node, ~model_variable, ~biological_variable, ~transformation,
  ~reference_dataset, ~scaling_base_variable,
  "tree_richness", "sc_logtr", "tree_r", "cap at 16; log2; z-scale",
  "network_plots", "log_tr",
  "parasitoid_richness", "par_rich", "par_rich", "identity", "pd_net",
  NA_character_,
  "focal_richness", "focal_rich", "focal_rich", "identity", "pd_net",
  NA_character_,
  "niche_overlap", "niche_overlap_hl", "niche_overlap_hl",
  "identity (proportion)", "pd_net", NA_character_,
  "host_richness", "sc_hr", "host_rich", "z-scale", "pd_net", "host_rich",
  "host_abundance", "sc_rtcells", "total_cells", "square root; z-scale",
  "pd_net", "rtcells",
  "parasitism", "cell_prate", "cell_prate", "identity (proportion)", "pd_net",
  NA_character_,
  "realized_parasitism", "par_abund", "par_abund", "identity", "pd_net",
  NA_character_,
  "forest_age", "sc_fa", "for_age", "z-scale", "pd_net", "for_age"
)
node_reference_table <-
  dplyr::bind_rows(lapply(seq_len(nrow(node_reference_spec)), function(i) {
    spec <- node_reference_spec[i, ]
    data <- node_reference_datasets[[spec$reference_dataset]]
    x <- data[[spec$model_variable]]
    biological <- data[[spec$biological_variable]]
    base <- spec$scaling_base_variable
    center <- if (is.na(base)) NA_real_ else mean(analysis_data$pd_all[[base]],
      na.rm = TRUE)
    divisor <- if (is.na(base)) NA_real_ else
      stats::sd(analysis_data$pd_all[[base]], na.rm = TRUE)
    stopifnot(!anyNA(x), all(is.finite(x)), stats::sd(x) > 0)
    # Check the stored full-data scaling without rescaling the network subset.
    if (!is.na(base)) stopifnot(isTRUE(all.equal(
      x, (data[[base]] - center) / divisor, check.attributes = FALSE
    )))
    dplyr::mutate(
      spec,
      observational_grain =
        if (reference_dataset == "network_plots") "unique plot" else
          "network plot-year",
      n_reference = length(x), mean_model = mean(x), sd_model = stats::sd(x),
      observed_low_model = min(x), observed_high_model = max(x),
      observed_low_biological = min(biological),
      observed_high_biological = max(biological),
      scaling_center = center, scaling_divisor = divisor
    )
  }))

relevant_range_spec <-
  dplyr::bind_rows(lapply(seq_len(nrow(node_reference_spec)), function(i) {
    spec <- node_reference_spec[i, ]
    values <- node_reference_datasets[[spec$reference_dataset]][[
      spec$biological_variable
    ]]
    integer_richness <- spec$node %in% c(
      "parasitoid_richness", "focal_richness", "host_richness"
    )
    limits <- stats::quantile(values, c(0.1, 0.9), names = FALSE,
      type = if (integer_richness) 1L else 7L)
    rationale <-
      if (integer_richness) "10th-90th observed quantiles (type 1)" else
        "10th-90th percentiles (type 7)"
    if (spec$node == "tree_richness") {
      limits <- c(1, 16)
      rationale <- "experimental limits"
    }
    if (spec$node == "forest_age") {
      limits <- c(4, 14)
      rationale <- "experimental limits"
    }
    tibble::tibble(node = spec$node, rr_low_biological = limits[1],
      rr_high_biological = limits[2], rr_rationale = rationale)
  }))
range_setup <- prepare_sem_relevant_ranges(node_reference_table,
  relevant_range_spec)
node_reference_table <- range_setup$node_reference
relevant_ranges_ready <- range_setup$ready
stopifnot(relevant_ranges_ready, all(node_reference_table$rr_supported))
reference_values <- stats::setNames(node_reference_table$mean_model,
  node_reference_table$model_variable)
reference_values[
  node_reference_table$model_variable[
    !is.na(node_reference_table$scaling_base_variable)]] <- 0

# ------------------------------------------------------------------------------
# 5. SEM terms and three standardized effect measures ####
# ------------------------------------------------------------------------------

if (!is.null(frich_bootstrap)) {

  # Inference stays on the fitted link/location scale. Standardized effects are
  # point estimates from equal-row averages with random effects set to zero.
  # Focal inference uses bootstrap draws; percentile intervals are not
  # converted to P values.
  focal_adapter <- sem_focal_dbi_inference_adapter(frich_sem, focal_data,
    frich_bootstrap)
  network_component_results <- stats::setNames(lapply(names(network_models),
    function(id) {
      model <- network_models[[id]]
      focal <- id == "frich_sem"
      response <- all.vars(stats::formula(model))[[1L]]
      result <- build_sem_component_results(
        model = model, component = network_components[[id]], model_id = id,
        data_id = network_data_ids[[id]], data_sets = node_reference_datasets,
        node_reference = node_reference_table,
        reference_values = reference_values,
        include_relevant_range = relevant_ranges_ready,
        coefficient_table = if (focal) focal_adapter$coefficients else NULL,
        inference_function = if (focal) focal_adapter$inference else
          sem_term_wald_test,
        predict_function = if (focal) sem_predict_focal_dbi else NULL,
        evaluation_function = sem_network_evaluation(response)
      )
      terms <- result$term_results
      terms$response_estimand <- if (focal) {
        "expected focal count from exact DBI probabilities"
      } else if (id == "niche_sem") {
        "overall expected overlap including exact zeros"
      } else if (id == "prt_sem") {
        "expected parasitism proportion"
      } else "expected parasitoid count"
      terms$oe_estimand <- if (focal) {
        "DBI location-logit OE extension using expected counts"
      } else if (id == "niche_sem") {
          paste0(
            "positive-overlap logit OE with overall response ",
            "correlation; constant zero probability"
          )
        } else "Menard observed-empirical link-scale coefficient"
      terms$effect_method <- ifelse(
        terms$term_role %in% c("age_moderation", "ecological_interaction"),
        "not calculated: interaction",
        "average numerical partial derivative; feasible one-sided limits"
      )
      if (focal) terms$effect_method[terms$predictor == "par_rich"] <-
        paste0(
          "integer-support secant per species; +/-1 internally, ",
          "one-sided at observed endpoints; bound updated"
        )
      terms$uncertainty_status <- if (focal) paste(
        focal_adapter$interval_status, ";", frich_bootstrap$summary$accepted,
        paste0(
          "accepted bootstrap draws; P values unavailable; ",
          "standardized effects are point estimates"
        )
      ) else paste0(
        "95% link-scale Wald intervals; standardized effects are ",
        "point estimates"
      )
      result$term_results <- terms
      result
    }), names(network_models))

  network_sem_term_results <- dplyr::bind_rows(lapply(network_component_results,
    `[[`, "term_results"))
  network_path_table <- network_sem_path_table(network_component_results)
  network_edge_spec <- dplyr::bind_rows(lapply(network_component_results, `[[`,
    "edge_spec"))
  network_standardized_effects <-
    dplyr::bind_rows(lapply(network_component_results, `[[`,
      "standardized_effects"))
  stopifnot(
    !anyDuplicated(network_edge_spec$edge_id),
    identical(network_path_table$edge_id, network_edge_spec$edge_id),
    all(network_sem_term_results$native_n == nrow(pd_net)),
    all(is.finite(network_standardized_effects$standardized_response_effect)),
    all(is.finite(network_standardized_effects$standardized_oe_effect)),
    all(is.finite(
      network_standardized_effects$standardized_relevant_range_effect
    ))
  )

  # The secondary maximal SEM differs only in its terminal response. Reuse the
  # three upstream network components and their focal-bootstrap inference, and
  # calculate only the realized-parasitism terminal results.
  realized_terminal <- build_sem_component_results(
    model = pa_sem,
    component = "realized_parasitism",
    model_id = "pa_sem",
    data_id = "pd_net",
    data_sets = node_reference_datasets,
    node_reference = node_reference_table,
    reference_values = reference_values,
    include_relevant_range = relevant_ranges_ready,
    evaluation_function = sem_network_evaluation("cell_prate")
  )
  realized_terminal$term_results$response_estimand <-
    "expected parasitoid abundance (realized parasitism)"
  realized_terminal$term_results$oe_estimand <-
    "Menard observed-empirical link-scale coefficient"
  realized_terminal$term_results$effect_method <- ifelse(
    realized_terminal$term_results$term_role %in%
      c("age_moderation", "ecological_interaction"),
    "not calculated: interaction",
    "average numerical partial derivative; feasible one-sided limits"
  )
  realized_terminal$term_results$uncertainty_status <- paste0(
    "95% link-scale Wald intervals; standardized effects are point estimates"
  )

  realized_network_component_results <- c(
    list(pa_sem = realized_terminal),
    network_component_results[c("niche_sem", "frich_sem", "pr_sem")]
  )
  realized_network_models <- c(
    realized_terminal$model,
    network_models[c("niche_sem", "frich_sem", "pr_sem")]
  )
  realized_network_sem_term_results <- dplyr::bind_rows(lapply(
    realized_network_component_results, `[[`, "term_results"
  ))
  realized_network_path_table <- network_sem_path_table(
    realized_network_component_results
  )
  realized_network_edge_spec <- dplyr::bind_rows(lapply(
    realized_network_component_results, `[[`, "edge_spec"
  ))
  realized_network_standardized_effects <- dplyr::bind_rows(lapply(
    realized_network_component_results, `[[`, "standardized_effects"
  ))
  realized_network_component_metadata <- dplyr::bind_rows(lapply(
    realized_network_component_results, `[[`, "component_metadata"
  ))
  realized_network_standardization_spec <- dplyr::bind_rows(lapply(
    realized_network_component_results, `[[`, "standardization_spec"
  ))

  stopifnot(
    identical(
      names(realized_network_models),
      realized_network_component_metadata$source_model
    ),
    !"parasitism" %in% realized_network_sem_term_results$component,
    "realized_parasitism" %in% realized_network_sem_term_results$component,
    !anyDuplicated(realized_network_edge_spec$edge_id),
    identical(
      realized_network_path_table$edge_id,
      realized_network_edge_spec$edge_id
    ),
    setequal(
      realized_network_standardization_spec$edge_id,
      realized_network_standardized_effects$edge_id
    ),
    all(realized_network_sem_term_results$native_n == nrow(pd_net)),
    all(is.finite(
      realized_network_standardized_effects$standardized_response_effect
    )),
    all(is.finite(
      realized_network_standardized_effects$standardized_oe_effect
    )),
    all(is.finite(
      realized_network_standardized_effects$standardized_relevant_range_effect
    ))
  )

  # ----------------------------------------------------------------------------
  # 6. Component R2 and reproducible handoff ####
  # ----------------------------------------------------------------------------

  # DBI and constant-zero-inflated beta components use exact response moments
  # and delta-method mean derivatives;
  # their link coordinates and approximation limits are recorded in the output.
  # No residual diagnostics, grouped CV, graph reduction or global-fit tests
  # here.
  network_component_r2 <- dplyr::bind_rows(lapply(names(network_models),
    function(id) {
      sem_network_component_r2(network_models[[id]],
        node_reference_datasets[[network_data_ids[[id]]]]) |>
        dplyr::mutate(component = network_components[[id]], model_id = id,
          data_id = network_data_ids[[id]], .before = 1)
    }))
  realized_terminal_r2 <- sem_network_component_r2(pa_sem, pd_net) |>
    dplyr::mutate(
      component = "realized_parasitism",
      model_id = "pa_sem",
      data_id = "pd_net",
      .before = 1
    )
  realized_network_component_r2 <- dplyr::bind_rows(
    realized_terminal_r2,
    network_component_r2 |>
      dplyr::filter(component != "parasitism")
  )
  stopifnot(
    identical(
      realized_network_component_r2$model_id,
      names(realized_network_models)
    ),
    all(is.finite(realized_network_component_r2$r2_marginal)),
    all(is.finite(realized_network_component_r2$r2_conditional))
  )
  network_source_paths <- c(
    analysis_data_path,
    "Network_SEM.R",
    "analysis_functions.R"
  )
  network_handoff <- list(
    sem_variant = "interaction_network", primary_sem_models = network_models,
    node_reference_datasets = node_reference_datasets,
    node_reference_table = node_reference_table,
    reference_values = reference_values,
    relevant_ranges_ready = relevant_ranges_ready,
    primary_component_metadata =
      dplyr::bind_rows(lapply(network_component_results, `[[`,
        "component_metadata")),
    primary_sem_term_results = network_sem_term_results,
    primary_edge_spec = network_edge_spec,
    primary_standardization_spec =
      dplyr::bind_rows(lapply(network_component_results, `[[`,
        "standardization_spec")),
    primary_standardized_effects = network_standardized_effects,
    path_table = network_path_table, component_r2_table = network_component_r2,
    focal_bootstrap = frich_bootstrap,
    interpretation = c(
      paste0(
        "Exploratory network pathway description conditional on ",
        "the network-eligible sample."
      ),
      paste0(
        "Focal effects at fixed total richness describe ",
        "replacement of nonfocal richness."
      ),
      paste0(
        "Total-to-focal count contrasts include changes in ",
        "min(4, total richness); the link slope describes ",
        "location only."
      ),
      paste0(
        "Relevant-range contrasts average only rows supporting ",
        "both endpoints; rr_n_excluded records restrictions."
      ),
      paste0(
        "Fixed-effect predictions set random effects to zero; ",
        "they do not integrate over random effects."
      ),
      paste0(
        "Link/location intervals do not quantify uncertainty in ",
        "standardized effects or nonlinear pathway products."
      ),
      paste0(
        "Focal bootstrap inference remains model-based and ",
        "conditional on the maximal specification; it does not ",
        "account for model selection."
      )
    ),
    source_checksums = data.frame(
      path = normalizePath(network_source_paths, mustWork = TRUE),
      md5 = unname(tools::md5sum(network_source_paths))
    ),
    session_info = utils::sessionInfo()
  )

  realized_network_handoff <- list(
    sem_variant = "interaction_network_realized_parasitism",
    primary_sem_models = realized_network_models,
    node_reference_datasets = node_reference_datasets,
    node_reference_table = node_reference_table,
    reference_values = reference_values,
    relevant_ranges_ready = relevant_ranges_ready,
    primary_component_metadata = realized_network_component_metadata,
    primary_sem_term_results = realized_network_sem_term_results,
    primary_edge_spec = realized_network_edge_spec,
    primary_standardization_spec =
      realized_network_standardization_spec,
    primary_standardized_effects = realized_network_standardized_effects,
    path_table = realized_network_path_table,
    component_r2_table = realized_network_component_r2,
    focal_bootstrap = frich_bootstrap,
    interpretation = c(
      paste0(
        "Secondary maximal network SEM differing from the primary variant ",
        "only in its terminal realized-parasitism response."
      ),
      paste0(
        "Realized parasitism is parasitoid abundance; its conditional paths ",
        "combine attack probability with realized counts while adjusting for ",
        "brood-cell abundance."
      ),
      paste0(
        "All niche-overlap, focal-richness and total-richness components, ",
        "including focal-bootstrap inference, are reused unchanged from the ",
        "primary maximal network SEM."
      ),
      paste0(
        "Focal effects at fixed total richness describe replacement of ",
        "nonfocal richness."
      ),
      paste0(
        "Relevant-range contrasts average only rows supporting both ",
        "endpoints; rr_n_excluded records restrictions."
      ),
      paste0(
        "Fixed-effect predictions set random effects to zero; they do not ",
        "integrate over random effects."
      ),
      paste0(
        "Link/location intervals do not quantify uncertainty in standardized ",
        "effects or nonlinear pathway products."
      ),
      paste0(
        "Focal bootstrap inference remains model-based and conditional on ",
        "the maximal specification; it does not account for model selection."
      )
    ),
    source_checksums = data.frame(
      path = normalizePath(network_source_paths, mustWork = TRUE),
      md5 = unname(tools::md5sum(network_source_paths))
    ),
    session_info = utils::sessionInfo()
  )

  # ----------------------------------------------------------------------------
  # 7. Export maximal network SEM results ####
  # ----------------------------------------------------------------------------

  # Keep a readable bootstrap export alongside the checkpoint and place each
  # network variant in its own subdirectory.
  export_focal_dbi_bootstrap(
    frich_bootstrap,
    file.path(network_output_dir, "focal_richness_bootstrap"),
    source_paths = network_source_paths
  )
  export_network_sem(
    network_handoff,
    file.path(network_output_dir, "maximal"),
    file_stem = "maximal_network_sem"
  )
  export_network_sem(
    realized_network_handoff,
    file.path(network_output_dir, "maximal_realized"),
    file_stem = "maximal_realized_network_sem"
  )

} # End bootstrap-dependent SEM reporting and exports.

# ------------------------------------------------------------------------------
# 8. Stepwise component reduction ####
# ------------------------------------------------------------------------------

## 8.1 Parasitism ##
parasitism_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_fa:sc_hr", "accept", FALSE,
  "LRT P = 0.899; delta AIC = -1.98; retained estimates changed by 
  <0.05 parent SE",
  2L, "par_rich:sc_fa", "accept", FALSE,
  "LRT P = 0.991; delta AIC = -2.00; maximum coefficient change =
  0.00479 parent SE.",
  3L, "sc_logtr:sc_fa", "accept", FALSE,
  "LRT P = 0.771; delta AIC = -1.92; maximum coefficient change =
  0.0383 parent SE.",
  4L, "focal_rich:sc_fa", "accept", FALSE,
  "LRT P = 0.587; delta AIC = -1.70. Maximum coefficient change = 
  0.189 parent SE but direction was unchanged
  and Wald test remained non-significant (P = 0.414 to 0.503)",
  5L, "sc_fa:sc_rtcells", "accept", FALSE,
  "LRT P = 0.503; delta AIC = -1.55. Maximum coefficient change = 
  0.181 parent SE but association remained positive
  and strongly supported. Other retained estimates changed little.",
  6L, "sc_logtr", "accept", FALSE,
  "LRT P = 0.561; delta AIC = -1.66. Maximum retained coefficient change =
  0.0482 parent SE; all directions and significance conclusions
  preserved.",
  7L, "sc_fa", "accept", TRUE,
  "LRT P = 0.596; delta AIC = -1.72. Maximum retained coefficient change =
  0.0939 parent SE; all directions preserved. Response-scale CV R2
  increased from 0.31949 to 0.32234 (delta = +0.00285)."
)

parasitism_reduction <- run_component_reduction_plan(
  component = "parasitism",
  maximal_model = network_models[["prt_sem"]],
  data = node_reference_datasets[[network_data_ids[["prt_sem"]]]],
  plan = parasitism_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

# Print the finalized reduction and selected model.
parasitism_reduction$removal_log
parasitism_reduction$comparison_history[["step_7"]]$coefficient_stability
summary(parasitism_reduction$current_model)

prt_sem_parsim <- parasitism_reduction$current_model
prt_sem_parsim_terms <- parasitism_reduction$current_term_results

## 8.1.1 Realized parasitism ##
real_parasitism_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_fa:sc_hr", "accept", FALSE,
  "LRT P = 0.832; delta AIC = -1.96;
  maximum retained-coefficient change = 0.0746 parent-model SE;
  all directions and inferential conclusions were preserved.",
  2L, "sc_logtr:sc_fa", "accept", FALSE,
  "LRT P = 0.565; delta AIC = -1.67; maximum retained
  reference-slope change = 0.0774 parent-model SE;
  all directions and inferential conclusions were preserved.",
  3L, "focal_rich:sc_fa", "accept", FALSE,
  "LRT P = 0.716; delta AIC = -1.87; maximum retained
  reference-slope change = 0.164 parent-model SE;
  all directions and inferential conclusions were preserved.",
  4L, "sc_fa:sc_rtcells", "accept", FALSE,
  "LRT P = 0.887; delta AIC = -1.98; maximum retained
  reference-effect change = 0.0700 parent-model SE.
  All directions and inferential conclusions were preserved.",
  5L, "par_rich:sc_fa", "accept", FALSE,
  "LRT P = 0.915; delta AIC = -1.99. The maximum retained-effect
  change was 0.0432 parent-model SE for parasitoid richness, which
  remained positive and strongly supported. All directions and
  inferential conclusions were preserved.",
  6L, "sc_fa", "accept", TRUE,
  "LRT P = 0.542; delta AIC = -1.63. Maximum retained-effect
  change was 0.252 parent-model SE; all directions and inferential
  conclusions were preserved. Grouped-CV R2 decreased from
  0.53727 to 0.53217 (delta = -0.00511).",
  7L, "sc_logtr", "accept", TRUE,
  "LRT P = 0.209; delta AIC = -0.420. Maximum retained-effect
  change was 0.171 parent-model SE, with all directions and
  inferential conclusions preserved. Grouped-CV R2 decreased from
  0.53217 to 0.53066 (delta = -0.00151)."
)

real_parasitism_reduction <- run_component_reduction_plan(
  component = "realized_parasitism",
  maximal_model = pa_sem,
  data = pd_net,
  plan = real_parasitism_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

real_parasitism_reduction$removal_log
real_parasitism_reduction$
  comparison_history[["step_7"]]$model_comparison
real_parasitism_reduction$
  comparison_history[["step_7"]]$likelihood_ratio_test
real_parasitism_reduction$
  comparison_history[["step_7"]]$coefficient_stability
summary(real_parasitism_reduction$current_model)

pa_sem_parsim <- real_parasitism_reduction$current_model
pa_sem_parsim_terms <- real_parasitism_reduction$current_term_results
# Retain sc_hr and sc_rtcells a priori to control for host-network breadth and
# host availability/sampling opportunity, irrespective of their Wald P values.

## 8.2 Niche overlap ##
niche_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "focal_rich:sc_fa", "accept", FALSE,
  "LRT P = 0.726; delta AIC = -1.88. Maximum coefficient change = 
  0.162 parent SE, but remained positive and supported",
  2L, "sc_fa:sc_hr", "accept", FALSE,
  "LRT P = 0.254; delta AIC = -0.70. Host richness changed by 0.438 parent SE
  but remained negative and supported. Host abundance × age changed
  by 0.428 parent SE, remaining negative and
  non-significant. Other associations remained stable.",
  3L, "sc_fa:sc_rtcells", "accept", FALSE,
  "LRT P = 0.244; delta AIC = -0.644. Parasitoid richness × age
  became more negative (change = 0.516 parent SE), with Wald P decreasing
  from 0.0467 to 0.00627. Direction was retained."
)

niche_reduction <- run_component_reduction_plan(
  component = "niche_overlap",
  maximal_model = network_models[["niche_sem"]],
  data = node_reference_datasets[[network_data_ids[["niche_sem"]]]],
  plan = niche_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

# Print the finalized reduction and selected model.
niche_reduction$removal_log
tail(niche_reduction$comparison_history, 1)[[1]]$coefficient_stability
summary(niche_reduction$current_model)

niche_sem_parsim <- niche_reduction$current_model
niche_sem_parsim_terms <- niche_reduction$current_term_results

## 8.3 Parasitoid richness ##
pr_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~run_cv, ~rationale,
  1L, "sc_fa:sc_hr", "accept", FALSE,
  "LRT P = 0.289; delta AIC = -0.876. Host abundance × age changed from -0.0859
  to -0.0942 (0.416 parent SE), remaining negative and strongly supported.
  Other retained estimates changed little."
)

pr_reduction <- run_component_reduction_plan(
  component = "parasitoid_richness",
  maximal_model = network_models[["pr_sem"]],
  data = node_reference_datasets[[network_data_ids[["pr_sem"]]]],
  plan = pr_reduction_plan,
  reference_values = reference_values,
  control = strict_control
)

pr_reduction$removal_log
tail(pr_reduction$comparison_history, 1)[[1]]$coefficient_stability
summary(pr_reduction$current_model)

pr_sem_parsim <- pr_reduction$current_model
pr_sem_parsim_terms <- pr_reduction$current_term_results

## 8.4 Focal parasitoid richness ##
# Candidate 1 removes the host-richness-by-age interaction while retaining the
# total-richness-by-age interaction for the next reduction step.
frich_no_hr_age <- update(
  frich_sem, . ~ . - sc_hr:sc_fa, data = focal_data
)
frich_comparison <- data.frame(
  model = c("maximal", "no_host_richness_age"),
  converged = c(frich_sem$converged, frich_no_hr_age$converged),
  n = c(frich_sem$noObs, frich_no_hr_age$noObs),
  AIC = c(AIC(frich_sem), AIC(frich_no_hr_age))
)
frich_comparison$delta_AIC <- frich_comparison$AIC - frich_comparison$AIC[1]
print(frich_comparison)
stopifnot(all(frich_comparison$converged),
  all(frich_comparison$n == nrow(focal_data)))

# GAMLSS summaries/SEs are approximate with random-effect smoothers.
frich_parent_summary <- summary(frich_sem)
frich_candidate_summary <- summary(frich_no_hr_age)
frich_retained <- coef(frich_no_hr_age, what = "mu")
frich_retained <- names(frich_retained)[!is.na(frich_retained)]
frich_coefficient_stability <- data.frame(
  term = frich_retained,
  parent = frich_parent_summary[frich_retained, "Estimate"],
  candidate = frich_candidate_summary[frich_retained, "Estimate"],
  parent_SE = frich_parent_summary[frich_retained, "Std. Error"],
  row.names = NULL
)
frich_coefficient_stability$change <- with(frich_coefficient_stability,
  candidate - parent)
frich_coefficient_stability$change_parent_SE <- with(
  frich_coefficient_stability, change / parent_SE
)
print(frich_coefficient_stability)

frich_reduction_plan <- tibble::tribble(
  ~step, ~removed_term, ~decision, ~rationale,
  1L, "sc_fa:sc_hr", "accept",
  "Delta AIC = +0.148; accepted for parsimony.
  Host abundance x age changed by 0.553 approximate parent SE,
  retaining its negative direction and approximate support.
  Other retained coefficients changed by <0.185 parent SE."
)

frich_sem_parsim <- frich_no_hr_age

summary(frich_sem_parsim)

# Candidate 2 removes the total-richness-by-age interaction from the accepted
# first candidate.
frich_no_par_age <- update(
  frich_sem_parsim, . ~ . - par_rich:sc_fa, data = focal_data
)
frich_comparison_2 <- data.frame(
  model = c("accepted_parent", "no_total_richness_age"),
  converged = c(frich_sem_parsim$converged, frich_no_par_age$converged),
  n = c(frich_sem_parsim$noObs, frich_no_par_age$noObs),
  AIC = c(AIC(frich_sem_parsim), AIC(frich_no_par_age))
)
frich_comparison_2$delta_AIC <-
  frich_comparison_2$AIC - frich_comparison_2$AIC[1]
print(frich_comparison_2)
stopifnot(all(frich_comparison_2$converged),
  all(frich_comparison_2$n == nrow(focal_data)))

# Approximate GAMLSS summaries and coefficient stability relative to the
# accepted parent.
frich_parent_summary_2 <- summary(frich_sem_parsim)
frich_candidate_summary_2 <- summary(frich_no_par_age)
frich_retained_2 <- coef(frich_no_par_age, what = "mu")
frich_retained_2 <- names(frich_retained_2)[!is.na(frich_retained_2)]
frich_coefficient_stability_2 <- data.frame(
  term = frich_retained_2,
  parent = frich_parent_summary_2[frich_retained_2, "Estimate"],
  candidate = frich_candidate_summary_2[frich_retained_2, "Estimate"],
  parent_SE = frich_parent_summary_2[frich_retained_2, "Std. Error"],
  row.names = NULL
)
frich_coefficient_stability_2$change <- with(frich_coefficient_stability_2,
  candidate - parent)
frich_coefficient_stability_2$change_parent_SE <- with(
  frich_coefficient_stability_2, change / parent_SE
)
print(frich_coefficient_stability_2)

frich_reduction_plan <- tibble::add_row(
  frich_reduction_plan,
  step = 2L,
  removed_term = "par_rich:sc_fa",
  decision = "accept",
  rationale = paste(
    "Delta AIC = -1.996;",
    "Maximum retained coefficient change = 0.036 approximate parent SE;",
    "all directions preserved. Approximate age SE decreased from",
    "0.1315 to 0.0518."
  )
)

frich_sem_parsim <- frich_no_par_age

# ------------------------------------------------------------------------------
# 9. Parsimonious focal-richness bootstrap ####
# ------------------------------------------------------------------------------

# Use the same replicate count and settings as the maximal reference, in a
# separate directory. Matching reruns only read this checkpoint; mismatches stop.
focal_parsim_bootstrap_dir <- file.path(
  network_output_dir, "focal_richness_bootstrap",
  paste0("parsimonious_", focal_bootstrap_settings$times)
)
focal_parsim_bootstrap_file <- file.path(focal_parsim_bootstrap_dir,
  "focal_dbi_bootstrap.rds")
frich_parsim_bootstrap <- if (file.exists(focal_parsim_bootstrap_file)) {
  saved <- readRDS(focal_parsim_bootstrap_file)
  matches <- tryCatch({
    stopifnot(identical(saved$cache_signature, focal_bootstrap_signature))
    sem_focal_dbi_inference_adapter(frich_sem_parsim, focal_data, saved)
    stopifnot(
      identical(frich_sem_parsim$family, saved$reference_model$family),
      identical(frich_sem_parsim$mu.link, saved$reference_model$mu.link),
      identical(frich_sem_parsim$sigma.link, saved$reference_model$sigma.link),
      isTRUE(frich_sem_parsim$converged),
      isTRUE(all.equal(frich_sem_parsim$mu.fv, saved$reference_model$mu.fv,
        tolerance = 1e-8)),
      isTRUE(all.equal(
        vapply(1:2, function(j) gamlss::getSmo(frich_sem_parsim, "mu",
          which = j)$sigb,
        numeric(1)),
        unname(c(saved$generating_parameters$plot_SD,
          saved$generating_parameters$year_SD)), tolerance = 1e-8
      ))
    )
    TRUE
  }, error = function(e) FALSE)
  if (!matches) stop(
    paste0(
      "Saved parsimonious focal bootstrap does not match the ",
      "current model, data, settings or implementation. "
    ),
    "Reference preserved at: ", focal_parsim_bootstrap_file,
    paste0(
      ". Inspect the mismatch; choose a new reference ",
      "directory for an intentional new run."
    )
  )
  message("Reusing parsimonious focal bootstrap: ", focal_parsim_bootstrap_file)
  saved
} else {
  result <- do.call(bootstrap_focal_dbi, c(
    list(model = frich_sem_parsim, data = focal_data), focal_bootstrap_settings
  ))
  result$cache_signature <- focal_bootstrap_signature
  export_focal_dbi_bootstrap(
    result, focal_parsim_bootstrap_dir,
    source_paths = c(
      analysis_data_path,
      "Network_SEM.R",
      "analysis_functions.R"
    )
  )
  result
}
if (!is.null(frich_parsim_bootstrap)) print(frich_parsim_bootstrap$summary)

# ------------------------------------------------------------------------------
# 10. Parsimonious network SEM ####
# ------------------------------------------------------------------------------

# Reuse the maximal node scales, ranges and adapters; intervals are conditional
# on the selected models. Standardized effects remain point estimates.
network_parsim_models <- list(
  prt_sem = prt_sem_parsim, niche_sem = niche_sem_parsim,
  frich_sem = frich_sem_parsim, pr_sem = pr_sem_parsim
)
focal_parsim_adapter <- sem_focal_dbi_inference_adapter(
  frich_sem_parsim, focal_data, frich_parsim_bootstrap
)
network_parsim_component_results <-
  stats::setNames(lapply(names(network_parsim_models), function(id) {
    model <- network_parsim_models[[id]]
    focal <- id == "frich_sem"
    result <- build_sem_component_results(
      model = model, component = network_components[[id]], model_id = id,
      data_id = network_data_ids[[id]], data_sets = node_reference_datasets,
      node_reference = node_reference_table,
      reference_values = reference_values,
      include_relevant_range = relevant_ranges_ready,
      coefficient_table =
        if (focal) focal_parsim_adapter$coefficients else NULL,
      inference_function = if (focal) focal_parsim_adapter$inference else
        sem_term_wald_test,
      predict_function = if (focal) sem_predict_focal_dbi else NULL,
      evaluation_function =
        sem_network_evaluation(all.vars(stats::formula(model))[[1L]])
    )
    terms <- result$term_results
    terms$response_estimand <- if (focal) {
      "expected focal count from exact DBI probabilities"
    } else if (id == "niche_sem") {
      "overall expected overlap including exact zeros"
    } else if (id == "prt_sem") {
      "expected parasitism proportion"
    } else "expected parasitoid count"
    terms$oe_estimand <- if (focal) {
      "DBI location-logit OE extension using expected counts"
    } else if (id == "niche_sem") {
        paste0(
          "positive-overlap logit OE with overall response ",
          "correlation; constant zero probability"
        )
      } else "Menard observed-empirical link-scale coefficient"
    terms$effect_method <- ifelse(
      terms$term_role %in% c("age_moderation", "ecological_interaction"),
      "not calculated: interaction",
      "average numerical partial derivative; feasible one-sided limits"
    )
    if (focal) terms$effect_method[terms$predictor == "par_rich"] <-
      paste0(
        "integer-support secant per species; +/-1 internally, ",
        "one-sided at observed endpoints; bound updated"
      )
    terms$uncertainty_status <- if (focal) paste(
      focal_parsim_adapter$interval_status, ";",
      frich_parsim_bootstrap$summary$accepted,
      paste0(
        "accepted bootstrap draws; P values unavailable; ",
        "standardized effects are point estimates"
      )
    ) else paste0(
      "95% link-scale Wald intervals; standardized effects are ",
      "point estimates"
    )
    terms$uncertainty_status <- paste(terms$uncertainty_status,
      "; conditional on selected specification")
    result$term_results <- terms
    result
  }), names(network_parsim_models))

network_parsim_terms <-
  dplyr::bind_rows(lapply(network_parsim_component_results, `[[`,
    "term_results"))
network_parsim_paths <- network_sem_path_table(network_parsim_component_results)
network_parsim_edges <-
  dplyr::bind_rows(lapply(network_parsim_component_results, `[[`, "edge_spec"))
network_parsim_effects <-
  dplyr::bind_rows(lapply(network_parsim_component_results, `[[`,
    "standardized_effects"))
network_parsim_r2 <- dplyr::bind_rows(lapply(names(network_parsim_models),
  function(id) {
    sem_network_component_r2(network_parsim_models[[id]],
      node_reference_datasets[[network_data_ids[[id]]]]) |>
      dplyr::mutate(component = network_components[[id]], model_id = id,
        data_id = network_data_ids[[id]], .before = 1)
  }))
stopifnot(
  !anyDuplicated(network_parsim_edges$edge_id),
  identical(network_parsim_paths$edge_id, network_parsim_edges$edge_id),
  all(network_parsim_terms$native_n == nrow(pd_net)),
  all(is.finite(network_parsim_effects$standardized_response_effect)),
  all(is.finite(network_parsim_effects$standardized_oe_effect)),
  all(is.finite(network_parsim_effects$standardized_relevant_range_effect)),
  all(is.finite(network_parsim_r2$r2_marginal)),
  all(is.finite(network_parsim_r2$r2_conditional))
)
network_parsim_reduction_log <- dplyr::bind_rows(
  parasitism_reduction$removal_log, niche_reduction$removal_log,
  pr_reduction$removal_log,
  dplyr::mutate(
    frich_reduction_plan, component = "focal_richness",
    converged = c(frich_no_hr_age$converged, frich_no_par_age$converged),
    parent_AIC = c(frich_comparison$AIC[1], frich_comparison_2$AIC[1]),
    candidate_AIC = c(frich_comparison$AIC[2], frich_comparison_2$AIC[2]),
    delta_AIC = candidate_AIC - parent_AIC,
    max_abs_change_in_parent_se = c(
      max(abs(frich_coefficient_stability$change_parent_SE)),
      max(abs(frich_coefficient_stability_2$change_parent_SE))
    )
  )
)

# The secondary parsimonious SEM differs only in its terminal response. Reuse
# the three selected upstream components and their selected-model focal
# bootstrap, and calculate only the realized-parasitism terminal results.
realized_parsim_terminal <- build_sem_component_results(
  model = pa_sem_parsim,
  component = "realized_parasitism",
  model_id = "pa_sem",
  data_id = "pd_net",
  data_sets = node_reference_datasets,
  node_reference = node_reference_table,
  reference_values = reference_values,
  include_relevant_range = relevant_ranges_ready,
  # The response name controls feasible predictor perturbations only; the
  # realized terminal uses the same network-richness constraints as parasitism.
  evaluation_function = sem_network_evaluation("cell_prate")
)
realized_parsim_terminal$term_results$response_estimand <-
  "expected parasitoid abundance (realized parasitism)"
realized_parsim_terminal$term_results$oe_estimand <-
  "Menard observed-empirical link-scale coefficient"
realized_parsim_terminal$term_results$effect_method <- ifelse(
  realized_parsim_terminal$term_results$term_role %in%
    c("age_moderation", "ecological_interaction"),
  "not calculated: interaction",
  "average numerical partial derivative; feasible one-sided limits"
)
realized_parsim_terminal$term_results$uncertainty_status <- paste0(
  "95% link-scale Wald intervals; standardized effects are point estimates; ",
  "conditional on selected specification"
)

realized_network_parsim_component_results <- c(
  list(pa_sem = realized_parsim_terminal),
  network_parsim_component_results[c("niche_sem", "frich_sem", "pr_sem")]
)
realized_network_parsim_models <- c(
  realized_parsim_terminal$model,
  network_parsim_models[c("niche_sem", "frich_sem", "pr_sem")]
)
realized_network_parsim_terms <- dplyr::bind_rows(lapply(
  realized_network_parsim_component_results, `[[`, "term_results"
))
realized_network_parsim_paths <- network_sem_path_table(
  realized_network_parsim_component_results
)
realized_network_parsim_edges <- dplyr::bind_rows(lapply(
  realized_network_parsim_component_results, `[[`, "edge_spec"
))
realized_network_parsim_effects <- dplyr::bind_rows(lapply(
  realized_network_parsim_component_results, `[[`, "standardized_effects"
))
realized_network_parsim_metadata <- dplyr::bind_rows(lapply(
  realized_network_parsim_component_results, `[[`, "component_metadata"
))
realized_network_parsim_standardization <- dplyr::bind_rows(lapply(
  realized_network_parsim_component_results, `[[`, "standardization_spec"
))
realized_network_parsim_r2 <- dplyr::bind_rows(
  sem_network_component_r2(pa_sem_parsim, pd_net) |>
    dplyr::mutate(
      component = "realized_parasitism",
      model_id = "pa_sem",
      data_id = "pd_net",
      .before = 1
    ),
  network_parsim_r2 |>
    dplyr::filter(component != "parasitism")
)
realized_network_parsim_reduction_log <- dplyr::bind_rows(
  real_parasitism_reduction$removal_log,
  network_parsim_reduction_log |>
    dplyr::filter(component != "parasitism")
)
stopifnot(
  identical(
    names(realized_network_parsim_models),
    realized_network_parsim_metadata$source_model
  ),
  !"parasitism" %in% realized_network_parsim_terms$component,
  "realized_parasitism" %in% realized_network_parsim_terms$component,
  !anyDuplicated(realized_network_parsim_edges$edge_id),
  identical(
    realized_network_parsim_paths$edge_id,
    realized_network_parsim_edges$edge_id
  ),
  setequal(
    realized_network_parsim_standardization$edge_id,
    realized_network_parsim_effects$edge_id
  ),
  all(realized_network_parsim_terms$native_n == nrow(pd_net)),
  all(is.finite(
    realized_network_parsim_effects$standardized_response_effect
  )),
  all(is.finite(realized_network_parsim_effects$standardized_oe_effect)),
  all(is.finite(
    realized_network_parsim_effects$standardized_relevant_range_effect
  )),
  identical(
    realized_network_parsim_r2$model_id,
    names(realized_network_parsim_models)
  ),
  all(is.finite(realized_network_parsim_r2$r2_marginal)),
  all(is.finite(realized_network_parsim_r2$r2_conditional))
)

# ------------------------------------------------------------------------------
# 11. D-separation assessment and parsimonious exports ####
# ------------------------------------------------------------------------------

# Include forest-age adjustment paths in the graph, alongside ecological paths.
# Exogenous predictors may covary freely; interactions are not separate claims.
network_maximal_graph_edges <- dplyr::bind_rows(lapply(names(network_models),
  function(id) {
    model <- network_models[[id]]
    predictors <- if (inherits(model, "gamlss")) {
      names(stats::coef(model, what = "mu"))
    } else names(glmmTMB::fixef(model)$cond)
    predictors <- intersect(predictors, node_reference_table$model_variable)
    tibble::tibble(
      source_node = node_reference_table$node[match(predictors,
        node_reference_table$model_variable)],
      destination_node = network_components[[id]],
      edge_id = paste(source_node, destination_node, sep = "__"),
      model_id = id, data_id = network_data_ids[[id]], predictor = predictors,
      response = all.vars(stats::formula(model))[[1L]]
    )
  }))
network_expected_removed_edges <- c("tree_richness__parasitism",
  "forest_age__parasitism")
network_parsim_graph <- assemble_reduced_sem_graph(
  models = network_parsim_models,
  maximal_edge_spec = network_maximal_graph_edges,
  expected_removed_edges = network_expected_removed_edges
)

# Both missing paths condition on the five retained parents of parasitism.
# These terms are already in the final model, so no extra conditioning
# is needed.
network_basis_conditioning <- tibble::tibble(
  edge_id = network_expected_removed_edges,
  conditioning_nodes = rep(list(c("parasitoid_richness", "focal_richness",
    "niche_overlap", "host_richness",
    "host_abundance")), 2L),
  conditioning_terms = rep(list(character()), 2L)
)
network_parsim_basis_set <- derive_sem_basis_set(
  maximal_edge_spec = network_maximal_graph_edges,
  removed_edges = network_parsim_graph$removed_edges,
  conditioning_spec = network_basis_conditioning, models = network_parsim_models
)
# Refit each one-df add-back against the final model, retaining the
# beta-binomial response, trial weights, dispersion and random structure
# on the same 448 rows.
network_final_basis_audit <- audit_sem_basis_set(
  basis_set = network_parsim_basis_set, models = network_parsim_models,
  data_sets = node_reference_datasets, model_data_ids = network_data_ids,
  component_metadata = dplyr::bind_rows(lapply(network_parsim_component_results,
    `[[`, "component_metadata")),
  reference_values = reference_values, control = strict_control
)
network_parsim_missing_paths <- network_final_basis_audit$results
network_parsim_global_fit <- sem_fisher_c(network_parsim_missing_paths) |>
  dplyr::mutate(
    n = nrow(pd_net),
    scope = paste0(
      "Exploratory post-selection Fisher's C from ",
      "native-response add-back LRTs; exogenous associations ",
      "unrestricted"
    )
  )

# Build the corresponding maximal edge registry for the secondary SEM.
realized_network_maximal_models <- list(
  pa_sem = pa_sem,
  niche_sem = niche_sem,
  frich_sem = frich_sem,
  pr_sem = pr_sem
)
realized_network_components <- c(
  pa_sem = "realized_parasitism",
  niche_sem = "niche_overlap",
  frich_sem = "focal_richness",
  pr_sem = "parasitoid_richness"
)
realized_network_data_ids <- c(
  pa_sem = "pd_net",
  niche_sem = "pd_net",
  frich_sem = "focal_data",
  pr_sem = "pd_net"
)
realized_network_maximal_graph_edges <- dplyr::bind_rows(lapply(
  names(realized_network_maximal_models),
  function(id) {
    model <- realized_network_maximal_models[[id]]
    predictors <- if (inherits(model, "gamlss")) {
      names(stats::coef(model, what = "mu"))
    } else names(glmmTMB::fixef(model)$cond)
    predictors <- intersect(predictors, node_reference_table$model_variable)
    tibble::tibble(
      source_node = node_reference_table$node[match(
        predictors, node_reference_table$model_variable
      )],
      destination_node = realized_network_components[[id]],
      edge_id = paste(source_node, destination_node, sep = "__"),
      model_id = id,
      data_id = realized_network_data_ids[[id]],
      predictor = predictors,
      response = all.vars(stats::formula(model))[[1L]]
    )
  }
))
realized_network_expected_removed_edges <- c(
  "tree_richness__realized_parasitism",
  "forest_age__realized_parasitism"
)
realized_network_parsim_graph <- assemble_reduced_sem_graph(
  models = realized_network_parsim_models,
  maximal_edge_spec = realized_network_maximal_graph_edges,
  expected_removed_edges = realized_network_expected_removed_edges
)

# Host richness and brood-cell abundance remain a priori controls for network
# breadth and host availability. Both omitted terminal paths condition on all
# five retained ecological parents of realized parasitism.
realized_network_basis_conditioning <- tibble::tibble(
  edge_id = realized_network_expected_removed_edges,
  conditioning_nodes = rep(list(c(
    "parasitoid_richness", "focal_richness", "niche_overlap",
    "host_richness", "host_abundance"
  )), 2L),
  conditioning_terms = rep(list(character()), 2L)
)
realized_network_parsim_basis_set <- derive_sem_basis_set(
  maximal_edge_spec = realized_network_maximal_graph_edges,
  removed_edges = realized_network_parsim_graph$removed_edges,
  conditioning_spec = realized_network_basis_conditioning,
  models = realized_network_parsim_models
)

# Refit each one-df add-back against the selected nbinom terminal model on the
# same 448 network-eligible plot-years.
realized_network_final_basis_audit <- audit_sem_basis_set(
  basis_set = realized_network_parsim_basis_set,
  models = realized_network_parsim_models,
  data_sets = node_reference_datasets,
  model_data_ids = realized_network_data_ids,
  component_metadata = realized_network_parsim_metadata,
  reference_values = reference_values,
  control = strict_control
)
realized_network_parsim_missing_paths <-
  realized_network_final_basis_audit$results
realized_network_parsim_global_fit <- sem_fisher_c(
  realized_network_parsim_missing_paths
) |>
  dplyr::mutate(
    n = nrow(pd_net),
    scope = paste0(
      "Exploratory post-selection Fisher's C from ",
      "native-response add-back LRTs; exogenous associations ",
      "unrestricted"
    )
  )

network_parsim_source_paths <- c(
  analysis_data_path,
  focal_parsim_bootstrap_file,
  "Network_SEM.R",
  "analysis_functions.R"
)
network_parsim_handoff <- list(
  sem_variant = "interaction_network_parsimonious",
  primary_sem_models = network_parsim_models,
  node_reference_datasets = node_reference_datasets,
  node_reference_table = node_reference_table,
  reference_values = reference_values,
  relevant_ranges_ready = relevant_ranges_ready,
  primary_component_metadata =
    dplyr::bind_rows(lapply(network_parsim_component_results, `[[`,
      "component_metadata")),
  primary_sem_term_results = network_parsim_terms,
  primary_edge_spec = network_parsim_edges,
  primary_standardization_spec =
    dplyr::bind_rows(lapply(network_parsim_component_results, `[[`,
      "standardization_spec")),
  primary_standardized_effects = network_parsim_effects,
  path_table = network_parsim_paths, component_r2_table = network_parsim_r2,
  focal_bootstrap = frich_parsim_bootstrap,
  reduction_log = network_parsim_reduction_log,
  graph = network_parsim_graph, basis_set = network_parsim_basis_set,
  missing_path_audits = network_final_basis_audit$audits,
  missing_path_results = network_parsim_missing_paths,
  global_fit = network_parsim_global_fit,
  interpretation = c(
    paste0(
      "Exploratory network SEM after interaction and ",
      "main-effect reduction; intervals exclude ",
      "model-selection uncertainty."
    ),
    paste0(
      "Fisher's C assesses the two omitted paths after ",
      "selection; exogenous associations are unrestricted."
    ),
    paste0(
      "Focal effects at fixed total richness describe ",
      "replacement of nonfocal richness."
    ),
    paste0(
      "Total-to-focal count contrasts update min(4, total ",
      "richness); the link slope describes location only."
    ),
    paste0(
      "Relevant-range contrasts use rows supporting both ",
      "endpoints; rr_n_excluded records restrictions."
    ),
    "Predictions set random effects to zero without integrating over them.",
    paste0(
      "Link/location intervals do not quantify uncertainty in ",
      "standardized effects or nonlinear pathway products."
    ),
    paste0(
      "Focal bootstrap intervals use accepted refits; failure ",
      "and warning audits are retained in focal_bootstrap."
    )
  ),
  source_checksums = data.frame(
    path = normalizePath(network_parsim_source_paths, mustWork = TRUE),
    md5 = unname(tools::md5sum(network_parsim_source_paths))
  ),
  session_info = utils::sessionInfo()
)

realized_network_parsim_handoff <- list(
  sem_variant = "interaction_network_realized_parasitism_parsimonious",
  primary_sem_models = realized_network_parsim_models,
  node_reference_datasets = node_reference_datasets,
  node_reference_table = node_reference_table,
  reference_values = reference_values,
  relevant_ranges_ready = relevant_ranges_ready,
  primary_component_metadata = realized_network_parsim_metadata,
  primary_sem_term_results = realized_network_parsim_terms,
  primary_edge_spec = realized_network_parsim_edges,
  primary_standardization_spec = realized_network_parsim_standardization,
  primary_standardized_effects = realized_network_parsim_effects,
  path_table = realized_network_parsim_paths,
  component_r2_table = realized_network_parsim_r2,
  focal_bootstrap = frich_parsim_bootstrap,
  reduction_log = realized_network_parsim_reduction_log,
  graph = realized_network_parsim_graph,
  basis_set = realized_network_parsim_basis_set,
  missing_path_audits = realized_network_final_basis_audit$audits,
  missing_path_results = realized_network_parsim_missing_paths,
  global_fit = realized_network_parsim_global_fit,
  interpretation = c(
    paste0(
      "Secondary selected network SEM differing from the primary variant ",
      "only in its terminal realized-parasitism response."
    ),
    paste0(
      "Realized parasitism is parasitoid abundance; host richness and ",
      "brood-cell abundance remain a priori controls for network breadth ",
      "and host availability."
    ),
    paste0(
      "Fisher's C assesses the two omitted terminal paths after selection; ",
      "exogenous associations are unrestricted."
    ),
    paste0(
      "The selected niche-overlap, focal-richness and total-richness ",
      "components, including focal-bootstrap inference, are reused unchanged ",
      "from the primary parsimonious network SEM."
    ),
    paste0(
      "Focal effects at fixed total richness describe replacement of ",
      "nonfocal richness."
    ),
    paste0(
      "Relevant-range contrasts use rows supporting both endpoints; ",
      "rr_n_excluded records restrictions."
    ),
    "Predictions set random effects to zero without integrating over them.",
    paste0(
      "Link/location intervals do not quantify uncertainty in standardized ",
      "effects, nonlinear pathway products or model selection."
    ),
    paste0(
      "Focal bootstrap intervals use accepted refits; failure and warning ",
      "audits are retained in focal_bootstrap."
    )
  ),
  source_checksums = data.frame(
    path = normalizePath(network_parsim_source_paths, mustWork = TRUE),
    md5 = unname(tools::md5sum(network_parsim_source_paths))
  ),
  session_info = utils::sessionInfo()
)

# Maximal results and both bootstrap checkpoints remain separate from the
# selected-model exports.
network_parsim_output_dir <- file.path(network_output_dir, "parsimonious")
export_network_sem(
  network_parsim_handoff,
  network_parsim_output_dir,
  file_stem = "parsimonious_network_sem"
)
utils::write.csv(network_parsim_reduction_log,
  file.path(network_parsim_output_dir,
    "parsimonious_network_sem_reduction_log.csv"), row.names = FALSE)

realized_network_parsim_output_dir <- file.path(
  network_output_dir, "parsimonious_realized"
)
export_network_sem(
  realized_network_parsim_handoff,
  realized_network_parsim_output_dir,
  file_stem = "parsimonious_realized_network_sem"
)
utils::write.csv(
  realized_network_parsim_reduction_log,
  file.path(
    realized_network_parsim_output_dir,
    "parsimonious_realized_network_sem_reduction_log.csv"
  ),
  row.names = FALSE
)

print(network_parsim_missing_paths, n = Inf, width = Inf)
print(network_parsim_global_fit, width = Inf)
print(realized_network_parsim_missing_paths, n = Inf, width = Inf)
print(realized_network_parsim_global_fit, width = Inf)

# ------------------------------------------------------------------------------
# 12. Full-dataset focal-richness sensitivity ####
# ------------------------------------------------------------------------------

# Compare both maximal parasitism responses with focal-richness adjustments on
# all 786 eligible plot-years, retaining capped richness and existing scaling.
# Attenuation of sc_fa:sc_pr10 is descriptive, not a mediated proportion.
# Coefficient tables, nested comparisons, and fitted-model RDS files are saved
# for both terminal responses.

analysis_data <- readRDS(analysis_data_path)
strict_control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000L, eval.max = 10000L)
)

# Full dataset for parasitism models
pd_parasitism <- analysis_data$pd_parasitism |>
  arrange(plot_id, year) |>
  droplevels()

prt_full <- glmmTMB(
  cell_prate ~
    sc_logtr + sc_fa + sc_rtfd + sc_rtsv +
      sc_rtcells + sc_hr + sc_pr10 +
      sc_logtr:sc_fa +
      sc_fa:sc_rtsv +
      sc_fa:sc_rtcells +
      sc_fa:sc_hr +
      sc_fa:sc_pr10 +
      sc_rtsv:sc_pr10 +
      sc_rtcells:sc_hr +
      (1 | site / plot_id) + (1 | year),
  data = pd_parasitism,
  weights = total_cells,
  family = betabinomial(link = "logit"),
  dispformula = ~ sc_cells,
  ziformula = ~ 0,
  REML = FALSE,
  na.action = na.fail,
  control = strict_control
)

prt_full_focal <- update(
  prt_full,
  . ~ . + sc_frich
)

prt_full_focal_age <- update(
  prt_full,
  . ~ . + sc_frich + sc_fa:sc_frich
)

print(summary(prt_full))
print(summary(prt_full_focal))
print(summary(prt_full_focal_age))
prt_full_comparison <- anova(prt_full, prt_full_focal, prt_full_focal_age)
print(prt_full_comparison)

# Prepare the full-data sensitivity records before writing any files.
prt_full_models <- list(
  prt_full = prt_full,
  prt_full_focal = prt_full_focal,
  prt_full_focal_age = prt_full_focal_age
)
prt_full_coefficients <- dplyr::bind_rows(lapply(
  names(prt_full_models), function(model_id) {
    model <- prt_full_models[[model_id]]
    coefficient_components <- summary(model)$coefficients
    dplyr::bind_rows(lapply(names(coefficient_components), function(component) {
      tab <- coefficient_components[[component]]
      if (is.null(tab) || nrow(tab) == 0L) return(NULL)
      data.frame(
        model_id = model_id,
        response = "parasitism",
        component = component,
        term = rownames(tab),
        estimate = tab[, "Estimate"],
        std_error = tab[, "Std. Error"],
        conf_low = tab[, "Estimate"] - qnorm(0.975) * tab[, "Std. Error"],
        conf_high = tab[, "Estimate"] + qnorm(0.975) * tab[, "Std. Error"],
        z_value = tab[, "z value"],
        p_value = tab[, "Pr(>|z|)"],
        n_obs = nobs(model),
        family = family(model)$family,
        coefficient_scale = if (component == "cond") "logit" else "log",
        interval_method = "95% Wald",
        p_adjustment = "none",
        row.names = NULL
      )
    }))
  }
))

# ANOVA rows retain the baseline and successive nested comparisons, including
# log likelihood, AIC, likelihood-ratio statistics, df and unadjusted P values.
prt_full_comparison_table <- as.data.frame(prt_full_comparison)
prt_full_comparison_table$model_id <- rownames(prt_full_comparison_table)
stopifnot(identical(
  prt_full_comparison_table$model_id, names(prt_full_models)
))
prt_full_comparison_table$compared_with <- c(
  NA_character_, "prt_full", "prt_full_focal"
)
prt_full_comparison_table$n_obs <- vapply(prt_full_models, nobs, numeric(1))
prt_full_comparison_table$p_adjustment <- "none"
rownames(prt_full_comparison_table) <- NULL

# Realized parasitism: use the saved maximal model, not
# the 448-row network model (which has a different negative-binomial family).
# Refit its unchanged specification on the original ordered data so every
# nested fit uses the same data expression as well as identical observations.
pa_reference_path <- file.path(
  model_output_dir,
  "maximal_models",
  "maximal_models.rds"
)
if (!file.exists(pa_reference_path)) {
  stop("Missing maximal_models.rds. Run Maximal_models.R first.")
}
pa_reference <- readRDS(pa_reference_path)$maximal_models$realized_parasitism
if (!inherits(pa_reference, "glmmTMB") || nobs(pa_reference) != 786L ||
    family(pa_reference)$family != "nbinom1" ||
    family(pa_reference)$link != "log") {
  stop("Expected the saved 786-row negative-binomial-1 realized-parasitism model.")
}
pa_full_data <- as.data.frame(analysis_data$pd_parasitism)
pa_reference_frame <- model.frame(pa_reference)
if (!all(c(names(pa_reference_frame), "sc_frich") %in% names(pa_full_data)) ||
    !isTRUE(all.equal(as.data.frame(pa_reference_frame),
                      pa_full_data[names(pa_reference_frame)],
                      check.attributes = FALSE)) ||
    any(!is.finite(pa_full_data$sc_frich))) {
  stop("The saved realized-parasitism model and prepared data do not match.")
}
pa_full <- update(
  pa_reference, data = pa_full_data, na.action = na.fail, control = strict_control
)
if (!isTRUE(all.equal(unlist(fixef(pa_full)), unlist(fixef(pa_reference)),
                      tolerance = 1e-5)) ||
    abs(as.numeric(logLik(pa_full) - logLik(pa_reference))) > 1e-5) {
  stop("The realized-parasitism baseline refit differs from its saved reference.")
}
pa_full_focal <- update(pa_full, . ~ . + sc_frich)
pa_full_focal_age <- update(pa_full, . ~ . + sc_frich + sc_fa:sc_frich)
pa_full_models <- list(
  pa_full = pa_full, pa_full_focal = pa_full_focal,
  pa_full_focal_age = pa_full_focal_age
)

# Do not export failed fits for either response.
for (model in c(prt_full_models, pa_full_models)) {
  if (nobs(model) != 786L || model$fit$convergence != 0L ||
      !isTRUE(model$sdr$pdHess) || !is.finite(as.numeric(logLik(model))) ||
      any(!is.finite(unlist(fixef(model)))) ||
      any(!is.finite(model$sdr$cov.fixed))) {
    stop("A full-data focal-sensitivity fit failed numerical or sample checks.")
  }
}
print(summary(pa_full))
print(summary(pa_full_focal))
print(summary(pa_full_focal_age))
pa_full_comparison <- anova(pa_full, pa_full_focal, pa_full_focal_age)
print(pa_full_comparison)
pa_full_coefficients <- dplyr::bind_rows(lapply(
  names(pa_full_models), function(model_id) {
    model <- pa_full_models[[model_id]]
    coefficient_components <- summary(model)$coefficients
    dplyr::bind_rows(lapply(names(coefficient_components), function(component) {
      tab <- coefficient_components[[component]]
      if (is.null(tab) || nrow(tab) == 0L) return(NULL)
      data.frame(
        model_id = model_id,
        response = "realized_parasitism",
        component = component,
        term = rownames(tab),
        estimate = tab[, "Estimate"],
        std_error = tab[, "Std. Error"],
        conf_low = tab[, "Estimate"] - qnorm(0.975) * tab[, "Std. Error"],
        conf_high = tab[, "Estimate"] + qnorm(0.975) * tab[, "Std. Error"],
        z_value = tab[, "z value"],
        p_value = tab[, "Pr(>|z|)"],
        n_obs = nobs(model),
        family = family(model)$family,
        coefficient_scale = "log",
        interval_method = "95% Wald",
        p_adjustment = "none",
        row.names = NULL
      )
    }))
  }
))
pa_full_comparison_table <- as.data.frame(pa_full_comparison)
pa_full_comparison_table$model_id <- rownames(pa_full_comparison_table)
stopifnot(identical(pa_full_comparison_table$model_id, names(pa_full_models)))
pa_full_comparison_table$compared_with <- c(
  NA_character_, "pa_full", "pa_full_focal"
)
pa_full_comparison_table$n_obs <- vapply(pa_full_models, nobs, numeric(1))
pa_full_comparison_table$p_adjustment <- "none"
rownames(pa_full_comparison_table) <- NULL

# Export the primary sensitivity and place the parallel realized-parasitism
# handoff in its own subdirectory.
prt_full_output_dir <- file.path(
  network_output_dir, "full_data_focal_sensitivity"
)
dir.create(prt_full_output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(prt_full_output_dir)) {
  stop("Could not create the full-data focal-sensitivity output directory.")
}
utils::write.csv(
  prt_full_coefficients,
  file.path(prt_full_output_dir, "coefficients.csv"), row.names = FALSE
)
utils::write.csv(
  prt_full_comparison_table,
  file.path(prt_full_output_dir, "model_comparisons.csv"), row.names = FALSE
)
saveRDS(
  list(
    models = prt_full_models,
    coefficients = prt_full_coefficients,
    comparisons = prt_full_comparison_table,
    interpretation = paste(
      "Full-data focal-richness sensitivity retaining capped richness and",
      "existing scaling. Attenuation of sc_fa:sc_pr10 is descriptive,",
      "not a mediated proportion."
    ),
    session_info = utils::sessionInfo()
  ),
  file.path(prt_full_output_dir, "full_data_focal_sensitivity.rds")
)

pa_full_output_dir <- file.path(prt_full_output_dir, "realized")
dir.create(pa_full_output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(pa_full_output_dir)) {
  stop("Could not create the realized-parasitism focal-sensitivity directory.")
}
utils::write.csv(
  pa_full_coefficients,
  file.path(pa_full_output_dir, "coefficients.csv"), row.names = FALSE
)
utils::write.csv(
  pa_full_comparison_table,
  file.path(pa_full_output_dir, "model_comparisons.csv"), row.names = FALSE
)
saveRDS(
  list(
    models = pa_full_models,
    coefficients = pa_full_coefficients,
    comparisons = pa_full_comparison_table,
    interpretation = paste(
      "Full-data realized-parasitism focal-richness sensitivity retaining",
      "the maximal negative-binomial-1 specification, capped richness and",
      "original scaling. Coefficient attenuation is descriptive, not mediation."
    ),
    baseline_source = pa_reference_path,
    baseline_source_md5 = unname(tools::md5sum(pa_reference_path)),
    session_info = utils::sessionInfo()
  ),
  file.path(pa_full_output_dir, "full_data_focal_sensitivity.rds")
)

message("Primary and realized-parasitism network SEMs ran successfully ✔️")
