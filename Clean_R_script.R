# ------------------------------------------------------------------------------
## Header ##
# Title: Forest growth strengthens diversity effects on multitrophic interactions
# Author: Massimo Martini
# Date: 06th May 2026
# R version: 4.5.1
# Description: Main analysis script. Reproduces all statistical models.

# ------------------------------------------------------------------------------

## Comments ##
# For simplicity, we leave out all diagnostic and plotting procedures

# Recommended workflow:
# 1. Open the .Rproj file
# 2. Recommended - restore package environment with:
  #  renv::restore()
# 3. Run analysis script

# clear workspace (only if not in fresh R session)
# rm(list = ls())
# graphics.off()
# cat("\014")

# ------------------------------------------------------------------------------
# 1. Load libraries ####
# ------------------------------------------------------------------------------

#load
library(dplyr)
library(readr)
library(lme4)
library(glmmTMB)
library(piecewiseSEM)
library(DHARMa)
library(performance)
library(MuMIn)
library(broom.mixed)
library(officer)
library(flextable)
library(tibble)
library(tidyr)

# ------------------------------------------------------------------------------
# 2. Read and prepare data ####
# ------------------------------------------------------------------------------

plot_data <- readr::read_csv("plot_data.csv")

source("final_functions.R") #upload helper functions

plot_data$tree_r[plot_data$tree_r>16] = 16 # we grouped the sparsely-replicated 24-tree-species-
# richness plots (n = 4) with the 16-species level to avoid unreliable model estimates at high tree
# richness levels and to have equal sample size among species mixture levels

plot_data$log_tr = log2(plot_data$tree_r) #log-transformation

plot_data$par_rich10 <- plot_data$par_rich
plot_data$par_rich10[plot_data$par_rich > 10] <- 10 # similarly as for tree richness, we binned
# sparsely-replicated, high-richness parasitoid communities (with more than 10 species, n = 4) with
# 10-parasitoid-species communities. Note that par_rich10 is only used as a predictor for parasitism
# models, to avoid overestimation of model predictions.

plot_data <- plot_data %>%
  arrange(plot_id, year) # arrange dataset by plot and year

plot_data$for_age <- ifelse(plot_data$site == "A",
                            plot_data$year - 2009,
                            plot_data$year - 2010) # calculate stand age based on year of site planting

plot_data$year = as.character(plot_data$year) # year as categorical factor

plot_data <- plot_data %>%
  mutate(cell_prate = (par_abund / total_cells)) #calculate parasitism rate as the fraction of 
# attacked vs total brood cells

plot_data$cell_prate[is.nan(plot_data$cell_prate)] <- 0 #change NAs to 0

plot_data$p_failure <- (plot_data$total_cells - plot_data$par_abund) #calculate parasitism failure

pd_climate <- plot_data %>%
  filter(!is.na(annual_temperature)) # subset dataframe with plot-years with climate data (n = 448)


#scale all predictors
#choose variables to scale
cols_to_scale <- c("log_tr", "stand_volume", "tree_fd", "for_age", "elevation", "eastness",
                   "northness", "slope","annual_temperature", "annual_humidity", "total_cells", 
                   "host_rich", "par_rich", "par_rich10","network_size", "n_links", 
                   "linkage_density", "h2", "mean_dprime_hl", "niche_overlap_hl", 
                   "interaction_evenness", "robustness_hl", "focal_rich")

#set new abbreviated scaled column names
name_map <- c(
  log_tr               = "sc_logtr",
  stand_volume         = "sc_sv",
  tree_fd              = "sc_fd",
  for_age              = "sc_fa",
  elevation            = "sc_elev",
  eastness             = "sc_east",
  northness            = "sc_north",
  slope                = "sc_slope",
  annual_temperature   = "sc_temp",
  annual_humidity      = "sc_humid",
  total_cells          = "sc_cells",
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


#scale all selected variables in both dataframes
for (old in names(name_map)) {
  new <- name_map[[old]]
  plot_data[[new]] <- as.numeric(scale(plot_data[[old]]))
}

for (old in names(name_map)) {
  new <- name_map[[old]]
  pd_climate[[new]] <- as.numeric(scale(pd_climate[[old]]))
}

#check correlation between plot and community variables
vars <- plot_data[, c("log_tr", "stand_volume", "tree_fd", "for_age", "elevation", "eastness",
                      "northness", "slope","annual_temperature", "annual_humidity", "total_cells", 
                      "host_rich","par_abund", "par_rich", "focal_rich")]

cortable <- cor(vars, use = "complete.obs", method = "pearson")

cortable %>% 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  filter(var1 < var2, abs(cor) > 0.60) %>% 
  arrange(desc(abs(cor))) 

#check correlation between network variables
vars2 <- plot_data[, c("network_size", "n_links", "linkage_density", "h2", "mean_dprime_hl", 
                      "niche_overlap_hl", "interaction_evenness", "robustness_hl", "focal_rich")]

cortable2 <- cor(vars2, use = "complete.obs", method = "pearson")

cortable2 %>% 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  filter(var1 < var2, abs(cor) > 0.60) %>% 
  arrange(desc(abs(cor))) # many network metrics correlated with one another

# ------------------------------------------------------------------------------
# 3. GENERALIZED LINEAR MIXED MODELS (GLMMs) ####
# ------------------------------------------------------------------------------

## Hosts ##
# Host species richness

hr_tree <- glmmTMB(host_rich ~ sc_logtr + sc_fa +
                     sc_logtr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois()) # basic model

hr_vol <- glmmTMB(host_rich ~ sc_logtr + sc_sv + sc_fa +
                    sc_logtr : sc_fa +
                    sc_sv : sc_fa +
                    (1 | site/plot_id) + (1 | year), 
                  data = plot_data,
                  family = genpois()) # with stand biomass

hr_full <- glmmTMB(host_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells +
                        sc_slope + sc_elev + sc_east + sc_north +
                        sc_logtr : sc_fa +
                        sc_sv : sc_fa +
                        (1 | site/plot_id) + (1 | year), 
                      data = plot_data,
                      family = genpois()) # full model

hr_parsim <- glmmTMB(host_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells +
                     sc_slope + sc_elev + sc_east + sc_north +
                     sc_sv : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois()) # parsimonious (non-sig. interactions removed)

# Host abundance (total brood cells constructed)
# Models include a dispersion parameter for temporal autocorrelation

ha_tree <- glmmTMB(total_cells ~ sc_logtr + sc_fa +
                     sc_logtr : sc_fa +
                     (1 | site/plot_id) + (1 | year),
                   data = plot_data,
                   family = genpois(),
                   dispformula = ~ar1(year + 0 | site/plot_id)) # basic model

ha_vol <- glmmTMB(total_cells ~ sc_logtr + sc_sv + sc_fa +
                    sc_logtr : sc_fa +
                    sc_sv : sc_fa +
                    (1 | site/plot_id) + (1 | year), 
                  data = plot_data,
                  family = genpois(),
                  dispformula = ~ar1(year + 0 | site/plot_id)) # with stand biomass

ha_full <- glmmTMB(total_cells ~ sc_logtr + sc_fd + sc_sv + sc_fa +
                     sc_slope + sc_elev + sc_east + sc_north +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   dispformula = ~ar1(year + 0 | site/plot_id)) # full model

ha_parsim <- glmmTMB(total_cells ~ sc_logtr + sc_fd + sc_sv + sc_fa +
                     sc_slope + sc_elev + sc_east + sc_north +
                     sc_sv : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   dispformula = ~ar1(year + 0 | site/plot_id)) # parsimonious
                                                                # (non-sig. interactions removed)

## Parasitoids ##
# Parasitoid species richness 

pr_tree <- glmmTMB(par_rich ~ sc_logtr + sc_fa +
                     sc_logtr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1) # basic model

pr_vol <- glmmTMB(par_rich ~ sc_logtr + sc_sv + sc_fa +
                    sc_logtr : sc_fa +
                    sc_sv : sc_fa +
                    (1 | site/plot_id) + (1 | year), 
                  data = plot_data,
                  family = genpois(),
                  ziformula = ~1) # with stand biomass

pr_full <- glmmTMB(par_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     sc_cells : sc_fa +
                     sc_hr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois()) # full model

pr_parsim <- glmmTMB(par_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells +
                     sc_cells : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois()) # parsimonious (non sig. interactions removed)

# Realized parasitism (or parasitoid abundance)

pa_tree <- glmmTMB(par_abund ~ sc_logtr + sc_fa +
                     sc_logtr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1) # basic model

pa_vol <- glmmTMB(par_abund ~ sc_logtr + sc_sv + sc_fa +
                    sc_logtr : sc_fa +
                    sc_sv : sc_fa +
                    (1 | site/plot_id) + (1 | year), 
                  data = plot_data,
                  family = genpois(),
                  ziformula = ~1) # with stand biomass

pa_full <- glmmTMB(par_abund ~  sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells + sc_pr10 +
                     sc_logtr : sc_fa + 
                     sc_sv : sc_fa +
                     sc_cells : sc_fa +
                     sc_hr : sc_fa +
                     sc_pr10 : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1) # full

pa_parsim <- glmmTMB(par_abund ~  sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells + sc_pr10 +
                     sc_sv : sc_fa +
                     sc_pr10 : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1) # parsimonious (non-sig. interactions removed)

## Parasitism ##

prt_tree <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_fa +
                      sc_logtr : sc_fa +
                      (1 | site/plot_id) + (1 | year),
                    data = plot_data,
                    family = betabinomial(link = "logit")) # basic model

prt_vol <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_sv + sc_fa +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa + 
                     (1 | site/plot_id) + (1 | year),
                   data = plot_data,
                   family = betabinomial(link = "logit")) # with stand volume

prt_pr <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_sv + sc_fa + sc_hr + sc_pr10 +
                    sc_logtr : sc_fa +
                    sc_sv : sc_fa + 
                    sc_hr : sc_fa +
                    sc_pr10 : sc_fa +
                    (1 | site/plot_id) + (1 | year),
                  data = plot_data,
                  family = betabinomial(link = "logit")) # with host and parasitoid richness

prt_full <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                      sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr10 +
                      sc_cells : sc_fa +
                      sc_hr : sc_fa +
                      sc_pr10 : sc_fa +
                      sc_sv : sc_fa +
                      sc_pr10 : sc_sv +
                      sc_hr : sc_cells,
                    data = plot_data,
                    family = betabinomial(link = "logit")) # full model

prt_parsim <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                      sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr10 +
                      sc_pr10 : sc_fa +
                      sc_pr10 : sc_sv +
                      sc_hr : sc_cells,
                    data = plot_data,
                    family = betabinomial(link = "logit")) # parsimonious

#put main model objects into a single list "models" 
models <- list(
  hr_tree   = hr_tree,
  hr_vol    = hr_vol,
  hr_full   = hr_full,
  hr_parsim = hr_parsim,
  ha_tree   = ha_tree,
  ha_vol    = ha_vol,
  ha_full   = ha_full,
  ha_parsim = ha_parsim,
  pr_tree   = pr_tree,
  pr_vol    = pr_vol,
  pr_full   = pr_full, 
  pr_parsim = pr_parsim,
  pa_tree   = pa_tree,
  pa_vol    = pa_vol,
  pa_full   = pa_full,
  pa_parsim = pa_parsim,
  prt_tree  = prt_tree,
  prt_vol   = prt_vol,
  prt_pr    = prt_pr,
  prt_full  = prt_full,
  prt_parsim = prt_parsim
)

rm(list = names(models))

output_dir <- "model_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# create new directory for model summaries

model_summaries <- export_glmmTMB_markdown(
  models = models,
  file = file.path(output_dir, "model_summaries.md"),
  digits = 3
) # export results as a clean .md file

model_summaries_docx <- export_glmmTMB_docx(
  models = models,
  file = file.path(output_dir, "model_summaries.docx"),
  digits = 3
) # export results as a .docx file

# ------------------------------------------------------------------------------
# 4. SENSITIVITY ANALYSES ####
# ------------------------------------------------------------------------------
# Robustness checks supporting our model specifications and findings.
# Output from these models is reported as supplementary information.

# Host species richness

hr_fam <- glmmTMB(host_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells +
                     sc_slope + sc_elev + sc_east + sc_north +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = poisson()) # change family

hr_quad <- glmmTMB(host_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + I(sc_fa^2) + 
                          sc_cells +
                          sc_slope + sc_elev + sc_east + sc_north +
                          sc_logtr : sc_fa +
                          sc_sv : sc_fa +
                          (1 | site/plot_id) + (1 | year), 
                        data = plot_data,
                        family = genpois()) # check for nonlinear effects of time

# Host abundance

ha_fam <- glmmTMB(total_cells ~ sc_logtr + sc_fd + sc_sv + sc_fa +
                     sc_slope + sc_elev + sc_east + sc_north +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = nbinom2(),
                   dispformula = ~ar1(year + 0 | site/plot_id)) # change family

ha_quad <- glmmTMB(total_cells ~ sc_logtr + sc_fd + sc_sv +
                         sc_fa + I(sc_fa^2) +
                         sc_slope + sc_elev + sc_east + sc_north +
                         sc_logtr : sc_fa +
                         sc_sv : sc_fa +
                         (1 | site/plot_id) + (1 | year), 
                       data = plot_data,
                       family = genpois(),
                       dispformula = ~ar1(year + 0 | site/plot_id))
                       # check for nonlinear effect of time

# Parasitoid species richness

pr_fam <- glmmTMB(par_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     sc_cells : sc_fa +
                     sc_hr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = poisson()) # change family

pr_quad <- glmmTMB(par_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + I(sc_fa^2) +
                     sc_hr + sc_cells +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa +
                     sc_cells : sc_fa +
                     sc_hr : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois()) # check for nonlinear effect of time

# Parasitoid abundance

pa_fam <- glmmTMB(par_abund ~  sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells + sc_pr10 +
                     sc_logtr : sc_fa + 
                     sc_pr10 : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = nbinom2(),
                   ziformula = ~1) # change family

pa_quad <- glmmTMB(par_abund ~  sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells + sc_pr10 +
                     I(sc_fa^2) +
                     sc_logtr : sc_fa + 
                     sc_pr10 : sc_fa +
                     (1 | site/plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1) # check for nonlinear effect of time

# Parasitism rate

prt_noz <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                     sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr10 +
                     sc_logtr : sc_fa +
                     sc_pr10 : sc_fa,
                   data =  subset(plot_data, par_rich != 0),
                   family = betabinomial(link = "logit")) # remove plot-years with no parasitism (n = 71)

prt_pr10 <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                     sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr +
                     sc_logtr : sc_fa +
                     sc_pr : sc_fa,
                   data =  plot_data,
                   family = betabinomial(link = "logit")) # use "uncapped" parasitoid richness

prt_nocell <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                        sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_pr10 +
                        sc_logtr : sc_fa +
                        sc_pr10 : sc_fa,
                      data = plot_data,
                      family = betabinomial(link = "logit")) # remove brood cells as predictor

prt_quad <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                      sc_logtr + sc_fd + sc_sv + sc_fa + I(sc_fa^2) +
                      sc_cells + sc_hr + sc_pr10 +
                      sc_logtr : sc_fa +
                      sc_pr10 : sc_fa,
                    data = plot_data,
                    family = betabinomial(link = "logit")) # consider non-linear temporal dynamics

prt_quad_noz <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                          sc_logtr + sc_fd + sc_sv + sc_fa + I(sc_fa^2) +
                          sc_cells + sc_hr + sc_pr10 +
                          sc_logtr : sc_fa +
                          sc_pr10 : sc_fa,
                        data =  subset(plot_data, par_rich != 0),
                        family = betabinomial(link = "logit")) # no zeroes + non-linear time effect

prt_frich <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_fa +
                       sc_logtr:sc_fa +
                       sc_frich +
                       (1 | site/plot_id) + (1 | year),
                     data = plot_data,
                     family = betabinomial(link = "logit")) # add co-occurrence of focal parasitoids
                                                            # to basic parasitism model  

prt_frich2 <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_fa + sc_pr10 +
                        sc_cells + sc_hr +
                       sc_logtr:sc_fa +
                       sc_pr10:sc_fa +
                       sc_frich +
                       (1 | site/plot_id) + (1 | year),
                     data = plot_data,
                     family = betabinomial(link = "logit")) # both interactions sensitive to focal
                                                            # parasitoid co-occurrence

# Niche overlap
niche_frich <- glmmTMB(niche_overlap_hl ~ sc_hr + sc_pr10 + sc_cells + sc_logtr + sc_fa +
                         sc_logtr : sc_fa +
                         sc_frich +                            
                         (1 | site/plot_id) + (1 | year),
                         family = gaussian, data = plot_data)

sensitivity_models <- list(
  hr_fam       = hr_fam,
  hr_quad      = hr_quad,
  ha_fam       = ha_fam,
  ha_quad      = ha_quad,
  pr_fam       = pr_fam,
  pr_quad      = pr_quad,
  pa_fam       = pa_fam,
  pa_quad      = pa_quad,
  prt_noz      = prt_noz,
  prt_pr10     = prt_pr10,
  prt_nocell   = prt_nocell,
  prt_quad     = prt_quad,
  prt_quad_noz = prt_quad_noz,
  prt_frich    = prt_frich,
  prt_frich2   = prt_frich2,
  niche_frich  = niche_frich
)

rm(list = names(sensitivity_models))

model_summaries <- export_glmmTMB_markdown(
  models = sensitivity_models,
  file = file.path(output_dir, "robustness_checks.md"),
  digits = 3
)

model_summaries_docx <- export_glmmTMB_docx(
  models = sensitivity_models,
  file = file.path(output_dir, "robustness_checks.docx"),
  digits = 3
)

# ------------------------------------------------------------------------------
# 5. PATH ANALYSES ####
# ------------------------------------------------------------------------------
# Path analyses were used to discern potential hypothesized mechanisms of temporal dynamics.
# Based on a series of linear GLM component models with normal distributions.
# Variables are transformed, when necessary, and scaled to approximate model assumptions
# Transformations used include empirical logit, log, and square-root

emp_logit <- function(y, m) {
  qlogis((y + 0.5) / (m + 1))
} # empirical logit function for transforming parasitism rate

plot_data$sem_cprate <- as.numeric(scale(emp_logit(plot_data$par_abund, plot_data$total_cells)))
plot_data$sem_totcells <- as.numeric(scale(sqrt(plot_data$total_cells)))
plot_data$sem_pa <- as.numeric(scale(sqrt(plot_data$par_abund)))
plot_data$sem_pr <- as.numeric(scale(sqrt(plot_data$par_rich10)))
plot_data$sem_fa <- plot_data$sc_fa
plot_data$sem_hr <- plot_data$sc_hr
plot_data$sem_tr <- plot_data$sc_logtr
plot_data$sem_sv <- as.numeric(scale(sqrt(plot_data$stand_volume)))
plot_data$sem_fd <- as.numeric(scale(log1p(plot_data$tree_fd)))
plot_data$sem_frich <- as.numeric(scale(plot_data$focal_rich))
plot_data$sem_niche <- as.numeric(scale(plot_data$niche_overlap_hl))


# Path analysis 1. Output visualized in Fig. 2 #
# Component models
prt_sem <- glmmTMB(sem_cprate ~ sem_pr + sem_hr + sem_totcells +
                     sem_fa + sem_tr +
                     sem_pr : sem_fa +
                     sem_tr : sem_fa +
                     (1|site/plot_id) + (1|year),
                   data = plot_data,
                   family = "gaussian")

pr_sem <- glmmTMB(sem_pr ~ sem_hr + sem_totcells + sem_sv + 
                    sem_fa +
                    sem_totcells : sem_fa +
                    (1|site/plot_id) + (1|year),
                  data = plot_data,
                  family = "gaussian")

hr_sem <- glmmTMB(sem_hr ~ sem_totcells + sem_sv + 
                    (1|site/plot_id) + (1|year),
                  data = plot_data,
                  family = "gaussian") 

ha_sem <- glmmTMB(sem_totcells ~ sem_sv + 
                    sem_fa + sem_fd +
                    sem_sv : sem_fa +
                    (1|site/plot_id) + (1|year),
                  data = plot_data,
                  family = "gaussian")

sv_sem <- glmmTMB(sem_sv ~ sem_tr + sem_fa +
                    sem_tr : sem_fa +
                    (1|site/plot_id) + (1|year),
                  data = plot_data,
                  family = "gaussian")

fd_sem <- glmmTMB(sem_fd ~ sem_tr +
                    (1|site),
                  data = plot_data,
                  family = "gaussian")

sem1_mods <- list(
  prt_sem = prt_sem,
  pr_sem  = pr_sem,
  hr_sem  = hr_sem,
  ha_sem = ha_sem,
  sv_sem  = sv_sem,
  fd_sem = fd_sem,
  sem_sv %~~% sem_cprate,
  sem_fd %~~% sem_sv
)

sem1_all <- do.call(piecewiseSEM::psem, c(sem1_mods, list(data = plot_data)))

summary(sem1_all)
sem1_coef <- as.data.frame(coefs(sem1_all))

# Path analysis 2. Output visualized in Fig. 3 #

sem_vars <- c("sem_cprate", "sem_niche", "sem_pr", 
              "sem_hr", "sem_totcells", "sem_sv", 
              "sem_fa", "sem_tr", "site", "plot_id", "year", "sem_frich") # select variables

sem_data <- plot_data %>%
  dplyr::select(all_of(sem_vars)) %>%
  na.omit() # new dataframe only for plot-years with valid networks

nrow(sem_data) # n = 448

# Component models
prt_sem2 <- glmmTMB(sem_cprate ~ sem_niche + sem_frich + sem_pr + sem_hr +
                      sem_totcells +
                      (1 | site/plot_id) + (1 | year),
                    data = sem_data, family = gaussian())

niche_sem2 <- glmmTMB(sem_niche ~ sem_frich + sem_hr +
                        sem_totcells + sem_pr +
                        (1 | site/plot_id) + (1 | year),
                      data = sem_data, family = gaussian())

frich_sem2 <- glmmTMB(sem_frich ~ sem_pr + sem_hr + sem_tr + sem_fa +
                        sem_totcells +
                        sem_pr : sem_fa +
                        sem_tr : sem_fa +
                        (1 | site/plot_id) + (1 | year),
                      data = sem_data, family = gaussian())

pr_sem2 <- glmmTMB(sem_pr ~ sem_tr + sem_fa +
                     sem_totcells + sem_hr +
                     (1 | site/plot_id) + (1 | year),
                   data = sem_data, family = gaussian())

sem2_mods <- list(
  prt_sem2 = prt_sem2,
  niche_sem2  = niche_sem2,
  frich_sem2 = frich_sem2,
  pr_sem2 = pr_sem2
)

sem2_all <- do.call(piecewiseSEM::psem, c(sem2_mods, list(data = sem_data)))
summary(sem2_all)



#clean global environment from SEM sub-model objects
rm_path_objects(sem1_mods, sem2_mods)

sems <- list(
  sem1_all = sem1_all,
  sem2_all = sem2_all
) 

sem_md <- export_piecewiseSEM_markdown(
  sems = sems,
  file = file.path(output_dir, "sem_results.md"),
  digits = 3
) # export results from path analyses as .md file

sem_docx <- export_piecewiseSEM_docx(
  sems = sems,
  file = file.path(output_dir, "sem_results.docx"),
  digits = 3
) # export as ,docx file



# Supplementary models #
# Support results from the path analysis

# Probability of occurrence of a focal parasitoid species
# Namely  Chrysis principalis, Amobia auriceps, Lycogaster violaceipennis, and Melittobia sosui
# focal_tree <- glmmTMB(cbind(focal_rich, 4 - focal_rich) ~ sc_logtr + sc_fa +
#                       sc_logtr : sc_fa +
#                       (1 | site/plot_id) + (1 | year),
#                     data = plot_data,
#                     family = betabinomial())
# 
# # Quantitative parasitoid community niche overlap
# niche_mod <- glmmTMB(niche_overlap_hl ~ sc_logtr*sc_fa +
#                         sc_cells + sc_hr + sc_pr10 +
#                      #  sc_frich +                             # add to test mediation
#                         (1 | site/plot_id) + (1 | year),
#                       family = gaussian, data = plot_data)


# ------------------------------------------------------------------------------
# 5. Anova Type I based on GLMMs ####
# ------------------------------------------------------------------------------

## Hosts ##
# Host richness full
HR_M <- glmmTMB(host_rich ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(HR_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)
Model_2 <- update(Model_1, . ~ . +  sc_slope + sc_elev + sc_east + sc_north) 
Model_3 <- update(Model_2, . ~ . + sc_cells) 
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)

make_named_list <- function(...) {
  lst <- list(...)
  names(lst) <- as.character(substitute(list(...)))[-1]
  lst
} # function

this <- anova(HR_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(HR_M, Model_1, Model_2, Model_3, Model_4, Model_5)
hr_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Host richness partial
HR_M <- glmmTMB(host_rich ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(HR_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . + sc_fa)
Model_3 <- update(Model_2, . ~ . +  sc_logtr : sc_fa) 
Model_4 <- update(Model_3, . ~ . + sc_sv) 
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)

this <- anova(HR_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(HR_M, Model_1, Model_2, Model_3, Model_4, Model_5)
hr_partial_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Host abundance full
HA_M <- glmmTMB(total_cells ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(), 
                dispformula = ~ar1(year + 0 | site/plot_id))

Model_1 <- update(HA_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)
Model_2 <- update(Model_1, . ~ . +  sc_slope + sc_elev + sc_east + sc_north) 
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa)
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)

this <- anova(HA_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(HA_M, Model_1, Model_2, Model_3, Model_4)
ha_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Host abundance partial
HA_M <- glmmTMB(total_cells ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(), 
                dispformula = ~ar1(year + 0 | site/plot_id))

Model_1 <- update(HA_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . + sc_fa)
Model_3 <- update(Model_2, . ~ . +  sc_logtr : sc_fa)
Model_4 <- update(Model_3, . ~ . + sc_sv)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)

this <- anova(HA_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(HA_M, Model_1, Model_2, Model_3, Model_4, Model_5)
ha_partial_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

## Parasitoids ##
# Parasitoid richness full
ER_M <- glmmTMB(par_rich ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(ER_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                      
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr) 
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa)
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)
Model_5 <- update(Model_4, . ~ . + sc_hr : sc_fa)
Model_6 <- update(Model_5, . ~ . + sc_cells : sc_fa)

this <- anova(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6, test = "Chisq")
mods <- make_named_list(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6)
pr_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Parasitoid richness partial
ER_M <- glmmTMB(par_rich ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(ER_M, . ~ . + sc_logtr)                      
Model_2 <- update(Model_1, . ~ . + sc_fa)                      
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa) 
Model_4 <- update(Model_3, . ~ . + sc_sv)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)

this <- anova(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5)
pr_partial_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

## Parasitism ##
# Total parasitized cells full
EA_M <- glmmTMB(par_abund ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(),  ziformula = ~1) # ziform. only sig. in Model_1
                
Model_1 <- update(EA_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                      
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr)
Model_3 <- update(Model_2, . ~ . + sc_pr10)
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)
Model_6 <- update(Model_5, . ~ . + sc_pr10 : sc_fa)

this <- anova(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6, test = "Chisq")
mods <- make_named_list(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6)
pa_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Total parasitized cells partial
EA_M <- glmmTMB(par_abund ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(),  ziformula = ~1)

Model_1 <- update(EA_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . + sc_fa)
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa) 
Model_4 <- update(Model_3, . ~ . + sc_sv)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)                                 

this <- anova(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5)
pa_partial_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

#Parasitism full
PRT_M <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_data,
                 family = betabinomial(link = "logit"))

Model_1 <- update(PRT_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                       
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr + sc_pr10)
Model_3 <- update(Model_2, . ~ . + sc_cells : sc_hr)
Model_4 <- update(Model_3, . ~ . + sc_hr : sc_fa)
Model_5 <- update(Model_4, . ~ . + sc_cells : sc_fa)
Model_6 <- update(Model_5, . ~ . + sc_pr10 : sc_fa)
Model_7 <- update(Model_6, . ~ . + sc_pr10 : sc_sv)

this <- anova(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6, Model_7, test = "Chisq")
mods <- make_named_list(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6, Model_7)
prt_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Parasitism partial
PRT_M <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_data,
                 family = betabinomial(link = "logit"))

Model_1 <- update(PRT_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . + sc_fa)
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa) 
Model_4 <- update(Model_3, . ~ . + sc_sv)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)  

this <- anova(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5)
prt_partial_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Parasitism partial with focal-parasitoid co-occurrence as covariable
PRT_M <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_data,
                 family = betabinomial(link = "logit"))

Model_1 <- update(PRT_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . + sc_fa)
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa) 
Model_4 <- update(Model_3, . ~ . + sc_sv)
Model_5 <- update(Model_4, . ~ . + sc_sv : sc_fa)  
Model_6 <- update(Model_5, . ~ . + sc_pr10)
Model_7 <- update(Model_6, . ~ . + sc_pr10 : sc_fa)
Model_8 <- update(Model_7, . ~ . + sc_frich)

this <- anova(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6,
              Model_7, Model_8, test = "Chisq")
mods <- make_named_list(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, Model_6, 
                        Model_7, Model_8)
prt_frich_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)


# Interaction-network parasitism
# Check if focal-parasiotid co-occurrence and niche overlap independently improve the model
# Only keep plot-years with minimum 2 species per trophic level and 2 interactions (n = 448)
plot_dat_net <- plot_data[!is.na(plot_data$network_size), ]

PRT_M2 <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_dat_net,
                 family = betabinomial(link = "logit"))

Model_1 <- update(PRT_M2, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                       
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr + sc_pr10)
Model_3 <- update(Model_2, . ~ . + sc_frich)
Model_4 <- update(Model_3, . ~ . + sc_niche)

this <- anova(PRT_M2, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(PRT_M2, Model_1, Model_2, Model_3, Model_4)
net_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

## Network structure ##
# Niche overlap
NO_M <- glmmTMB(niche_overlap_hl ~ 1 + (1 | site/plot_id) + (1 | year),
                data = plot_dat_net, family = gaussian)

Model_1 <- update(NO_M, . ~ . + sc_hr + sc_pr10)
Model_2 <- update(Model_1, . ~ . + sc_cells)
Model_3 <- update(Model_2, . ~ . + sc_logtr + sc_fa)  
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa)  
Model_5 <- update(Model_4, . ~ . + sc_frich)  

this <- anova(NO_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(NO_M, Model_1, Model_2, Model_3, Model_4, Model_5)
niche_tbl <- make_type1_table(mods, this, pretty_map, digits = 3)

# Export results
anova_tables <- list(
  hr_tbl          = hr_tbl,
  hr_partial_tbl  = hr_partial_tbl,
  ha_tbl          = ha_tbl,
  ha_partial_tbl  = ha_partial_tbl,
  pr_tbl          = pr_tbl,
  pr_partial_tbl  = pr_partial_tbl,
  pa_tbl          = pa_tbl,
  pa_partial_tbl  = pa_partial_tbl,
  prt_tbl         = prt_tbl,
  prt_partial_tbl = prt_partial_tbl,
  net_tbl         = net_tbl,
  niche_tbl       = niche_tbl
  )

rm(list = names(anova_tables))
rm(list = ls(pattern = "^Model_")) # clean environment

anova_ns <- c(
  hr_tbl = nobs(HR_M),
  hr_partial_tbl = nobs(HR_M),
  ha_tbl = nobs(HA_M),
  ha_partial_tbl = nobs(HA_M),
  pr_tbl = nobs(ER_M),
  pr_partial_tbl = nobs(ER_M),
  pa_tbl = nobs(EA_M),
  pa_partial_tbl = nobs(EA_M),
  prt_tbl = nobs(PRT_M),
  prt_partial_tbl = nobs(PRT_M),
  net_tbl = nobs(PRT_M2),
  niche_tbl = nobs(NO_M)
)

anova_md <- export_table_list_markdown(
  tables = anova_tables,
  file = file.path(output_dir, "anova_type1_tables.md"),
  digits = 3,
  n = anova_ns
) # export anova tables as .md file

anova_docx <- export_table_list_docx(
  tables = anova_tables,
  file = file.path(output_dir, "anova_type1_tables.docx"),
  digits = 3,
  n = anova_ns
) # export anova tables as .docx file

message("Script ran successfully ✔️")
