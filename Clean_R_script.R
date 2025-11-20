# clearing workspace
rm(list = ls())
graphics.off()
cat("\014") 

# ------------------------------------------------------------------------------
## Header ##
# Title: Time amplifies multitrophic diversity-functioning relationships in forests
# Author: Massimo Martini
# Date: 20th November 2025

# ------------------------------------------------------------------------------

## Comments ##
# For simplicity, I leave out all diagnostic and plotting procedures

# ------------------------------------------------------------------------------
# 1. Load libraries and set working directory
# ------------------------------------------------------------------------------

# Set working directory
# setwd() 


#load
library(vegan)
library(janitor)
library(dplyr)
library(tidyverse)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(performance)
library(DHARMa)
library(emmeans)
library(piecewiseSEM)
library(patchwork)
library(broom)
library(broom.mixed)
library(kableExtra)
library(MuMIn)
library(parameters)
library(partR2)

# ------------------------------------------------------------------------------
# 2. Read and prepare data
# ------------------------------------------------------------------------------

#upload data
plot_data <- readr::read_csv("plot_data.csv")

str(plot_data)

source("final_functions.R") #upload helper functions

plot_data$tree_r[plot_data$tree_r>16] = 16 # reduce the 24-tree species richness plots to 16 tree species (2 plots per site)
plot_data$log_tr = log2(plot_data$tree_r) #log-transformation

plot_data$par_rich10 <- plot_data$par_rich
plot_data$par_rich10[plot_data$par_rich > 10] <- 10 #reduce >10 enemy species communities to 10 (n = 4) for parasitism rate models

plot_data <- plot_data %>%
  arrange(plot_id, year) #arrange dataset by plot and year

plot_data$for_age <- ifelse(plot_data$site == "A", 
                            plot_data$year - 2009, 
                            plot_data$year - 2010) #calculate stand age based on year of site planting

plot_data$year = as.character(plot_data$year) # year as categorical factor

plot_data <- plot_data %>%
  mutate(cell_prate = (par_abund / total_cells)) #calculate parasitism rate as the fraction of attacked vs total brood cells

plot_data$cell_prate[is.nan(plot_data$cell_prate)] <- 0 #change NAs to 0

plot_data <- plot_data %>%
  group_by(plot_id) %>%
  mutate(cells_lag = lag(total_cells, n = 1)) %>%
  ungroup() #add column with lagged host abundance (abundance from the year before)

plot_data$sqrt_sv <- sqrt(plot_data$stand_volume) #transform stand volume for path analyses
plot_data$sqrt_fd <- sqrt(plot_data$tree_fd) #transform tree_fd for path analyses

#create pd_climate by only keeping plot-years with climate data (n = 448, no data for year 2014, 56 plots remain)
pd_climate <- plot_data %>%
  filter(!is.na(annual_temperature))

#scale all predictors
#choose variables to scale
cols_to_scale <- c("log_tr", "stand_volume", "tree_fd", "sqrt_sv", "sqrt_fd", "for_age", "elevation", "eastness",
                   "northness", "slope","annual_temperature", "annual_humidity", "total_cells", "host_rich",
                   "par_abund", "par_rich", "par_rich10","cells_lag")

#set new abbreviated scaled column names
name_map <- c(
  log_tr       = "sc_logtr",
  stand_volume = "sc_sv",
  tree_fd      = "sc_fd",
  for_age      = "sc_fa",
  sqrt_sv      = "sc_rtsv",
  sqrt_fd      = "sc_rtfd",
  elevation    = "sc_elev",
  eastness     = "sc_east",
  northness    = "sc_north",
  slope        = "sc_slope",
  annual_temperature   = "sc_temp",
  annual_humidity      = "sc_humid",
  total_cells  = "sc_cells",
  host_rich    = "sc_hr",
  par_abund    = "sc_pa",
  par_rich     = "sc_pr",
  par_rich10   = "sc_pr10",
  cells_lag    = "sc_cellag"
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


#check correlation between variables
vars <- plot_data[, c("log_tr", "stand_volume", "tree_fd", "for_age", "elevation", "eastness",
                          "northness", "slope","annual_temperature", "annual_humidity", "total_cells", "host_rich",
                          "par_abund", "par_rich", "cells_lag")]

cortable <- cor(vars, use = "complete.obs", method = "pearson")

cortable %>% 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  filter(var1 < var2, abs(cor) > 0.60) %>% 
  arrange(desc(abs(cor))) #only enemy abundance and richness reach Pearson threshold (0.7)


# ------------------------------------------------------------------------------
# 3. GLMM analyses
# ------------------------------------------------------------------------------

####
#Parasitism
####

plot_data$p_failure <- (plot_data$total_cells - plot_data$par_abund) #calculate parasitism failure

#using sc_pr10 to limit model extrapolation at highest enemy richness levels
#parasitism full
prt_full <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                        sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr10 +
                        sc_cells : sc_fa +
                        sc_hr : sc_fa +
                        sc_pr10 : sc_fa +
                        sc_pr10 : sc_sv +
                        # sc_pr10 : sc_logtr +
                        # sc_pr10 : sc_fd +
                        sc_cells : sc_hr,
                      data = plot_data,
                      family = betabinomial(link = "logit"))

#parasitism full without zeroes (n=721)
prt_full_nz <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | site/plot_id) + (1 | year) + 
                      sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr + sc_pr10 +
                      sc_cells : sc_fa +
                      sc_hr : sc_fa +
                      sc_pr10 : sc_fa +
                      sc_pr10 : sc_sv +
                      # sc_pr10 : sc_logtr +
                      # sc_pr10 : sc_fd +
                      sc_cells : sc_hr,
                      data = subset(plot_data, par_rich != 0),
                    family = betabinomial(link = "logit"))

compare_models(prt_full, prt_full_nz) #no qualitative difference in models
#Figure 4 slope change also remains qualitatively unchanged (increases at high-richness enemy communities)
#zeroes are not driving the relationship


#parasitism full with enemy abundance (pa)
#model is under-dispersed due to enemy abundance circularity and enemy richness - abundance correlation (0.69)
#interpret with caution - only useful to compare temporal patterns of enemy richness vs enemy abundance for relative slope calculations (Figure 2)
prt_full_wpa <- glmmTMB(cbind(par_abund, p_failure) ~ (1 | year) + 
                            sc_logtr + sc_fd + sc_sv + sc_fa +
                            sc_cells + sc_hr + sc_pa + sc_pr10 +
                            sc_cells : sc_fa +
                            sc_hr : sc_fa +
                            sc_pa : sc_fa +
                            sc_pr10 : sc_fa +
                            sc_pr10 : sc_sv +
                            # sc_pr10 : sc_logtr +
                            # sc_pr10 : sc_fd +
                            sc_cells : sc_hr,
                          data = plot_data,
                          family = betabinomial(link = "logit"))

compare_models(prt_full, prt_full_wpa) #sanity check

#parasitism partial with only tree richness and time
prt_tree <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_fa +
                      sc_logtr : sc_fa +
                      (1 | site/plot_id) + (1 | year),
                    data = plot_data,
                    family = betabinomial(link = "logit"))

#parasitism partial with forest variables and forest age
prt_vol <- glmmTMB(cbind(par_abund, p_failure) ~ sc_logtr + sc_sv + sc_fa +
                     sc_logtr : sc_fa +
                     sc_sv : sc_fa + 
                     (1 | site/plot_id) + (1 | year),
                   data = plot_data,
                   family = betabinomial(link = "logit"))

#### 
#Natural enemies
#### 

#enemy richness full
pr_full <- glmmTMB(par_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_hr + sc_cells + sc_pa +
                        sc_logtr : sc_fa +
                        sc_sv : sc_fa +
                        sc_cells : sc_fa +
                        sc_hr : sc_fa +
                        sc_hr : sc_sv +
                        sc_hr : sc_cells +
                        (1 | year), 
                      data = plot_data,
                      family = genpois())

#enemy richness partial with tree rich and time
pr_tree <- glmmTMB(par_rich ~ sc_logtr + sc_fa +
                        sc_logtr : sc_fa +
                        (1 | plot_id) + (1 | year), 
                      data = plot_data,
                      family = genpois(),
                      ziformula = ~1)

#enemy richness partial with forest variables and time
pr_vol <- glmmTMB(par_rich ~ sc_logtr + sc_sv + sc_fa +
                       sc_logtr : sc_fa +
                       sc_sv : sc_fa +
                       (1 | plot_id) + (1 | year), 
                     data = plot_data,
                     family = genpois(),
                     ziformula = ~1)


#enemy abundance full
pa_full <- glmmTMB(par_abund ~  sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells + sc_hr +
                        sc_logtr : sc_fa +
                        sc_sv : sc_fa +
                        sc_cells : sc_fa +
                        sc_hr : sc_fa +
                        sc_hr : sc_sv +
                        sc_hr : sc_cells +
                        (1 | site/plot_id) + (1 | year), 
                      data = plot_data,
                      family = genpois())

#enemy abundance partial with tree richness and time
pa_tree <- glmmTMB(par_abund ~ sc_logtr + sc_fa +
                     sc_logtr : sc_fa +
                     (1 | plot_id) + (1 | year), 
                   data = plot_data,
                   family = genpois(),
                   ziformula = ~1)

#enemy richness partial with forest variables and time
pa_vol <- glmmTMB(par_abund ~ sc_logtr + sc_sv + sc_fa +
                       sc_logtr : sc_fa +
                       sc_sv : sc_fa +
                       (1 | plot_id) + (1 | year), 
                     data = plot_data,
                     family = genpois(),
                     ziformula = ~1)


#### 
#Hosts
#### 

#host richness full
hr_full<- glmmTMB(host_rich ~ sc_logtr + sc_fd + sc_sv + sc_fa + sc_cells +
                        sc_slope + sc_elev + sc_east + sc_north +
                        sc_logtr : sc_fa +
                        sc_sv : sc_fa +
                        (1 | year), 
                      data = plot_data,
                      family = genpois())

#host richness partial with tree richness and time
hr_tree <- glmmTMB(host_rich ~ sc_logtr + sc_fa +
                        sc_logtr : sc_fa +
                        (1 | year), 
                      data = plot_data,
                      family = genpois())

#host richness partial with forest variables and time
hr_vol <- glmmTMB(host_rich ~ sc_logtr + sc_sv + sc_fa +
                       sc_logtr : sc_fa +
                       sc_sv : sc_fa +
                       (1 | year), 
                     data = plot_data,
                     family = genpois())

#host abundance full
ha_full <- glmmTMB(total_cells ~ sc_logtr + sc_fd + sc_sv + sc_fa +
                        sc_slope + sc_elev + sc_east + sc_north +
                        sc_logtr : sc_fa +
                        sc_sv : sc_fa +
                        (1 | site/plot_id) + (1 | year), 
                      data = plot_data,
                      family = genpois(),
                      dispformula = ~ar1(year + 0 | site/plot_id), #add dispersion parameter for temporal autocorrelation
                      ziformula = ~1)

#host abundance partial with tree richness and time
ha_tree <- glmmTMB(total_cells ~ sc_logtr + sc_fa +
                        sc_logtr : sc_fa +
                        (1 | site/plot_id) + (1 | year), 
                      data = plot_data,
                      family = genpois(),
                      dispformula = ~ar1(year + 0 | site), #add dispersion parameter for temporal autocorrelation
                      ziformula = ~1)

#host abundance partial with forest variables and time
ha_vol <- glmmTMB(total_cells ~ sc_logtr + sc_sv + sc_fa +
                       sc_logtr : sc_fa +
                       sc_sv : sc_fa +
                       (1 | site/plot_id) + (1 | year), 
                     data = plot_data,
                     family = genpois(),
                     dispformula = ~ar1(year + 0 | site/plot_id), #add dispersion parameter for temporal autocorrelation
                     ziformula = ~1)


#put all model objects into a single list "models" 
models <- list(
  ha_full = ha_full,
  ha_vol = ha_vol,
  ha_tree = ha_tree,
  hr_full = hr_full,
  hr_vol = hr_vol,
  hr_tree = hr_tree,
  pa_full = pa_full,
  pa_vol = pa_vol,
  pa_tree = pa_tree,
  pr_full = pr_full,
  pr_vol = pr_vol,
  pr_tree = pr_tree,
  prt_full = prt_full,
  prt_full_wpa = prt_full_wpa,
  prt_vol = prt_vol,
  prt_tree = prt_tree)

rm(list = names(models))

model_summaries <- pretty_glmmTMB_summaries(models = models,
                                file   = "model_summaries.md",
                                digits = 3)

#export a table showing model diagnostics (dispersion, AIC, mean and max VIF cvalues)
disp_tbl <- imap_dfr(models, get_disp_row) %>%
  relocate(Family, Dispersion, P_value, .after = Model)

disp_tbl <- disp_tbl %>%
  mutate(
    AIC = map_dbl(Model, ~ tryCatch(AIC(models[[.x]]), error = function(e) NA_real_))
  ) #add AIC 

vif_summary <- imap_dfr(models, function(mod, name) {
  vifs <- try(performance::check_collinearity(mod), silent = TRUE)
  if (inherits(vifs, "try-error") || is.null(vifs)) {
    tibble(Model = name, Mean_VIF = NA_real_, Max_VIF = NA_real_)
  } else {
    tibble(Model = name,
           Mean_VIF = mean(vifs$VIF, na.rm = TRUE),
           Max_VIF = max(vifs$VIF, na.rm = TRUE))
  }
})

disp_tbl <- left_join(disp_tbl, vif_summary, by = "Model") #add mean and maximum predictor variance inflation factor

write.csv(disp_tbl,"model_diagnostics.csv")


# ------------------------------------------------------------------------------
# 4. Path analyses
# ------------------------------------------------------------------------------

# We used sqrt-transformed stand biomass (sc_rtsv) and tree FD (sc_rtfd) only as outcomes 
# in their respective submodels to ensure residual normality and meet model assumptions.
# Elsewhere, predictors and outcomes remain untransformed for interpretability.
# Predictors are not transformed so that coefficients remain directly interpretable 
# on their original scale throughout the SEM.

#Main SEM: All core models and relationships (untransformed outcomes)
prt_sem <- glmmTMB(cell_prate ~ 
                     sc_pa + sc_cells +
                     sc_pr10 + sc_hr + sc_fa +
                     (1|year),
                   data = plot_data,
                   family = "gaussian")

pr_sem <- glmmTMB(sc_pr10 ~ sc_pa + sc_hr + sc_cells + sc_fa + sc_logtr +
                    (1|plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")

pa_sem <- glmmTMB(sc_pa ~ sc_cells + sc_logtr + sc_sv +
                    (1|site/plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")

hr_sem <- glmmTMB(sc_hr ~ sc_cells + sc_logtr + sc_sv +
                    (1|site/plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")

ha_sem <- glmmTMB(sc_cells ~ sc_fd + sc_sv +
                     (1|plot_id)+(1|year),
                   data = plot_data,
                   family = "gaussian")

sv_sem <- glmmTMB(sc_sv ~ sc_logtr + sc_fa +
                    (1|site/plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")


fd_sem <- glmmTMB(sc_fd ~ sc_logtr +
                    (1|site),
                  data = plot_data,
                  family = "gaussian")

sem_mods <- list(
  prt_sem = prt_sem,
  pr_sem  = pr_sem,
  pa_sem  = pa_sem,
  hr_sem  = hr_sem,
  ha_sem  = ha_sem,
  sv_sem  = sv_sem,
  fd_sem  = fd_sem,
  sc_sv %~~% sc_fd,
  sc_pr10 %~~% sc_sv #correlation term to remove model artifact
)

sem_all <- do.call(piecewiseSEM::psem, c(sem_mods, list(data = plot_data)))
summary(sem_all)


# Submodel SEM: Only for forest-level sqrt-transformed outcomes, for improved diagnostics
# This is used for model summaries for stand biomass and tree FD
# whose residuals improve with normalizing transformations. Relationships for all other
# outcomes and predictors should be interpreted from the main SEM above.

sv_sem2 <- glmmTMB(sc_rtsv ~ sc_logtr + sc_fa + (1|site/plot_id)+(1|year), 
                  data = plot_data, family = "gaussian") # sqrt-transformed stand biomass

fd_sem2 <- glmmTMB(sc_rtfd ~ sc_logtr + (1|site),
                   data = plot_data, family = "gaussian") # sqrt-transformed tree FD


sem_mods2 <- list(
  sv_sem2  = sv_sem2,
  fd_sem2  = fd_sem2,
  sc_rtsv %~~% sc_rtfd
)

sem_all2 <- do.call(piecewiseSEM::psem, c(sem_mods2, list(data = plot_data)))
summary(sem_all2)


#####
# Path analysis with only host and enemy richness pathways (Figure S2)
# Useful for improved legibility and easier interpretation on bottom-up effects
# Same as above, sqrt-transformed forest sub-models are used for final figure


prt_semR <- glmmTMB(cell_prate ~ 
                     sc_pr10 + sc_hr + sc_sv + sc_fd +
                     (1|plot_id) + (1|year),
                   data = plot_data,
                   family = "gaussian")

pr_semR <- glmmTMB(sc_pr10 ~ sc_hr + sc_sv + sc_fa + sc_logtr +
                    (1|site/plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")

hr_semR <- glmmTMB(sc_hr ~ sc_logtr + sc_sv +
                    (1|site/plot_id)+(1|year), 
                  data = plot_data,
                  family = "gaussian")

sem_mods_rich <- list(
  prt_semR = prt_semR,
  pr_semR  = pr_semR,
  hr_semR  = hr_semR,
  sv_sem  = sv_sem,
  fd_sem  = fd_sem,
  sc_sv %~~% sc_fd
)

sem_rich <- do.call(piecewiseSEM::psem, c(sem_mods_rich, list(data = plot_data)))
summary(sem_rich)


#####
# Path analysis with microclimate data (Figure S6)
# This uses a subset of 56 plots and years 2015 - 2023, n = 448
# Used to determine whether microclimate mediates effects of stand biomass on hosts and natural enemies
# Results show only partial mediation - we can assume that remaining direct effects of stand biomass are mediated by
# structural complexity and resource availability

prt_clim <- glmmTMB(cell_prate ~ sc_pr10 + sc_pa + sc_hr + sc_cells + sc_fa +
                      sc_temp +
                      (1|year),
                    data = pd_climate,
                    family = "gaussian")

pr_clim <- glmmTMB(sc_pr10 ~ sc_pa + sc_hr +
                     sc_humid +
                     (1|year), 
                   data = pd_climate,
                   family = "gaussian")

pa_clim <- glmmTMB(sc_pa ~ sc_hr + sc_cells + sc_sv +
                     (1|site/plot_id)+(1|year), 
                   data = pd_climate,
                   family = "gaussian")

hr_clim <- glmmTMB(sc_hr ~ sc_cells + sc_temp +
                     (1|site/plot_id)+(1|year), 
                   data = pd_climate,
                   family = "gaussian")

ha_clim <- glmmTMB(sc_cells ~ sc_sv + sc_humid +
                     (1|site/plot_id)+(1|year), 
                   data = pd_climate,
                   family = "gaussian")

temp_clim <- glmmTMB(sc_temp ~ sc_sv + sc_fd +
                       (1|plot_id)+(1|year), 
                     data = pd_climate,
                     family = "gaussian")

humid_clim <- glmmTMB(sc_humid ~ sc_sv +
                        (1|site/plot_id)+(1|year), 
                      data = pd_climate,
                      family = "gaussian")

sv_clim <- glmmTMB(sc_sv ~ sc_logtr + sc_fa +
                     (1|site/plot_id), 
                   data = pd_climate,
                   family = "gaussian")

fd_clim <- glmmTMB(sc_fd ~ sc_logtr +
                     (1|site), 
                   data = pd_climate,
                   family = "gaussian")

clim_mods <- list(
  prt_clim = prt_clim,
  pr_clim  = pr_clim,
  pa_clim  = pa_clim,
  hr_clim  = hr_clim,
  ha_clim  = ha_clim,
  temp_clim = temp_clim,
  humid_clim = humid_clim,
  sv_clim  = sv_clim,
  fd_clim  = fd_clim,
  sc_temp %~~% sc_humid,
  sc_fd %~~% sc_sv,
  sc_temp %~~% sc_fa, # correlation term to remove model artifact
  sc_humid %~~% sc_fa # correlation term to remove model artifact
)

sem_clim <- do.call(piecewiseSEM::psem, c(clim_mods, list(data = pd_climate)))
summary(sem_clim)

#Submodel climate SEM: Only for forest-level sqrt-transformed outcomes

sv_clim2 <- glmmTMB(sc_rtsv ~ sc_logtr + sc_fa +
                                        (1|site/plot_id), 
                                      data = pd_climate,
                                      family = "gaussian")

fd_clim2 <- glmmTMB(sc_rtfd ~ sc_logtr +
                      (1|site), 
                    data = pd_climate,
                    family = "gaussian")

clim_mods2 <- list(
  sv_clim2  = sv_clim2,
  fd_clim2  = fd_clim2,
  sc_rtsv %~~% sc_rtfd
)

sem_clim2 <- do.call(piecewiseSEM::psem, c(clim_mods2, list(data = pd_climate)))
summary(sem_clim2)


#exporting standardized coefficients from all path analyses as a single .csv
sems <- list(
  sem_all = sem_all,
  sem_all2 = sem_all2,
  sem_rich = sem_rich,
  sem_clim = sem_clim,
  sem_clim2 = sem_clim2
)

#clean global environment from SEM sub-model objects
rm_path_objects(sem_mods, sem_mods2, sem_mods_rich, clim_mods, clim_mods2)


#extract standardized coefficients and export as .csv file
coef_combined <- map_df(names(sems), function(nm) {
  
  df <- coefs(sems[[nm]], standardize = "scale")
  names(df)[names(df) == "" | is.na(names(df))] <- "sig."
  df %>%
    as_tibble() %>%       
    mutate(SEM = nm)
})


write.csv(coef_combined, "all_sem_coefficients.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# 5. Anova Type I based on GLMMs
# ------------------------------------------------------------------------------

#Parasitism full
PRT_M <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_data,
                 family = betabinomial(link = "logit"))

# add terms in your chosen order
Model_1 <- update(PRT_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                       # forest 
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr + sc_pr10)                           # host-enemy 
Model_3 <- update(Model_2, . ~ . + sc_cells : sc_hr)                                     # interdependence  
Model_4 <- update(Model_3, . ~ . + sc_cells : sc_fa + sc_hr : sc_fa + sc_pr10 : sc_fa)   # temporal dynamics
Model_5 <- update(Model_4, . ~ . + sc_pr10 : sc_sv)                                      # structural complexity

# Type I (sequential) LR tests:
this <- anova(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")

make_named_list <- function(...) {
  lst <- list(...)
  names(lst) <- as.character(substitute(list(...)))[-1]
  lst
}

mods <- make_named_list(PRT_M, Model_1, Model_2, Model_3, Model_4, Model_5)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl,
          file = "PRT_full_anova.csv",
          row.names = FALSE,
          na = "")

#Parasitism partial
PRT_M <- glmmTMB(cbind(par_abund, p_failure) ~ 1 + (1 | site/plot_id) + (1 | year),
                 data = plot_data,
                 family = betabinomial(link = "logit"))

# add terms in your chosen order
Model_1 <- update(PRT_M, . ~ . + sc_logtr)                     
Model_2 <- update(Model_1, . ~ . + sc_logtr : sc_fa)                           
Model_3 <- update(Model_2, . ~ . + sc_sv)                                      
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)  

this <- anova(PRT_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(PRT_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "PRT_partial_anova.csv",  row.names = FALSE, na = "")


#Natural enemy richness full
ER_M <- glmmTMB(par_rich ~ 1 + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(ER_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                      
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr + sc_pa) 
Model_3 <- update(Model_2, . ~ . + sc_cells : sc_hr)
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa + sc_sv : sc_fa + sc_cells : sc_fa + sc_hr : sc_fa)                                   
Model_5 <- update(Model_4, . ~ . + sc_hr : sc_sv)

this <- anova(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(ER_M, Model_1, Model_2, Model_3, Model_4, Model_5)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "ER_full_anova.csv",  row.names = FALSE, na = "")


#Natural enemy richness partial
ER_M <- glmmTMB(par_rich ~ 1 + (1 | plot_id) + (1 | year), , data = plot_data, family = genpois())

Model_1 <- update(ER_M, . ~ . + sc_logtr)                      
Model_2 <- update(Model_1, . ~ . + sc_logtr : sc_fa) 
Model_3 <- update(Model_2, . ~ . + sc_sv)
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)

this <- anova(ER_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(ER_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "ER_partial_anova.csv",  row.names = FALSE, na = "")

#Natural enemy abundance full
EA_M <- glmmTMB(par_abund ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(EA_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)                      
Model_2 <- update(Model_1, . ~ . + sc_cells + sc_hr) 
Model_3 <- update(Model_2, . ~ . + sc_cells : sc_hr)
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa + sc_sv : sc_fa + sc_cells : sc_fa + sc_hr : sc_fa)                                 
Model_5 <- update(Model_4, . ~ . + sc_hr : sc_sv)

this <- anova(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5, test = "Chisq")
mods <- make_named_list(EA_M, Model_1, Model_2, Model_3, Model_4, Model_5)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "EA_full_anova.csv",  row.names = FALSE, na = "")

#Natural enemy abundance partial
EA_M <- glmmTMB(par_abund ~ 1 + (1 | plot_id) + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(EA_M, . ~ . + sc_logtr)                      
Model_2 <- update(Model_1, . ~ . + sc_logtr : sc_fa) 
Model_3 <- update(Model_2, . ~ . + sc_sv)
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)                                 

this <- anova(EA_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(EA_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "EA_partial_anova.csv",  row.names = FALSE, na = "")


#Host richness full
HR_M <- glmmTMB(host_rich ~ 1 + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(HR_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)
Model_2 <- update(Model_1, . ~ . +  sc_slope + sc_elev + sc_east + sc_north) 
Model_3 <- update(Model_2, . ~ . + sc_cells) 
Model_4 <- update(Model_3, . ~ . + sc_logtr : sc_fa + sc_sv : sc_fa)

this <- anova(HR_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(HR_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "HR_full_anova.csv",  row.names = FALSE, na = "")

#Host richness partial
HR_M <- glmmTMB(host_rich ~ 1 + (1 | year), data = plot_data, family = genpois())

Model_1 <- update(HR_M, . ~ . + sc_logtr)
Model_2 <- update(Model_1, . ~ . +  sc_logtr : sc_fa) 
Model_3 <- update(Model_2, . ~ . + sc_sv) 
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)

this <- anova(HR_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(HR_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "HR_partial_anova.csv",  row.names = FALSE, na = "")

#Host abundance full
HA_M <- glmmTMB(total_cells ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(), 
                dispformula = ~ar1(year + 0 | site/plot_id), 
                ziformula = ~1)

Model_1 <- update(HA_M, . ~ . + sc_logtr + sc_fd + sc_sv + sc_fa)
Model_2 <- update(Model_1, . ~ . +  sc_slope + sc_elev + sc_east + sc_north) 
Model_3 <- update(Model_2, . ~ . + sc_logtr : sc_fa + sc_sv : sc_fa)

this <- anova(HA_M, Model_1, Model_2, Model_3, test = "Chisq")
mods <- make_named_list(HA_M, Model_1, Model_2, Model_3)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "HA_full_anova.csv",  row.names = FALSE, na = "")

#Host abundance partial
HA_M <- glmmTMB(total_cells ~ 1 + (1 | site/plot_id) + (1 | year), data = plot_data, 
                family = genpois(), 
                ziformula = ~1)

Model_1 <- update(HA_M, . ~ . + sc_logtr, dispformula = ~ar1(year + 0 | site/plot_id))
Model_2 <- update(Model_1, . ~ . +  sc_logtr : sc_fa, dispformula = ~ 1 )                #Model_2 does not converge with temp. autocorr.
Model_3 <- update(Model_2, . ~ . + sc_sv, dispformula = ~ar1(year + 0 | site/plot_id))
Model_4 <- update(Model_3, . ~ . + sc_sv : sc_fa)

this <- anova(HA_M, Model_1, Model_2, Model_3, Model_4, test = "Chisq")
mods <- make_named_list(HA_M, Model_1, Model_2, Model_3, Model_4)
tbl <- make_type1_table(mods, this, pretty_map, digits = 3)
write.csv(tbl, file = "HA_partial_anova.csv",  row.names = FALSE, na = "")

message("Script ran successfully ✔️")
