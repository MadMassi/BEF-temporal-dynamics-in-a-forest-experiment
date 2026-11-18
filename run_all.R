## Header ##
# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Complete analysis workflow
# Author: Massimo Martini

# Run this script from the project root. Restoring renv first is recommended
# for exact package versions, but the workflow can use the active R library.

# Keep sourced expressions hidden for a concise console log.
# Set to TRUE to echo the code; messages and printed results appear either way.
echo_code <- FALSE

analysis_scripts <- c(
  "Design_models.R",
  "Maximal_models.R",
  "Maximal_SEM.R",
  "Parsimonious_SEM.R",
  "Network_SEM.R"
)

message(
  "For exact reproducibility, restore renv before running this workflow. ",
  "Continuing with the active R library."
)

for (script in analysis_scripts) {
  message("\nRunning ", script, " ...")
  source(script, echo = echo_code)
}
