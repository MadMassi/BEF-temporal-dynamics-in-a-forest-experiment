# Forest growth strengthens diversity effects on host-parasitoid interactions

This repository contains the plot-level data and R code used for the study
*"Forest growth strengthens diversity effects on host-parasitoid interactions"*.

The study examines cavity-nesting Hymenoptera and their natural enemies across
88 experimental forest plots at the BEF-China platform. Trap nests were sampled
in nine years between 2014 and 2023, yielding 792 plot-year records. The analyses
test how tree diversity, stand age, and stand biomass
relate to host and parasitoid diversity, parasitism, and interaction
networks.

This repository reproduces the core statistical analyses.

## Data

- `plot_data.csv`: analysis dataset with 792 rows and 40 columns. Each row is one
  plot in one sampling year.
- `plot_data_metadata.csv`: variable descriptions, units, data types, permitted
  or observed values, derivations, missing-value conventions, and non-missing
  observation counts.

The dataset includes the experimental design, forest structure and topography,
below-canopy microclimate, trap-nest responses, sampling coverage, interaction-network
metrics, and occurrence of four focal natural-enemy species. Empty CSV fields
represent missing values and are read as `NA`.

The data are stored separately on Figshare.

## Analysis files

- `run_all.R` runs the five analysis scripts in sequence.
- `Design_models.R` prepares the analysis datasets, selects the repeated-measures
  and dispersion structures, and fits the five design models.
- `Maximal_models.R` fits and validates the five maximal GLMMs and runs the
  prespecified sensitivity analyses.
- `Maximal_SEM.R` assembles the maximal parasitism and realized-parasitism
  structural equation models (SEMs).
- `Parsimonious_SEM.R` applies the recorded component reductions and evaluates
  the reduced SEMs using d-separation and Fisher's C.
- `Network_SEM.R` fits the maximal and parsimonious interaction-network SEMs,
  including focal-richness bootstraps and full-data sensitivity analyses.
- `analysis_functions.R` contains the shared functions sourced by the five
  analysis scripts; it should not be run separately.
- `renv.lock` records the package versions used for the analyses.

## Software

The analyses were run with **R 4.6.1**. From the repository root, restore the R
package environment with:

```r
install.packages("renv") # omit if renv is already installed
renv::restore(retry = FALSE)
renv::status()
```

## Running the analyses

After restoring the package environment, run the complete workflow from the
repository root with:

```r
source("run_all.R")
```

We recommend restoring `renv` before running. The script does not install or restore packages automatically. Set
`echo_code <- TRUE` near the top of `run_all.R` to print each sourced R
expression to the console; messages and printed results appear either way.

Alternatively, run individual stages in the following order:

```r
source("Design_models.R")
source("Maximal_models.R")
source("Maximal_SEM.R")
source("Parsimonious_SEM.R")
source("Network_SEM.R")
```

The scripts exchange saved files through `model_outputs/` and do not depend on
objects left in the R workspace. When run individually, they can therefore be
started in fresh R sessions.

Approximate runtimes on the development computer, a Lenovo ThinkPad E14 Gen 2 on Linux Ubuntu
were:

| Script | Approximate runtime |
|---|---:|
| `Design_models.R` | 7-8 minutes |
| `Maximal_models.R` | 18-20 minutes |
| `Maximal_SEM.R` | 1 minute |
| `Parsimonious_SEM.R` | 5 minutes |
| `Network_SEM.R`, first run | allow 4-5 hours or 30 minutes if quick_review <- TRUE |

By default, the first `Network_SEM.R` run performs two sequential
1,000-replicate bootstraps. For a quicker review of the workflow, set
`quick_review <- TRUE` near the beginning of the script; this uses 100
replicates for each bootstrap and saves all results separately under
`model_outputs/network_sem/quick_review/`. Keep the switch at `FALSE` for the
publication analysis. Each checkpoint is saved only after all requested
replicates finish, so an interrupted partial bootstrap cannot be resumed.
Completed checkpoints are validated and reused on later runs.

## Outputs and interpretation

Each script writes to its own subdirectory under `model_outputs/`. RDS files
retain the fitted models, analysis data, validation records, and session
information; CSV files provide readable result tables.

The design models provide the primary test of experimental-design effects.
The maximal models are the primary inferential reference for conditional associations. Parsimonious SEMs are
exploratory pathway summaries.
Network analyses apply only to plot-years meeting the stated network-eligibility
criteria; see `plot_data_metadata.csv` for details.

## Citation and archives

When using these data or scripts, cite both the associated article and the
archived release.

- Article citation and DOI:

Massimo Martini, Felix Fornoff, Finn Rehling et al. Forest growth strengthens diversity effects on multitrophic interactions, 04 June 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-9904409/v1]

- Release version: 3.0.0
- GitHub project URL: https://github.com/MadMassi/BEF-temporal-dynamics-in-a-forest-experiment.git
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.17661629
- Code licence: GNU General Public License v3.0 or later

## Contact

For questions, contact:

Massimo Martini  
massimo.martini@nature.uni-freiburg.de
