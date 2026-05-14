# README for data "plot_data.csv" and associated R scripts

This repository contains the R analysis code for plot-year-level data on host-parasitoid diversity from the BEF-China platform. The data were collected from 2014 to 2023 using trap nests for cavity-nesting Hymenoptera and their parasitoids.

This is version 2 of the analysis code.

## Data

The main analysis uses `plot_data.csv`, a plot-year-level dataset containing host and parasitoid abundance and richness, parasitism-rate variables, and plot-level environmental predictors such as tree richness, stand volume, tree functional diversity, stand age, and climate variables.

The data are stored separately on Figshare:

DOI:
URL:

Included with the dataset release, `plot_data_metadata.csv` describes the variables in `plot_data.csv`.

## Files

- `Clean_R_script.R`: main analysis script
- `final_functions.R`: helper functions sourced by `Clean_R_script.R`
- `renv.lock`: package-version lockfile for reproducibility
- `.Rprofile` and `renv/activate.R`: files used by `renv` to activate the project environment

Data are stored on Figshare, while the analysis scripts are stored on GitHub and archived on Zenodo.

## Software

The analysis was run in R version 4.5.1.

Package versions are recorded in `renv.lock`. 
After downloading or cloning the repository, restore the package environment with:

```r
renv::restore()
```
Note: This project was tested with R 4.5.1. For the most reliable restore, 
use R 4.5.1 before running `renv::restore()`.

Then run the analysis with:

```r
source("Clean_R_script.R", echo = TRUE)
```

## Code Archive

The analysis code is archived on Zenodo:

DOI:
URL: 

## Contact

For questions, please contact:

Massimo Martini  
massimo.martini@nature.uni-freiburg.de
