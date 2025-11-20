# BEF-temporal-dynamics-in-a-forest-experiment
clean R script for data analysis of the manuscript: Time amplifies multitrophic diversity-functioning relationships in forests. Files include the main R script, and a secondary script with helper functions sourced into the main analysis

# README for data "plot_data.csv" and associated R scripts

[Access this dataset on Dryad](Dataset DOI)

Data and analysis of host-natural enemy diversity from the 'BEF-China' platform. Collected from 2014-2023 using trapnests for cavity-nesting Hymenoptera and their parasitoids.


## Description of the data and file structure

Plot-year-level data on the abundance and species richness of cavity-nesting Hymenoptera and their natural enemies (parasites, parasitoids, kleptoparasites).
Includes plot-level environmental variables such as tree richness (the only manipulated variable in this large-scale forest experiment), stand volume, tree functional diversity, mean annual temperature, etc. As well as response variables of host (solitary bee and wasp) abundance and richness, and natural enemy abundance and richness. The ecosystem function of parasitism rate can be calculated and analysed as a response variable as well.

Main data file is "plot_data.csv"
Metadata file is "plot_data_metadata.csv"
Main data analysis R script is "Clean_R_script.R"
Secondary R script with helper functions is "final_functions.R". This secondary script is sourced into the main script for data analysis.

Data is stored in Dryad, scripts are stored on Github + Zenodo

For questions please contact:
Massimo Martini - massimo.martini@nature.uni-freiburg.de


## Code/Software

All code was written and runs in R. 
R version 4.5.1 was used for the data analysis.
