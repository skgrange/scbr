# **scbr**

## Introduction

**scbr** is an R package designed for O~2~-CO~2~ atmospheric partitioning and global carbon budget calculations. It provides a robust, reproducible framework to allocate carbon and oxygen changes from multi-site atmospheric monitoring networks across three primary reservoirs: the atmospheric reservoir, the net land sink, and the ocean sink, utilizing stoichiometric mass-balance principles. 

## Installation

To install the development version from GitHub, ensure you have the [**remotes**](https://github.com/r-lib/remotes) package installed:

```R
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install scbr
remotes::install_github("skgrange/scbr")
```

## Quick start

To get going, here is a code snippet: 

```R
# Load packages
library(scbr)
library(dplyr)

# Apply the O2-CO2 partitioning method
budget_results <- apply_o2_co2_partitioning(
  df = data_obs_wide,
  conversion_factor_gigatonnes_carbon_to_ppm = 2.124, 
  conversion_factor_gigatonnes_oxygen_to_ppm = 4.76, 
  oxidation_ratio_fossil_fuels = 1.37,
  oxidation_ratio_biosphere = 1.1,
  o2_ocean_outgasing_year = 0.0
)

# View partitioned data
head(budget_results)
```

