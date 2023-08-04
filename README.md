This repository contains the code used for creating all the results in the paper "Fast spatial simulation of extreme 
high-resolution radar precipitation data using INLA", available as a preprint on
[https://arxiv.org/abs/2307.11390](https://arxiv.org/abs/2307.11390).

# Setup

Start by calling 
```r
renv::load()
```
and follow the instructions that appear on the screen to download/install the correct libraries in
order to ensure reproducible results. You might have to install `renv` first, if you have not
already done so. Then you need to call
```r
devtools::install()
```
to install the `R` package. Then,
you can load the package using
```r
library(extremePrecipitation)
```
Finally, you need to call
```r
make_cgeneric("all")
```
in order to compile and link the necessary `cgeneric` models. Before calling this function, you might
need to change the C compiler `GCC` (by default equal to `gcc-12`) in the makefile
`cgeneric/Makefile` to a compiler that is available on your system.
In order to run the script for downloading all data from the online repositories you must first
install the Climate Data Operators (CDO) program. This is freely available at
[https://code.mpimet.mpg.de/projects/cdo](https://code.mpimet.mpg.de/projects/cdo).

# Running the code

All scripts for reproducing the results found in the paper are available in the `exec/`
folder. These scripts are numbered in the order they should be run. An overview follows:

- `1-download-data.R`  
  This script downloads everything needed for modelling the extreme radar precipitation data.
  The results are saved as `raw_data/downloads/dem.tif` and `raw_data/downloads/radar.nc`.
- `2-process-data.R`  
  This script processes the data from `download_data.R` to create the file
  `raw_data/downloads/radar.rds`.
- `3-data_plotting.R`  
  This script creates some nice plots based on the downloaded data, and saves them in
  the folder `raw_data/images/`.
- `4-gamma_marginals.R`  
  This script performs modelling of the bulk of the marginal distributions of nonzero precipitation
  intensity over our spatial domain, using a latent Gaussian model with a gamma likelihood.
  The results are saved in the file `raw_data/results/gamma_model.rds`,
  and some plots based on the exploratory analysis and model evaluation are saved in the folder
  `raw_data/images/`
- `5-gp_marginals.R`  
  This script performs modelling of the tails of the marginal distributions of nonzero precipitation
  intensity over our spatial domain using a latent Gaussian model with a generalised Pareto likelihood.
  The results are saved in the file `raw_data/results/gp_model.rds`,
  and some plots based on the exploratory analysis and model evaluation are saved in the folder
  `raw_data/images/`.
- `6-intensity_process.R`  
  This script performs modelling of the extremal dependence structure of the nonzero precipitation intensities
  using the spatial conditional extremes model.
  The results are saved in the file `raw_data/results/intensity_process.rds`,
  and some plots based on the exploratory analysis and model evaluation are saved in the folder
  `raw_data/images/`.
- `7-occurrence_process.R`  
  This script performs modelling of the precipitation occurrence process, using four different competing models.
  The results are saved in the file `raw_data/results/occurrence_process.rds`,
  and some plots based on the exploratory analysis and model evaluation are saved in the folder
  `raw_data/images/`.
- `8-create_simulations.R`  
  This script combines the model fits from the four previous scripts, and use them for simulating
  extreme precipitation over the Stordalselva catchment. Plots are created for describing the simulations
  and for comparing them to each other and to the observed data. These are saved in the folder
  `raw_data/images/`.