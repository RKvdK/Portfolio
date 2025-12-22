# This code defines helper functions to create file paths for different project directories

# In case the "here" package is not installed yet
# install.packages("here") 

# Load the here package

library(here) 

# Build function to operationalize the here function

path_R <- function(...) here::here("R", ...) 

# Script to run the entire simulation pipeline 

# Load path definitions

source(path_R("00_paths.R")) 

# Load required packages

source(path_scripts("01_packages.R"))

# Load utility functions

source(path_R("02_utils.R")) 

# Load parameter grid definitions

source(path_R("03_parametergrid.R")) 

# Load SIM function to generate populations

source(path_R("04_generate_population.R")) 

# Load MC function to perform Monte Carlo simulations

source(path_R("05_perform_mc.R")) 

# Load RUN function to run simulations

source(path_R("06_run_simulation.R")) 

# Load results processing functions

source(path_scripts("07_results.R")) 

# Load visualization functions

source(path_scripts("08_visualizations.R")) 
