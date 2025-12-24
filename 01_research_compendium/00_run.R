# This script serves as master script to run all other scripts in the project

# Set the working directory to this script location
# In this way, relative paths work correctly and is reproducible

setwd(dirname(normalizePath(sys.frame(1)$ofile))) 

source("01_packages.R")
source("02_utils.R")
source("03_parameter_grid.R")
source("04_selectivity_patterns.R")
source("05_generate_population.R")
source("06_monte_carlo.R")
source("07_run_grid.R")
source("08_output_data.R")
