# This script serves as master script to run all other scripts in the project

# Set the working directory to this script location
# In this way, relative paths work correctly and is reproducible

setwd(dirname(normalizePath(sys.frame(1)$ofile))) 

# In my thesis project, I intent to perform 1000 Monte Carlo iterations per parameter combination
# However, for this course, I lowered this number to 10 to reduce computation time
# Now the computation takes less than 30 seconds (instead of multiple days)

MCnum <- 10

source("01_packages.R")
source("02_utils.R")
source("03_parameter_grid.R")
source("04_selectivity_patterns.R")

# In my thesis project, I intend to analyze all possible parameter combinations
# However, again, I reduced the number of combinations here to reduce computation time

simgrid <- simgrid[1:5, ]  # Select first 5 parameter combinations only

source("05_generate_population.R")
source("06_monte_carlo.R")
source("07_run_grid.R")
source("08_output_data.R")
source("09_visualization_data.R")
