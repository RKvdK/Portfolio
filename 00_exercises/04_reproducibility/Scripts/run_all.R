# This script runs the entire simulation study by sourcing the necessary scripts

if (!requireNamespace("here", quietly = TRUE)) install.packages("here") # Install the "here" package if not already installed
library(here) # Load the package

source(here("R", "packages.R")) # Load the packages
source(here("R", "data.R")) # Read the data
source(here("Scripts", "simulation.R")) # Run the simulation and create the visualization
