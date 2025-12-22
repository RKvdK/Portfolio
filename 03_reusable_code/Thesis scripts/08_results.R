# In this script, all specified functions and objects are combined to retrieve the actual results
# The SIM, MC and RUN functions are required for this script to run adequately

# Run the simulation with 1000 Monte Carlo iterations per parameter combination

# simdat <- RUN(simgrid, MCnum = 1000) 

# Save the simulation results

# saveRDS(simdat, file = path_data("simdat.rds")) 

# Load the simulation results

simdat <- readRDS(path_data("simdat.rds")) 

# Display the first rows of the simulation results

head(simdat) 

