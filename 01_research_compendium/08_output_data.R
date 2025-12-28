# In this script, all specified functions and objects are combined to retrieve the actual results
# The SIM, MC and RUN functions are required for this script to run adequately

# Run the simulation with 1000 Monte Carlo iterations per parameter combination

mlarpis_dat <- RUN(simgrid, MCnum = MCnum) 

# Save the simulation results

saveRDS(mlarpis_dat, file = "mlarpis_dat.rds")
